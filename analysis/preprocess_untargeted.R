library(tidyverse)
library(broom)
library(magrittr)
library(made4)
library(ggridges)
library(modelr)
library(splines)
library(quantreg)
library(gbm3)
library(patchwork)

positive_mode_data <- read_csv("~/course/ND_Metabolomics/data/got-ms/measurement_ESI+.csv")

negative_mode_data <- read_csv("~/course/ND_Metabolomics/data/got-ms/measurement_ESI-.csv")

format_data <- function(df, type="Unknown") {



    meta_data_columns <- c("Compound", "Neutral mass (Da)", "m/z",
                                  "Charge", "Retention time (min)",
                                  "Chromatographic peak width (min)",
                                  "Identifications", "Isotope Distribution",
                                  "Maximum Abundance", "Minimum CV%",
                                  "Accepted Compound ID",
                                  "Accepted Description",
                                  "Adducts", "Compound Link",
                                  "Formula", "Score", "Fragmentation Score",
                                  "Isotope Similarity", "Mass Error (ppm)",
                                  "Retention Time Error (mins)",
                                  "Compound Link")

    cnms <- setdiff(colnames(df), meta_data_columns)
    RunIndex <- 1:length(cnms)
    names(RunIndex) <- cnms
    
    df <- df %>% gather(key = Name, value = Abundance,
                        -one_of(meta_data_columns)) 

    df <- df %>% mutate(Abundance = ifelse(Abundance == 0, NA, Abundance))
    df <- df %>% mutate(Abundance = log(as.numeric(Abundance)))

    df %<>% mutate(Mode = ifelse(grepl("QC", Name),
                                 sprintf("QC_%s", type),
                                 type))

    df %<>% mutate(SampleId = as.numeric(gsub(toupper(type), "", Name)))
    df %<>% mutate(RunIndex = recode(Name, !!!RunIndex))
    
    df %<>% rename(Metabolite = Compound)
    
}

negative_mode_data <- format_data(negative_mode_data, "neg")
positive_mode_data <- format_data(positive_mode_data, "pos")

## combine data into one long table
csf_long <- bind_rows(negative_mode_data,
                      positive_mode_data)

## Split QC and sample data
QC_long <- csf_long %>% filter(grepl("QC", Mode))

nd_data <- csf_long %>% filter(!grepl("QC", Mode))

## 0 abundace peps are NA for now
nd_data$Abundance[nd_data$Abundance == -Inf] <- NA
nd_data$Abundance[nd_data$Abundance == 0] <- NA

## Load Tracking data
tracking_data <- read_csv("~/course/ND_Metabolomics/data/NDTracking.csv")

## incorporate supplementary PD info provided by Cyrus

pd_supplement <- read_csv("../data/PANUC/panuc.csv")
pd_supplement <- pd_supplement %>%
    dplyr::select(summary_id, ApoE, GBAStatus, GBA_T369M, cognitive_status) %>%
    rename(Id = summary_id)
pd_supplement <- pd_supplement %>%
    rename(APOE = ApoE) %>%
    mutate(APOE = sub(".", "", APOE, fixed=TRUE))

tracking_data <- tracking_data %>%
    mutate(Id = gsub("IPD.+;\\s", "", Id)) %>%
    mutate(APOE = as.character(APOE))

## join parkinsons supplement with tracking
tracking_data <- left_join(tracking_data, pd_supplement, by="Id") %>%
    mutate(APOE = coalesce(APOE.x, APOE.y)) %>%
    dplyr::select(-c(APOE.x, APOE.y))

## Join Tracking data and metabolomics
subject_data <- inner_join(tracking_data, nd_data, by = c("Index" = "SampleId"))
subject_data %<>% mutate_at(vars(c("Type", "Gender", "APOE")), as.factor)
subject_data %<>% mutate_at(vars(c("Index")), as.integer)

subject_data %<>% rename(Raw = Abundance)

subject_data %<>% group_by(Mode, RunIndex) %>%
    mutate(Abundance = Raw - median(Raw, na.rm=TRUE)) %>% ungroup

fit_boosted_model <- function(df, ntrees=10000, cv.folds=10, shrinkage=0.001,
                              min_node_size=10, fig_dir=NULL) {

    print(df$Metabolite[1])
    if( sum( !is.na( df$Abundance) ) > 50 ){
        
        df$Abundance[!is.finite(df$Abundance)] <- NA
        not_na_indices <- which(!is.na(df$Abundance))
        df.omitted <- df %>% drop_na(Abundance)
        
        boost_fit_cv <- gbm(Abundance ~ RunIndex + Age + Type + Gender + APOE,
                            data = df.omitted, distribution = "laplace",
                            shrinkage = shrinkage,
                            n.trees = ntrees, cv.folds = 10,
                            n.minobsinnode = min_node_size, verbose = FALSE)

        opt_iters <- gbm.perf(boost_fit_cv, method = "cv")
        print(opt_iters)

        ## Note: this isn't ordered by run index 
        trend <- plot(boost_fit_cv, 1, num_trees=opt_iters,
                      grid_levels=df.omitted$RunIndex, return_grid=TRUE)$y
        
        abund <- df$Abundance
        trend_fit <- rep(NA, length(abund))

        abund[not_na_indices] <- df.omitted$Abundance - trend
        trend_fit[not_na_indices] <- trend

        df$RawScaled <- df$Abundance
        
        df$Abundance <- abund
        df$Trend <- trend_fit
        
        
    } else {

        df$Abundance[!is.finite(df$Abundance)] <- NA
        not_na_indices <- which(!is.na(df$Abundance))
        df.omitted <- df %>% drop_na(Abundance)
        
        med_value <- median(df.omitted$Abundance)
        trend <- rep(med_value, nrow(df.omitted))
        
        abund <- df$Abundance
        trend_fit <- rep(NA, length(abund))

        abund[not_na_indices] <- df.omitted$Abundance - trend
        trend_fit[not_na_indices] <- trend

        df$RawScaled <- df$Abundance
        
        df$Abundance <- abund
        df$Trend <- trend_fit
    }

    ## save figures
    if(!is.null(fig_dir)) {

        trend_plot <- ggplot(df[not_na_indices,]) +
            geom_point(aes(x=RunIndex, y=RawScaled, shape=Type, col=Type)) +
            geom_line(aes(x=RunIndex, y=Trend), size=2, col="red")
        residual_plot <- ggplot(df[not_na_indices,]) +
            geom_point(aes(x=RunIndex, y=Abundance, shape=Type, col=Type))
        combo_plot <- trend_plot + residual_plot
        ggsave(paste0(fig_dir, make.names(df$Metabolite[1]), ".pdf"),
               combo_plot, width=14)
    }
    
    df 

}

## boosted model 
detrended <- subject_data %>%
    group_by(Metabolite, Mode) %>%
    do(fit_boosted_model(., fig_dir = "../results/untargeted_trend_data/"))

subject_data <- detrended

## Type2 collapse (CO, CM, CY) to C
subject_data <-  subject_data %>% 
    mutate(Type2 = ifelse(Type %in% c("CY", "CM", "CO"),
                          "C",
                          as.character(.$Type))) %>%
    mutate(Type2 = factor(Type2))

## APOE NA's to "Unknown"
subject_data <- subject_data %>%
    mutate(APOE=ifelse(is.na(APOE),
                       "Unknown",
                       as.character(APOE))) %>%
    mutate(APOE=factor(APOE)) %>%
    ungroup

save(subject_data, QC_long,
     file = sprintf("preprocessed_untargeted_data-%s.RData", Sys.Date()))


