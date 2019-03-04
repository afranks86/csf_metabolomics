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

negative_mode_data <- read_csv("~/course/ND_Metabolomics/data/got-ms/NEG-AUTO.csv")
positive_mode_data_01 <- read_csv("~/course/ND_Metabolomics/data/got-ms/POS-AUTO-01.csv")
positive_mode_data_02 <- read_csv("~/course/ND_Metabolomics/data/got-ms/POS-AUTO-02.csv")
positive_mode_data_03 <- read_csv("~/course/ND_Metabolomics/data/got-ms/POS-AUTO-03.csv")

format_data <- function(df, type="Unknown") {

    df %<>%  mutate_at(vars(-one_of(c("Name", "Data File"))), funs(log(.)))
    df %<>% mutate(Mode = ifelse(grepl("QC", Name),
                                 sprintf("QC_%s", type),
                                 type))
    df$RunIndex <- 1:nrow(df)
    df %<>% mutate(SampleId = as.numeric(gsub("(Sample)|(CSF)", "", Name)))

    df
    
}

negative_mode_data <- format_data(negative_mode_data, "neg")
positive_mode_data_01 <- format_data(positive_mode_data_01, "pos01")
positive_mode_data_02 <- format_data(positive_mode_data_02, "pos02")
positive_mode_data_03 <- format_data(positive_mode_data_03, "pos03")

## Gather mode data in long table format
neg_long <- negative_mode_data %>%
    gather(key = Metabolite, value = Abundance,
           -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))

pos_long_01 <- positive_mode_data_01 %>%
    gather(key = Metabolite, value = Abundance,
           -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))

pos_long_02 <- positive_mode_data_02 %>%
    gather(key = Metabolite, value = Abundance,
           -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))

pos_long_03 <- positive_mode_data_03 %>%
    gather(key = Metabolite, value = Abundance,
           -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))

## combine data into one long table
csf_long <- bind_rows(neg_long,
                     pos_long_01,
                     pos_long_02,
                     pos_long_03)

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
                              min_node_size=5, fig_dir=NULL) {

    print(df$Metabolite[1])

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

    best_fit <- predict(boost_fit_cv, df.omitted, opt_iters)

    ## Note: this isn't ordered by run index 
    trend <- plot(boost_fit_cv, 1, num_trees=opt_iters,
                  grid_levels=df.omitted$RunIndex, return_grid=TRUE)$y
    
    abund <- df$Abundance
    bf <- rep(NA, length(abund))

    abund[not_na_indices] <- df.omitted$Abundance - trend
    bf[not_na_indices] <- best_fit

    df$RawScaled <- df$Abundance
    
    df$Abundance <- abund
    df$Trend <- trend

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
    do(fit_boosted_model(., fig_dir = "../results/gotms_trend_data/"))

## detrended <- subject_data %>%
##     group_by(Metabolite, Mode) %>%
##     filter(Metabolite == "1391-7.063 Results") %>% 
##     do(fit_boosted_model(., fig_dir = "../results/gotms_trend_data/"))

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
     file = sprintf("preprocessed_gotms_data-%s.RData", Sys.Date()))


