library(tidyverse)
library(broom)
library(magrittr)
library(made4)
library(ggridges)
library(modelr)
library(splines)
library(quantreg)
library(patchwork)
library(gbm3)

lipid_data <- read_csv("~/course/ND_Metabolomics/data/lipidomics/lipid_concentrations.csv",
                       na=c("", ".", "na"),
                       col_types=paste0(c("c", rep("d", 1070)), collapse=""))

## Add Index columns, and take log of concentrations
format_data <- function(df) {

    df %<>%  mutate_at(vars(-c("Name")), funs(log(.)))
    df %<>% mutate(Mode = ifelse(grepl("QC", Name),
                                 "QC","Sample" ))
    df$RunIndex <- 1:nrow(df)
    df$SampleId <- as.numeric(sub("CSF-", replacement="", df$Name))

    df
    
}

lipid_data <-  format_data(lipid_data)

csf_long <- lipid_data %>%
    gather(key = Lipid, value = Abundance,
           -one_of("Name", "Mode", "SampleId", "RunIndex"))


## Split QC and sample data
QC_long <- csf_long %>% filter(grepl("QC", Mode))

lipid_long <- csf_long %>% filter(!grepl("QC", Mode))

## 0 abundace peps are NA for now
lipid_long$Abundance[lipid_long$Abundance == -Inf] <- NA
lipid_long$Abundance[lipid_long$Abundance == 0] <- NA

## Load Tracking data
tracking_data <- read_csv("~/course/ND_Metabolomics/data/NDTracking.csv")

## incorporate supplementary PD info provided by Cyrus
pd_supplement <- read_csv("../data/PANUC/panuc.csv")
pd_supplement <- pd_supplement %>%
    dplyr::select(summary_id, ApoE, GBAStatus, GBA_T369M, cognitive_status) %>%
    rename(Id = summary_id)

## fix APOE format
pd_supplement <- pd_supplement %>% rename(APOE = ApoE) %>%
    mutate(APOE = sub(".", "", APOE, fixed=TRUE))

## grab IPD id
tracking_data <- tracking_data %>%
    mutate(Id = gsub("IPD.+;\\s", "", Id)) %>%
    mutate(APOE = as.character(APOE))

## merge PANUC with AD tracking
tracking_data <- left_join(tracking_data, pd_supplement, by="Id") %>%
    mutate(APOE = coalesce(APOE.x, APOE.y)) %>%
    dplyr::select(-c(APOE.x, APOE.y))

## Join Tracking data and metabolomics
subject_data <- inner_join(tracking_data, lipid_long,
                           by = c("Index" = "SampleId"))
subject_data %<>% mutate_at(vars(c("Type", "Gender", "APOE")), as.factor)
subject_data %<>% mutate_at(vars(c("Index")), as.integer)

subject_data %<>% rename(Raw = Abundance)

subject_data %<>% group_by(Mode, RunIndex) %>%
    mutate(Abundance = Raw - median(Raw, na.rm=TRUE)) %>% ungroup

fit_spline_model <- function(df, knots, Boundary.knots, fig_dir = NULL) {

    ## At least 20 observations
    if( sum( !is.na( df$Abundance) ) > 40 ){

        train <- df %>% select(RunIndex, Gender, Type) %>%
            mutate_at(vars(Gender, Type), as.factor) %>%
            model_matrix(~ RunIndex) %>% as.matrix

        not_na_indices <- !is.na(df$Abundance)
        df.omitted <- df %>% filter(not_na_indices)

        rq_fit <- rq(Abundance ~ bs(RunIndex, knots=knots, Boundary.knots),
                     data = df.omitted, tau=0.5)

        abund <- df.omitted$Abundance

        trend <- predict(rq_fit)
        
        abund <- df.omitted$Abundance - trend

        df$RawScaled <- df$Abundance
        
        df$Abundance[not_na_indices] <- abund
        df$Trend <- df$Abundance
        df$Trend[not_na_indices] <- trend

        ## save figures
        if(!is.null(fig_dir)) {
            print(df$Lipid[1])
            trend_plot <- ggplot(df[not_na_indices,]) +
                geom_point(aes(x=RunIndex, y=RawScaled, shape=Type, col=Type)) +
                geom_line(aes(x=RunIndex, y=trend), size=2, col="red")
            residual_plot <- ggplot(df[not_na_indices,]) +
                geom_point(aes(x=RunIndex, y=Abundance, shape=Type, col=Type))
            combo_plot <- trend_plot + residual_plot
            ggsave(paste0(fig_dir, make.names(df$Lipid[1]), ".pdf"),
                   combo_plot, width=14)
        }
    }

    df
}

fit_boosted_model <- function(df, ntrees=20000, cv.folds=10,
                              min_node_size=10, fig_dir=NULL) {
  print(df$Lipid[1])
  if( sum( !is.na( df$Abundance) ) > 50 ){
      
      df$Abundance[!is.finite(df$Abundance)] <- NA
      not_na_indices <- which(!is.na(df$Abundance))
      df.omitted <- df %>% drop_na(Abundance)

      ntrees <- ntrees/2
      opt_iters <- ntrees
      while(opt_iters == ntrees) {
          ntrees <- 2*ntrees
          boost_fit_cv <- gbm(Abundance ~ RunIndex,
                              data = df.omitted, distribution = "laplace",
                              n.trees = ntrees, cv.folds = 10,
                              n.minobsinnode = min_node_size, verbose = FALSE)
          opt_iters <- gbm.perf(boost_fit_cv, method = "cv")
          print(opt_iters)
      }

      best_fit <- predict(boost_fit_cv, df.omitted, opt_iters)

      abund <- df$Abundance
      bf <- rep(NA, length(abund))

      abund[not_na_indices] <- df.omitted$Abundance - best_fit
      bf[not_na_indices] <- best_fit

      df$RawScaled <- df$Abundance
      
      df$Abundance <- abund
      df$Trend <- bf

      ## save figures
      if(!is.null(fig_dir)) {
          print(df$Lipid[1])
          trend_plot <- ggplot(df[not_na_indices,]) +
              geom_point(aes(x=RunIndex, y=RawScaled, shape=Type, col=Type)) +
              geom_line(aes(x=RunIndex, y=Trend), size=2, col="red")
          residual_plot <- ggplot(df[not_na_indices,]) +
              geom_point(aes(x=RunIndex, y=Abundance, shape=Type, col=Type))
          combo_plot <- trend_plot + residual_plot
          ggsave(paste0(fig_dir, make.names(df$Lipid[1]), ".pdf"),
                 combo_plot, width=14)
      }
  } 
  df 

}





total <- subject_data %>%
    group_by(RunIndex) %>%
    summarise(mean_abundance = mean(Abundance, na.rm=TRUE),
              Gender=Gender[[1]], Type = Type[[1]])

total_res <- fit_boosted_model(rename(total, Abundance=mean_abundance))

total_res %>%
    ggplot() +
    geom_line(aes(x=RunIndex, y=RawScaled)) +
    geom_line(aes(x=RunIndex, y=Trend), size=2, col="red")

subject_data %<>%
    group_by(Lipid) %>%
    mutate(Abundance = Abundance - total$mean_abundance) %>% ungroup

detrended <-
    subject_data %>%
    group_by(Lipid, Mode) %>%
    do(fit_boosted_model(., fig_dir = "../results/lipid_trend_data/"))

subject_data <- detrended

## Type2 collapse (CO, CM, CY) to C
subject_data %<>%
    mutate(Type2 = ifelse(Type %in% c("CY", "CM", "CO"),
                          "C",
                          as.character(.$Type))) %>%
    mutate(Type2 = factor(Type2))

## APOE NA's to "Unknown"
subject_data <-
    subject_data %>%
    mutate(APOE=ifelse(is.na(APOE), "Unknown", as.character(APOE))) %>%
    mutate(APOE=factor(APOE)) %>% 
    ungroup

save(subject_data, QC_long,
     file = sprintf("preprocessed_lipid_data-%s.RData", Sys.Date()))



