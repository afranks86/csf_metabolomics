library(tidyverse)
library(broom)
library(magrittr)
library(made4)
library(ggjoy)
library(gbm3)

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

neg_long <- negative_mode_data %>% gather(key = Metabolite, value = Abundance, -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))
pos_long_01 <- positive_mode_data_01 %>% gather(key = Metabolite, value = Abundance, -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))
pos_long_02 <- positive_mode_data_02 %>% gather(key = Metabolite, value = Abundance, -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))
pos_long_03 <- positive_mode_data_03 %>% gather(key = Metabolite, value = Abundance, -one_of("Name", "RunIndex", "Mode", "Data File", "SampleId"))

csf_long <- bind_rows(neg_long,
                     pos_long_01,
                     pos_long_02,
                     pos_long_03)

################################
## Detrend data
################################

## Split QC and sample data
QC_long <- csf_long %>% filter(grepl("QC", Mode))

nd_data <- csf_long %>% filter(!grepl("QC", Mode))

## 0 abundace peps are NA for now
nd_data$Abundance[nd_data$Abundance == -Inf] <- NA
nd_data$Abundance[nd_data$Abundance == 0] <- NA

## Join with Tracking data
tracking_data <- read_csv("~/course/ND_Metabolomics/data/NDTracking.csv")
subject_data <- inner_join(tracking_data, nd_data, by = c("Index" = "SampleId"))
subject_data %<>% mutate_at(vars(c("Type", "Gender", "APOE")), as.factor)
subject_data %<>% mutate_at(vars(c("Index")), as.integer)

subject_data %<>% rename(Raw = Abundance)

subject_data %<>% group_by(Mode, RunIndex) %>% mutate(Abundance = Raw - median(Raw, na.rm=TRUE)) %>% ungroup


fit_boosted_model <- function(df, ntrees=20000, cv.folds=10, min_node_size=10, shrinkage=0.001) {
    
  print(df$Metabolite[1])

  not_na_indices <- !is.na(df$Abundance)
  df.omitted <- df %>% filter(not_na_indices)

  boost_fit_cv <- gbm(
    Abundance ~ RunIndex + Gender + Type, data = df.omitted, distribution = "laplace",
    shrinkage=shrinkage, n.trees = ntrees, cv.folds = 10, n.minobsinnode = min_node_size, verbose = FALSE
  )
  opt_iters <- gbm.perf(boost_fit_cv, method = "cv")
  print(sprintf("Opt iters: %g, err: %f", opt_iters, boost_fit_cv$cv_error[opt_iters]))

  abund <- df.omitted$Abundance

  trend <- plot(boost_fit_cv, 1, num_trees=opt_iters, grid_levels=df.omitted$RunIndex, return_grid=TRUE)$y
  
    abund <- df.omitted$Abundance - trend

    df$RawScaled <- df$Abundance
    
    df$Abundance[not_na_indices] <- abund
    df$Trend <- df$Abundance
    df$Trend[not_na_indices] <- trend

    df
}

detrended <- subject_data %>%  group_by(Metabolite, Mode) %>% do(fit_boosted_model(., shrinkage=0.001))

subject_data <- detrended

## Type2 collapse (CO, CM, CY) to C
subject_data %<>% mutate(Type2 = ifelse(Type %in% c("CY", "CM", "CO"), "C", as.character(.$Type))) %>% mutate(Type2 = factor(Type2))

## APOE NA's to "Unknown"
subject_data <- subject_data %>% mutate(APOE=ifelse(is.na(APOE), "Unknown", as.character(APOE))) %>% mutate(APOE=factor(APOE))
subject_data %<>% ungroup

save(subject_data, QC_long, file = sprintf("preprocessed_gotms_data-%s.RData", Sys.Date()))


## Experiment with deterending
if(FALSE) {
    met <- sample(subject_data$Metabolite, 1)
    tmp <- subject_data %>%  filter(Metabolite == met) %>% group_by(Metabolite, Mode) %>% do(fit_boosted_model(., shrinkage=0.01))
    p1 <- tmp %>% filter(Metabolite == met) %>% ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) + geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=RawScaled))

    tmp <- subject_data %>%  filter(Metabolite == met) %>% group_by(Metabolite, Mode) %>% do(fit_boosted_model(., shrinkage=0.001))
    p2 <- tmp %>% filter(Metabolite == met) %>% ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) + geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=RawScaled))
    p1 + p2
}
