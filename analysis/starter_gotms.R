library(tidyverse)
library(ggExtra)
library(magrittr)
library(microbenchmark)
library(mvtnorm)
library(mvnfast)
library(rstiefel)
library(mgCov)
library(Amelia)
library(modelr)
library(robust)
library(ggridges)
library(Matrix)
library(car)
library(patchwork)
library(glmnet)

source("utility.R")

processed_files <- dir(pattern="^preprocessed_gotms_data*")
## Most recent file
load(max(processed_files[grep("-20+", processed_files)]))
load("../data/got-ms/identification_map.RData")

wide_data <- subject_data %>%     
  filter(!(Type %in% c("Other"))) %>%
  unite("Metabolite", c("Metabolite", "Mode")) %>% 
  mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
  dplyr::select(-one_of("Raw", "RawScaled", "Trend",
                        "RunIndex", "Name","Data File")) %>%
  spread(key=Metabolite, value=Abundance)
dim(wide_data)

## Impute missing values
Y <- wide_data %>%
  dplyr::select(-one_of("Type", "Type2", "Gender", "Age", "APOE", "Batch",
                        "Data File", "Index", "GBAStatus", "Id",
                        "GBA_T369M", "cognitive_status")) %>%
  as.matrix()
Y[Y==0] <- NA
dim(Y)

Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
Y <- t(Yt) %>% as.matrix
