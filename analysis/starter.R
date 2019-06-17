library(tidyverse)
library(ggridges)
library(gbm3)
library(patchwork)
library(scales)
library(quantreg)
library(glmnet)
library(mice)

source("utility.R")

processed_files <- dir(pattern="^preprocessed_lipid_data*")


processed_files <- dir(pattern="^preprocessed_gotms_data*")


processed_files <- dir(pattern="^preprocessed_untargeted_data*")

## load("preprocessed_gotms_data.RData")

## Most recent file
load(max(processed_files[grep("-20+", processed_files)]))


wide_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Batch", "Name",
                          "Id", "Index", "GBAStatus", "GBA_T369M",
                          "cognitive_status")) %>%
    spread(key=Lipid, value=Abundance)


Yage <- wide_data %>% dplyr::select(Age) %>% as.matrix

Ytype <- wide_data %>% dplyr::select(Type) %>% as.matrix

X <- wide_data %>% dplyr::select(-one_of(c("Age", "Gender", "Type", "APOE",
                                           "Mode", "RunIndex", "Type2"))) %>%
    as.matrix

apply(X, 2, function(x) mean(is.na(x))) %>% summary
X[, 100]


X <- X[, colMeans(is.na(X)) < 0.2]
colnames(X) <- make.names(colnames(X))
mice_X <- mice(X)

Xcomp <- complete(mice_X) %>% as.matrix

res <- cv.glmnet(x=Xcomp, y=Yage, family="gaussian", alpha=0.5)
coef(res, res$lambda.min)
plot(res)
