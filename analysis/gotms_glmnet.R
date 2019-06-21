library(tidyverse)
library(ggExtra)
library(magrittr)
library(microbenchmark)
#library(mvtnorm)
#library(mvnfast)
#library(rstiefel)
#library(mgCov)
library(Amelia)
library(modelr)
#library(robust)
library(ggridges)
library(Matrix)
#library(car)
#library(patchwork)
library(glmnet)
library(ROCR)

source("utility.R")

data_path <- '~/course/ND_Metabolomics/'
processed_files <- dir(path = data_path, pattern="^preprocessed_gotms_data*")
## Most recent file
load(max(file.path(data_path, processed_files[grep("-20+", processed_files)])))
load(file.path(data_path, 'data', 'got-ms',"identification_map.RData"))

wide_data <- subject_data %>%     
  filter(!(Type %in% c("Other"))) %>%
  unite("Metabolite", c("Metabolite", "Mode")) %>% 
  mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
  dplyr::select(-one_of("Raw", "RawScaled", "Trend",
                        "RunIndex", "Name","Data File")) %>%
  spread(key=Metabolite, value=Abundance)
dim(wide_data)

#' filter type and impute using amelia
#' col is a vector of strings, one of the levels of type
filter_and_impute <- function(types){
  ## Impute missing values
  filtered <- wide_data %>%
    filter(Type %in% types)
  
  type <- filtered$Type %>%
    droplevels
  
  Y <- filtered %>% 
    dplyr::select(-one_of("Type", "Type2", "Gender", "Age", "APOE", "Batch",
                        #"Data File",  (not found in dataset, so removed)
                        "Index", "GBAStatus", "Id",
                        "GBA_T369M", "cognitive_status")) %>%
    as.matrix()
  Y[Y==0] <- NA
  dim(Y)
  
  Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
  Y <- t(Yt) %>% as.matrix
  
  
  return(list(Y, type))
}


all_types <- wide_data$Type %>% 
  unique %>%
  as.character

imputed_all <- filter_and_impute(all_types)
imputed_all_Y <- imputed_all[[1]]
imputed_all_labels <- imputed_all[[2]] 




### PD vs CO analysis ###

#filter dataset to only include PD, CO
imputed_pd_co <- filter_and_impute(c('PD', 'CO'))
imputed_pd_co_y <- imputed_pd_co[[1]]
imputed_pd_co_labels <- imputed_pd_co[[2]]

#create foldid so we can test different alphas on same sets
set.seed(1)
foldid <- sample(nrow(imputed_pd_co_y))

#does loocv for logistic regression, using misclass error loss (maybe better to use deviance?), eval at different elastic net alphas
fit_pd_co_lasso <- cv.glmnet(imputed_pd_co_y, imputed_pd_co_labels, family = 'binomial', 
                       type.measure = 'class', nfolds = nrow(imputed_pd_co_y),
                       foldid = foldid, alpha = 1, grouped = FALSE, standardize = FALSE)

fit_pd_co_ridge <- cv.glmnet(imputed_pd_co_y, imputed_pd_co_labels, family = 'binomial', 
                             type.measure = 'class', nfolds = nrow(imputed_pd_co_y),
                             foldid = foldid, alpha = 0, grouped = FALSE, standardize = FALSE)

fit_pd_co_half <- cv.glmnet(imputed_pd_co_y, imputed_pd_co_labels, family = 'binomial', 
                             type.measure = 'class', nfolds = nrow(imputed_pd_co_y),
                             foldid = foldid, alpha = 0.5, grouped = FALSE, standardize = FALSE)


## predict using response (for prob) or class (.5 threshold) using min lambda (complexity determined by cv)
pred_pd_co_lasso <- predict(fit_pd_co_lasso, newx = imputed_pd_co_y, type = 'response', s = 'lambda.min')
pred_pd_co_ridge <- predict(fit_pd_co_ridge, newx = imputed_pd_co_y, type = 'class', s = 'lambda.min')
pred_pd_co_half <- predict(fit_pd_co_half, newx = imputed_pd_co_y, type = 'class', s = 'lambda.min')


## get idea of variable importance. note that this relies on our variables being standardized before 
# also glmnet doesn't give us a good idea of standard error because it fits using coord descent (algorithmic, not statistical)
  # to fix, maybe refit model in glm?
# numbers are not on any kind of scale, but give some idea of relative importance
importance <- function(fit){
  coefficients <- coef(fit, s = 'lambda.min') %>% 
    abs %>%
    as.matrix
  coefficients_sorted <- coefficients[order(coefficients, decreasing = TRUE) & coefficients > 0,]
  return(coefficients_sorted)
  
}

pd_co_lasso_importance <- importance(fit_pd_co_lasso)
pd_co_ridge_importance <- importance(fit_pd_co_ridge)
pd_co_half_importance <- importance(fit_pd_co_half)


## basic roc plot. perfect prediction for lasso/.5
rocpred_pd_co_lasso <- ROCR::prediction(pred_pd_co_lasso, imputed_pd_co_labels)
rocperf_pd_co_lasso <- ROCR::performance(rocpred_pd_co_lasso, measure = 'tpr', x.measure = 'fpr')

## todo: look into https://www.ggplot2-exts.org/plotROC.html
plot(rocperf_pd_co_lasso)
abline(a = 0, b = 1, lty= 2)


