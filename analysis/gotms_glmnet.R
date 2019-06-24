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

#create foldid so we can test different alphas on same sets
set.seed(1)
foldid <- sample(nrow(imputed_pd_co_y))


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

#does loocv for logistic regression, using misclass error loss (maybe better to use deviance?), eval at different elastic net alphas
fit_pd_co <- function(alpha, measure = 'class'){
  set.seed(1)
  #set foldid so that same folds every time (it's loocv so it's just an ordering)
  foldid <- sample(nrow(imputed_pd_co_y))
  
  fit <- cv.glmnet(imputed_pd_co_y, imputed_pd_co_labels, family = 'binomial', 
                   type.measure = measure, nfolds = nrow(imputed_pd_co_y),
                   foldid = foldid, alpha = alpha, grouped = FALSE, standardize = FALSE)
}


#function to return fpr, trp given prediction and true label
fpr_tpr <- function(pred, label){
  rocpred <- ROCR::prediction(pred, label)
  rocperf <- ROCR::performance(rocpred, measure = 'tpr', x.measure = 'fpr')
  return(tibble(fpr = deframe(rocperf@x.values), 
                tpr = deframe(rocperf@y.values)))
}

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


#get the types in our dataset (ie AD, PD, CO, ..)
all_types <- wide_data$Type %>% 
  unique %>%
  as.character

#impute all types using amelia for an example
  #Y holds features (keeping with prior notation), labels holds type
imputed_all <- filter_and_impute(all_types)
imputed_all_Y <- imputed_all[[1]]
imputed_all_labels <- imputed_all[[2]] 


### START PD vs CO analysis ###

#filter dataset to only include PD, CO
imputed_pd_co <- filter_and_impute(c('PD', 'CO'))
imputed_pd_co_y <- imputed_pd_co[[1]]
imputed_pd_co_labels <- imputed_pd_co[[2]]


## START sample analysis with a few extreme alphas ##
# fit
fit_pd_co_lasso <- fit_pd_co(alpha = 1, measure = 'class')
fit_pd_co_ridge <- fit_pd_co(alpha = 0, measure = 'class')
fit_pd_co_half <- fit_pd_co(alpha = 0.5, measure = 'class')

# predict using response (for prob) or class (.5 threshold) using min lambda (complexity determined by cv)
pred_pd_co_lasso <- predict(fit_pd_co_lasso, newx = imputed_pd_co_y, type = 'response', s = 'lambda.min')
pred_pd_co_ridge <- predict(fit_pd_co_ridge, newx = imputed_pd_co_y, type = 'class', s = 'lambda.min')
pred_pd_co_half <- predict(fit_pd_co_half, newx = imputed_pd_co_y, type = 'class', s = 'lambda.min')

# some measure of variable importance
importance_pd_co_lasso <- importance(fit_pd_co_lasso)
importance_pd_co_ridge <- importance(fit_pd_co_ridge)
importance_pd_co_half <- importance(fit_pd_co_half)

# roc plot for lasso only
rocpred_pd_co_lasso <- ROCR::prediction(pred_pd_co_lasso, imputed_pd_co_labels)
rocperf_pd_co_lasso <- ROCR::performance(rocpred_pd_co_lasso, measure = 'tpr', x.measure = 'fpr')

# todo: look into https://www.ggplot2-exts.org/plotROC.html
plot(rocperf_pd_co_lasso)
abline(a = 0, b = 1, lty= 2)

## END sample analysis with a few extreme alphas##

## START try 10 alphas between 0 and 1
fit_pd_co_list <- lapply(seq(0, 1, .1), function(x) fit_pd_co(alpha = x, measure = 'class'))

#fit models with each of the alphas
pred_pd_co_list <- lapply(fit_pd_co_list, 
                          function(x) predict(x, newx = imputed_pd_co_y, 
                                              type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_co_list <- lapply(fit_pd_co_list, function(x) importance(x))

#roc for each of the alphas
roc_pd_co_list <- lapply(pred_pd_co_list, function(x) fpr_tpr(x, imputed_pd_co_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha = seq(0,1,.1) %>%    #match with the actual alpha value
           rep(each = length(imputed_pd_co_labels)+1) %>%  #+1 for the 0,0 point
           as.factor)

#plot for all alphas
ggplot(roc_pd_co_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line()

### END PD vs CO analysis ###



### START {AD,PD} vs CO analysis ###
imputed_adpd_co <- filter_and_impute(c('AD', 'PD', 'CO'))
imputed_adpd_co_y <- imputed_adpd_co[[1]]
#Group AD, PD into D (for diseased)
imputed_adpd_co_labels <- imputed_adpd_co[[2]] %>% 
  fct_collapse(D = c('AD', 'PD'))
  



