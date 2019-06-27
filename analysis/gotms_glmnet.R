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
  
  #drop unused levels
  type <- filtered$Type %>%
    droplevels
  
  Y <- filtered %>% 
    dplyr::select(-one_of("Type", "Type2", "Gender", "Age", "APOE", "Batch",
                        #"Data File",  (not found in dataset, so removed)
                        "Index", "GBAStatus", "Id",
                        "GBA_T369M", "cognitive_status")) %>%
    as.matrix()
  Y[Y==0] <- NA
  
  #do imputation
  Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
  #our functions take matrices, and we add age/gender (1 = F, 2 = M)
  Y_tmp <- t(Yt) %>% 
    as_tibble %>%
    mutate(Age = filtered$Age,
           Gender = filtered$Gender)
  #convert gender to a dummy var (1 if male, 0 if female)
  Y_tmp <- model.matrix(~., Y_tmp)
  #remove intercept column created by model.matrix
  Y <- Y_tmp[,-1]
    
    
  
  
  return(list(Y, type))
}

#does loocv for logistic regression, using misclass error loss (maybe better to use deviance?), eval at different elastic net alphas
fit_glmnet <- function(features, labels, alpha, measure = 'class'){
  set.seed(1)
  #set foldid so that same folds every time (it's loocv so it's just an ordering)
  foldid <- sample(nrow(features))
  
  fit <- cv.glmnet(features, labels, family = 'binomial', 
                   type.measure = measure, nfolds = nrow(features),
                   foldid = foldid, alpha = alpha, grouped = FALSE, standardize = FALSE)
  return(fit)
}


#function to return fpr, tpr given prediction and true label
fpr_tpr <- function(pred, label){
  rocpred <- ROCR::prediction(pred, label)
  rocfpr_tpr <- ROCR::performance(rocpred, measure = 'tpr', x.measure = 'fpr')
  rocauc <- ROCR::performance(rocpred, measure = 'auc')
  return(tibble(fpr = deframe(rocfpr_tpr@x.values), 
                tpr = deframe(rocfpr_tpr@y.values),
                auc = as.numeric(rocauc@y.values)))
}

## get idea of variable importance. note that this relies on our variables being standardized before 
# also glmnet doesn't give us a good idea of standard error because it fits using coord descent (algorithmic, not statistical)
# to fix, maybe refit model in glm?
# numbers are not on any kind of scale, but give some idea of relative importance
importance <- function(fit){
  coefficients <- coef(fit, s = 'lambda.min') %>% 
    abs %>%
    as.matrix
  #order the coefficients for weak measure of importance, and remove coefs with 0 coeff
  coefficients_sorted <- coefficients[order(coefficients, decreasing = TRUE) & coefficients > 0,]
  #map names to the metabolite name
  names(coefficients_sorted) <- str_replace_all(names(coefficients_sorted), ' Result.*', '') %>%
    sapply(function(x) all_matches[match(x, all_matches$Name), 'Metabolite'] %>% deframe)
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

# #write new file as csv
# imputed_pd_co_y %>% 
#   as.data.frame() %>%
#   mutate(Type = imputed_pd_co_labels) %>%
#   write_csv(path = 'gotms_imputed_pd_co.csv')

## START sample analysis with a few extreme alphas ##
# fit
fit_pd_co_lasso <- fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = 1, measure = 'class')
fit_pd_co_ridge <- fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = 0, measure = 'class')
fit_pd_co_half <- fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = 0.5, measure = 'class')

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
fit_pd_co_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = x, measure = 'class'))

#fit models with each of the alphas
pred_pd_co_list <- lapply(fit_pd_co_list, 
                          function(x) predict(x, newx = imputed_pd_co_y, 
                                              type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_co_list <- lapply(fit_pd_co_list, function(x) importance(x))

#write one to csv
#choice of 7th element (ie alpha = .6) is mostly arbitrary.
  #normally, between 30-50 predictors are significant, this one has 40.
importance_pd_co_list[[7]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite)) %>%  #na metabolite is the intercept
  arrange(desc(coefficient)) %>% 
  write_csv(path = 'gotms_glmnet.6_predictors.csv')

#roc for each of the alphas
roc_pd_co_list <- lapply(pred_pd_co_list, function(x) fpr_tpr(x, imputed_pd_co_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha = seq(0,1,.1) %>%    #match with the actual alpha value
           rep(each = length(imputed_pd_co_labels)+1) %>%  #+1 for the 0,0 point
           as.factor)

#plot for all alphas
ggplot(roc_pd_co_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC: PD vs CO')

#look at auc's for each alpha
roc_pd_co_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)

### END PD vs CO analysis ###



### START {AD,PD} vs CO analysis ###

#imputation
imputed_adpd_co <- filter_and_impute(c('AD', 'PD', 'CO'))
imputed_adpd_co_y <- imputed_adpd_co[[1]]
#Group AD, PD into D (for diseased)
imputed_adpd_co_labels <- imputed_adpd_co[[2]] %>% 
  fct_collapse(D = c('AD', 'PD'))

imputed_adpd_co_y %>% 
  as.data.frame() %>%
  mutate(Type = imputed_adpd_co_labels) %>%
  write_csv(path = 'gotms_imputed_adpd_co.csv')


#fit models with each of the alphas
fit_adpd_co_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_adpd_co_y, imputed_adpd_co_labels, alpha = x, measure = 'class'))

pred_adpd_co_list <- lapply(fit_adpd_co_list, 
                          function(x) predict(x, newx = imputed_adpd_co_y, 
                                              type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_adpd_co_list <- lapply(fit_adpd_co_list, function(x) importance(x))

#roc for each of the alphas
#note: both alpha = 0 and alpha = 1 give constant prediction probability for all observations
  #something is probably wrong with the code, need to look at it.
roc_adpd_co_list <- lapply(pred_adpd_co_list, function(x) fpr_tpr(x, imputed_adpd_co_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha = c(0,0, seq(0.1,.9,.1) %>%    #match with the actual alpha value
          rep(each = length(imputed_adpd_co_labels)+1),  #+1 for the 0,0 point
          1,1) %>% as.factor)

#plot for all alphas
  #note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_adpd_co_list, mapping = aes(tpr, fpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {AD,PD} vs CO')

#look at auc for each alpha. need to do (1- auc) to show the flip
roc_adpd_co_list %>% 
  mutate(auc = 1 - auc) %>%
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)


## look at the ridge regression in particular, since it gave constant results.
  # same thing happened in the lasso case. not sure why yet
fit_adpd_co_ridge <- fit_glmnet(imputed_adpd_co_y, imputed_adpd_co_labels, alpha = 0, measure = 'class')
pred_adpd_co_ridge <- predict(fit_adpd_co_ridge, newx = imputed_adpd_co_y, s = 'lambda.min', type = 'response')
rocpred_adpd_co_ridge <- ROCR::prediction(pred_adpd_co_ridge, imputed_adpd_co_labels)
rocperf_adpd_co_ridge <- ROCR::performance(rocpred_adpd_co_ridge, measure = 'tpr', x.measure = 'fpr')

# todo: look into https://www.ggplot2-exts.org/plotROC.html
plot(rocperf_adpd_co_ridge)
abline(a = 0, b = 1, lty= 2)




