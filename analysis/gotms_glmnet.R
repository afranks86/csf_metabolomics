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
library(latex2exp)
library(here)

source(here("analysis/utility.R"))

set.seed(1)
# data_path <- '~/course/ND_Metabolomics/'
# processed_files <- dir(path = data_path, pattern="^preprocessed_gotms_data*")
# ## Most recent file
# load(max(file.path(data_path, processed_files[grep("-20+", processed_files)])))
# load(file.path(data_path, 'data', 'got-ms',"identification_map.RData"))

data_path <- file.path('E:', 'Projects', 'metabolomics', 'ND_Metabolomics')
processed_files <- dir(path = file.path(data_path, 'analysis'), pattern="^preprocessed_gotms_data*")
## Most recent file
load(max(file.path(data_path, 'analysis', processed_files[grep("-20+", processed_files)])))
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
filter_and_impute <- function(data, types){
  ## Impute missing values
  filtered <- data %>%
    filter(Type %in% types) %>%
    #remove columns that are ALL NA
    select_if(function(x) any(!is.na(x)))
  
  #drop unused levels
  type <- filtered$Type %>%
    droplevels
  
  #keep gba for potential gba analysis
  gba <- filtered$GBAStatus %>%
    as.factor
    
  
  
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
    
    
  
  
  return(list(Y, type, gba))
}

#does loocv for logistic regression, using deviance loss by default (check?), eval at different elastic net alphas
#if penalize_age_gender is false, set penalty coeff to 0
fit_glmnet <- function(features, labels, alpha, measure = 'deviance', penalize_age_gender = TRUE){
  set.seed(1)
  #set foldid so that same folds every time (it's loocv so it's just an ordering)
  foldid <- sample(nrow(features))
  
  #set penalty factors. (1 for each is default)
  p_factors <- rep(1, ncol(features))
  #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
  if(penalize_age_gender == FALSE){
    age_gender_index <- which(colnames(features) %in% c('Age', 'GenderM'))
    p_factors[age_gender_index] <- 0
  }
  
  fit <- cv.glmnet(features, labels, family = 'binomial', 
                   type.measure = measure, nfolds = nrow(features),
                   foldid = foldid, alpha = alpha, grouped = FALSE, standardize = FALSE, penalty.factor = p_factors)
  
  return(fit)
}

# function to fit glmnet on n-1 observations and predict on 1, and do n times.
# lambda can be min lambda from fit_glmnet
# index is the one observation to remove
# only pass in binomial or gaussian
#https://stats.stackexchange.com/questions/304440/building-final-model-in-glmnet-after-cross-validation
loo_pred_glmnet <- function(lambda, index, features, labels, alpha, penalize_age_gender = TRUE, family = 'binomial'){
  #features and label, leaving out one observation for training
  loo_features <- features[-index,]
  loo_labels <- labels[-index]
  
  #features and labels on the held out observation
  
  new_features <- features[index,] %>% matrix(nrow = 1, dimnames = list('', colnames(loo_features) ))
  new_label <- labels[index]
  
  #set penalty factors. (1 for each is default)
  p_factors <- rep(1, ncol(loo_features))
  #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
  if(penalize_age_gender == FALSE){
    age_gender_index <- which(colnames(loo_features) %in% c('Age', 'GenderM'))
    p_factors[age_gender_index] <- 0
  }
  
  #documentation says we should avoid passing in single value of lambda. how else to do?
  #fit on n-1
  fit <- glmnet(loo_features, loo_labels, family = family, lambda = lambda, alpha = alpha, standardize = FALSE, penalty.factor = p_factors)
  
  #predict on 1
  pred <- predict(fit, newx = new_features, type = 'response', s = lambda)
  return(pred)
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
# numbers are not on any kind of scale, but give some idea of relative importance\
#metabolites is a flag for determining whether we need to map names to their metabolite name.
importance <- function(fit, metabolites = TRUE){
  coefficients <- coef(fit, s = 'lambda.min') %>% 
    as.matrix
  #order the coefficients for weak measure of importance, and remove coefs with 0 coeff
  #coefficients_sorted <- coefficients[order(abs(coefficients), decreasing = TRUE) & abs(coefficients) > 0,]
  coefficients_sorted_with_zeroes <- coefficients[order(abs(coefficients), decreasing = TRUE),]
  coefficients_sorted <- coefficients_sorted_with_zeroes[coefficients_sorted_with_zeroes != 0]
  
  
  if(metabolites == TRUE){
    #map metabolites to their names. keep the names of gender and age, since they aren't metabolites
    names(coefficients_sorted) <- if_else(names(coefficients_sorted) %in% c('Age', 'GenderM', 'TypeCM', 'TypeCO', 'TypeCY', 'TypePD'),
                                          names(coefficients_sorted), 
                                          str_replace_all(names(coefficients_sorted), ' Result.*', "") %>%
                                            str_replace_all('`', '') %>%
                                            str_replace_all('\\\\', '') %>%
                                            sapply(function(x) all_matches[match(x, all_matches$Name), 'Metabolite'] %>% deframe))
  }
  
  return(coefficients_sorted)
  
}


#get the types in our dataset (ie AD, PD, CO, ..)
all_types <- wide_data$Type %>% 
  unique %>%
  as.character

#impute all types using amelia for an example
  #Y holds features (keeping with prior notation), labels holds type
imputed_all <- filter_and_impute(wide_data,all_types)
imputed_all_Y <- imputed_all[[1]]
imputed_all_labels <- imputed_all[[2]] 


### START PD vs CO analysis ###

#filter dataset to only include PD, CO
imputed_pd_co <- filter_and_impute(wide_data,c('PD', 'CO'))
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
pred_pd_co_half <- predict(fit_pd_co_half, newx = imputed_pd_co_y, type = 'response', s = 'lambda.min')

# some measure of variable importance
importance_pd_co_lasso <- importance(fit_pd_co_lasso)
importance_pd_co_ridge <- importance(fit_pd_co_ridge)
importance_pd_co_half <- importance(fit_pd_co_half)

# roc plot for half only
rocpred_pd_co_half <- ROCR::prediction(pred_pd_co_half, imputed_pd_co_labels)
rocperf_pd_co_half <- ROCR::performance(rocpred_pd_co_half, measure = 'tpr', x.measure = 'fpr')

# todo: look into https://www.ggplot2-exts.org/plotROC.html
plot(rocperf_pd_co_half)
abline(a = 0, b = 1, lty= 2)

#alt: use ggplot
roc_pd_co_half <- fpr_tpr(pred_pd_co_half, imputed_pd_co_labels)
ggplot(roc_pd_co_half) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: PD vs CO",
       subtitle = TeX('Metabolites,$\\alpha = 0.5$'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', roc_pd_co_half$auc[1])) #auc is the same for every row in the df
ggsave(filename = 'gotms_roc_pdco.png')

## END sample analysis with a few extreme alphas##


## START try 10 alphas between 0 and 1

fit_pd_co_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = x, measure = 'deviance'))

#fit models with each of the alphas
pred_pd_co_list <- lapply(fit_pd_co_list, 
                          function(x) predict(x, newx = imputed_pd_co_y, 
                                              type = 'response', s = 'lambda.min'))
#some measure of variable importance
  #positive means higher abundance => more likely to have pd
  #negative means higher abundance => more likely to be control
importance_pd_co_list <- lapply(fit_pd_co_list, function(x) importance(x))

#write one to csv
#choice of 6th element (ie alpha = .5) is mostly arbitrary.
  #normally, between 30-50 predictors are significant, this one has 40.
importance_pd_co_list[[6]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite)) %T>%  #na metabolite is the intercept 
  write_csv(path = 'gotms_glmnet.5_predictors.csv')

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
imputed_adpd_co <- filter_and_impute(wide_data,c('AD', 'PD', 'CO'))
#features
imputed_adpd_co_y <- imputed_adpd_co[[1]]
#Group AD, PD into D (for diseased)
imputed_adpd_co_labels <- imputed_adpd_co[[2]] %>% 
  fct_collapse(D = c('AD', 'PD'))

# #write the dataset used to csv if imputation takes a long time
# imputed_adpd_co_y %>% 
#   as.data.frame() %>%
#   mutate(Type = imputed_adpd_co_labels) %>%
#   write_csv(path = 'gotms_imputed_adpd_co.csv')


#fit models with each of the alphas
fit_adpd_co_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_adpd_co_y, imputed_adpd_co_labels, alpha = x, measure = 'deviance'))

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
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_adpd_co_labels) + 1) %>%
           as.factor)
  
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


## look at the ridge regression in particular, since it gave weird results.
fit_adpd_co_ridge <- fit_glmnet(imputed_adpd_co_y, imputed_adpd_co_labels, alpha = 0, measure = 'class')
pred_adpd_co_ridge <- predict(fit_adpd_co_ridge, newx = imputed_adpd_co_y, s = 'lambda.min', type = 'response')
rocpred_adpd_co_ridge <- ROCR::prediction(pred_adpd_co_ridge, imputed_adpd_co_labels)
rocperf_adpd_co_ridge <- ROCR::performance(rocpred_adpd_co_ridge, measure = 'tpr', x.measure = 'fpr')

# todo: look into https://www.ggplot2-exts.org/plotROC.html
plot(rocperf_adpd_co_ridge)
abline(a = 0, b = 1, lty= 2)



### End {AD,PD} vs CO analysis ###


### Start AD vs PD vs CO analysis ###

#start with same imputation as the {ad,pd} vs CO analysis, just without grouping ad and pd
imputed_ad_pd_co_y <- imputed_adpd_co[[1]]
imputed_ad_pd_co_labels <- imputed_adpd_co[[2]]

set.seed(1)
#set foldid so that same folds every time (it's loocv so it's just an ordering)
foldid <- sample(nrow(imputed_ad_pd_co_y))

fit_ad_pd_co_ridge <- cv.glmnet(imputed_ad_pd_co_y, imputed_ad_pd_co_labels, family = 'multinomial', 
                 type.measure = 'deviance', nfolds = nrow(imputed_ad_pd_co_y),
                 foldid = foldid, alpha = 0, grouped = FALSE, standardize = FALSE)
pred_ad_pd_co_ridge <- predict(fit_ad_pd_co_ridge, newx = imputed_ad_pd_co_y, 'lambda.min', 'class')

accuracy_ad_pd_co_ridge <- (data.frame(predicted = pred_ad_pd_co_ridge, truth = imputed_ad_pd_co_labels) %>%
  table %>%
  diag %>%
  sum)/length(pred)


confusion_ad_pd_co_ridge <- (data.frame(predicted = pred_ad_pd_co_ridge, truth = imputed_ad_pd_co_labels) %>%
                              table %>%
                              diag %>%
                              sum)/length(pred)


importance_ad_pd_co_ridge <- importance(fit_ad_pd_co_ridge)


#now do the list with multiple alphas
fit_ad_pd_co_list <- lapply(seq(0, 1, .1), function(x) cv.glmnet(imputed_ad_pd_co_y, imputed_ad_pd_co_labels, family = 'multinomial', 
                                                                 type.measure = 'deviance', nfolds = nrow(imputed_ad_pd_co_y),
                                                                 foldid = foldid, alpha = x, grouped = FALSE, standardize = FALSE))
pred_ad_pd_co_list <- lapply(fit_ad_pd_co_list, 
                            function(x) predict(x, newx = imputed_ad_pd_co_y, 
                                                type = 'class', s = 'lambda.min'))
#some measure of variable importance
importance_ad_pd_co_list <- lapply(fit_ad_pd_co_list, function(x) importance(x))


accuracy_ad_pd_co_list <- lapply(pred_ad_pd_co_list, function(x) (data.frame(predicted = x, truth = imputed_ad_pd_co_labels) %>%
                                                                                              table %>%
                                                                                              diag %>%
                                                                                              sum)/length(pred)
)

confusion_ad_pd_co_list <- lapply(pred_ad_pd_co_list, function(x) data.frame(predicted = x, truth = imputed_ad_pd_co_labels) %>%
                                                                     table)


### End AD vs PD vs CO analysis ###




### Start PD vs {CO, CM, CY} analysis ###

imputed_pd_c <- filter_and_impute(wide_data, c('PD', 'CO', 'CM', 'CY'))
imputed_pd_c_y <- imputed_pd_c[[1]]
#Group AD, PD into D (for diseased)
imputed_pd_c_labels <- imputed_pd_c[[2]] %>% 
  fct_collapse(C = c('CO', 'CM', 'CY'))

# #write the dataset used to csv if imputation takes a long time
# imputed_adpd_co_y %>% 
#   as.data.frame() %>%
#   mutate(Type = imputed_adpd_co_labels) %>%
#   write_csv(path = 'gotms_imputed_adpd_co.csv')


## First try with age/gender penalty

#fit models with each of the alphas
fit_pd_c_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_y, imputed_pd_c_labels, alpha = x, measure = 'deviance'))

pred_pd_c_list <- lapply(fit_pd_c_list, 
                            function(x) predict(x, newx = imputed_pd_c_y, 
                                                type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list <- lapply(fit_pd_c_list, function(x) importance(x))

#roc for each of the alphas
roc_pd_c_list <- lapply(pred_pd_c_list, function(x) fpr_tpr(x, imputed_pd_c_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_pd_c_labels) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_pd_c_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, PD vs {CO, CM, CY}',
       subtitle = 'With Age, Gender Penalty')
ggsave('roc_pd_c_with_penalty.png')


#look at auc for each alpha.
roc_pd_c_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)


## try same analysis, but without penalizing age/gender

#fit models with each of the alphas
fit_pd_c_list_no_penalty <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_y, imputed_pd_c_labels, alpha = x, measure = 'deviance', penalize_age_gender = FALSE))

pred_pd_c_list_no_penalty <- lapply(fit_pd_c_list_no_penalty, 
                         function(x) predict(x, newx = imputed_pd_c_y, 
                                             type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list_no_penalty <- lapply(fit_pd_c_list_no_penalty, function(x) importance(x))

#roc for each of the alphas
roc_pd_c_list_no_penalty <- lapply(pred_pd_c_list_no_penalty, function(x) fpr_tpr(x, imputed_pd_c_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_pd_c_labels) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_pd_c_list_no_penalty, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, PD vs {CO, CM, CY}',
       subtitle = 'No Age, Gender Penalty')
ggsave('roc_pd_c_no_penalty.png')

#look at auc for each alpha.
roc_pd_c_list_no_penalty %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)




## Do same analysis, but fitting n models on n-1 observations

#let's test on the alpha = .5
min_lambda_5 <- fit_pd_c_list_no_penalty[[6]]$lambda.min
id <- sample(nrow(imputed_pd_c_y))

loo_pred <- sapply(id, function(x) (loo_pred_glmnet(lambda = min_lambda_5, index = x, features = imputed_pd_c_y, labels = imputed_pd_c_labels, alpha = 0.5, penalize_age_gender = FALSE)))

roc_loo_pd_c_half <- fpr_tpr(loo_pred, imputed_pd_c_labels)
ggplot(roc_loo_pd_c_half) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('Metabolites,$\\alpha = 0.5$, loo'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', roc_loo_pd_c_half$auc[1])) #auc is the same for every row in the df
ggsave('roc_pd_c_loo.png')



# to try the list, can do below. need to figure out how to get each of the respective lambdas into the result though.
# loo_pred_list <- lapply(seq(0,1,.1), function(y)
#   sapply(id, function(x) (loo_pred_glmnet(lambda = min_lambda_5, index = x, features = imputed_pd_c_y, labels = imputed_pd_c_labels, alpha = y, penalize_age_gender = FALSE))))





### End PD vs {CO, CM, CY} analysis ###




