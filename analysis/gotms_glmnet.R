source('analysis/starter.R')


############################

### PD vs CO ###

############################


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
fit_pd_co_lasso <- fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = 1)
fit_pd_co_ridge <- fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = 0)
fit_pd_co_half <- fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = 0.5)

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

fit_pd_co_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_co_y, imputed_pd_co_labels, alpha = x))

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
importance_pd_co_list[[8]] %>% 
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
  labs(title = 'ROC: PD vs CO',
       subtitle = 'GOT')

#look at auc's for each alpha
roc_pd_co_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)





############################

### {AD, PD} vs CO ###
## With Age/Gender Penalty ##

############################

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
fit_adpd_co_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_adpd_co_y, imputed_adpd_co_labels, alpha = x))

pred_adpd_co_list <- lapply(fit_adpd_co_list, 
                          function(x) predict(x, newx = imputed_adpd_co_y, 
                                              type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_adpd_co_list <- lapply(fit_adpd_co_list, function(x) importance(x))

importance_adpd_co_list[[7]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 
  

#roc for each of the alphas
#TODO: look into why we have flipped positive and negative classes here.
roc_adpd_co_list <- lapply(pred_adpd_co_list, function(x) fpr_tpr(x, imputed_adpd_co_labels, label_ordering = c('D', 'CO'))) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_adpd_co_labels) + 1) %>%
           as.factor)
  
#plot for all alphas
  #note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_adpd_co_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {AD,PD} vs CO',
       subtitle = 'GOT, with age/gender penalty')

#look at auc for each alpha. need to do (1- auc) to show the flip
roc_adpd_co_list %>% 
  #mutate(auc = 1 - auc) %>%
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)


## look at the ridge regression in particular, since it gave weird results.
fit_adpd_co_ridge <- fit_glmnet(imputed_adpd_co_y, imputed_adpd_co_labels, alpha = 0)
pred_adpd_co_ridge <- predict(fit_adpd_co_ridge, newx = imputed_adpd_co_y, s = 'lambda.min', type = 'response')
rocpred_adpd_co_ridge <- ROCR::prediction(pred_adpd_co_ridge, imputed_adpd_co_labels)
rocperf_adpd_co_ridge <- ROCR::performance(rocpred_adpd_co_ridge, measure = 'tpr', x.measure = 'fpr')

# todo: look into https://www.ggplot2-exts.org/plotROC.html
plot(rocperf_adpd_co_ridge)
abline(a = 0, b = 1, lty= 2)




############################

### {AD, PD} vs CO ###
## no Age/Gender Penalty ##

############################

#fit models with each of the alphas
fit_adpd_co_no_penalty_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_adpd_co_y, imputed_adpd_co_labels, alpha = x, penalize_age_gender = FALSE))

pred_adpd_co_no_penalty_list <- lapply(fit_adpd_co_no_penalty_list, 
                            function(x) predict(x, newx = imputed_adpd_co_y, 
                                                type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_adpd_co_no_penalty_list <- lapply(fit_adpd_co_no_penalty_list, function(x) importance(x))

importance_adpd_co_no_penalty_list[[8]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 


#roc for each of the alphas
#TODO: look into why we have flipped positive and negative classes here.
roc_adpd_co_no_penalty_list <- lapply(pred_adpd_co_no_penalty_list, function(x) fpr_tpr(x, imputed_adpd_co_labels, label_ordering = c('D', 'CO'))) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_adpd_co_labels) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_adpd_co_no_penalty_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {AD,PD} vs CO',
       subtitle = 'GOT, without age/gender penalty')

#look at auc for each alpha. need to do (1- auc) to show the flip
roc_adpd_co_no_penalty_list %>% 
  #mutate(auc = 1 - auc) %>%
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)



############################

### {AD, PD} vs CO ###
## no Age/Gender Penalty, rawScaled ##

############################

#imputation
imputed_adpd_co_rawScaled <- filter_and_impute(wide_data_rawscaled,c('AD', 'PD', 'CO'))
#features
imputed_adpd_co_rawScaled_y <- imputed_adpd_co_rawScaled[[1]]
#Group AD, PD into D (for diseased)
imputed_adpd_co_rawScaled_labels <- imputed_adpd_co_rawScaled[[2]] %>% 
  fct_collapse(D = c('AD', 'PD'))


#fit models with each of the alphas
fit_adpd_co_no_penalty_rawScaled_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_adpd_co_rawScaled_y, imputed_adpd_co_rawScaled_labels, alpha = x, penalize_age_gender = FALSE))

pred_adpd_co_no_penalty_rawScaled_list <- lapply(fit_adpd_co_no_penalty_rawScaled_list, 
                                       function(x) predict(x, newx = imputed_adpd_co_rawScaled_y, 
                                                           type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_adpd_co_no_penalty_rawScaled_list <- lapply(fit_adpd_co_no_penalty_rawScaled_list, function(x) importance(x))

importance_adpd_co_no_penalty_rawScaled_list[[5]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 


#roc for each of the alphas
#TODO: look into why we have flipped positive and negative classes here.
roc_adpd_co_no_penalty_rawScaled_list <- lapply(pred_adpd_co_no_penalty_rawScaled_list, function(x) fpr_tpr(x, imputed_adpd_co_rawScaled_labels, label_ordering = c('D', 'CO'))) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_adpd_co_rawScaled_labels) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_adpd_co_no_penalty_rawScaled_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {AD,PD} vs CO',
       subtitle = 'GOT, without age/gender penalty, rawScaled')

#look at auc for each alpha. need to do (1- auc) to show the flip
roc_adpd_co_no_penalty_rawScaled_list %>% 
  #mutate(auc = 1 - auc) %>%
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)




############################

### {AD, PD} vs CO ###
## no Age/Gender Penalty, Raw ##

############################

#imputation
imputed_adpd_co_raw <- filter_and_impute(wide_data_raw,c('AD', 'PD', 'CO'))
#features
imputed_adpd_co_raw_y <- imputed_adpd_co_raw[[1]]
#Group AD, PD into D (for diseased)
imputed_adpd_co_raw_labels <- imputed_adpd_co_raw[[2]] %>% 
  fct_collapse(D = c('AD', 'PD'))


#fit models with each of the alphas
fit_adpd_co_no_penalty_raw_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_adpd_co_raw_y, imputed_adpd_co_raw_labels, alpha = x, penalize_age_gender = FALSE))

pred_adpd_co_no_penalty_raw_list <- lapply(fit_adpd_co_no_penalty_raw_list, 
                                                 function(x) predict(x, newx = imputed_adpd_co_raw_y, 
                                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_adpd_co_no_penalty_raw_list <- lapply(fit_adpd_co_no_penalty_raw_list, function(x) importance(x))

importance_adpd_co_no_penalty_raw_list[[5]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 


#roc for each of the alphas
#TODO: look into why we have flipped positive and negative classes here.
roc_adpd_co_no_penalty_raw_list <- lapply(pred_adpd_co_no_penalty_raw_list, function(x) fpr_tpr(x, imputed_adpd_co_raw_labels, label_ordering = c('D', 'CO'))) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_adpd_co_raw_labels) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_adpd_co_no_penalty_raw_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {AD,PD} vs CO',
       subtitle = 'GOT, without age/gender penalty, raw')

#look at auc for each alpha. need to do (1- auc) to show the flip
roc_adpd_co_no_penalty_raw_list %>% 
  #mutate(auc = 1 - auc) %>%
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)





############################

### AD vs PD vs CO ###

############################


#start with same imputation as the {ad,pd} vs CO analysis, just without grouping ad and pd
imputed_ad_pd_co_y <- imputed_adpd_co[[1]]
imputed_ad_pd_co_labels <- imputed_adpd_co[[2]]

set.seed(1)
#set foldid so that same folds every time (it's loocv so it's just an ordering)
foldid <- sample(nrow(imputed_ad_pd_co_y))

fit_ad_pd_co_ridge <- cv.glmnet(imputed_ad_pd_co_y, imputed_ad_pd_co_labels, family = 'multinomial', 
                 type.measure = 'deviance', nfolds = nrow(imputed_ad_pd_co_y),
                 foldid = foldid, alpha = 0,  standardize = TRUE)
pred_ad_pd_co_ridge <- predict(fit_ad_pd_co_ridge, newx = imputed_ad_pd_co_y, 'lambda.min', 'class')

accuracy_ad_pd_co_ridge <- (data.frame(predicted = pred_ad_pd_co_ridge, truth = imputed_ad_pd_co_labels) %>%
  table %>%
  diag %>%
  sum)/length(pred)


confusion_ad_pd_co_ridge <- (data.frame(predicted = pred_ad_pd_co_ridge, truth = imputed_ad_pd_co_labels) %>%
                              table %>%
                              diag %>%
                              sum)/length(pred)





#now do the list with multiple alphas
fit_ad_pd_co_list <- lapply(seq(0, 1, .1), function(x) cv.glmnet(imputed_ad_pd_co_y, imputed_ad_pd_co_labels, family = 'multinomial', 
                                                                 type.measure = 'deviance', nfolds = nrow(imputed_ad_pd_co_y),
                                                                 foldid = foldid, alpha = x, standardize = TRUE))
pred_ad_pd_co_list <- lapply(fit_ad_pd_co_list, 
                            function(x) predict(x, newx = imputed_ad_pd_co_y, 
                                                type = 'class', s = 'lambda.min'))

# Need to do extra work to get coefficients for the multinomial case, since each class has its own set of coefficients
multinomial_importance <- function(fit, metabolites = FALSE){
  coefficients <- coef(fit, s = 'lambda.min') %>% 
    lapply(as.matrix)
  
  #order the coefficients for weak measure of importance, and remove coefs with 0 coeff
  #coefficients_sorted <- coefficients[order(abs(coefficients), decreasing = TRUE) & abs(coefficients) > 0,]
  coefficients_sorted_with_zeroes <- lapply(coefficients, function(x) x[order(abs(x), decreasing = TRUE),])
  coefficients_sorted <- lapply(coefficients_sorted_with_zeroes, function(x) x[x != 0])
  if(metabolites == TRUE){
    #map metabolites to their names. keep the names of gender and age, since they aren't metabolites
    coefficients_sorted <- lapply(coefficients_sorted, function(x) setNames(x, if_else(names(x) %in% c('Age', 'GenderM', 'TypeCM', 'TypeCO', 'TypeCY', 'TypePD'),
                                                                names(x), 
                                                                str_replace_all(names(x), ' Result.*', "") %>%
                                                                  str_replace_all('`', '') %>%
                                                                  str_replace_all('\\\\', '') %>%
                                                                  sapply(function(y) all_matches[match(y, all_matches$Name), 'Metabolite'] %>% deframe))))
    
  }
  return(coefficients_sorted)
}
importance_ad_pd_co_list <- lapply(fit_ad_pd_co_list, function(x) multinomial_importance(x, metabolites = TRUE))



accuracy_ad_pd_co_list <- lapply(pred_ad_pd_co_list, function(x) (data.frame(predicted = x, truth = imputed_ad_pd_co_labels) %>%
                                                                                              table %>%
                                                                                              diag %>%
                                                                                              sum)/length(pred)
)

confusion_ad_pd_co_list <- lapply(pred_ad_pd_co_list, function(x) data.frame(predicted = x, truth = imputed_ad_pd_co_labels) %>%
                                                                     table)





############################

### PD vs {CO, CM, CY} ###

############################


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
fit_pd_c_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_y, imputed_pd_c_labels, alpha = x))

pred_pd_c_list <- lapply(fit_pd_c_list, 
                            function(x) predict(x, newx = imputed_pd_c_y, 
                                                type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list <- lapply(fit_pd_c_list, function(x) importance(x))

importance_pd_c_list[[8]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 


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
       subtitle = 'GOT, With Age, Gender Penalty')
ggsave('roc_pd_c_with_penalty.png')


#look at auc for each alpha.
roc_pd_c_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)







############################

### PD vs {CO, CM, CY} ###
 ## No Gender/Age Penalty ##

############################


#fit models with each of the alphas
fit_pd_c_list_no_penalty <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_y, imputed_pd_c_labels, alpha = x, penalize_age_gender = FALSE))

pred_pd_c_list_no_penalty <- lapply(fit_pd_c_list_no_penalty, 
                         function(x) predict(x, newx = imputed_pd_c_y, 
                                             type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list_no_penalty <- lapply(fit_pd_c_list_no_penalty, function(x) importance(x))


importance_pd_c_list_no_penalty[[7]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 


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


############################

### PD vs {CO, CM, CY} ###
## No Gender/Age Penalty ##
## using loocv models to predict 1 observation at a time (fit on n to get lambda, fit on n-1 using that lambda, predict on 1, repeat n times) ## 


############################


## Do same analysis, but fitting n models on n-1 observations, using the min lambda from fitting on full data ## 

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
       subtitle = TeX('GOT,$\\alpha = 0.5$, loo, using min lambda'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', roc_loo_pd_c_half$auc[1])) #auc is the same for every row in the df
ggsave('roc_pd_c_loo.png')



# to try the list, can do below. need to figure out how to get each of the respective lambdas into the result though.
# loo_pred_list <- lapply(seq(0,1,.1), function(y)
#   sapply(id, function(x) (loo_pred_glmnet(lambda = min_lambda_5, index = x, features = imputed_pd_c_y, labels = imputed_pd_c_labels, alpha = y, penalize_age_gender = FALSE))))





############################

### PD vs {CO, CM, CY} ###
## No Gender/Age Penalty ##
## using loocv models to predict 1 observation at a time (fit on n-1 with different lamda, predict on 1, repeat n times) ## 

############################


  #try alpha = 0.5
id <- sample(nrow(imputed_pd_c_y))

loo_fitcv_half <- lapply(id, function(x) loo_cvfit_glmnet(index = x, features = imputed_pd_c_y, labels = imputed_pd_c_labels, alpha = 0.5, penalize_age_gender = FALSE, family = 'binomial'))

loo_predcv_half <- lapply(loo_fitcv_half, function(x) x[[2]]) %>%
  unlist

roc_loocv_pd_c_half <- fpr_tpr(loo_predcv_half, imputed_pd_c_labels[id])
ggplot(roc_loocv_pd_c_half) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('GOT,$\\alpha = 0.5$, loo fits'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', roc_loocv_pd_c_half$auc[1])) #auc is the same for every row in the df
ggsave('roc_pd_c_loo_fits.png')







############################

### PD vs {CO, CM, CY} ###
## No Gender/Age Penalty ##
## the results of the loo are messed up, so let's try a simple 80/20 train/test ##

############################


set.seed(1)
random_index <- sample(nrow(imputed_pd_c_y), size = floor(nrow(imputed_pd_c_y)*.8))
pd_c_y_train <- imputed_pd_c_y[random_index,]
pd_c_label_train <- imputed_pd_c_labels[random_index]
pd_c_y_test <- imputed_pd_c_y[-random_index,]
pd_c_label_test <- imputed_pd_c_labels[-random_index]

fit_pd_c_list_no_penalty_split <- lapply(seq(0, 1, .1), function(x) fit_glmnet(pd_c_y_train, pd_c_label_train, alpha = x, penalize_age_gender = FALSE))

#go straight to test set
pred_pd_c_list_no_penalty_split <- lapply(fit_pd_c_list_no_penalty_split, 
                                    function(x) predict(x, newx = pd_c_y_test, 
                                                        type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list_no_penalty_split <- lapply(fit_pd_c_list_no_penalty_split, function(x) importance(x))

#roc for each of the alphas
roc_pd_c_list_no_penalty_split <- lapply(pred_pd_c_list_no_penalty_split, function(x) fpr_tpr(x, pd_c_label_test)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = as.factor((as.numeric(alpha) - 1)*.1))

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_pd_c_list_no_penalty_split, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, PD vs {CO, CM, CY}',
       subtitle = 'GOT, No Age, Gender Penalty, 80/20 train/test')
ggsave('roc_pd_c_no_penalty_split.png')

#look at auc for each alpha.
roc_pd_c_list_no_penalty_split %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)

#(results are still good.)






############################

### PD vs {CO, CM, CY} ###
## No Gender/Age Penalty ##
## With raw scaled instead of abundances ##

############################

imputed_pd_c_rawscaled <- filter_and_impute(wide_data_rawscaled, c('PD', 'CO', 'CM', 'CY'))
imputed_pd_c_y_rawscaled <- imputed_pd_c_rawscaled[[1]]
#Group AD, PD into D (for diseased)
imputed_pd_c_labels_rawscaled <- imputed_pd_c_rawscaled[[2]] %>% 
  fct_collapse(C = c('CO', 'CM', 'CY'))


#note that we now set standardize to be true, since we aren't using pre-standardized data
# but cv.glmnet takes us back to original scale after its done
fit_pd_c_list_no_penalty_rawscaled <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_y_rawscaled, imputed_pd_c_labels_rawscaled, alpha = x, penalize_age_gender = FALSE))

pred_pd_c_list_no_penalty_rawscaled <- lapply(fit_pd_c_list_no_penalty_rawscaled, 
                                    function(x) predict(x, newx = imputed_pd_c_y_rawscaled, 
                                                        type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list_no_penalty_rawscaled <- lapply(fit_pd_c_list_no_penalty_rawscaled, function(x) importance(x))

importance_pd_c_list_no_penalty_rawscaled[[8]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 


#roc for each of the alphas
roc_pd_c_list_no_penalty_rawscaled <- lapply(pred_pd_c_list_no_penalty_rawscaled, function(x) fpr_tpr(x, imputed_pd_c_labels_rawscaled)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_pd_c_labels_rawscaled) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_pd_c_list_no_penalty_rawscaled, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, PD vs {CO, CM, CY}',
       subtitle = 'GOT, rawScaled')
ggsave('roc_pd_c_no_penalty_rawscaled.png')

#look at auc for each alpha.
roc_pd_c_list_no_penalty_rawscaled %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)



############################

### PD vs {CO, CM, CY} ###
## No Gender/Age Penalty ##
## With raw instead of abundances ##

############################


imputed_pd_c_raw <- filter_and_impute(wide_data_raw, c('PD', 'CO', 'CM', 'CY'))
imputed_pd_c_y_raw <- imputed_pd_c_raw[[1]]
#Group AD, PD into D (for diseased)
imputed_pd_c_labels_raw <- imputed_pd_c_raw[[2]] %>% 
  fct_collapse(C = c('CO', 'CM', 'CY'))


fit_pd_c_list_no_penalty_raw <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_y, imputed_pd_c_labels, alpha = x, penalize_age_gender = FALSE))

pred_pd_c_list_no_penalty_raw <- lapply(fit_pd_c_list_no_penalty_raw, 
                                              function(x) predict(x, newx = imputed_pd_c_y, 
                                                                  type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list_no_penalty_raw <- lapply(fit_pd_c_list_no_penalty_raw, function(x) importance(x))

importance_pd_c_list_no_penalty_rawscaled[[8]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 

 
#roc for each of the alphas
roc_pd_c_list_no_penalty_raw <- lapply(pred_pd_c_list_no_penalty_raw, function(x) fpr_tpr(x, imputed_pd_c_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_pd_c_labels) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_pd_c_list_no_penalty_raw, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, PD vs {CO, CM, CY}',
       subtitle = 'GOT, No Age, Gender Penalty, Raw')
ggsave('roc_pd_c_no_penalty_raw.png')

#look at auc for each alpha.
roc_pd_c_list_no_penalty_raw %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)






