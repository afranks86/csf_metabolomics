source('analysis/starter.R')


#### START LIPIDS ONLY ANALYSIS ####





############################

### PD vs CO ###

############################

set.seed(1)

#filter dataset to only include PD, CO
imputed_pd_co_lipids <- filter_and_impute(wide_data_lipids,c('PD', 'CO'))
imputed_pd_co_y_lipids <- imputed_pd_co_lipids[[1]]
imputed_pd_co_labels_lipids <- imputed_pd_co_lipids[[2]]



## START sample analysis with a few extreme alphas ##
# fit
fit_pd_co_half_lipids <- fit_glmnet(imputed_pd_co_y_lipids, imputed_pd_co_labels_lipids, alpha = 0.5)

# predict using response (for prob) or class (.5 threshold) using min lambda (complexity determined by cv)
pred_pd_co_half_lipids <- predict(fit_pd_co_half_lipids, newx = imputed_pd_co_y_lipids, type = 'response', s = 'lambda.min')

# some measure of variable importance
importance_pd_co_half_lipids <- importance(fit_pd_co_half_lipids, metabolites = FALSE)

# roc plot for half only
rocpred_pd_co_half_lipids <- ROCR::prediction(pred_pd_co_half_lipids, imputed_pd_co_labels_lipids)
rocperf_pd_co_half_lipids <- ROCR::performance(rocpred_pd_co_half_lipids, measure = 'tpr', x.measure = 'fpr')

# todo: look into https://www.ggplot2-exts.org/plotROC.html
plot(rocperf_pd_co_half_lipids)
abline(a = 0, b = 1, lty= 2)

#alt: use ggplot
roc_pd_co_half_lipids <- fpr_tpr(pred_pd_co_half_lipids, imputed_pd_co_labels_lipids)
ggplot(roc_pd_co_half_lipids) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: PD vs CO",
       subtitle = TeX('Lipids, $\\alpha = 0.5$'),
       x = 'False Positive Rate',
       y = 'True Positive Rate')
ggsave(filename = 'gotms_roc_pdco.png')



## Fitting on the list of 0 -> 1 ##


fit_pd_co_list_lipids <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_co_y_lipids, imputed_pd_co_labels_lipids, alpha = x))

#fit models with each of the alphas
pred_pd_co_list_lipids <- lapply(fit_pd_co_list_lipids, 
                                 function(x) predict(x, newx = imputed_pd_co_y_lipids, 
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
#positive means higher abundance => more likely to have pd
#negative means higher abundance => more likely to be control
importance_pd_co_list_lipids <- lapply(fit_pd_co_list_lipids, function(x) importance(x, metabolites = FALSE))

#write one to csv
#choice of 8th element (ie alpha = .7) is mostly arbitrary.
#normally, between 15-30 predictors are significant.
importance_pd_co_list_lipids[[6]] %>% 
  enframe(name = 'lipid', value= 'coefficient') %>% 
  filter(lipid != '(Intercept)') %T>%
  #arrange(desc(coefficient)) %T>% 
  write_csv(path = 'lipids_glmnet.7_predictors.csv')

#roc for each of the alphas
roc_pd_co_list_lipids <- lapply(pred_pd_co_list_lipids, function(x) fpr_tpr(x, imputed_pd_co_labels_lipids)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha = seq(0,1,.1) %>%    #match with the actual alpha value
           rep(each = length(imputed_pd_co_labels_lipids)+1) %>%  #+1 for the 0,0 point
           as.factor)

#plot for all alphas
ggplot(roc_pd_co_list_lipids, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  theme_minimal() + 
  labs(title = 'ROC: PD vs CO',
       subtitle = 'Lipids')




#look at auc's for each alpha
auc_lipids <- roc_pd_co_list_lipids %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc) 
auc_lipids


## Example plot for just alpha = 0.7
ggplot(roc_pd_co_list_lipids %>% filter(alpha == 0.7), mapping = aes(fpr, tpr))+ 
  geom_line() + 
  theme_minimal() + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: PD vs CO",
       subtitle = TeX('Lipids, $\\alpha = 0.7$'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', (auc_lipids %>% filter(alpha == 0.7))$auc %>% round(digits = 3)))
ggsave(filename = 'lipids_roc_pdco.png')




## Fitting loo for better prediction ##


fitpred_pd_co_loo_lipids <- lapply(1:nrow(imputed_pd_co_y_lipids), function(x) loo_cvfit_glmnet(index = x, features = imputed_pd_co_y_lipids, labels = imputed_pd_co_labels_lipids, 
                                                                                alpha = 0.5, penalize_age_gender = FALSE, family = 'binomial'))
fit_pd_co_loo_lipids <- lapply(fitpred_pd_co_loo_lipids, function(x) x[[1]])
pred_pd_co_loo_lipids <- lapply(fitpred_pd_co_loo_lipids, function(x) x[[2]]) %>% unlist

roc_pd_co_loo_lipids <- fpr_tpr(pred_pd_co_loo_lipids, imputed_pd_co_labels_lipids)
ggplot(roc_pd_co_loo_lipids) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('Lipids,$\\alpha = 0.5$, loo'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', roc_pd_co_loo_lipids$auc[1])) #auc is the same for every row in the df
ggsave(filename = 'roc_lipids_pd_co_loo_5.png')





############################

### GBA vs lipids ###

############################


imputed_pd <- filter_and_impute(wide_data_lipids,c('PD'))
imputed_pd_features <- imputed_pd[[1]]
imputed_pd_gba <- imputed_pd[[4]] %>%
  fct_collapse(Carrier = c('E326K Carrier', 'Pathogenic Carrier', 'CT'))

fit_carrier_pd_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_features, imputed_pd_gba, penalize_age_gender = FALSE, alpha = x))

pred_carrier_pd_list <- lapply(fit_carrier_pd_list, 
                            function(x) predict(x, newx = imputed_pd_features, 
                                                type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_carrier_pd_list <- lapply(fit_carrier_pd_list, function(x) importance(x, metabolites = FALSE))
importance_carrier_pd_list[[7]] %>% 
  enframe(name = 'lipid', value= 'coefficient') %>% 
  filter(lipid != '(Intercept)')

#roc for each of the alphas
#note: both alpha = 0 and alpha = 1 give constant prediction probability for all observations
#something is probably wrong with the code, need to look at it.
roc_carrier_pd_list <- lapply(pred_carrier_pd_list, function(x) fpr_tpr(x, imputed_pd_gba)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = as.factor((as.numeric(alpha) - 1)*.1))


#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_carrier_pd_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {T369M_CT, E326K, Pathogenic Carrier} vs Non-Carrier',
       subtitle = 'lipids')
ggsave(filename = 'lipids_roc_gba.png')

#look at auc for each alpha. need to do (1- auc) to show the flip
roc_carrier_pd_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)


#confusion matrix for one of them
data.frame(pred = ifelse(pred_carrier_pd_list[[5]] > .5, 'Non-Carrier', 'Carrier'), 
           truth = imputed_pd_gba) %>% 
  table





############################

### GBA vs {Lipids, GOT} ###

############################

wide_data_combined <- wide_data %>%
  #remove duplicate columns
  select(-c(Age, Type, Gender, Batch, Index, GBAStatus, GBA_T369M, cognitive_status, APOE, Type2)) %>%
  left_join(wide_data_lipids, by = 'Id') 


imputed_pd_comb <- filter_and_impute(wide_data_combined,c('PD'))
imputed_pd_comb_features <- imputed_pd_comb[[1]]
imputed_pd_comb_gba <- imputed_pd_comb[[4]] %>%
  fct_collapse(Carrier = c('E326K Carrier', 'Pathogenic Carrier', 'CT'))

fit_carrier_pd_comb_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_comb_features, imputed_pd_comb_gba, penalize_age_gender = FALSE, alpha = x))

pred_carrier_pd_comb_list <- lapply(fit_carrier_pd_comb_list, 
                               function(x) predict(x, newx = imputed_pd_comb_features, 
                                                   type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_carrier_pd_comb_list <- lapply(fit_carrier_pd_comb_list, function(x) importance(x, metabolites = FALSE))

importance_carrier_pd_comb_list[[5]] %>% 
  enframe(name = 'lipid', value= 'coefficient') %>% 
  filter(lipid != '(Intercept)')


#roc for each of the alphas
#note: both alpha = 0 and alpha = 1 give constant prediction probability for all observations
#something is probably wrong with the code, need to look at it.
roc_carrier_pd_comb_list <- lapply(pred_carrier_pd_comb_list, function(x) fpr_tpr(x, imputed_pd_comb_gba)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = as.factor((as.numeric(alpha) - 1)*.1))


#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_carrier_pd_comb_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {T369M_CT, E326K, Pathogenic Carrier} vs Non-Carrier',
       subtitle = 'lipids and got')
ggsave(filename = 'lipids_and_got_roc_gba.png')

#look at auc for each alpha. need to do (1- auc) to show the flip
roc_carrier_pd_comb_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)


#confusion matrix for one of them
data.frame(pred = ifelse(pred_carrier_pd_comb_list[[5]] > .5, 'Non-Carrier', 'Carrier'), 
           truth = imputed_pd_comb_gba) %>% 
  table




############################

### PD vs C ###
## {Lipids, GOT} ##

############################

imputed_pd_c_comb <- filter_and_impute(wide_data_combined, c('PD', 'CO', 'CM', 'CY'))
imputed_pd_c_comb_y <- imputed_pd_c_comb[[1]]
imputed_pd_c_comb_labels <- imputed_pd_c_comb[[2]] %>% 
  fct_collapse(C = c('CO', 'CM', 'CY'))

# #write the dataset used to csv if imputation takes a long time
# imputed_adpd_co_y %>% 
#   as.data.frame() %>%
#   mutate(Type = imputed_adpd_co_labels) %>%
#   write_csv(path = 'gotms_imputed_adpd_co.csv')


## First try with age/gender penalty

#fit models with each of the alphas
fit_pd_c_comb_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_comb_y, imputed_pd_c_comb_labels, alpha = x, penalize_age_gender = FALSE))

pred_pd_c_comb_list <- lapply(fit_pd_c_comb_list, 
                         function(x) predict(x, newx = imputed_pd_c_comb_y, 
                                             type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_comb_list <- lapply(fit_pd_c_comb_list, function(x) importance(x))

importance_pd_c_comb_list[[7]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 


#roc for each of the alphas
roc_pd_c_comb_list <- lapply(pred_pd_c_comb_list, function(x) fpr_tpr(x, imputed_pd_c_comb_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_pd_c_comb_labels) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_pd_c_comb_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, PD vs {CO, CM, CY}',
       subtitle = 'GOT + Lipids, no Age, Gender Penalty')
ggsave('roc_pd_c_comb.png')


#look at auc for each alpha.
roc_pd_c_comb_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)




############################

### AD vs C ###
## {Lipids, GOT} ##

############################

imputed_ad_c_comb <- filter_and_impute(wide_data_combined, c('AD', 'CO', 'CM', 'CY'))
imputed_ad_c_comb_y <- imputed_ad_c_comb[[1]]
imputed_ad_c_comb_labels <- imputed_ad_c_comb[[2]] %>% 
  fct_collapse(C = c('CO', 'CM', 'CY'))

#fit models with each of the alphas
fit_ad_c_comb_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_ad_c_comb_y, imputed_ad_c_comb_labels, alpha = x, penalize_age_gender = FALSE))
#in sample prediction
pred_ad_c_comb_list <- lapply(fit_ad_c_comb_list, 
                              function(x) predict(x, newx = imputed_ad_c_comb_y, 
                                                  type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_ad_c_comb_list <- lapply(fit_ad_c_comb_list, function(x) importance(x))

# importance_ad_c_comb_list[[7]] %>% 
#   enframe(name = 'metabolite', value= 'coefficient') %>% 
#   filter(!is.na(metabolite))  #na metabolite is the intercept 


#roc for each of the alphas
roc_ad_c_comb_list <- lapply(pred_ad_c_comb_list, function(x) fpr_tpr(x, imputed_ad_c_comb_labels)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = as.factor((as.numeric(alpha) - 1)*.1))

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_ad_c_comb_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, AD vs {CO, CM, CY}',
       subtitle = 'GOT + Lipids, no Age, Gender Penalty')
ggsave('roc_ad_c_comb.png')


#look at auc for each alpha.
roc_ad_c_comb_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)


## out of sample prediction using the loo fits
## alpha = 0.3 had best roc

fitpred_ad_c_loo <- lapply(1:nrow(imputed_ad_c_comb_y), function(x) loo_cvfit_glmnet(index = x, features = imputed_ad_c_comb_y, labels = imputed_ad_c_comb_labels, 
                                                                                                alpha = 0.3, penalize_age_gender = FALSE, family = 'binomial'))
fit_ad_c_loo <- lapply(fitpred_ad_c_loo, function(x) x[[1]])
pred_ad_c_loo <- lapply(fitpred_ad_c_loo, function(x) x[[2]]) %>% unlist

roc_ad_c_loo <- fpr_tpr(pred_ad_c_loo, imputed_ad_c_comb_labels)
ggplot(roc_ad_c_loo) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: AD vs C",
       subtitle = TeX('GOT + Lipids ,$\\alpha = 0.3$, loo'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', round(roc_ad_c_loo$auc[1], 3)))
ggsave(filename = 'roc_comb_ad_c_loo_3.png')






############################

### PD vs C ###
## Lipids, matched ##

############################



imputed_pd_c_lipids_matched <- filter_and_impute(wide_data_lipids_matched_pd_c, c('PD', 'CO', 'CM', 'CY'))
imputed_pd_c_y_lipids_matched <- imputed_pd_c_lipids_matched[[1]]
imputed_pd_c_labels_lipids_matched <- imputed_pd_c_lipids_matched[[2]] %>% 
  fct_collapse(C = c('CO', 'CM', 'CY'))


fit_pd_c_list_lipids_matched <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_c_y_lipids_matched, imputed_pd_c_labels_lipids_matched, alpha = x, penalize_age_gender = FALSE))
#in sample prediction
pred_pd_c_list_lipids_matched <- lapply(fit_pd_c_list_lipids_matched, 
                                 function(x) predict(x, newx = imputed_pd_c_y_lipids_matched, 
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_pd_c_list_lipids_matched <- lapply(fit_pd_c_list_lipids_matched, function(x) importance(x))

# importance_pd_c_list_lipids_matchedscaled[[8]] %>% 
#   enframe(name = 'metabolite', value= 'coefficient') %>% 
#   filter(!is.na(metabolite))  #na metabolite is the intercept 


#roc for each of the alphas
roc_pd_c_list_lipids_matched <- lapply(pred_pd_c_list_lipids_matched, function(x) fpr_tpr(x, imputed_pd_c_labels_lipids_matched)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = as.factor((as.numeric(alpha) - 1)*.1))


#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_pd_c_list_lipids_matched, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, PD vs {CO, CM, CY} (lipids_matched)',
       subtitle = 'lipids, No Age, Gender Penalty, in sample')
ggsave('roc_pd_c_lipids_matched.png')

#look at auc for each alpha.
roc_pd_c_list_lipids_matched %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)



#out of sample prediction using the loo method
#choice of alpha is kinda arb. using the number of nonzero coefficients as a kinda benchmark
fitpred_pd_c_lipids_matched_loo <- lapply(1:nrow(imputed_pd_c_y_lipids_matched), function(x) loo_cvfit_glmnet(index = x, features = imputed_pd_c_y_lipids_matched, labels = imputed_pd_c_labels_lipids_matched, 
                                                                                                alpha = 0.5, penalize_age_gender = FALSE, family = 'binomial'))
fit_pd_c_lipids_matched_loo <- lapply(fitpred_pd_c_lipids_matched_loo, function(x) x[[1]])
pred_pd_c_lipids_matched_loo <- lapply(fitpred_pd_c_lipids_matched_loo, function(x) x[[2]]) %>% unlist

roc_pd_c_lipids_matched_loo <- fpr_tpr(pred_pd_c_lipids_matched_loo, imputed_pd_c_labels_lipids_matched)
ggplot(roc_pd_c_lipids_matched_loo) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: PD vs C (lipids_matched)",
       subtitle = TeX('lipids,$\\alpha = 0.5$, loo'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', round(roc_pd_c_lipids_matched_loo$auc[1],3))) #auc is the same for every row in the df
ggsave(filename = 'roc_lipids_pd_c_loo_matched5.png')







############################

### Gender vs C ###
## {Lipids, GOT} ##

############################

imputed_c_comb <- filter_and_impute(wide_data_combined, c('CO', 'CM', 'CY'))
imputed_c_gender_comb_y <- imputed_c_comb[[1]][, -which(colnames(imputed_c_comb[[1]]) == 'GenderM')]
imputed_c_gender_comb_type <- imputed_c_comb[[2]] %>% 
  fct_collapse(C = c('CO', 'CM', 'CY'))
imputed_c_comb_gender <- imputed_c_comb[[1]][, 'GenderM'] %>%
  as.factor %>%
  fct_recode(M = '1', F = '0')

#fit models with each of the alphas
fit_c_gender_comb_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_c_gender_comb_y, imputed_c_comb_gender, alpha = x, penalize_age_gender = FALSE))
#in sample prediction
pred_c_gender_comb_list <- lapply(fit_c_gender_comb_list, 
                              function(x) predict(x, newx = imputed_c_gender_comb_y, 
                                                  type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_gender_comb_list <- lapply(fit_c_gender_comb_list, function(x) importance(x))

importance_c_gender_comb_list[[6]] %>%
  enframe(name = 'name', value= 'coefficient') 
  

importance_c_gender_comb_list%>% 
  lapply(function(x) names(x) %>% 
           na.omit) %>%
  unlist %>%
  table %>% enframe(name = 'name', value = 'num_appearances') %>%
  arrange(desc(num_appearances))


#roc for each of the alphas
roc_c_gender_comb_list <- lapply(pred_c_gender_comb_list, function(x) fpr_tpr(x, imputed_c_comb_gender)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = as.factor((as.numeric(alpha) - 1)*.1))

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_c_gender_comb_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, Gender vs {CO, CM, CY}',
       subtitle = 'GOT + Lipids, no Age, Gender Penalty')
ggsave('roc_gender_c_comb_insample.png')


#look at auc for each alpha.
roc_c_gender_comb_list %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc)


## out of sample prediction using the loo fits
## alpha = 0.5 had bestish

fitpred_c_gender_loo <- lapply(1:nrow(imputed_c_gender_comb_y), function(x) loo_cvfit_glmnet(index = x, features = imputed_c_gender_comb_y, labels = imputed_c_comb_gender, 
                                                                                     alpha = 0.5, penalize_age_gender = FALSE, family = 'binomial'))
fit_c_gender_loo <- lapply(fitpred_c_gender_loo, function(x) x[[1]])
pred_c_gender_loo <- lapply(fitpred_c_gender_loo, function(x) x[[2]]) %>% unlist
importance_c_gender_loo_list <- lapply(fit_c_gender_loo, function(x) importance(x))
importance_c_gender_loo_list%>% 
  lapply(function(x) names(x) %>% 
           na.omit) %>%
  unlist %>%
  table %>%
  sort() %>% names




roc_c_gender_loo <- fpr_tpr(pred_c_gender_loo, imputed_c_comb_gender)
ggplot(roc_c_gender_loo) + 
  geom_line(mapping = aes(fpr, tpr)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  theme_minimal() + 
  labs(title = "ROC: Gender vs C",
       subtitle = TeX('GOT + Lipids ,$\\alpha = 0.5$, loo'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = 0.9, y = 0,label = paste0('AUC:', round(roc_c_gender_loo$auc[1], 3)))
ggsave(filename = 'roc_comb_gender_c_loo_5.png')



