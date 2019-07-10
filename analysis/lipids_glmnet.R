source('analysis/gotms_glmnet.R')


#### START LIPIDS ONLY ANALYSIS ####

# processed_files_lipids <- dir(path = data_path, pattern="^preprocessed_lipid_data*")
# load(max(file.path(data_path, processed_files_lipids[grep("-20+", processed_files_lipids)])))


processed_files_lipids <- dir(path = file.path(data_path, 'analysis'), pattern="^preprocessed_lipid_data*")
## Most recent file
load(max(file.path(data_path, 'analysis', processed_files_lipids[grep("-20+", processed_files_lipids)])))



wide_data_lipids <- subject_data %>%     
  filter(!(Type %in% c("Other"))) %>%
  mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
  dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Name", "Mode", "RunIndex")) %>%
  spread(key=Lipid, value=Abundance)


### Start lipids PD vs CO ###



set.seed(1)

#filter dataset to only include PD, CO
imputed_pd_co_lipids <- filter_and_impute(wide_data_lipids,c('PD', 'CO'))
imputed_pd_co_y_lipids <- imputed_pd_co_lipids[[1]]
imputed_pd_co_labels_lipids <- imputed_pd_co_lipids[[2]]

foldid <- sample(nrow(imputed_pd_co_y_lipids))


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
       subtitle = TeX('$\\alpha = 0.5$'),
       x = 'False Positive Rate',
       y = 'True Positive Rate')
ggsave(filename = 'gotms_roc_pdco.png')


#try the list
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
importance_pd_co_list_lipids[[8]] %>% 
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
  labs(title = 'ROC: PD vs CO')




#look at auc's for each alpha
auc_lipids <- roc_pd_co_list_lipids %>% 
  group_by(alpha) %>%
  slice(1) %>%
  select(alpha, auc) 
auc_lipids


#example plot for just alpha = 0.7
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


### start predict GBA vs metabolites for pd ###

imputed_pd <- filter_and_impute(wide_data_lipids,c('PD'))
imputed_pd_features <- imputed_pd[[1]]
imputed_pd_gba <- imputed_pd[[3]] %>%
  fct_collapse(Carrier = c('E326K Carrier', 'Pathogenic Carrier', 'CT'))


fit_carrier_pd_list <- lapply(seq(0, 1, .1), function(x) fit_glmnet(imputed_pd_features, imputed_pd_gba, penalize_age_gender = FALSE, alpha = x))

pred_carrier_pd_list <- lapply(fit_carrier_pd_list, 
                            function(x) predict(x, newx = imputed_pd_features, 
                                                type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_carrier_pd_list <- lapply(fit_carrier_pd_list, function(x) importance(x, metabolites = FALSE))

#roc for each of the alphas
#note: both alpha = 0 and alpha = 1 give constant prediction probability for all observations
#something is probably wrong with the code, need to look at it.
roc_carrier_pd_list <- lapply(pred_carrier_pd_list, function(x) fpr_tpr(x, imputed_pd_gba)) %>%
  bind_rows(.id = 'alpha') %>%      #convert to long format with new id column alpha
  mutate(alpha  = seq(0,1,.1) %>%
           rep(each = length(imputed_pd_gba) + 1) %>%
           as.factor)

#plot for all alphas
#note: tpr and fpr are switched because the roc curve was going the wrong way
ggplot(roc_carrier_pd_list, mapping = aes(fpr, tpr, color = alpha))+ 
  geom_line() + 
  labs(title = 'ROC, {E326K, Pathogenic Carrier} vs Non-Carrier')
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

### end predict GBA vs metabolites for pd ###



#### END LIPIDS ONLY ANALYSIS ####
