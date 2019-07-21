source('analysis/starter.R')

#' helper function for scatterplot
age_apoe_plot <- function(data, x, y){
  ggplot(data) + 
    geom_point(aes(!!sym(x), !!sym(y), color = apoe)) + 
    scale_color_brewer(type = 'qual', palette = 'Set1') +
    labs(title = paste0('Control: ', x, ' vs ', y),
         x = x,
         y = y)
}




## use previously created imputation for total dataset
#impute all types using amelia for an example
  #Y holds features (keeping with prior notation), labels holds type
imputed_all <- filter_and_impute(wide_data,c('CO', 'CY', 'CM'))
imputed_all_Y <- imputed_all[[1]]
imputed_all_labels <- imputed_all[[2]]
imputed_all_apoe <- imputed_all[[3]]



imputed_all_age <- imputed_all_Y[,'Age']
# readd type as a feature in this analysis
imputed_all_features_age_tmp <- imputed_all_Y %>% 
  as_tibble %>%
  #mutate(Type = imputed_all_labels) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_all_features_age_tmp<- model.matrix(~., imputed_all_features_age_tmp)
#remove intercept column created by model.matrix
imputed_all_features_age <- imputed_all_features_age_tmp[,-1]



#######################################

### GOT age analysis ###
### Single alpha case: Ridge Regression ###

#########################################

foldid <- sample(nrow(imputed_all_features_age))
fit_age_ridge <- cv.glmnet(imputed_all_features_age, imputed_all_age, family = 'gaussian', 
                           type.measure = 'mse', nfolds = nrow(imputed_all_features_age),
                           foldid = foldid, alpha = 0, standardize = TRUE)

pred_age_ridge <- predict(fit_age_ridge, newx = imputed_all_features_age, s= 'lambda.min')

mse_age_ridge <- (pred_age_ridge - imputed_all_age)^2 %>% mean
mae_age_ridge <- (pred_age_ridge - imputed_all_age) %>% 
  abs %>%
  mean

importance_age_ridge <- importance(fit_age_ridge)




#######################################

### GOT age analysis ###
### multiple alphas ###

#########################################


fit_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_all_features_age, imputed_all_age, family = 'gaussian', alpha = x, penalize_age_gender = FALSE))
# fit_age_list <- lapply(seq(0,1,.1), function(x) cv.glmnet(imputed_all_features_age, imputed_all_age, family = 'gaussian', 
#                                                           type.measure = 'mse', nfolds = nrow(imputed_all_features_age),
#                                                           foldid = foldid, alpha = x, standardize = TRUE))
pred_age_list <- lapply(fit_age_list, 
                        function(x) predict(x, newx = imputed_all_features_age, 
                                            type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_age_list <- lapply(fit_age_list, function(x) importance(x))


importance_adpd_co_list[[6]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 



mse_age_list <- lapply(pred_age_list, function(x) mean((x - imputed_all_age)^2))

#which has lowest mse?
  #it's alpha = 0.8, but same for .7, .6
which.min(mse_age_list)

residuals_age_list <- lapply(pred_age_list, function(x) x - imputed_all_age)


### Normality analysis ###

#shapiro-wilkes test tests H0: data is normal
sw_age_list <- lapply(residuals_age_list, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#all plots look pretty normal!
resid_9 <- residuals_age_list[[9]]

plot(resid_9)
abline(h = 0)

hist(resid_9)

qqnorm(resid_9)
qqline(resid_9)


##########




### Plots for a particular alpha (0.3) ###

age_table_3 <- tibble(truth = imputed_all_age, 
                               pred = as.numeric(pred_age_list[[4]]),
                               resid = truth - pred,
                               apoe = imputed_all_apoe,
                               type = imputed_all_labels
)

ggplot(age_table_3) + 
  geom_point(aes(truth, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'GOT, Residuals = Truth - Pred, alpha = 0.7',
       x = 'True Age',
       y = 'Predicted Age') #+ 
#stat_ellipse(data = filter(age_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')





res_7 <- residuals_age_list[[8]]




## try with loo method to generate standard error. try with the 9th again (alpha = 0.8)
loo_pred_8 <- sapply(foldid, function(x) loo_pred_glmnet(lambda = fit_age_list[[9]]$lambda.min, index = x, features =imputed_all_features_age, labels = imputed_all_age, alpha = 0.8, penalize_age_gender = FALSE, family= 'gaussian'))

#see (standard error of esimate) = sqrt(SS residuals / df), where df = n - number of coefficients in the model (including intercept)
  #how to do standard error if df > n? i do rmse instead 
rmse_loo_pred_8 <- (sum((loo_pred_8 - imputed_all_age)^2) / length(imputed_all_age)) %>% sqrt








#effective degrees





















##########################

### GOT and Lipids Age analysis for controls ###
### In sample ###

###########################

#impute data and split into parts
set.seed(1)
imputed_all_combined <- filter_and_impute(wide_data_combined,c('CO', 'CM', 'CY'))
imputed_all_combined_Y <- imputed_all_combined[[1]]
imputed_all_combined_labels <- imputed_all_combined[[2]]
imputed_all_combined_apoe <- imputed_all_combined[[3]]
imputed_all_combined_age <- imputed_all_combined_Y[,'Age']

# we include apoe in this anlaysis
imputed_all_features_combined_age_tmp <- imputed_all_combined_Y %>% 
  as_tibble %>%
  mutate(APOE = imputed_all_combined_apoe) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_all_features_combined_age_tmp<- model.matrix(~., imputed_all_features_combined_age_tmp)
#remove intercept column created by model.matrix
imputed_all_features_combined_age <- imputed_all_features_combined_age_tmp[,-1]


#now try with list of multiple alphas
fit_combined_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_all_features_combined_age, imputed_all_combined_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

# fit_age_list <- lapply(seq(0,1,.1), function(x) cv.glmnet(imputed_all_features_age, imputed_all_age, family = 'gaussian', 
#                                                           type.measure = 'mse', nfolds = nrow(imputed_all_features_age),
#                                                           foldid = foldid, alpha = x, standardize = TRUE))

#in sample prediction
pred_combined_age_list <- lapply(fit_combined_age_list, 
                        function(x) predict(x, newx = imputed_all_features_combined_age, 
                                            type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_combined_age_list <- lapply(fit_combined_age_list, function(x) importance(x))

mse_combined_age_list <- lapply(pred_combined_age_list, function(x) mean((x - imputed_all_combined_age)^2))

#get the residuals
residuals_combined_age_list <- lapply(pred_combined_age_list, function(x) x - imputed_all_combined_age)

#shapiro-wilkes test tests H0: data is normal
sw_combined_age_list <- lapply(residuals_combined_age_list, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#which has lowest mse?
#it's alpha = 0, but next losest is alph - 0.4
which.min(mse_combined_age_list)


#all plots look pretty normal!
resid_combined_4 <- residuals_combined_age_list[[5]]

qplot(resid_combined_4, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
ggsave('got_lipids_age_control_resid_hist.png')


qqnorm(resid_combined_4)
qqline(resid_combined_4)





### look at alpha = 0.4
age_combined_table_4 <- tibble(truth = imputed_all_combined_age, 
                               pred = as.numeric(pred_combined_age_list[[5]]),
                               resid = truth - pred,
                               apoe = imputed_all_combined_apoe,
                               type = imputed_all_combined_labels
                               )


#pred vs resid
ggplot(age_combined_table_4) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Predicted Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)')
ggsave('age_pred_resid_c_insample_4.png', width = 7.26, height = 7.26, units = 'in')


#truth vs pred
ggplot(age_combined_table_4) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted (in-sample) Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4',
       x = 'True Age',
       y = 'Predicted Age (in-sample)')
ggsave('age_true_pred_c_insample_4.png', width = 7.26, height = 7.26, units = 'in')


#looking at distribution of apoe across the control group
age_combined_table_4 %>% 
  mutate(type = fct_relevel(type, c('CY', 'CM', 'CO'))) %>%
  select(apoe, type) %>%
  table


#quick lm to see if there's a linear relationship between apoe and pred
fit <- lm(pred~apoe+ truth, data = age_combined_table_3)
summary(fit)

fit2 <- lm(pred ~ apoe, data = age_combined_table_3)
summary(fit2)







##########################

### GOT and Lipids Age analysis for controls ###
### Leave one out ###

###########################

### Note: It's really expensive to do a list of alphas like we do with in sample prediction, so I'm just pulling the single alpha based on the insample list


#now try with alpha = 0.4
fitpred_combined_loo_age <- lapply(1:nrow(imputed_all_features_combined_age), function(x) loo_cvfit_glmnet(x, imputed_all_features_combined_age, imputed_all_combined_age, 
                                                           alpha = 0.4, family = 'gaussian', penalize_age_gender = FALSE))

fit_combined_loo_age <- lapply(fitpred_combined_loo_age, function(x) x[[1]])
pred_combined_loo_age <- lapply(fitpred_combined_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_combined_loo_age <- lapply(fit_combined_loo_age, function(x) importance(x))
mse_combined_loo_age <- mean((pred_combined_loo_age - imputed_all_combined_age)^2)
resid_combined_loo_age <- pred_combined_loo_age - imputed_all_combined_age


shapiro.test(resid_combined_loo_age)

qplot(resid_combined_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(resid_combined_loo_age)
qqline(resid_combined_loo_age)


### look at alpha = 0.4
loo_age_table <- tibble(truth = imputed_all_combined_age, 
                               pred = pred_combined_loo_age,
                               resid = truth - pred,
                               apoe = imputed_all_combined_apoe,
                               type = imputed_all_combined_labels
)

#pred vs rsiduals
ggplot(loo_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_control_loo_5.png')


ggplot(loo_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('pred_truth_control_loo_5.png')


#####




