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

imputed_comb_all <- filter_and_impute(wide_data_combined, c('AD', 'PD', 'CY', 'CM', 'CO'))



## use previously created imputation for total dataset
#impute controls using amelia for an example
  #Y holds features (keeping with prior notation), labels holds type
  #this is for got
imputed_c_got <- filter_and_impute(wide_data,c('CO', 'CY', 'CM'))
imputed_c_got_Y <- imputed_c_got[[1]]
imputed_c_got_labels <- imputed_c_got[[2]]
imputed_c_got_apoe <- imputed_c_got[[3]]
imputed_c_got_gender <- imputed_c_got_Y[,'GenderM']


imputed_c_got_age <- imputed_c_got_Y[,'Age']
# readd type as a feature in this analysis
imputed_c_got_features_age_tmp <- imputed_c_got_Y %>% 
  as_tibble %>%
  #mutate(Type = imputed_c_got_labels) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_c_got_features_age <- model.matrix(~., imputed_c_got_features_age_tmp)



#######################################

### GOT age analysis ###
### Single alpha case: Ridge Regression ###

#########################################

foldid <- sample(nrow(imputed_c_got_features_age))
fit_age_ridge <- cv.glmnet(imputed_c_got_features_age, imputed_c_got_age, family = 'gaussian', 
                           type.measure = 'mse', nfolds = nrow(imputed_c_got_features_age),
                           foldid = foldid, alpha = 0, standardize = TRUE)

pred_age_ridge <- predict(fit_age_ridge, newx = imputed_c_got_features_age, s= 'lambda.min')

mse_age_ridge <- (pred_age_ridge - imputed_c_got_age)^2 %>% mean
mae_age_ridge <- (pred_age_ridge - imputed_c_got_age) %>% 
  abs %>%
  mean

importance_age_ridge <- importance(fit_age_ridge)




#######################################

### GOT age analysis for control ###
### multiple alphas ###

#########################################


fit_c_got_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_got_features_age, imputed_c_got_age, family = 'gaussian', alpha = x, penalize_age_gender = FALSE))
# fit_c_got_age_list <- lapply(seq(0,1,.1), function(x) cv.glmnet(imputed_c_got_features_age, imputed_c_got_age, family = 'gaussian', 
#                                                           type.measure = 'mse', nfolds = nrow(imputed_c_got_features_age),
#                                                           foldid = foldid, alpha = x, standardize = TRUE))
pred_c_got_age_list <- lapply(fit_c_got_age_list, 
                        function(x) predict(x, newx = imputed_c_got_features_age, 
                                            type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_got_age_list <- lapply(fit_c_got_age_list, function(x) importance(x, metabolites = FALSE))


importance_c_got_age_list[[6]] %>% 
  enframe(name = 'metabolite', value= 'coefficient') %>% 
  filter(!is.na(metabolite))  #na metabolite is the intercept 



mse_c_age_list <- lapply(pred_c_got_age_list, function(x) mean((x - imputed_c_got_age)^2))

#which has lowest mse?
  #it's alpha = 0.8, but same for .7, .6
which.min(mse_c_age_list)

residuals_c_age_list <- lapply(pred_c_got_age_list, function(x) x - imputed_c_got_age)


### Normality analysis ###

#shapiro-wilkes test tests H0: data is normal
sw_c_age_list <- lapply(residuals_c_age_list, function(x) (shapiro.test(x))$p.value)

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

age_c_table_3 <- tibble(truth = imputed_c_got_age, 
                               pred = as.numeric(pred_c_got_age_list[[4]]),
                               resid = truth - pred,
                               apoe = imputed_c_got_apoe,
                               type = imputed_c_got_labels
)

ggplot(age_c_table_3) + 
  geom_point(aes(truth, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'GOT, Residuals = Truth - Pred, alpha = 0.7',
       x = 'True Age',
       y = 'Predicted Age') #+ 
#stat_ellipse(data = filter(age_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')





### leave one out using chosen alpha from above

fitpred_c_got_loo_age <- lapply(1:nrow(imputed_c_got_features_age), function(x) loo_cvfit_glmnet(x, imputed_c_got_features_age, imputed_c_got_age, 
                                                                                                     alpha = 0.4, family = 'gaussian', penalize_age_gender = FALSE))

fit_c_got_loo_age <- lapply(fitpred_c_got_loo_age, function(x) x[[1]])
pred_c_got_loo_age <- lapply(fitpred_c_got_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_c_got_loo_age <- lapply(fit_c_got_loo_age, function(x) importance(x))
mse_c_got_loo_age <- mean((pred_c_got_loo_age - imputed_c_got_age)^2)
resid_c_got_loo_age <- pred_c_got_loo_age - imputed_c_got_age


shapiro.test(resid_c_got_loo_age)

qplot(resid_c_got_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_got_age_control_resid_hist.png')
qqnorm(resid_c_got_loo_age)
qqline(resid_c_got_loo_age)


### look at alpha = 0.4
got_c_loo_age_table <- tibble(truth = imputed_c_got_age, 
                               pred = pred_c_got_loo_age,
                               resid = truth - pred,
                               apoe = imputed_c_got_apoe,
                               type = imputed_c_got_labels,
                               gender = imputed_c_got_gender,
                               apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)


#### colored by APOE  ####
#pred vs residuals
ggplot(got_c_loo_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'GOT, alpha = 0.4, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_got_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_got_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_control_got_loo_5.png')


ggplot(got_c_loo_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'GOT, alpha = 0.4, loo',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(intercept = 0, slope = 1)
ggsave('pred_truth_control_got_loo_5.png')













# res_7 <- residuals_age_list[[8]]
# 
# 
# 
# 
# ## try with loo method to generate standard error. try with the 9th again (alpha = 0.8)
# #loo_pred_8 <- sapply(foldid, function(x) loo_pred_glmnet(lambda = fit_c_got_age_list[[9]]$lambda.min, index = x, features =imputed_c_got_features_age, labels = imputed_c_got_age, alpha = 0.8, penalize_age_gender = FALSE, family= 'gaussian'))
# 
# #see (standard error of esimate) = sqrt(SS residuals / df), where df = n - number of coefficients in the model (including intercept)
#   #how to do standard error if df > n? i do rmse instead 
# rmse_loo_pred_8 <- (sum((loo_pred_8 - imputed_c_got_age)^2) / length(imputed_c_got_age)) %>% sqrt








#effective degrees







#######################################

### LIPIDS age analysis for control ###

#########################################


imputed_c_lipids <- filter_and_impute(wide_data_lipids,c('CO', 'CY', 'CM'))
imputed_c_lipids_Y <- imputed_c_lipids[[1]]
imputed_c_lipids_labels <- imputed_c_lipids[[2]]
imputed_c_lipids_apoe <- imputed_c_lipids[[3]]
imputed_c_lipids_gender <- imputed_c_lipids_Y[,'GenderM'] %>% as.factor %>%
  fct_recode(M = '1', F = '0')





imputed_c_lipids_age <- imputed_c_lipids_Y[,'Age']
# readd type as a feature in this analysis
imputed_c_features_lipids_age_tmp <- imputed_c_lipids_Y %>% 
  as_tibble %>%
  #mutate(Type = imputed_c_labels) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_c_features_lipids_age <- model.matrix(~., imputed_c_features_lipids_age_tmp)


### list to determine alphas
fit_c_lipids_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_features_lipids_age, imputed_c_lipids_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))


#in sample prediction
pred_c_lipids_age_list <- lapply(fit_c_lipids_age_list, 
                                 function(x) predict(x, newx = imputed_c_features_lipids_age, 
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_lipids_age_list <- lapply(fit_c_lipids_age_list, function(x) importance(x, metabolites = FALSE))

mse_c_lipids_age_list <- lapply(pred_c_lipids_age_list, function(x) mean((x - imputed_c_lipids_age)^2))




### leave one out using chosen alpha from above

fitpred_c_lipids_loo_age <- lapply(1:nrow(imputed_c_features_lipids_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_lipids_age, imputed_c_lipids_age, 
                                                                                                         alpha = 0.4, family = 'gaussian', penalize_age_gender = FALSE))

fit_c_lipids_loo_age <- lapply(fitpred_c_lipids_loo_age, function(x) x[[1]])
pred_c_lipids_loo_age <- lapply(fitpred_c_lipids_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_c_lipids_loo_age <- lapply(fit_c_lipids_loo_age, function(x) importance(x))
mse_c_lipids_loo_age <- mean((pred_c_lipids_loo_age - imputed_c_lipids_age)^2)
resid_c_lipids_loo_age <- pred_c_lipids_loo_age - imputed_c_lipids_age


shapiro.test(resid_c_lipids_loo_age)

qplot(resid_c_lipids_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(resid_c_lipids_loo_age)
qqline(resid_c_lipids_loo_age)


### look at alpha = 0.4
lipids_c_loo_age_table <- tibble(truth = imputed_c_lipids_age, 
                        pred = pred_c_lipids_loo_age,
                        resid = truth - pred,
                        apoe = imputed_c_lipids_apoe,
                        type = imputed_c_lipids_labels,
                        gender = imputed_c_lipids_gender,
                        apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)

#### colored by APOE  ####
#pred vs residuals
ggplot(lipids_c_loo_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Lipid, alpha = 0.4, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_lipids_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_lipids_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_control_lipids_loo_5.png')


ggplot(lipids_c_loo_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Lipids, alpha = 0.4, loo',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(intercept = 0, slope = 1)
ggsave('pred_truth_control_lipids_loo_5.png')





##### spearman correlation between got and lipids #####
cor(got_c_loo_age_table$resid, lipids_c_loo_age_table$resid, method = 'spearman')
cor.test(got_c_loo_age_table$resid, lipids_c_loo_age_table$resid, method = 'spearman')

got_lipids_resid <- tibble(got = got_c_loo_age_table$resid,
                           lipids = lipids_c_loo_age_table$resid)

ggplot(got_lipids_resid) + 
  geom_point(aes(got, lipids))



##### regularized canonical correlation between got and lipids #####
# First, get data without the common variables age, gender, intercept (first column and last 2)
got_only_c_Y <- imputed_c_got_Y[,-c(1, ncol(imputed_c_got_Y)-1, ncol(imputed_c_got_Y))]
lipids_only_c_Y <- imputed_c_lipids_Y[,-c(1, ncol(imputed_c_lipids_Y)-1, ncol(imputed_c_lipids_Y))]

#grid for the loocv to find ridge param. defaults to 5x5 gridsearch
# params <- CCA::estim.regul(got_only_c_Y, lipids_only_c_Y)
# rcca <- CCA::rcc(got_only_c_Y, lipids_only_c_Y, params$lambda1, params$lamda2)
# barplot(rcca$cor, xlab = "Dimension",
#         + ylab = "Canonical correlations", ylim = c(0,1))


##########################

### GOT and Lipids Age analysis for controls ###
### In sample ###
### Without APOE as predictor ###

###########################

# #impute data and split into parts
set.seed(1)
imputed_c_combined <- filter_and_impute(wide_data_combined,c('CO', 'CM', 'CY'))
imputed_c_combined_Y <- imputed_c_combined[[1]]
imputed_c_combined_labels <- imputed_c_combined[[2]]
imputed_c_combined_apoe <- imputed_c_combined[[3]]
imputed_c_combined_age <- imputed_c_combined_Y[,'Age']

# c_index <- which(imputed_comb_all[[2]] %in% c('CM', 'CY', 'CO'))
# imputed_c_combined_Y <- (imputed_comb_all[[1]])[c_index,]
# imputed_c_combined_labels <- (imputed_comb_all[[2]])[c_index]
# imputed_c_combined_apoe <- (imputed_comb_all[[3]])[c_index]
# imputed_c_combined_age <- imputed_c_combined_Y[, 'Age']
imputed_c_combined_gender <- imputed_c_combined_Y[,'GenderM'] %>% as.factor %>%
  fct_recode(M = '1', F = '0')


## without APOE as a predictor ##
imputed_c_features_combined_age_tmp <- imputed_c_combined_Y %>% 
  as_tibble(.name_repair = 'minimal') %>%
  select(-c(Age))

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_c_features_combined_age <- model.matrix(~., imputed_c_features_combined_age_tmp)

###

#now try with list of multiple alphas
fit_combined_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_features_combined_age, imputed_c_combined_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

# fit_age_list <- lapply(seq(0,1,.1), function(x) cv.glmnet(imputed_c_features_age, imputed_c_age, family = 'gaussian', 
#                                                           type.measure = 'mse', nfolds = nrow(imputed_c_features_age),
#                                                           foldid = foldid, alpha = x, standardize = TRUE))

#in sample prediction
pred_combined_age_list <- lapply(fit_combined_age_list, 
                        function(x) predict(x, newx = imputed_c_features_combined_age, 
                                            type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_combined_age_list <- lapply(fit_combined_age_list, function(x) importance(x))

mse_combined_age_list <- lapply(pred_combined_age_list, function(x) mean((x - imputed_c_combined_age)^2))

#get the residuals
residuals_combined_age_list <- lapply(pred_combined_age_list, function(x) x - imputed_c_combined_age)

#shapiro-wilkes test tests H0: data is normal
sw_combined_age_list <- lapply(residuals_combined_age_list, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#which has lowest mse?
#it's alpha = 0, but next losest is alph - 0.4
which.min(mse_combined_age_list)


#c plots look pretty normal!
resid_combined_4 <- residuals_combined_age_list[[5]]

qplot(resid_combined_4, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
ggsave('got_lipids_age_control_resid_hist.png')


qqnorm(resid_combined_4)
qqline(resid_combined_4)





### look at alpha = 0.4
age_combined_table_4 <- tibble(truth = imputed_c_combined_age, 
                               pred = as.numeric(pred_combined_age_list[[5]]),
                               resid = truth - pred,
                               apoe = imputed_c_combined_apoe,
                               type = imputed_c_combined_labels
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

age_combined_table_3 <- age_combined_table_3 %>%
  mutate(apoe = ifelse(apoe %in% c(24, 34, 44), 4, 0))

#quick lm to see if there's a linear relationship between apoe and pred
fit <- lm(pred~apoe+ truth, data = age_combined_table_3)
summary(fit)

fit2 <- lm(pred ~ apoe, data = age_combined_table_3)
summary(fit2)

fit3 <- lm(resid ~ apoe, data = age_combined_table_3)







##########################

### GOT and Lipids Age analysis for controls ###
### Leave one out  (without apoe as predictor) ###

###########################

### Note: It's really expensive to do a list of alphas like we do with in sample prediction, so I'm just pulling the single alpha based on the insample list


#now try with alpha = 0.4
fitpred_combined_loo_age <- lapply(1:nrow(imputed_c_features_combined_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_combined_age, imputed_c_combined_age, 
                                                           alpha = 0.4, family = 'gaussian', penalize_age_gender = FALSE))

fit_combined_loo_age <- lapply(fitpred_combined_loo_age, function(x) x[[1]])
pred_combined_loo_age <- lapply(fitpred_combined_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_combined_loo_age <- lapply(fit_combined_loo_age, function(x) importance(x))
mse_combined_loo_age <- mean((pred_combined_loo_age - imputed_c_combined_age)^2)
resid_combined_loo_age <- pred_combined_loo_age - imputed_c_combined_age


shapiro.test(resid_combined_loo_age)

qplot(resid_combined_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(resid_combined_loo_age)
qqline(resid_combined_loo_age)


### look at alpha = 0.4
loo_age_table <- tibble(truth = imputed_c_combined_age, 
                        pred = pred_combined_loo_age,
                        resid = truth - pred,
                        apoe = imputed_c_combined_apoe,
                        type = imputed_c_combined_labels,
                        gender = imputed_c_combined_gender,
                        apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)

#### colored by APOE  ####
#pred vs residuals
ggplot(loo_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_control_loo_4.png')


ggplot(loo_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4, loo',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('pred_truth_control_loo_4.png')


#####




#### colored by Type ####
#pred vs residuals
ggplot(loo_age_table) + 
  geom_point(aes(pred, resid, color = type)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_control_loo_5.png')




ggplot(loo_age_table) + 
  geom_point(aes(truth, pred, color = type)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('pred_truth_control_loo_5.png')

ggplot(loo_age_table) + 
  geom_histogram(aes(resid, fill = apoe4), position = 'identity', alpha = 0.5, bins = 10) + 
  labs(title = 'Residuals Histogram , split by APOE4',
       subtitle = 'GOT + Lipids (limited to Controls), alpha = 0.5, loo',
       x = 'Residuals (Truth - Pred)')
ggsave('age_resid_c_apoe4_5.png')





### colored by gender
ggplot(loo_age_table) + 
  geom_point(aes(pred, resid, color = gender)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') +
  geom_hline(yintercept = 0)
ggsave('pred_age_residuals_control_loo_gender5.png')

ggplot(loo_age_table) + 
  geom_point(aes(truth, pred, color = gender)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(intercept = 0, slope = 1)
ggsave('pred_truth_control_loo_gender5.png')



### Predict on ad/pd #####
### Combined ###

# imputed_adpd_combined <- filter_and_impute(wide_data_combined,c('AD', 'PD'))
# #features
# imputed_adpd_combined_Y <- imputed_adpd_combined[[1]]
# imputed_adpd_combined_labels <- imputed_adpd_combined[[2]]
# imputed_adpd_combined_apoe <- imputed_adpd_combined[[3]]
# imputed_adpd_combined_age <- imputed_adpd_combined_Y[,'Age']

adpd_index <- which(imputed_comb_all[[2]] %in% c('AD', 'PD'))
imputed_adpd_combined_Y <- (imputed_comb_all[[1]])[adpd_index,]
imputed_adpd_combined_labels <- (imputed_comb_all[[2]])[adpd_index] %>%
  droplevels
imputed_adpd_combined_apoe <- (imputed_comb_all[[3]])[adpd_index]
imputed_adpd_combined_age <- imputed_adpd_combined_Y[, 'Age']
imputed_adpd_combined_gender <- imputed_adpd_combined_Y[,'GenderM'] %>% as.factor %>%
  fct_recode(M = '1', F = '0')

# without apoe
imputed_adpd_features_combined_age_tmp <- imputed_adpd_combined_Y %>% 
  as_tibble %>%
  #mutate(APOE = imputed_c_combined_apoe) %>%
  select(-Age)

#turn factors into dummy vars.
imputed_adpd_features_combined_age <- model.matrix(~., imputed_adpd_features_combined_age_tmp)



pred_age_adpd_4 <- predict(fit_combined_age_list[[5]], newx = imputed_adpd_features_combined_age, s = 'lambda.min')
adpd_age_table_4 <- tibble(truth = imputed_adpd_combined_age, 
                           pred = pred_age_adpd_4,
                           resid = truth - pred,
                           apoe = imputed_adpd_combined_apoe,
                           type = imputed_adpd_combined_labels,
                           apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
                                  )


#pred vs residuals
ggplot(adpd_age_table_4) + 
  geom_point(aes(pred, resid, color = type)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'AD/PD: Age vs Residuals (trained on Control)',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_adpd_4.png')


ggplot(adpd_age_table_4) + 
  geom_point(aes(truth, pred, color = type)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'AD/PD: True vs Predicted Age (trained on control)',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(slope = 1, intercept = 0)
ggsave('pred_age_truth_adpd_4.png')

ggplot(adpd_age_table_4) + 
  geom_point(aes(truth, pred, color = apoe4)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'AD/PD: True vs Predicted Age (trained on control)',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('pred_truth_adpd_apoe4_4.png')

ggplot(adpd_age_table_4) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age (trained on control)',
       subtitle = 'Combined GOT and Lipid, alpha = 0.4',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('pred_truth_adpd_apoe_4.png')


#side by side residuals hist
residuals_comparison <- tibble(resid = adpd_age_table_4$resid,
                               id = "adpd") %>%
  bind_rows(tibble(resid = age_combined_table_4$resid,
                   id = 'control')) %>%
  mutate(id = as.factor(id))

ggplot(residuals_comparison) +
  geom_histogram(aes(resid, fill = id), position = 'identity', alpha = 0.5, bins = 10) + 
  labs(title = 'Residuals Histogram, split by type',
       subtitle = 'GOT + Lipids, alpha = 0.4',
       x = 'Residuals (Truth - Pred)')
ggsave('age_resid_type_4.png')



#################



##### GOT and Lipids AGE analysis for controls ###
##### With apoe as predictor #####

## with apoe as a predictor ##
imputed_c_features_combined_age_apoe_tmp <- imputed_c_combined_Y %>% 
  as_tibble %>%
  mutate(APOE = imputed_c_combined_apoe) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_c_features_combined_age_apoe <- model.matrix(~., imputed_c_features_combined_age_apoe_tmp)


fit_c_combined_age_apoe_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_features_combined_age_apoe, imputed_c_combined_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

#in sample prediction
pred_c_combined_age_apoe_list <- lapply(fit_c_combined_age_apoe_list, 
                                 function(x) predict(x, newx = imputed_c_features_combined_age_apoe, 
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_combined_age_apoe_list <- lapply(fit_c_combined_age_apoe_list, function(x) importance(x, metabolites = TRUE))

importance_c_combined_age_apoe_list[[5]] %>% 
  enframe(name = 'name', value= 'coefficient') %>%
  filter(name != '(Intercept)')


mse_c_combined_age_apoe_list <- lapply(pred_c_combined_age_apoe_list, function(x) mean((x - imputed_c_combined_age)^2))

#get the residuals
residuals_c_combined_age_apoe_list <- lapply(pred_c_combined_age_apoe_list, function(x) x - imputed_c_combined_age)

#shapiro-wilkes test tests H0: data is normal
sw_c_combined_age_apoe_list <- lapply(residuals_c_combined_age_apoe_list, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#which has lowest mse?
#it's alpha = 0, but next losest is alph - 0.4
which.min(mse_c_combined_age_apoe_list)


#c plots look pretty normal!
resid_c_combined_apoe4 <- residuals_c_combined_age_apoe_list[[5]]

qplot(resid_c_combined_apoe4, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')


qqnorm(resid_c_combined_apoe4)
qqline(resid_c_combined_apoe4)





### look at alpha = 0.4
age_apoe_c_combined_table_4 <- tibble(truth = imputed_c_combined_age, 
                               pred = as.numeric(pred_c_combined_age_apoe_list[[5]]),
                               resid = truth - pred,
                               apoe = imputed_c_combined_apoe,
                               type = imputed_c_combined_labels
)


#pred vs resid
ggplot(age_apoe_c_combined_table_4) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Predicted Age vs Residuals',
       subtitle = 'Combined GOT and Lipid (including apoe), alpha = 0.4',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)')
ggsave('age_pred_resid_c_apoe_insample_4.png', width = 7.26, height = 7.26, units = 'in')


#truth vs pred
ggplot(age_apoe_c_combined_table_4) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted (in-sample) Age',
       subtitle = 'Combined GOT and Lipid (including apoe), alpha = 0.4',
       x = 'True Age',
       y = 'Predicted Age (in-sample)')
ggsave('age_true_pred_c_apoe_insample_4.png', width = 7.26, height = 7.26, units = 'in')




###loo

























##########################

### GOT and Lipids Age analysis for controls ###
### Random Forest ###

###########################
library(randomForest)
set.seed(1)
#need univeral name repair to fit the random forest
imputed_c_combined_df <- imputed_c_combined_Y %>% as_tibble(.name_repair = 'universal')
imputed_adpd_combined_df <- imputed_adpd_combined_Y %>% as_tibble(.name_repair = 'universal')

#best performance with somewhere between 22 and 44 predictors (22 is better). lets say 30
cvrf_age_c <- randomForest::rfcv(imputed_c_features_combined_age, imputed_c_combined_age, cv.fold = nrow(imputed_c_features_combined_age))
#look for best mtry
best_mtry <- randomForest:::tuneRF(imputed_c_features_combined_age, imputed_c_combined_age)


fitrf_age_c <- randomForest::randomForest(Age ~ ., data = imputed_c_combined_df)
#oob prediction on training set (control)2
predrf_age_c <- predict(fitrf_age_c)
#out of sample prediction on adpd
predrf_age_adpd <- predict(fitrf_age_c, newdata = select(imputed_adpd_combined_df, -Age))

rf_age_c_table <- tibble(truth = imputed_c_combined_age, 
                        pred = predrf_age_c,
                        resid = truth - pred,
                        apoe = imputed_c_combined_apoe,
                        type = imputed_c_combined_labels,
                        gender = imputed_c_combined_gender,
                        apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)


ggplot(rf_age_c_table) + 
  geom_point(aes(truth, pred, color = gender)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age (trained on control)',
       subtitle = 'Combined GOT and Lipid, Random Forest, oob',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(slope = 1, intercept = 0)


rf_age_adpd_table <- tibble(truth = imputed_adpd_combined_age, 
                         pred = predrf_age_adpd,
                         resid = truth - pred,
                         apoe = imputed_adpd_combined_apoe,
                         type = imputed_adpd_combined_labels,
                         gender = imputed_adpd_combined_gender,
                         apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)

ggplot(rf_age_adpd_table) + 
  geom_point(aes(truth, pred, color = gender)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'AD/PD: True vs Predicted Age (trained on control)',
       subtitle = 'Combined GOT and Lipid, Random Forest',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(slope = 1, intercept = 0)




##########################

### GOT and Lipids Age analysis for all ###

###########################

#impute data and split into parts
set.seed(1)
imputed_all_combined <- filter_and_impute(wide_data_combined,c('CO', 'CM', 'CY', 'PD', 'AD'))
imputed_all_combined_Y <- imputed_all_combined[[1]]
imputed_all_combined_labels <- imputed_all_combined[[2]] %>%
  fct_collapse(C = c('CO', 'CM', 'CY'))
imputed_all_combined_apoe <- imputed_all_combined[[3]]
imputed_all_combined_age <- imputed_all_combined_Y[,'Age']

# we include apoe in this anlaysis
imputed_all_features_combined_age_tmp <- imputed_all_combined_Y %>%
  as_tibble %>%
  mutate(APOE = imputed_all_combined_apoe) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_all_features_combined_age <- model.matrix(~., imputed_all_features_combined_age_tmp)


#now try with list of multiple alphas
fit_combined_all_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_all_features_combined_age, imputed_all_combined_age,
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

# fit_age_list <- lapply(seq(0,1,.1), function(x) cv.glmnet(imputed_c_got_features_age, imputed_all_age, family = 'gaussian',
#                                                           type.measure = 'mse', nfolds = nrow(imputed_c_got_features_age),
#                                                           foldid = foldid, alpha = x, standardize = TRUE))

#in sample prediction
pred_combined_all_age_list <- lapply(fit_combined_all_age_list,
                                 function(x) predict(x, newx = imputed_c_got_features_combined_age,
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_combined_all_age_list <- lapply(fit_combined_all_age_list, function(x) importance(x))

mse_combined_all_age_list <- lapply(pred_combined_all_age_list, function(x) mean((x - imputed_all_combined_age)^2))

#get the residuals
residuals_combined_all_age_list <- lapply(pred_combined_all_age_list, function(x) x - imputed_all_combined_age)

#shapiro-wilkes test tests H0: data is normal
sw_combined_all_age_list <- lapply(residuals_combined_all_age_list, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#which has lowest mse?
#it's alpha = 0, but next losest is alph - 0.4
which.min(mse_combined_all_age_list)



### look at alpha = 0.7, which has the lowest mse
age_combined_all_table_7 <- tibble(truth = imputed_all_combined_age,
                               pred = as.numeric(pred_combined_all_age_list[[8]]),
                               resid = truth - pred,
                               apoe = imputed_all_combined_apoe,
                               type = imputed_all_combined_labels
)


#pred vs resid
ggplot(age_combined_all_table_7) +
  geom_point(aes(pred, resid, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'All: Predicted Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.7',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)')
ggsave('age_pred_resid_all_insample_7.png', width = 7.26, height = 7.26, units = 'in')


#truth vs pred
ggplot(age_combined_all_table_7) +
  geom_point(aes(truth, pred, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'All: True vs Predicted (in-sample) Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.7',
       x = 'True Age',
       y = 'Predicted Age (in-sample)')
ggsave('age_true_pred_all_insample_7.png', width = 7.26, height = 7.26, units = 'in')




##########################

### GOT and Lipids Age analysis for all ###
### Leave one out ###

###########################

### Note: It's recy expensive to do a list of alphas like we do with in sample prediction, so I'm just pulling the single alpha based on the insample list


#now try with alpha = 0.4
fitpred_combined_all_loo_age <- lapply(1:nrow(imputed_all_features_combined_age), function(x) loo_cvfit_glmnet(x, imputed_all_features_combined_age, imputed_all_combined_age,
                                                                                                         alpha = 0.7, family = 'gaussian', penalize_age_gender = FALSE))

fit_combined_all_loo_age <- lapply(fitpred_combined_all_loo_age, function(x) x[[1]])
pred_combined_all_loo_age <- lapply(fitpred_combined_all_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_combined_all_loo_age <- lapply(fit_combined_all_loo_age, function(x) importance(x))
mse_combined_all_loo_age <- mean((pred_combined_all_loo_age - imputed_c_combined_age)^2)
resid_combined_all_loo_age <- pred_combined_all_loo_age - imputed_c_combined_age


shapiro.test(resid_combined_all_loo_age)

qplot(resid_combined_all_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(resid_combined_all_loo_age)
qqline(resid_combined_all_loo_age)


### look at alpha = 0.4
loo_age_all_table <- tibble(truth = imputed_all_combined_age,
                        pred = pred_combined_all_loo_age,
                        resid = truth - pred,
                        apoe = imputed_all_combined_apoe,
                        type = imputed_all_combined_labels
)

#### colored by APOE  ####
#pred vs residuals
ggplot(loo_age_all_table) +
  geom_point(aes(pred, resid, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'All: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") +
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_all_loo_5.png')


ggplot(loo_age_all_table) +
  geom_point(aes(truth, pred, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'All: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('pred_truth_all_loo_5.png')


#####


##########################

### GOT and Lipids Age analysis for controls ###
### Split by Gender ###

###########################

male_index <- which(imputed_c_combined_gender == 'M')
imputed_c_male_combined_Y <- imputed_c_combined_Y[male_index,]
imputed_c_male_combined_labels <- imputed_c_combined_labels[male_index]
imputed_c_male_combined_apoe <- imputed_c_combined_apoe[male_index]
imputed_c_male_combined_age <- imputed_c_combined_age[male_index]

imputed_c_female_combined_Y <- imputed_c_combined_Y[-male_index,]
imputed_c_female_combined_labels <- imputed_c_combined_labels[-male_index]
imputed_c_female_combined_apoe <- imputed_c_combined_apoe[-male_index]
imputed_c_female_combined_age <- imputed_c_combined_age[-male_index]


### Male first

## without APOE as a predictor ##
imputed_c_male_features_combined_age_tmp <- imputed_c_male_combined_Y %>% 
  as_tibble(.name_repair = 'minimal') %>%
  select(-c(Age))

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_c_male_features_combined_age <- model.matrix(~., imputed_c_male_features_combined_age_tmp)

#now try with list of multiple alphas
fit_c_male_combined_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_male_features_combined_age, imputed_c_male_combined_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

#in sample prediction
pred_c_male_combined_age_list <- lapply(fit_c_male_combined_age_list, 
                                 function(x) predict(x, newx = imputed_c_male_features_combined_age, 
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_male_combined_age_list <- lapply(fit_c_male_combined_age_list, function(x) importance(x, metabolites = FALSE))

mse_c_male_combined_age_list <- lapply(pred_c_male_combined_age_list, function(x) mean((x - imputed_c_male_combined_age)^2))


## leave one out (alpha = 0.7)
fitpred_c_male_combined_loo_age <- lapply(1:nrow(imputed_c_male_features_combined_age), function(x) loo_cvfit_glmnet(x, imputed_c_male_features_combined_age, imputed_c_male_combined_age,
                                                                                                                 alpha = 0.7, family = 'gaussian', penalize_age_gender = FALSE))

fit_c_male_combined_loo_age <- lapply(fitpred_c_male_combined_loo_age, function(x) x[[1]])
pred_c_male_combined_loo_age <- lapply(fitpred_c_male_combined_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_c_male_combined_loo_age <- lapply(fit_c_male_combined_loo_age, function(x) importance(x))
mse_c_male_combined_loo_age <- mean((pred_c_male_combined_loo_age - imputed_c_male_combined_age)^2)
resid_c_male_combined_loo_age <- pred_c_male_combined_loo_age - imputed_c_male_combined_age


shapiro.test(resid_c_male_combined_loo_age)

qplot(resid_c_male_combined_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(resid_c_male_combined_loo_age)
qqline(resid_c_male_combined_loo_age)


### look at alpha = 0.4
loo_c_male_age_table <- tibble(truth = imputed_c_male_combined_age,
                            pred = pred_c_male_combined_loo_age,
                            resid = truth - pred,
                            apoe = imputed_c_male_combined_apoe,
                            apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33')),
                            type = imputed_c_male_combined_labels
)

#### colored by APOE  ####
#pred vs residuals
ggplot(loo_c_male_age_table) +
  geom_point(aes(pred, resid, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Male Controls: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.7, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") +
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_c_male_loo_7.png')


ggplot(loo_c_male_age_table) +
  geom_point(aes(truth, pred, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Male Controls: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.7, loo',
       x = 'True Age',
       y = 'Predicted Age') +
  geom_abline(intercept = 0, slope=1)
ggsave('pred_truth_c_male_loo_7.png')





### Females

## without APOE as a predictor ##
imputed_c_female_features_combined_age_tmp <- imputed_c_female_combined_Y %>% 
  as_tibble(.name_repair = 'minimal') %>%
  select(-c(Age))

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_c_female_features_combined_age <- model.matrix(~., imputed_c_female_features_combined_age_tmp)

#now try with list of multiple alphas
fit_c_female_combined_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_female_features_combined_age, imputed_c_female_combined_age, 
                                                                           family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

#in sample prediction
pred_c_female_combined_age_list <- lapply(fit_c_female_combined_age_list, 
                                        function(x) predict(x, newx = imputed_c_female_features_combined_age, 
                                                            type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_female_combined_age_list <- lapply(fit_c_female_combined_age_list, function(x) importance(x, metabolites = FALSE))

mse_c_female_combined_age_list <- lapply(pred_c_female_combined_age_list, function(x) mean((x - imputed_c_female_combined_age)^2))


## leave one out (alpha = 0.3)
fitpred_c_female_combined_loo_age <- lapply(1:nrow(imputed_c_female_features_combined_age), function(x) loo_cvfit_glmnet(x, imputed_c_female_features_combined_age, imputed_c_female_combined_age,
                                                                                                                     alpha = 0.3, family = 'gaussian', penalize_age_gender = FALSE))

fit_c_female_combined_loo_age <- lapply(fitpred_c_female_combined_loo_age, function(x) x[[1]])
pred_c_female_combined_loo_age <- lapply(fitpred_c_female_combined_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_c_female_combined_loo_age <- lapply(fit_c_female_combined_loo_age, function(x) importance(x))
mse_c_female_combined_loo_age <- mean((pred_c_female_combined_loo_age - imputed_c_female_combined_age)^2)
resid_c_female_combined_loo_age <- pred_c_female_combined_loo_age - imputed_c_female_combined_age


shapiro.test(resid_c_female_combined_loo_age)

qplot(resid_c_female_combined_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(resid_c_female_combined_loo_age)
qqline(resid_c_female_combined_loo_age)


### look at alpha = 0.4
loo_c_female_age_table <- tibble(truth = imputed_c_female_combined_age,
                               pred = pred_c_female_combined_loo_age,
                               resid = truth - pred,
                               apoe = imputed_c_female_combined_apoe,
                               apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33')),
                               type = imputed_c_female_combined_labels
)

#### colored by APOE  ####
#pred vs residuals
ggplot(loo_c_female_age_table) +
  geom_point(aes(pred, resid, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Female Controls: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.3, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") +
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_c_female_loo_7.png')


ggplot(loo_c_female_age_table) +
  geom_point(aes(truth, pred, color = apoe)) +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Female Controls: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.3, loo',
       x = 'True Age',
       y = 'Predicted Age') +
  geom_abline(intercept = 0, slope=1)
ggsave('pred_truth_c_female_loo_7.png')



##########################

### GOT and Lipids Age analysis for controls ###
### With Gender interaction term ###

###########################

## without APOE as a predictor ##
imputed_c_features_combined_age_tmp <- imputed_c_combined_Y %>% 
  as_tibble(.name_repair = 'universal') %>%
  select(-c(Age))

#this gives us a formula with every column being column*gender + column2*gender + ..
# this takes advantage of r adding in the individual terms when given *
gender_interaction_formula <- imputed_c_features_combined_age_tmp %>%
  names %>% 
  #add colname*GenderM to each of the terms
  map_chr(function(x) paste(x, 'GenderM', sep = '*')) %>%
  #removing the last term, which is GenderM * GenderM
  head(-1) %>%
  paste(collapse = '+') %>%
  #model.matrix doens't need anything except ~ on the lhs
  paste0('~',.) %>%
  formula

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_c_features_interaction_combined_age <- model.matrix(gender_interaction_formula, imputed_c_features_combined_age_tmp)

###

#try with list of multiple alphas
fit_c_interaction_combined_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_features_interaction_combined_age, imputed_c_combined_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

#in sample prediction
pred_c_interaction_combined_age_list <- lapply(fit_c_interaction_combined_age_list, 
                                 function(x) predict(x, newx = imputed_c_features_interaction_combined_age, 
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_interaction_combined_age_list <- lapply(fit_c_interaction_combined_age_list, function(x) importance(x, metabolites = FALSE))

mse_c_interaction_combined_age_list <- lapply(pred_c_interaction_combined_age_list, function(x) mean((x - imputed_c_combined_age)^2))

#get the residuals
residuals_c_interaction_combined_age_list <- lapply(pred_c_interaction_combined_age_list, function(x) x - imputed_c_combined_age)

#shapiro-wilkes test tests H0: data is normal
sw_c_interaction_combined_age_list <- lapply(residuals_c_interaction_combined_age_list, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#which has lowest mse?
which.min(mse_c_interaction_combined_age_list)


#c plots look pretty normal!
resid_c_interaction_combined_4 <- residuals_c_interaction_combined_age_list[[5]]
qplot(resid_c_interaction_combined_4, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')


qqnorm(resid_c_interaction_combined_4)
qqline(resid_c_interaction_combined_4)





### look at alpha = 0.4
age_c_interaction_combined_table_4 <- tibble(truth = imputed_c_combined_age, 
                               pred = as.numeric(pred_c_interaction_combined_age_list[[5]]),
                               resid = truth - pred,
                               apoe = imputed_c_combined_apoe,
                               type = imputed_c_combined_labels
)


#pred vs resid
ggplot(age_c_interaction_combined_table_4) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Predicted Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, with gender interaction, alpha = 0.4',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)')
ggsave('age_pred_resid_c_interaction_insample_4.png', width = 7.26, height = 7.26, units = 'in')


#truth vs pred
ggplot(age_c_interaction_combined_table_4) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted (in-sample) Age',
       subtitle = 'Combined GOT and Lipid, with gender interaction, alpha = 0.4',
       x = 'True Age',
       y = 'Predicted Age (in-sample)')
ggsave('age_true_pred_c_interaction_insample_4.png', width = 7.26, height = 7.26, units = 'in')



### Leave one out

#now try with alpha = 0.4
fitpred_c_interaction_combined_loo_age <- lapply(1:nrow(imputed_c_features_interaction_combined_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_interaction_combined_age, imputed_c_combined_age, 
                                                                                                         alpha = 0.4, family = 'gaussian', penalize_age_gender = FALSE))

fit_c_interaction_combined_loo_age <- lapply(fitpred_c_interaction_combined_loo_age, function(x) x[[1]])
pred_c_interaction_combined_loo_age <- lapply(fitpred_c_interaction_combined_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_c_interaction_combined_loo_age <- lapply(fit_c_interaction_combined_loo_age, function(x) importance(x))


#look at ones with gender
importance_c_interaction_combined_loo_age%>% 
  lapply(function(x) names(x) %>% 
           str_extract('GenderM.*') %>%
           na.omit) %>%
  unlist %>%
  table
  
  



mse_c_interaction_combined_loo_age <- mean((pred_c_interaction_combined_loo_age - imputed_c_combined_age)^2)
resid_c_interaction_combined_loo_age <- pred_c_interaction_combined_loo_age - imputed_c_combined_age


shapiro.test(resid_c_interaction_combined_loo_age)

qplot(resid_c_interaction_combined_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.4')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(resid_c_interaction_combined_loo_age)
qqline(resid_c_interaction_combined_loo_age)


### look at alpha = 0.4
loo_c_interaction_age_table <- tibble(truth = imputed_c_combined_age, 
                        pred = pred_c_interaction_combined_loo_age,
                        resid = truth - pred,
                        apoe = imputed_c_combined_apoe,
                        type = imputed_c_combined_labels,
                        gender = imputed_c_combined_gender,
                        apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)

#### colored by APOE  ####
#pred vs residuals
ggplot(loo_c_interaction_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid with gender interaction, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_c_interaction_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_c_interaction_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_interation_control_loo_5.png')


ggplot(loo_c_interaction_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid with gender interaction, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('pred_truth_interaction_control_loo_5.png')


ggplot(loo_c_interaction_age_table) + 
  geom_point(aes(truth, pred, color = gender)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid with gender interaction, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age')

#### #maybe look into group lasso?
#https://stats.stackexchange.com/questions/296581/group-elastic-net



###################

### PCA on combined GOT + Lipids ####

###################

library(ggbiplot)

#change dataset to include co/cy/cm, apoe
  #drop intercept column
imputed_c_combined_full <- do.call(cbind, list(imputed_c_combined_Y, type = imputed_c_combined_labels, apoe = imputed_c_combined_apoe))[,-1]

combined_c_pca <- prcomp(imputed_c_combined_full, scale = T, center = T)
ggbiplot::ggscreeplot(combined_c_pca)
#ggbiplot::ggbiplot(combined_c_pca, obs.scale = 1, choices = 1:2)

combined_c_pca$x %>%
  as_tibble %>%
  ggplot(aes(PC1, PC2)) + 
  geom_point(aes(color = imputed_c_combined_full[,'Age']), size = 3) + 
  scale_color_distiller(type = 'seq') + 
  labs(title = 'First Two Principal Components',
       subtitle = 'Combined GOT + Lipids, age + gender + type + apoe',
       color = 'Age')

combined_c_pca$x %>%
  as_tibble %>%
  ggplot(aes(PC1, PC2)) + 
  geom_point(aes(color = as.factor(imputed_c_combined_full[,'GenderM'])), size = 3) + 
  scale_color_brewer(palette = 'Set2') + 
  labs(title = 'First Two Principal Components',
       subtitle = 'Combined GOT + Lipids, age + gender + type + apoe',
       color = 'Male')

combined_c_pca$x %>%
  as_tibble %>%
  ggplot(aes(PC1, PC2)) + 
  geom_point(aes(color = imputed_c_combined_apoe), size = 3) + 
  scale_color_brewer(palette = 'Set2') + 
  labs(title = 'First Two Principal Components',
       subtitle = 'Combined GOT + Lipids, age + gender + type + apoe',
       color = 'APOE')

combined_c_pca$x %>%
  as_tibble %>%
  ggplot(aes(PC1, PC2)) + 
  geom_point(aes(color = imputed_c_combined_apoe%>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))), size = 3) + 
  scale_color_brewer(palette = 'Set2') + 
  labs(title = 'First Two Principal Components',
       subtitle = 'Combined GOT + Lipids, age + gender + type + apoe',
       color = 'APOE4')







# ggplot() + 
#   geom_line(aes(1:length(combined_c_pca$sdev), combined_c_pca$sdev^2/sum(combined_c_pca$sdev^2))) + 
#   geom_point() + 
#   labs(title = 'Proportion of Variance Explained',
#        subtitle = 'Combined GOT + Lipids + type + apoe + age + gender, all controls',
#        x = 'Principal Component #',
#        y = '% Total Variance Explained')


###################

### Looking at missingness ####

###################


#get sample size by type
n_by_type_comb <- wide_data_combined %>%
  group_by(Type) %>% 
  tally() %>%
  spread(key = 'Type', value = 'n') %>%
  rename_all(function(x) paste('n', x, sep = '_'))

missingness_by_type_comb_counts <- wide_data_combined %>%
  group_by(Type) %>%
  group_map(~ map_int(.x, function(y) sum(is.na(y)))) %>%
  set_names(wide_data_combined$Type %>% levels) %>%
  lapply(function(x) enframe(x, name = 'name', value = 'num_missing')) 



missingness_by_type_all <- reduce(missingness_by_type_comb_counts, inner_join, by = 'name') %>%
  #set the names to show type
  set_names(c('name', paste('num_missing',levels(wide_data_combined$Type), sep = '_'))) %>%
  filter(!(name %in% c('GBAStatus', 'cognitive_status', 'GBA_T369M'))) %>%
  cbind(n_by_type_comb) %>%
  mutate(pct_missing_AD = num_missing_AD/n_AD,
         pct_missing_CM = num_missing_CM/n_CM,
         pct_missing_CO = num_missing_CO/n_CO,
         pct_missing_CY = num_missing_CY/n_CY,
         pct_missing_PD = num_missing_PD/n_PD
         ) %>%
  #filter(reduce(list(pct_missing_AD, pct_missing_CM,pct_missing_CO,pct_missing_CY,pct_missing_PD), `==`)) %>%
  rowwise() %>%
  #mutate(p_value = (prop.test(x =  str_subset(names(.), 'num_missing'), n = str_subset(names(.), 'n_')))$p.value)
  mutate(p_value = (prop.test(x = c(num_missing_AD, num_missing_CM, num_missing_CO, num_missing_CY, num_missing_PD), n = c(n_AD,n_CM,n_CO,n_CY,n_PD)))$p.value) %>%
  cbind('bh_q_value' = p.adjust(.$p_value, method = 'BH')) %>%
  filter(bh_q_value < 0.05) %>%
  gather('Type', 'pct_missing', contains('pct_missing'))

ggplot(missingness_by_type_all, aes(pct_missing, name)) +
  geom_point(aes(color = Type), size = 3, position = position_jitter(width = 0.01, height = 0,seed = 1)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  theme(axis.text.x = element_text(angle= 90, hjust =1)) + 
  labs(title = 'Percent Missingness by Type',
       subtitle = 'GOT + Lipids, filtered using 2-tailed Pearson Chi squared test, BH q< 0.05')
ggsave('plots/combined_pct_missing.png')
  



cy_co_comb_missingness <- inner_join(missingness_by_type_comb_counts$CY, missingness_by_type_comb_counts$CO, by = 'name') %>%
  rename(num_missing_cy = num_missing.x, num_missing_co = num_missing.y) %>%
  filter(num_missing_cy != num_missing_co) %>%
  mutate(n_cy = n_by_type_comb$CY,
         n_co = n_by_type_comb$CO,
         pct_missing_cy = num_missing_cy / n_cy,
         pct_missing_co = num_missing_co / n_co
         ) %>%
  #go rowwise to get p values
  rowwise() %>%
  mutate(p_value = (prop.test(x = c(num_missing_cy, num_missing_co), n= c(n_cy, n_co)))$p.value) %>%
  #using 0.05 p value cutouff
  filter(p_value <= 0.05) %>%
  gather('Type', 'pct_missing', 'cy' = pct_missing_cy, 'co' = pct_missing_co)
  


ggplot(cy_co_comb_missingness) +
  geom_col(aes(name,pct_missing, fill = Type)) + 
  labs(title = 'Percent Missing in GOT/Lipid Data',
       subtitle = 'CO and CY, filtered using ')
  coord_flip() 






missingness_by_type_lipids <- wide_data_lipids %>% 
  group_by(Type) %>%
  group_map(~ map_dbl(.x, function(y) round(sum(is.na(y))/length(y), digits = 3))) %>%
  set_names(wide_data_combined$Type %>% levels) %>%
  lapply(function(x) enframe(x, name = 'lipid', value = 'pct_missing'))

wide_data_lipids %>% group_by(Type) %>% summarize(n = n())

co_cy_lipids_missingness <- inner_join(missingness_by_type_lipids$CO, missingness_by_type_lipids$CY, by = 'lipid') %>%
  rename(pct_missing_co = pct_missing.x, pct_missing_cy = pct_missing.y) %>%
  mutate(co_minus_cy = pct_missing_co - pct_missing_cy)%>%
  arrange(desc(abs(co_minus_cy)))

pd_co_comb_missingness <- inner_join(missingness_by_type_comb$PD, missingness_by_type_comb$CO, by = 'name') %>%
  rename(pct_missing_pd = pct_missing.x, pct_missing_co = pct_missing.y) %>%
  mutate(pd_minus_co = pct_missing_pd - pct_missing_co)%>%
  arrange(desc(abs(pd_minus_co)))

ad_co_comb_missingness <- inner_join(missingness_by_type_comb$AD, missingness_by_type_comb$CO, by = 'name') %>%
  rename(pct_missing_ad = pct_missing.x, pct_missing_co = pct_missing.y) %>%
  mutate(ad_minus_co = pct_missing_ad - pct_missing_co)%>%
  arrange(desc(abs(ad_minus_co)))




missingness_by_type_got <- wide_data %>% 
  group_by(Type) %>%
  group_map(~ map_dbl(.x, function(y) round(sum(is.na(y))/length(y), digits = 3))) %>%
  set_names(wide_data$Type %>% levels) %>%
  lapply(function(x) enframe(x, name = 'metabolite', value = 'pct_missing'))

co_cy_got_missingness <- inner_join(missingness_by_type_got$CO, missingness_by_type_got$CY, by = 'metabolite') %>%
  rename(pct_missing_co = pct_missing.x, pct_missing_cy = pct_missing.y) %>%
  mutate(co_minus_cy = pct_missing_co - pct_missing_cy)%>%
  arrange(desc(abs(co_minus_cy)))







##########################

### Untargeted Age analysis for controls ###
###  ###

###########################

imputed_c_untargeted <- filter_and_impute(wide_data_untargeted_dropped, types = c('CO', 'CY', 'CM'))

imputed_c_untargeted_Y <- imputed_c_untargeted[[1]]
imputed_c_untargeted_labels <- imputed_c_untargeted[[2]]
imputed_c_untargeted_apoe <- imputed_c_untargeted[[3]]
imputed_c_untargeted_gender <- imputed_c_untargeted_Y[,'GenderM'] %>% as.factor %>%
  fct_recode(M = '1', F = '0')





imputed_c_untargeted_age <- imputed_c_untargeted_Y[,'Age']
# no type as a feature in this analysis. it's just metabolites + gender
imputed_c_features_untargeted_age_tmp <- imputed_c_untargeted_Y %>% 
  as_tibble %>%
  #mutate(Type = imputed_c_labels) %>%
  select(-Age)

#transform factor vars into dummies
imputed_c_features_untargeted_age <- model.matrix(~., imputed_c_features_untargeted_age_tmp)


### list to determine alphas
fit_c_untargeted_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_features_untargeted_age, imputed_c_untargeted_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))


#in sample prediction
pred_c_untargeted_age_list <- lapply(fit_c_untargeted_age_list, 
                                 function(x) predict(x, newx = imputed_c_features_untargeted_age, 
                                                     type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_untargeted_age_list <- lapply(fit_c_untargeted_age_list, function(x) importance(x, metabolites = FALSE))

mse_c_untargeted_age_list <- lapply(pred_c_untargeted_age_list, function(x) mean((x - imputed_c_untargeted_age)^2))




### leave one out using chosen alpha from above

fitpred_c_untargeted_loo_age <- lapply(1:nrow(imputed_c_features_untargeted_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_untargeted_age, imputed_c_untargeted_age, 
                                                                                                       alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE))

fit_c_untargeted_loo_age <- lapply(fitpred_c_untargeted_loo_age, function(x) x[[1]])
pred_c_untargeted_loo_age <- lapply(fitpred_c_untargeted_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_c_untargeted_loo_age <- lapply(fit_c_untargeted_loo_age, function(x) importance(x))
mse_c_untargeted_loo_age <- mean((pred_c_untargeted_loo_age - imputed_c_untargeted_age)^2)
resid_c_untargeted_loo_age <- pred_c_untargeted_loo_age - imputed_c_untargeted_age


shapiro.test(resid_c_untargeted_loo_age)

qplot(resid_c_untargeted_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.5')
#ggsave('got_untargeted_age_control_resid_hist.png')
qqnorm(resid_c_untargeted_loo_age)
qqline(resid_c_untargeted_loo_age)


### look at alpha = 0.5
untargeted_c_loo_age_table <- tibble(truth = imputed_c_untargeted_age, 
                                 pred = pred_c_untargeted_loo_age,
                                 resid = truth - pred,
                                 apoe = imputed_c_untargeted_apoe,
                                 type = imputed_c_untargeted_labels,
                                 gender = imputed_c_untargeted_gender,
                                 apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)

#### colored by APOE  ####
#pred vs residuals
ggplot(untargeted_c_loo_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Untargeted, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_untargeted_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_untargeted_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_control_untargeted_loo_5.png')


ggplot(untargeted_c_loo_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'untargeted, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(intercept = 0, slope = 1)
ggsave('pred_truth_control_untargeted_loo_5.png')



##########################

### targeted Age analysis for controls ###
###  ###

###########################

imputed_c_targeted <- filter_and_impute(wide_data_targeted, types = c('CO', 'CY', 'CM'))

imputed_c_targeted_Y <- imputed_c_targeted[[1]]
imputed_c_targeted_labels <- imputed_c_targeted[[2]]
imputed_c_targeted_apoe <- imputed_c_targeted[[3]]
imputed_c_targeted_gender <- imputed_c_targeted_Y[,'GenderM'] %>% as.factor %>%
  fct_recode(M = '1', F = '0')





imputed_c_targeted_age <- imputed_c_targeted_Y[,'Age']
# no type as a feature in this analysis. it's just metabolites + gender
imputed_c_features_targeted_age_tmp <- imputed_c_targeted_Y %>% 
  as_tibble %>%
  #mutate(Type = imputed_c_labels) %>%
  select(-Age)

#transform factor vars into dummies
imputed_c_features_targeted_age <- model.matrix(~., imputed_c_features_targeted_age_tmp)


### list to determine alphas
fit_c_targeted_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_c_features_targeted_age, imputed_c_targeted_age, 
                                                                        family = 'gaussian', alpha = x, penalize_age_gender = FALSE))


#in sample prediction
pred_c_targeted_age_list <- lapply(fit_c_targeted_age_list, 
                                     function(x) predict(x, newx = imputed_c_features_targeted_age, 
                                                         type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_c_targeted_age_list <- lapply(fit_c_targeted_age_list, function(x) importance(x, metabolites = FALSE))

importance_c_targeted_age_names <- importance_c_targeted_age_list[[6]] %>% 
  enframe(name = 'name', value= 'coefficient') %>%
  filter(name != '(Intercept)') %>%
  select(name) %>% deframe


mse_c_targeted_age_list <- lapply(pred_c_targeted_age_list, function(x) mean((x - imputed_c_targeted_age)^2))




### leave one out using chosen alpha from above

fitpred_c_targeted_loo_age <- lapply(1:nrow(imputed_c_features_targeted_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_targeted_age, imputed_c_targeted_age, 
                                                                                                               alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE))

fit_c_targeted_loo_age <- lapply(fitpred_c_targeted_loo_age, function(x) x[[1]])
pred_c_targeted_loo_age <- lapply(fitpred_c_targeted_loo_age, function(x) x[[2]]) %>%
  unlist

#some measure of variable importance
importance_c_targeted_loo_age <- lapply(fit_c_targeted_loo_age, function(x) importance(x))
mse_c_targeted_loo_age <- mean((pred_c_targeted_loo_age - imputed_c_targeted_age)^2)
resid_c_targeted_loo_age <- pred_c_targeted_loo_age - imputed_c_targeted_age


shapiro.test(resid_c_targeted_loo_age)

qplot(resid_c_targeted_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.5')
#ggsave('got_targeted_age_control_resid_hist.png')
qqnorm(resid_c_targeted_loo_age)
qqline(resid_c_targeted_loo_age)


### look at alpha = 0.5
targeted_c_loo_age_table <- tibble(truth = imputed_c_targeted_age, 
                                     pred = pred_c_targeted_loo_age,
                                     resid = truth - pred,
                                     apoe = imputed_c_targeted_apoe,
                                     type = imputed_c_targeted_labels,
                                     gender = imputed_c_targeted_gender,
                                     apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
)

#### colored by APOE  ####
#pred vs residuals
ggplot(targeted_c_loo_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Targeted, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_targeted_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_targeted_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')
ggsave('pred_age_residuals_control_targeted_loo_5.png')


ggplot(targeted_c_loo_age_table) + 
  geom_point(aes(truth, pred, color = gender)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Targeted, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_text(x = 80, y = 20, label = paste0('MAE:', abs(targeted_c_loo_age_table$truth - targeted_c_loo_age_table$pred) %>% mean %>% round(3))) 
ggsave('pred_truth_control_targeted_loo_5.png')









#########################

### univariate regression on MISSINGNESS ~ Age + Gender for each metabolite ###

########################

#' Returns list of overall p value, p value for gender, p value for age only.
#' @param data is must have some NA
#' @param form is a formula for the logistic regression
#' @param metabolite is lipid/metabolite, as a string
missingness_metabolite_p <- function(data, form, metabolite){
  metabolite <- sym(metabolite)
  df <- data %>% 
    mutate(missing = ifelse(is.na(eval(metabolite)), 1, 0))
  
  fit <- glm(form, data = df, family = 'binomial', maxit = 100)
  null_fit <- glm(missing ~ 1, data = df, family = 'binomial', maxit = 100)
  
  #need to manually get an overall p value from the logistic regression, by comparing dispersion to null model
  compare <- anova(fit, null_fit, test = 'Chisq')
  
  overall_p <- compare$`Pr(>Chi)`[2]
  # -1 removes intercept
  individual_p <- summary(fit)$coefficients[-1,'Pr(>|z|)'] %>% as.list
  c('overall' = overall_p, individual_p)
  
}



missingness_some_na_names <- wide_data_combined %>% 
  #leaving only metabolites/lipids
  dplyr::select(-one_of("Type2", "Type", "Gender", "Age", "APOE", "Batch",
                        "Index", "GBAStatus", "Id", 
                        "GBA_T369M", "cognitive_status")) %>%
  #need at least one NA, but not all NA
  dplyr::select_if(function(x) any(is.na(x)) & !all(is.na(x))) %>%
  names

form_missing_genderAge <- formula(missing~Gender + Age)
missingness_overall_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined,form_missing_genderAge, .x)[[1]]) %>%
  unlist
missingness_gender_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined,form_missing_genderAge, .x)[[2]]) %>%
  unlist
missingness_age_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined,form_missing_genderAge, .x)[[3]]) %>%
  unlist


#bh corrected controls for false discovery rate.
missingness_p_table <- bind_cols('name' = missingness_some_na_names, 
                                 'og_overall_p_value' = missingness_overall_p_values,
                                 'og_gender_p_value' = missingness_gender_p_values,
                                 'og_age_p_value' = missingness_age_p_values,
                                 'bh_overall_q_value' = p.adjust(missingness_overall_p_values, method = 'BH'),
                                 'bh_gender_q_value' = p.adjust(missingness_gender_p_values, method = 'BH'),
                                 'bh_age_q_value' = p.adjust(missingness_age_p_values, method = 'BH')
                                 ) 

missingness_p_table$name <- if_else(!str_detect(missingness_p_table$name, 'Result'), #take advantage of fact that all metabolites have "results" in name
  missingness_p_table$name, 
  str_replace_all(missingness_p_table$name, 'Result.*', "") %>%
    str_replace_all('\\.$', '') %>%
    str_replace_all('^\\.+', '') %>%
    str_trim() %>%
    sapply(function(x) all_matches[match(x, all_matches$Name), 'Metabolite'] %>% deframe))


  
missingness_p_table %>% 
  arrange(bh_q_value) %>%
  #filter(bh_q_value < 0.05) %>%
  write_csv('combined_missingness_p.csv')


missingness_p_table %>%
  arrange(bh_gender_q_value) 

bh_missingness_age <- missingness_p_table %>%
  arrange(bh_age_q_value) %>%
  filter(bh_age_q_value < 0.01)

bh_missingness_overall <- missingness_p_table %>%
  arrange(bh_overall_q_value) %>%
  filter(bh_overall_q_value < 0.01)

intersect(bh_missingness_age$name, bh_missingness_overall$name)


#########################

### univariate regression on MISSINGNESS ~ Age +Type for each metabolite ###

########################


# indicators for AD, PD
wide_data_combined_adpd <- wide_data_combined %>%
  mutate(typeAD = ifelse(Type == 'AD', 1,0),
         typePD = ifelse(Type == 'PD', 1, 0))


formula_missing_ageType <- formula(missing ~ Age + typeAD + typePD)
missingness_overall_ageType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_ageType, .x)[[1]]) %>%
  unlist
missingness_age_ageType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_ageType, .x)[[2]]) %>%
  unlist
missingness_ad_ageType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_ageType, .x)[[3]]) %>%
  unlist
missingness_pd_ageType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_ageType, .x)[[4]]) %>%
  unlist



#bh corrected controls for false discovery rate.
missingness_ageType_p_table <- bind_cols('name' = missingness_some_na_names, 
                                 'og_overall_p_value' = missingness_overall_ageType_p_values,
                                 'og_age_p_value' = missingness_age_ageType_p_values,
                                 'og_ad_p_value' = missingness_ad_ageType_p_values,
                                 'og_pd_p_value' = missingness_pd_ageType_p_values,
                                 'bh_overall_q_value' = p.adjust(missingness_overall_ageType_p_values, method = 'BH'),
                                 'bh_age_q_value' = p.adjust(missingness_age_ageType_p_values, method = 'BH'),
                                 'bh_ad_q_value' = p.adjust(missingness_ad_ageType_p_values, method = 'BH'),
                                 'bh_pd_q_value' = p.adjust(missingness_pd_ageType_p_values, method = 'BH')
) 

missingness_ageType_p_table$name <- if_else(!str_detect(missingness_ageType_p_table$name, 'Result'), #take advantage of fact that all metabolites have "results" in name
                                    missingness_ageType_p_table$name, 
                                    str_replace_all(missingness_ageType_p_table$name, 'Result.*', "") %>%
                                      str_replace_all('\\.$', '') %>%
                                      str_replace_all('^\\.+', '') %>%
                                      str_trim() %>%
                                      sapply(function(x) all_matches[match(x, all_matches$Name), 'Metabolite'] %>% deframe))



missingness_ageType_p_table %>% 
  arrange(bh_overall_q_value) %>%
  #filter(bh_q_value < 0.05) %>%
  write_csv('combined_missingness_ageType_p.csv')

missingness_ageType_p_table %>%
  arrange(bh_pd_q_value)

missingness_ageType_p_table %>%
  arrange(bh_ad_q_value)

missingness_ageType_p_table %>%
  filter(bh_ad_q_value < 0.05)

a <- missingness_ageType_p_table %>%
  filter(bh_age_q_value < 0.01)

b <- missingness_ageType_p_table %>%
  filter(bh_overall_q_value < 0.01)



  


#########################

### univariate regression on MISSINGNESS ~ Gender +Type for each metabolite ###

########################


# indicators for AD, PD

formula_missing_genderType <- formula(missing ~ Gender + typeAD + typePD)
missingness_overall_genderType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_genderType, .x)[[1]]) %>%
  unlist
missingness_gender_genderType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_genderType, .x)[[2]]) %>%
  unlist
missingness_ad_genderType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_genderType, .x)[[3]]) %>%
  unlist
missingness_pd_genderType_p_values <- purrr::map(missingness_some_na_names, ~missingness_metabolite_p(wide_data_combined_adpd, formula_missing_genderType, .x)[[4]]) %>%
  unlist



#bh corrected controls for false discovery rate.
missingness_genderType_p_table <- bind_cols('name' = missingness_some_na_names, 
                                         'og_overall_p_value' = missingness_overall_genderType_p_values,
                                         'og_gender_p_value' = missingness_gender_genderType_p_values,
                                         'og_ad_p_value' = missingness_ad_genderType_p_values,
                                         'og_pd_p_value' = missingness_pd_genderType_p_values,
                                         'bh_overall_q_value' = p.adjust(missingness_overall_genderType_p_values, method = 'BH'),
                                         'bh_gender_q_value' = p.adjust(missingness_gender_genderType_p_values, method = 'BH'),
                                         'bh_ad_q_value' = p.adjust(missingness_ad_genderType_p_values, method = 'BH'),
                                         'bh_pd_q_value' = p.adjust(missingness_pd_genderType_p_values, method = 'BH')
) 

missingness_genderType_p_table$name <- if_else(!str_detect(missingness_genderType_p_table$name, 'Result'), #take advantgender of fact that all metabolites have "results" in name
                                            missingness_genderType_p_table$name, 
                                            str_replace_all(missingness_genderType_p_table$name, 'Result.*', "") %>%
                                              str_replace_all('\\.$', '') %>%
                                              str_replace_all('^\\.+', '') %>%
                                              str_trim() %>%
                                              sapply(function(x) all_matches[match(x, all_matches$Name), 'Metabolite'] %>% deframe))



missingness_genderType_p_table %>% 
  arrange(bh_overall_q_value) %>%
  #filter(bh_q_value < 0.05) %>%
  write_csv('combined_missingness_genderType_p.csv')

missingness_genderType_p_table %>% lapply(min)

missingness_genderType_p_table %>%
  arrange(bh_overall_q_value) %>%
  filter(bh_overall_q_value < 0.1)
  




#########################

### univariate regression on MISSINGNESS ~ Age + Gender for each metabolite ###
## on Untargeted ##

########################


missingness_some_na_names_untargeted <- wide_data_untargeted %>% 
  #leaving only metabolites/lipids
  dplyr::select(-one_of("Type", "Gender", "Age", "APOE", "GBAStatus", "Id", 
                        "GBA_T369M")) %>%
  #need at least one NA, but not all NA
  dplyr::select_if(function(x) any(is.na(x)) & !all(is.na(x))) %>%
  names

missingness_overall_p_values_untargeted <- map(missingness_some_na_names_untargeted, ~missingness_metabolite_p(wide_data_untargeted, form, .x)[[1]]) %>%
  unlist
missingness_age_p_values_untargeted <- map(missingness_some_na_names_untargeted, ~missingness_metabolite_p(wide_data_untargeted, form, .x)[[2]]) %>%
  unlist
missingness_gender_p_values_untargeted <- map(missingness_some_na_names_untargeted, ~missingness_metabolite_p(wide_data_untargeted, form, .x)[[3]]) %>%
  unlist


#bh corrected controls for false discovery rate.
missingness_p_table_untargeted <- bind_cols('name' = missingness_some_na_names_untargeted, 
                                 'og_p_value' = missingness_p_values_untargeted,
                                 'bh_q_value' = p.adjust(missingness_p_values_untargeted, method = 'BH'))



missingness_p_table_untargeted %>% 
  arrange(bh_q_value) %>%
  #filter(bh_q_value < 0.05) %>%
  write_csv('untargeted_missingness_p.csv')




#########################

### univariate regression on Age ~ Metabolite for each metabolite ###

########################

age_metabolite_p <- function(data, metabolite){
  metabolite <- sym(metabolite)
  df <- data %>%
    select(Age, !!metabolite)
  fit <- glm(Age ~ ., data = df)
  
  #get the t-value and p score (excluding intercept)
  tryCatch({
    enframe(summary(fit)$coefficients[-1,3:4]) %>% spread(key = name, value = value) %>%
    cbind('name' = rlang::as_string(metabolite))
  },
  error = function(w) cat('metabolite is ', metabolite)
  )

}

imputed_all_combined_df <- imputed_all_combined[[1]] %>% 
  as_tibble %>%
  select(-c('(Intercept)', 'GenderM'))

 age_metabolite_p_values <- imputed_all_combined_df %>% 
  names %>%
   setdiff('Age') %>%
   purrr::map(~age_metabolite_p(imputed_all_combined_df, .x)) %>%
   purrr::reduce(rbind) %>%
   as_tibble() %>%
   rename('og_p_value' = `Pr(>|t|)`)
 
 age_metabolite_p_table <- cbind(age_metabolite_p_values, 
                                 'bh_p_value' = p.adjust(age_metabolite_p_values$og_p_value, method = 'BH'))
 
 
age_metabolite_p_table %>%
  write_csv('combined_age_metabolite_p.csv')


#### Now do the same for controls only
imputed_c_combined_df <- imputed_c_combined_Y %>%
  as_tibble() %>%
  select(-c('(Intercept)', GenderM))

age_c_combined_p_values <- imputed_c_combined_df %>%
  names %>%
  setdiff('Age') %>%
  purrr::map(~age_metabolite_p(imputed_c_combined_df, .x)) %>%
  purrr::reduce(rbind) %>%
  as_tibble() %>%
  dplyr::rename('og_p_value' =  `Pr(>|t|)`)


age_c_combined_p_table <- cbind(age_c_combined_p_values, 
                              'bh_p_value' = p.adjust(age_c_combined_p_values$og_p_value, method = 'BH'))

age_combined_bh_cut <- age_combined_p_table %>%
  mutate(name = str_replace_all(name, '`', ''))
  filter(bh_p_value < 0.01)
  
  
## for mummichog
#age_combined_
















#### now do the same for untargeted
imputed_c_untargeted_df <- imputed_c_untargeted_Y %>%
  as_tibble() %>%
  select(-c('(Intercept)', GenderM))

age_untargeted_p_values <- imputed_c_untargeted_df %>%
  names %>%
  setdiff('Age') %>%
  purrr::map(~age_metabolite_p(imputed_c_untargeted_df, .x)) %>%
  purrr::reduce(rbind) %>%
  as_tibble() %>%
  dplyr::rename('og_p_value' =  `Pr(>|t|)`) %>%
  dplyr::mutate(name = str_replace_all(name, "`", "")) %>%
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2)

age_untargeted_p_table <- cbind(age_untargeted_p_values, 
                                'bh_p_value' = p.adjust(age_untargeted_p_values$og_p_value, method = 'BH'))


mz_retention_untargeted <- raw_data_untargeted %>%
  group_by(Metabolite, Mode) %>%
  slice(1) %>%
  select(Metabolite, Mode, `m/z`, `Retention time (min)`)

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_untargeted <- age_untargeted_p_table %>% 
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  filter(bh_p_value < 0.05) %>%
  select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 't-score' = `t value`, Metabolite, Mode)

mummichog_untargeted %>% 
  filter(Mode == 'neg') %>%
  write_tsv('mummichog_untargeted_neg.txt')

mummichog_untargeted %>% 
  filter(Mode == 'pos') %>%
  write_tsv('mummichog_untargeted_pos.txt')



##### now do the same for lipids only
imputed_c_lipids_df <- imputed_c_lipids_Y %>%
 as_tibble() %>%
  select(-c('(Intercept)', GenderM))
age_lipids_p_values <- imputed_c_lipids_df %>%
  names %>%
  map(~age_metabolite_p(imputed_c_lipids_df, .x)) %>%
  unlist %>%
  enframe(value = 'og_p_value')

age_lipids_p_table <- cbind(age_lipids_p_values, 
                                'bh_p_value' = p.adjust(age_lipids_p_values$og_p_value, method = 'BH'))
age_lipids_p_table %>% 
  filter(bh_p_value < 0.05)


#### now do the same for got only
imputed_c_got_df <- imputed_c_got_Y %>%
  as_tibble() %>%
  select(-c('(Intercept)', GenderM))
age_got_p_values <- imputed_c_got_df %>%
  names %>%
  map(~age_metabolite_p(imputed_c_got_df, .x)) %>%
  unlist %>%
  enframe(value = 'og_p_value')

age_got_p_table <- cbind(age_got_p_values, 
                            'bh_p_value' = p.adjust(age_got_p_values$og_p_value, method = 'BH'))
age_got_p_table %>% 
  filter(bh_p_value < 0.05)



#### now do the same for targeted
imputed_c_targeted_df <- imputed_c_targeted_Y %>%
  as_tibble() %>%
  select(-c('(Intercept)', GenderM))

age_targeted_p_values <- imputed_c_targeted_df %>%
  names %>%
  setdiff('Age') %>%
  purrr::map(~age_metabolite_p(imputed_c_targeted_df, .x)) %>%
  purrr::reduce(rbind) %>%
  as_tibble() %>%
  dplyr::rename('og_p_value' =  `Pr(>|t|)`) %>%
  dplyr::mutate(name = str_replace_all(name, "`", "")) #%>%
  #tidyr::separate(col = name, into = c('Metabolite', 'Mode'), sep = '_')
  

age_targeted_p_table <- cbind(age_targeted_p_values, 
                                'bh_p_value' = p.adjust(age_targeted_p_values$og_p_value, method = 'BH'))

univar_targeted_names <- age_targeted_p_table %>%
  filter(bh_p_value < 0.01) %>% select(name) %>% deframe

setdiff(importance_c_targeted_age_names, univar_targeted_names)
setdiff(univar_targeted_names, importance_c_targeted_age_names)
intersect(univar_targeted_names, importance_c_targeted_age_names)

# look at name overlap between here and the .5 glmnet







#########################

### univariate regression on Gender ~ Metabolite for each metabolite/lipid ###

########################


imputed_c_combined_gender_df <- imputed_c_combined_Y %>%
  as_tibble() %>%
  rename_all(function(x) str_replace_all(x, "`", "")) %>%
  select(-c('(Intercept)', Age))

gender_metabolite_p <- function(data, metabolite){
  metabolite <- sym(metabolite)
  formula <- formula(paste0("GenderM ~ `", metabolite, '`'))
  fit <- glm(formula, data = data, family = 'binomial')
  
  summary(fit)$coefficients[,'Pr(>|z|)'][2]
}


metabolite_lipid_names <- names(imputed_c_combined_gender_df) %>% setdiff('GenderM')
gender_p_values <- purrr::map(metabolite_lipid_names, ~gender_metabolite_p(imputed_c_combined_gender_df, .x)) %>%
  unlist

#bh corrected controls for false discovery rate.
gender_p_table <- bind_cols('name' = metabolite_lipid_names, 
                            'og_p_value' = gender_p_values,
                            'bh_q_value' = p.adjust(gender_p_values, method = 'BH')) 

gender_p_table$name <- if_else(!str_detect(gender_p_table$name, 'Result'), #take advantage of fact that all metabolites have "results" in name
                                    gender_p_table$name, 
                                    str_replace_all(gender_p_table$name, 'Result.*', "") %>%
                                      str_replace_all('\\.$', '') %>%
                                      str_replace_all('^\\.+', '') %>%
                                      str_trim() %>%
                                      sapply(function(x) all_matches[match(x, all_matches$Name), 'Metabolite'] %>% deframe))



gender_p_table %>% 
  arrange(bh_q_value) %>%
  #filter(bh_q_value < 0.05) %>%
  write_csv('combined_gender_p.csv')





### maybe play around with metaboanalyst package
