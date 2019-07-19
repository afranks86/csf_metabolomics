source('analysis/starter.R')

### Start AGE analysis ###

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


#now try with list of multiple alphas
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




### look at alpha = 0.3
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








### End AGE analysis ###
















### Combined analysis
set.seed(1)
imputed_all_combined <- filter_and_impute(wide_data_combined,c('CO', 'CM', 'CY'))
imputed_all_combined_Y <- imputed_all_combined[[1]]
imputed_all_combined_labels <- imputed_all_combined[[2]]
imputed_all_combined_apoe <- imputed_all_combined[[3]]



imputed_all_combined_age <- imputed_all_combined_Y[,'Age']
# readd type as a feature in this analysis
imputed_all_features_combined_age_tmp <- imputed_all_combined_Y %>% 
  as_tibble %>%
  mutate(APOE = imputed_all_combined_apoe) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_all_features_combined_age_tmp<- model.matrix(~., imputed_all_features_combined_age_tmp)
#remove intercept column created by model.matrix
imputed_all_features_combined_age <- imputed_all_features_combined_age_tmp[,-1]

foldid <- sample(nrow(imputed_all_features_combined_age))

#now try with list of multiple alphas
fit_combined_age_list <- lapply(seq(0,1,.1), function(x) fit_glmnet(imputed_all_features_combined_age, imputed_all_combined_age, 
                                                                    family = 'gaussian', alpha = x, penalize_age_gender = FALSE))

fit_combined_age_list
# fit_age_list <- lapply(seq(0,1,.1), function(x) cv.glmnet(imputed_all_features_age, imputed_all_age, family = 'gaussian', 
#                                                           type.measure = 'mse', nfolds = nrow(imputed_all_features_age),
#                                                           foldid = foldid, alpha = x, standardize = TRUE))
pred_combined_age_list <- lapply(fit_combined_age_list, 
                        function(x) predict(x, newx = imputed_all_features_combined_age, 
                                            type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_combined_age_list <- lapply(fit_combined_age_list, function(x) importance(x))




mse_combined_age_list <- lapply(pred_combined_age_list, function(x) mean((x - imputed_all_combined_age)^2))

#which has lowest mse?
#it's alpha = 0.8, but same for .7, .6
which.min(mse_combined_age_list)

residuals_combined_age_list <- lapply(pred_combined_age_list, function(x) x - imputed_all_combined_age)


#shapiro-wilkes test tests H0: data is normal
sw_combined_age_list <- lapply(residuals_combined_age_list, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#all plots look pretty normal!
resid_combined_3 <- residuals_combined_age_list[[4]]

qplot(resid_combined_3, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.3')
ggsave('got_lipids_age_control_resid_hist.png')


abline(h = 0)

hist(resid_combined_3)

qqnorm(resid_combined_3)
qqline(resid_combined_3)





### look at alpha = 0.4
age_combined_table_3 <- tibble(truth = imputed_all_combined_age, 
                               pred = as.numeric(pred_combined_age_list[[4]]),
                               resid = truth - pred,
                               apoe = imputed_all_combined_apoe,
                               type = imputed_all_combined_labels
                               )

#Truth vs rsiduals
ggplot(age_combined_table_3) + 
  geom_point(aes(truth, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.3',
       x = 'True Age',
       y = 'Residuals (Truth - Pred)') #+ 
  #stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
  #stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')

ggsave('got_lipids_age_control_resid.png')



#pred vs resid
ggplot(age_combined_table_3) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Predicted Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.3',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)')
ggsave('pred_age_residuals_control_3.png')


ggplot(age_combined_table_3) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.3',
       x = 'True Age',
       y = 'Predicted Age')
ggsave('age_true_pred_controls.png')


age_combined_table_3 %>% 
  mutate(type = fct_relevel(type, c('CY', 'CM', 'CO'))) %>%
  select(apoe, type) %>%
  table %>% View


fit <- lm(pred~apoe+ truth, data = age_combined_table_3)
summary(fit)

fit2 <- lm(pred ~ apoe, data = age_combined_table_3)
summary(fit2)







#### Try leave one out #####
imputed_all_combined <- filter_and_impute(wide_data_combined,c('CO', 'CM', 'CY'))
imputed_all_combined_Y <- imputed_all_combined[[1]]
imputed_all_combined_labels <- imputed_all_combined[[2]]
imputed_all_combined_apoe <- imputed_all_combined[[3]]



imputed_all_combined_age <- imputed_all_combined_Y[,'Age']
# readd type as a feature in this analysis
imputed_all_features_combined_age_tmp <- imputed_all_combined_Y %>% 
  as_tibble %>%
  mutate(APOE = imputed_all_combined_apoe) %>%
  select(-Age)

#turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
imputed_all_features_combined_age_tmp<- model.matrix(~., imputed_all_features_combined_age_tmp)
#remove intercept column created by model.matrix
imputed_all_features_combined_age <- imputed_all_features_combined_age_tmp[,-1]

foldid <- sample(nrow(imputed_all_features_combined_age))

#now try with alpha = 0.5
fitpred_loo_age <- lapply(foldid, function(x) loo_cvfit_glmnet(x, imputed_all_features_combined_age, imputed_all_combined_age, 
                                                           alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE))

fit_loo_age <- lapply(fit_loo_age, function(x) x[[1]])
pred_loo_age <- lapply(fit_loo_age, function(x) x[[2]]) %>%
  unlist



#some measure of variable importance
importance_combined_age_list <- lapply(fit_combined_age_list, function(x) importance(x))
mse_loo_age <- mean((pred_loo_age - imputed_all_combined_age)^2)
residuals_loo_age <- pred_loo_age - imputed_all_combined_age


shapiro.test(residuals_loo_age)

qplot(residuals_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.3')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(residuals_loo_age)
qqline(residuals_loo_age)


### look at alpha = 0.4
loo_age_table <- tibble(truth = imputed_all_combined_age, 
                               pred = pred_loo_age,
                               resid = truth - pred,
                               apoe = imputed_all_combined_apoe,
                               type = imputed_all_combined_labels
)

#Truth vs rsiduals
ggplot(loo_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')


ggplot(loo_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age')

#####






##### Results are weird, so try loo with fixed lambda

#take a min lambda from one of the fits using the whole data (alpha = .6)
lambda <- fit_combined_age_list[[7]]$lambda.min
pred_loo_age_onelambda <- lapply(foldid, function(x) loo_pred_glmnet(lambda, x, imputed_all_features_combined_age, imputed_all_combined_age, alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE)) %>%
  unlist


#some measure of variable importance
mse_loo_age_onelambda <- mean((pred_loo_age_onelambda - imputed_all_combined_age)^2)
residuals_loo_age_onelambda <- pred_loo_age_onelambda - imputed_all_combined_age


shapiro.test(residuals_loo_age_onelambda)

qplot(residuals_loo_age_onelambda, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.5')
#ggsave('got_lipids_age_onelambda_control_resid_hist.png')
qqnorm(residuals_loo_age_onelambda)
qqline(residuals_loo_age_onelambda)


### look at alpha = 0.4
loo_age_onelambda_table <- tibble(truth = imputed_all_combined_age, 
                        pred = pred_loo_age_onelambda,
                        resid = truth - pred,
                        apoe = imputed_all_combined_apoe,
                        type = imputed_all_combined_labels
)

#Truth vs rsiduals
ggplot(loo_age_onelambda_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo, lambda = 2.4',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_onelambda_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_onelambda_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')


ggplot(loo_age_onelambda_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age')











#### results are still weird so try 80/20 
set.seed(1)
age_random_test_index <- sample(nrow(imputed_all_features_combined_age), size = floor(nrow(imputed_all_features_combined_age)*.2))
age_y_train <- imputed_all_features_combined_age[-age_random_test_index,]
age_label_train <- imputed_all_combined_age[-age_random_test_index]
age_apoe_train <- imputed_all_combined_apoe[-age_random_test_index]


age_y_test <- imputed_all_features_combined_age[age_random_test_index,]
age_label_test <- imputed_all_combined_age[age_random_test_index]
age_apoe_test <- imputed_all_combined_apoe[age_random_test_index]


fit_list_8020  <- lapply(seq(0,1,.1), function(x) fit_glmnet(age_y_train, age_label_train, family = 'gaussian', alpha = x, penalize_age_gender = FALSE))
pred_list_80 <- lapply(fit_list_8020, 
                                 function(x) predict(x, newx = age_y_train, 
                                                     type = 'response', s = 'lambda.min'))
pred_list_20 <- lapply(fit_list_8020, 
                       function(x) predict(x, newx = age_y_test, 
                                           type = 'response', s = 'lambda.min'))

importance_list_8020 <- lapply(fit_list_8020, function(x) importance(x))
mse_list_80 <- lapply(pred_list_80, function(x) mean((x - age_label_train)^2))
mse_list_20 <- lapply(pred_list_20, function(x) mean((x - age_label_test)^2))

#which has lowest mse?
#it's alpha = 0.8, but same for .7, .6
which.min(mse_list_80)

residuals_list_80 <- lapply(pred_list_80, function(x) x - age_label_train)
residuals_list_20 <- lapply(pred_list_20, function(x) x - age_label_test)

#shapiro-wilkes test tests H0: data is normal
sw_list_80 <- lapply(residuals_list_80, function(x) (shapiro.test(x))$p.value)
sw_list_20 <- lapply(residuals_list_20, function(x) (shapiro.test(x))$p.value)

#performance is similar across the alphas, so let's just pick one pretty arb
#all plots look pretty normal!
resid_3_80 <- residuals_list_80[[4]]
qplot(resid_3_80, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.3')
qqnorm(resid_3_80)
qqline(resid_3_80)





### look at alpha = 0.3
age_table_3_80 <- tibble(truth = age_label_train, 
                               pred = as.numeric(pred_list_80[[4]]),
                               resid = truth - pred,
                               apoe = age_apoe_train#,
                               #type = imputed_all_combined_labels
)

age_table_3_20 <- tibble(truth = age_label_test,
                         pred = as.numeric(pred_list_20[[4]]),
                         resid = truth - pred,
                         apoe = age_apoe_test)

#Truth vs rsiduals
ggplot(age_table_3_80) + 
  geom_point(aes(truth, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.3',
       x = 'True Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')

age_apoe_plot(age_table_3_80, 'pred', 'resid')
age_apoe_plot(age_table_3_20, 'pred', 'resid')

age_apoe_plot(age_table_3_80, 'truth', 'pred') + 
  geom_abline(intercept = 0, slope = 1)
age_apoe_plot(age_table_3_20, 'truth', 'pred') + 
  geom_abline(intercept = 0, slope = 1)








#predict on 1
pred <- predict(fit, newx = new_features, type = 'response', s = 'lambda.min')
return(list(fit, pred))


#with the function
fitpred_age_8020 <- loo_cvfit_glmnet(age_random_test_index, imputed_all_features_combined_age, imputed_all_combined_age, 
                                     alpha = 0.5, penalize_age_gender = FALSE, family = 'gaussian')
fit_8020_age <- fitpred_age_8020[[1]]
#this is out of sample prediction
pred_8020_age <- fitpred_age_8020[[2]]



mse_8020_age <- mean((pred_8020_age- age_label_test)^2)
residuals_8020_age <- pred_8020_age - imputed_all_combined_age


shapiro.test(residuals_8020_age)

qplot(residuals_8020_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.5')
#ggsave('got_lipids_age_control_resid_hist.png')
qqnorm(residuals_8020_age)
qqline(residuals_8020_age)


### 8020k at alpha = 0.4
8020_age_table <- tibble(truth = imputed_all_combined_age, 
                                  pred = pred_8020_age,
                                  resid = truth - pred,
                                  apoe = imputed_all_combined_apoe,
                                  type = imputed_all_combined_labels
)

#Truth vs rsiduals
ggplot(8020_age_table) + 
  geom_point(aes(pred, resid, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: Age vs Residuals',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, 8020, lambda = 2.4',
       x = 'Predicted Age',
       y = 'Residuals (Truth - Pred)') #+ 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CO'), aes(truth, resid), size=1, colour="red") + 
#stat_ellipse(data = filter(age_combined_table_4, type == 'CM'), aes(truth, resid), size = 1, color = 'blue')


ggplot(8020_age_table) + 
  geom_point(aes(truth, pred, color = apoe)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'Combined GOT and Lipid, alpha = 0.5, 8020',
       x = 'True Age',
       y = 'Predicted Age')











pred_8020_age_test
