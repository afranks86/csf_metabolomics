source('analysis/starter.R')

### Start AGE analysis ###
## use previously created imputation for total dataset
#impute all types using amelia for an example
  #Y holds features (keeping with prior notation), labels holds type
imputed_all <- filter_and_impute(wide_data,all_types)
imputed_all_Y <- imputed_all[[1]]
imputed_all_labels <- imputed_all[[2]]



imputed_all_age <- imputed_all_Y[,'Age']
# readd type as a feature in this analysis
imputed_all_features_age_tmp <- imputed_all_Y %>% 
  as_tibble %>%
  mutate(Type = imputed_all_labels) %>%
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
  #it's lasso, which is not what we want.
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




## try with loo method to generate standard error. try with the 9th again (alpha = 0.8)
loo_pred_8 <- sapply(foldid, function(x) loo_pred_glmnet(lambda = fit_age_list[[9]]$lambda.min, index = x, features =imputed_all_features_age, labels = imputed_all_age, alpha = 0.8, penalize_age_gender = FALSE, family= 'gaussian'))

#see (standard error of esimate) = sqrt(SS residuals / df), where df = n - number of coefficients in the model (including intercept)
  #how to do standard error if df > n? i do rmse instead 
rmse_loo_pred_8 <- (sum((loo_pred_8 - imputed_all_age)^2) / length(imputed_all_age)) %>% sqrt








#effective degrees








### End AGE analysis ###






