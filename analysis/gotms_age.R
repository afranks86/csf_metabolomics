source('analysis/gotms_glmnet.R')

### Start AGE analysis ###
## use previously created imputation for total dataset
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
                           foldid = foldid, alpha = 0, standardize = FALSE, grouped = FALSE)

pred_age_ridge <- predict(fit_age_ridge, newx = imputed_all_features_age, s= 'lambda.min')

mse_age_ridge <- (pred_age_ridge - imputed_all_age)^2 %>% mean
mae_age_ridge <- (pred_age_ridge - imputed_all_age) %>% 
  abs %>%
  mean

importance_age_ridge <- importance(fit_age_ridge)


#now try with list of multiple alphas
fit_age_list <- lapply(seq(0,1,.1), function(x) cv.glmnet(imputed_all_features_age, imputed_all_age, family = 'gaussian', 
                                                          type.measure = 'mse', nfolds = nrow(imputed_all_features_age),
                                                          foldid = foldid, alpha = x, standardize = FALSE, grouped = FALSE))
pred_age_list <- lapply(fit_age_list, 
                        function(x) predict(x, newx = imputed_all_features_age, 
                                            type = 'response', s = 'lambda.min'))
#some measure of variable importance
importance_age_list <- lapply(fit_age_list, function(x) importance(x))


### End AGE analysis ###

