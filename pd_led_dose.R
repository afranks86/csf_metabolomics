panuc_data <- readxl::read_xlsx(file.path(data_path, 'data', 'PANUC', 'panuc-0133-2019_07_09.xlsx'))

wide_data_untargeted_panuc <- panuc_data %>% 
  select(Id = subject_id, led, no_meds_reported) %>%
  inner_join(wide_data_untargeted, by = "Id") %>%
  #drop na for now  (8 obs)
  filter(!is.na(led))


imputed_c_untargeted5_panuc <- filter_and_impute_multi(wide_data_untargeted_panuc, c('PD'), empri = 8, include_led = TRUE)
# make copy and modify led to be binary reponse
imputed_c_untargeted5_binary_led <- imputed_c_untargeted5_panuc
for(i in 1:5){
  imputed_c_untargeted5_binary_led[[i]][[1]][,"led"] <- ifelse(imputed_c_untargeted5_panuc[[i]][[1]][,"led"] == 0, 0, 1) 
}

#### --- logistic

led_logistic_fits <- list()
for(i in 1:5){
  set.seed(1)
  
  labels <- imputed_c_untargeted5_binary_led[[i]][[1]][,"led"]
  features <- imputed_c_untargeted5_binary_led[[i]][[1]][,-which(colnames(imputed_c_untargeted5_binary_led[[i]][[1]]) == "led")]
  #set penalty factors. (1 for each is default)
  p_factors_panuc <- rep(1, ncol(features))
  #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
  no_penalty_index <- which(colnames(features) %in% c('Age', 'GenderM'))
  p_factors[no_penalty_index] <- 0
  
  #grouped = false is already enforced for small folds. I'm just writing it explicitly to avoid the warning.
  binary_led_fit <- glmnet(features, labels, family = 'binomial', nlambda = nlambda, 
                           alpha = alpha, standardize = TRUE, penalty.factor = p_factors)
  
  pred_led <- predict(binary_led_fit, newx = features, type = "class")

  mse_led <- apply(pred_led, 2, function(x) (x - labels)^2 %>% mean)
  
  importance_led <- importance(binary_led_fit)
  
  pred_led[which(labels == 0), ]
  
  plot(binary_led_fit, xvar = "dev", label = T)
}







untargeted_led_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted5_panuc, target = "led", name = "untargeted", color = NULL, imp_num = .x, nlambda = 200))
untargeted_led_pred_df <- imputation_df(untargeted_led_analysis)
(untargeted_led_predtruth <- predtruth_plot(untargeted_led_pred_df, name = "Untargeted led")) 

untargeted_avg_retained_led <- untargeted_led_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_in_all_led <- untargeted_led_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

#plot colored by different things
predtruth_plot(untargeted_led_pred_df, color = "apoe", errorbar = FALSE, name = "Untargeted")
predtruth_plot(untargeted_led_pred_df, color = "apoe4", errorbar = FALSE, name = "Untargeted")


#-------------------------------------------
# targeted
wide_data_targeted_panuc <- panuc_data %>% 
  select(Id = subject_id, led, no_meds_reported) %>%
  inner_join(wide_data_targeted, by = "Id") %>%
  #drop na for now  (8 obs)
  filter(!is.na(led))


imputed_c_targeted5_panuc <- filter_and_impute_multi(wide_data_targeted_panuc, c('PD'), empri = 8, include_led = TRUE)


targeted_led_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_targeted5_panuc, target = "led", name = "targeted", color = NULL, imp_num = .x, nlambda = 200))
targeted_led_pred_df <- imputation_df(targeted_led_analysis)
(targeted_led_predtruth <- predtruth_plot(targeted_led_pred_df, name = "targeted led")) 

targeted_avg_retained_led <- targeted_led_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

targeted_in_all_led <- targeted_led_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

