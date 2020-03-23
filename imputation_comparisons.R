lipids_age <- function(imp, method = "amelia"){
  num_controls_lipids <- wide_data_lipids %>%
    filter(Type %in% c("CY", "CM", "CO")) %>% 
    nrow
  
  
  fitpred_lipids_separate_loo_age <- lapply(1:num_controls_lipids, function(x) impute_c_loo_cvfit_glmnet(x, wide_data_lipids, alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE, imp_num = imp, method = method))
  
  fit_lipids_separate_loo_age <- lapply(fitpred_lipids_separate_loo_age, function(x) x[[1]])
  pred_lipids_separate_loo_age <- lapply(fitpred_lipids_separate_loo_age, function(x) x[[2]]) %>%
    unlist
  truth_lipids_separate_loo_age <- lapply(fitpred_lipids_separate_loo_age, function(x) x[[3]]) %>% unlist
  type_lipids_separate_loo_age <- lapply(fitpred_lipids_separate_loo_age, function(x) x[[4]]) %>% unlist
  apoe_lipids_separate_loo_age <- lapply(fitpred_lipids_separate_loo_age, function(x) x[[5]]) %>% unlist
  gender_lipids_separate_loo_age <- lapply(fitpred_lipids_separate_loo_age, function(x) x[[6]]) %>% unlist
  
  #some measure of variable importance
  importance_lipids_separate_loo_age <- lapply(fit_lipids_separate_loo_age, function(x) importance(x))
  
  importance_lipids_separate_loo_median_age <- importance_lipids_separate_loo_age %>% 
    purrr::map(~importance_consolidated_loo(.x)) %>%
    bind_rows() %>%
    #select only the variables that were present in >95% of fits
    select_if(function(x) sum(is.na(x))/length(x) < .05) %>%
    map_dbl(~median(.x, na.rm = T))
  
  
  
  mse_lipids_separate_loo_age <- mean((pred_lipids_separate_loo_age - truth_lipids_separate_loo_age)^2)
  resid_lipids_separate_loo_age <- pred_lipids_separate_loo_age - truth_lipids_separate_loo_age
  
  
  shapiro.test(resid_lipids_separate_loo_age)
  
  qplot(resid_lipids_separate_loo_age, bins = 10, xlab = 'Residuals', main = 'Histogram of Age Residuals, alpha = 0.5')
  #ggsave('lipids_lipids_age_control_resid_hist.png')
  qqnorm(resid_lipids_separate_loo_age)
  qqline(resid_lipids_separate_loo_age)
  
  
  ### look at alpha = 0.5
  loo_lipids_separate_age_table <- tibble(truth = truth_lipids_separate_loo_age, 
                                          pred = pred_lipids_separate_loo_age,
                                          resid = truth - pred,
                                          apoe = apoe_lipids_separate_loo_age,
                                          type = type_lipids_separate_loo_age,
                                          gender = gender_lipids_separate_loo_age,
                                          apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
  )
  
  
  pred_truth_c_lipids_separate <- ggplot(loo_lipids_separate_age_table) + 
      geom_point(aes(truth, pred, color = apoe)) + 
      scale_color_brewer(type = 'qual', palette = 'Set1') +
      labs(title = 'Control: True vs Predicted Age',
           subtitle = 'lipids (imputed separately), alpha = 0.5, loo',
           x = 'True Age',
           y = 'Predicted Age') + 
      geom_abline(intercept = 0, slope = 1) + 
      geom_text(aes(x = Inf, y = -Inf, vjust = -1, hjust = 1.25, label = paste0("R: ", cor(truth, pred, method = "pearson") %>% round(2), 
                                                                                "\nRMSE: ", (truth - pred)^2 %>% mean %>% sqrt %>% round(2))))
  
  
  
  return(list(loo_lipids_separate_age_table, pred_truth_c_lipids_separate))
}


lipids_list <- list(lipids_age(1), lipids_age(2), lipids_age(3))

lipids_df <- lipids_list %>%
  purrr::map(~.x[[1]]$pred) %>%
  enframe %>%
  mutate(name = str_replace(name, "\\d", paste0("imp", name))) %>%
  spread(name, value) %>%
  unnest %>%
  # add truth (all of the truth columns are the same) 
  mutate(truth = lipids_list[[1]][[1]]$truth,
         gender = lipids_list[[1]][[1]]$gender) %>%
  rowwise() %>%
  mutate(imp_avg = mean(c(imp1, imp2, imp3)))

ggplot(lipids_df) + 
  geom_point(aes(truth, imp_avg, color = gender)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = 'lipids (imputed separately), alpha = 0.5, loo',
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_text(aes(x = Inf, y = -Inf, vjust = -1, hjust = 1.25, label = paste0("R: ", cor(truth, imp_avg, method = "pearson") %>% round(2), 
                                                                            "\nRMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2))))


  
  



lipids_list_mice <- list(lipids_age(1, method = "mice"), lipids_age(2, method = "mice"), lipids_age(3, method = "mice"))

