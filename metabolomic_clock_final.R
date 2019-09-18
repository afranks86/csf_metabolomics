
#########################################################


### Age analysis Function ####


############################################################


#' does full glmnet age analysis
#' @param imputations is the output of filter_and_impute_multi()
#' @param name is a string describing the dataset (eg GOT, Lipids, Combined)
#' @param color is a string, the variable what we want the plots to be colored by. options are gender, type, apoe, apoe4 
#' @param imp_num imputation number. 1-5
#' 
#' @return median importance, shapiro test, results table, predtruth plot
age_control_analysis <- function(imputations, name, color = NULL, imp_num = 1){
  
  if(!is.null(color)){
    color <- sym(color)
  }
  
  imputed_c <- imputations[[imp_num]]
  imputed_c_Y <- imputed_c[[1]]
  imputed_c_type <- imputed_c[[2]]
  imputed_c_apoe <- imputed_c[[3]]
  imputed_c_gender <- imputed_c_Y[,'GenderM']
  
  
  imputed_c_age <- imputed_c_Y[,'Age']
  # readd type as a feature in this analysis
  imputed_c_features_age_tmp <- imputed_c_Y %>% 
    as_tibble %>%
    #mutate(Type = imputed_c_type) %>%
    select(-Age)
  
  #turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
  imputed_c_features_age <- model.matrix(~., imputed_c_features_age_tmp)
  
  
  fitpred_c_loo_age <- lapply(1:nrow(imputed_c_features_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_age, imputed_c_age, 
                                                                                           alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE))
  
  fit_c_loo_age <- lapply(fitpred_c_loo_age, function(x) x[[1]])
  pred_c_loo_age <- lapply(fitpred_c_loo_age, function(x) x[[2]]) %>%
    unlist
  
  #some measure of variable importance
  importance_c_loo_age <- lapply(fit_c_loo_age, function(x) importance(x))
  
  importance_c_loo_median_age <- importance_c_loo_age %>% 
    purrr::map(~importance_consolidated_loo(.x)) %>%
    bind_rows() %>%
    #select only the variables that were present in >95% of fits
    select_if(function(x) sum(is.na(x))/length(x) < .05) %>%
    map_dbl(~median(.x, na.rm = T))
  
  
  resid_c_loo_age <- pred_c_loo_age - imputed_c_age
  shapiro <- shapiro.test(resid_c_loo_age)
  
  
  ### look at alpha = 0.5
  c_loo_age_table <- tibble(truth = imputed_c_age, 
                            pred = pred_c_loo_age,
                            resid = truth - pred,
                            apoe = imputed_c_apoe,
                            type = imputed_c_type,
                            gender = imputed_c_gender,
                            apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
  )
  
  
  pred_truth_c <- ggplot(c_loo_age_table) + 
    geom_point(aes(truth, pred, color = !!color)) + 
    scale_color_brewer(type = 'qual', palette = 'Set1') +
    labs(title = 'Control: True vs Predicted Age',
         subtitle = paste0(name, ', alpha = 0.5, loo'),
         x = 'True Age',
         y = 'Predicted Age') + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0("R: ", cor(truth, pred, method = "pearson") %>% round(2), 
                                                                           "\nRMSE: ", (truth - pred)^2 %>% mean %>% sqrt %>% round(2), 
                                                                           "\nMAE: ", (truth - pred) %>% abs %>% mean %>% round(2))))
  
  
  
  return(list(importance_c_loo_median_age, shapiro, c_loo_age_table, pred_truth_c)) 
}




########

### GOT

########


imputed_c_got5 <- filter_and_impute_multi(wide_data, c('CO', 'CY', 'CM'))
got_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_got5, name = "GOT", color = NULL, imp_num = .x))

got_pred_df <- got_age_analysis %>% 
  purrr::map(~.x[[3]]$pred) %>%
  enframe %>%
  mutate(name = str_replace(name, "\\d", paste0("imp", name))) %>%
  spread(name, value) %>%
  unnest %>%
  # add truth (all of the truth columns are the same) 
  mutate(truth = got_age_analysis[[1]][[3]]$truth,
         gender = got_age_analysis[[1]][[3]]$gender,
         apoe = got_age_analysis[[1]][[3]]$apoe,
         apoe4 = got_age_analysis[[1]][[3]]$apoe4,
         type = got_age_analysis[[1]][[3]]$type
         ) %>%
  rowwise() %>%
  mutate(imp_avg = mean(c(imp1, imp2, imp3, imp4, imp5)),
         imp_min = min(c(imp1, imp2, imp3, imp4, imp5)),
         imp_max = max(c(imp1, imp2, imp3, imp4, imp5)))

# final plot
(got_age_predtruth <- ggplot(got_pred_df) + 
  geom_point(aes(truth, imp_avg, color = NULL)) + 
    geom_errorbar(aes(x = truth, ymin = imp_min, ymax = imp_max)) + 
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'Control: True vs Predicted Age',
       subtitle = "GOT averaged over 5 imputations, alpha = 0.5, loo",
       x = 'True Age',
       y = 'Predicted Age') + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0("R: ", cor(truth, imp_avg, method = "pearson") %>% round(2), 
                                                                         "\nRMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2), 
                                                                         "\nMAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2))))
)


########

### Lipids

########

imputed_c_lipids5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM'))
lipids_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_lipids5, name = "Lipids", color = NULL, imp_num = .x))

lipids_pred_df <- lipids_age_analysis %>% 
  purrr::map(~.x[[3]]$pred) %>%
  enframe %>%
  mutate(name = str_replace(name, "\\d", paste0("imp", name))) %>%
  spread(name, value) %>%
  unnest %>%
  # add truth (all of the truth columns are the same) 
  mutate(truth = lipids_age_analysis[[1]][[3]]$truth,
         gender = lipids_age_analysis[[1]][[3]]$gender,
         apoe = lipids_age_analysis[[1]][[3]]$apoe,
         apoe4 = lipids_age_analysis[[1]][[3]]$apoe4,
         type = lipids_age_analysis[[1]][[3]]$type
  ) %>%
  rowwise() %>%
  mutate(imp_avg = mean(c(imp1, imp2, imp3, imp4, imp5)),
         imp_min = min(c(imp1, imp2, imp3, imp4, imp5)),
         imp_max = max(c(imp1, imp2, imp3, imp4, imp5)))

# final plot
(lipids_age_predtruth <- ggplot(lipids_pred_df) + 
    geom_point(aes(truth, imp_avg, color = NULL)) + 
    geom_errorbar(aes(x = truth, ymin = imp_min, ymax = imp_max)) + 
    scale_color_brewer(type = 'qual', palette = 'Set1') +
    labs(title = 'Control: True vs Predicted Age',
         subtitle = "Lipids averaged over 5 imputations, alpha = 0.5, loo",
         x = 'True Age',
         y = 'Predicted Age') + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0("R: ", cor(truth, imp_avg, method = "pearson") %>% round(2), 
                                                                           "\nRMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2), 
                                                                           "\nMAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2))))
)






########

### Combined GOT + lipids

########


#' Helper functions to combine the imputed lipids/GOT data
#' @param got_imp is a single element of filter_and_impute_multi(data = wide_data)
#' @param lipid_imp is a single element of filter_and_impute_multi(data = wide_data_lipids)
combine_got_lipids <- function(got_imp, lipid_imp){
  # add the id columns to make the join easier
  imputed_c_got_separate_withid <- got_imp[[1]] %>% 
    as_tibble() %>%
    mutate(id = got_imp[[4]])
  
  # we only need to track apoe/gender/age once
  imputed_c_lipids_separate_withid <- lipid_imp[[1]] %>%
    as_tibble() %>%
    dplyr::select(-c(Age, GenderM, `(Intercept)`)) %>%
    dplyr::mutate(id = lipid_imp[[4]], 
                  apoe = lipid_imp[[3]],
                  type = lipid_imp[[2]])
  
  # join
  imputed_c_combined_separate_df <- imputed_c_got_separate_withid %>%
    inner_join(imputed_c_lipids_separate_withid, by = 'id') %>%
    dplyr::select(-c(id))
  
  # transform into matrix
  Y <- imputed_c_combined_separate_df %>%
    dplyr::select(-c(apoe, type)) %>%
    as.matrix()
  
  type <- imputed_c_combined_separate_df$type %>% as.factor()
  apoe <- imputed_c_combined_separate_df$apoe %>% as.factor()
  id <- imputed_c_got_separate_withid$id
  
  message("list order is Y, type, apoe, id")
  return(list(Y, type, apoe, id))
}


imputed_c_combined5 <- purrr::map2(imputed_c_got5, imputed_c_lipids5, ~combine_got_lipids(.x, .y))

combined_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_combined5, name = "Combined GOT + Lipids", color = NULL, imp_num = .x))

combined_pred_df <- combined_age_analysis %>% 
  purrr::map(~.x[[3]]$pred) %>%
  enframe %>%
  mutate(name = str_replace(name, "\\d", paste0("imp", name))) %>%
  spread(name, value) %>%
  unnest %>%
  # add truth (all of the truth columns are the same) 
  mutate(truth = combined_age_analysis[[1]][[3]]$truth,
         gender = combined_age_analysis[[1]][[3]]$gender,
         apoe = combined_age_analysis[[1]][[3]]$apoe,
         apoe4 = combined_age_analysis[[1]][[3]]$apoe4,
         type = combined_age_analysis[[1]][[3]]$type
  ) %>%
  rowwise() %>%
  mutate(imp_avg = mean(c(imp1, imp2, imp3, imp4, imp5)),
         imp_min = min(c(imp1, imp2, imp3, imp4, imp5)),
         imp_max = max(c(imp1, imp2, imp3, imp4, imp5)))

# final plot
(combined_age_predtruth <- ggplot(combined_pred_df) + 
    geom_point(aes(truth, imp_avg, color = NULL)) + 
    geom_errorbar(aes(x = truth, ymin = imp_min, ymax = imp_max)) + 
    scale_color_brewer(type = 'qual', palette = 'Set1') +
    labs(title = 'Control: True vs Predicted Age',
         subtitle = "Combined GOT + Lipids averaged over 5 imputations, alpha = 0.5, loo",
         x = 'True Age',
         y = 'Predicted Age') + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0("R: ", cor(truth, imp_avg, method = "pearson") %>% round(2), 
                                                                           "\nRMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2), 
                                                                           "\nMAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2))))
)


