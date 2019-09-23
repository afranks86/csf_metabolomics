
#########################################################


### Analysis Functions ####


############################################################

gender_control_analysis <- function(imputations, imp_num = 1){
  
  imputed_c <- imputations[[imp_num]]
  imputed_c_Y <- imputed_c[[1]]
  imputed_c_type <- imputed_c[[2]]
  imputed_c_apoe <- imputed_c[[3]]
  imputed_c_gender <- imputed_c_Y[,'GenderM']
  
  
  imputed_c_age <- imputed_c_Y[,'Age']
  # readd type as a feature in this analysis
  imputed_c_features_gender_tmp <- imputed_c_Y %>% 
    as_tibble %>%
    #mutate(Type = imputed_c_type) %>%
    select(-GenderM)
  
  #turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
  imputed_c_features_gender <- model.matrix(~., imputed_c_features_gender_tmp)
  
  
  fitpred_c_loo_gender <- lapply(1:nrow(imputed_c_features_gender), function(x) loo_cvfit_glmnet(x, imputed_c_features_gender, imputed_c_gender, 
                                                                                           alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE))
  
  fit_c_loo_gender <- lapply(fitpred_c_loo_gender, function(x) x[[1]])
  pred_c_loo_gender <- lapply(fitpred_c_loo_gender, function(x) x[[2]]) %>%
    unlist
  
  #some measure of variable importance
  importance_c_loo_gender <- lapply(fit_c_loo_gender, function(x) importance(x))
  
  importance_c_loo_median_gender <- importance_c_loo_gender %>% 
    purrr::map(~importance_consolidated_loo(.x)) %>%
    bind_rows() %>%
    #select only the variables that were present in >95% of fits
    select_if(function(x) sum(is.na(x))/length(x) < .05) %>%
    map_dbl(~median(.x, na.rm = T))


  roc_c_gender_loo <- fpr_tpr(pred_c_loo_gender, imputed_c_gender)
  roc_gender_plot <- ggplot(roc_c_gender_loo) + 
    geom_line(mapping = aes(fpr, tpr)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    theme_minimal() + 
    labs(title = "ROC: Gender vs C",
         subtitle = TeX('GOT + Lipids ,$\\alpha = 0.5$, loo'),
         x = 'False Positive Rate',
         y = 'True Positive Rate') + 
    geom_label(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0('AUC:', round(roc_c_gender_loo$auc[1], 3)))

  return(list(importance_c_loo_median_gender, roc_gender_plot)) 
}




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
  
  #turn factors into dummy variables
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
  
  
  
  return(list(c_loo_age_table, importance_c_loo_median_age, shapiro, pred_truth_c)) 
}


#' create df to combine/average imputations
#' @param analysis is a list created by multiple calls to age_control_analysis (eg. got_age_analysis)
#' 
imputation_df <- function(analysis){
  analysis %>% 
    purrr::map(~.x[[1]]$pred) %>%
    enframe %>%
    mutate(name = str_replace(name, "\\d", paste0("imp", name))) %>%
    spread(name, value) %>%
    unnest %>%
    # add truth (all of the truth columns are the same) 
    mutate(truth = analysis[[1]][[1]]$truth,
           gender = analysis[[1]][[1]]$gender %>%
             as.factor %>%
             fct_recode(M = '1', F = '0'),
           apoe = analysis[[1]][[1]]$apoe,
           apoe4 = analysis[[1]][[1]]$apoe4,
           type = analysis[[1]][[1]]$type
    ) %>%
    rowwise() %>%
    mutate(imp_avg = mean(c(imp1, imp2, imp3, imp4, imp5)),
           imp_min = min(c(imp1, imp2, imp3, imp4, imp5)),
           imp_max = max(c(imp1, imp2, imp3, imp4, imp5)))
  
}


#' create ggplot for pred vs truth in these grouped imputation tables
#' @param df is a dataframe with columns truth, imp_avg, imp_min, imp_max, apoe, gender, type
#' @param name is what we want to call the plots in the title (a string)
#' @param color is a string, representing the variable to color points by, one of apoe, gender, type
#' @param errorbar is boolean for whether to plot an errorbar. defaults to true
predtruth_plot <- function(df, name, color = NULL, errorbar = TRUE, data_name = "Control"){
  if(!is.null(color)){
    color <- sym(color)
  }
  if(errorbar){
    ggplot(df) + 
      geom_point(aes(truth, imp_avg, color = !!color), size = 2.5) + 
      geom_errorbar(aes(x = truth, ymin = imp_min, ymax = imp_max), alpha = 0.5) + 
      scale_color_brewer(type = 'qual', palette = 'Set1') +
      labs(title = paste0(data_name, ': True vs Predicted Age'),
           subtitle = paste0(name, " averaged over 5 imputations, alpha = 0.5, loo"),
           x = 'True Age',
           y = 'Predicted Age') + 
      geom_abline(intercept = 0, slope = 1) + 
      geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0("R: ", cor(truth, imp_avg, method = "pearson") %>% round(2), 
                                                                             "\nRMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2), 
                                                                             "\nMAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2))))
  } else {
    # same as above, but just without errorbar
    ggplot(df) + 
      geom_point(aes(truth, imp_avg, color = !!color), size = 2.5) + 
      scale_color_brewer(type = 'qual', palette = 'Set1') +
      labs(title = paste0(data_name, ': True vs Predicted Age'),
           subtitle = paste0(name, " averaged over 5 imputations, alpha = 0.5, loo"),
           x = 'True Age',
           y = 'Predicted Age') + 
      geom_abline(intercept = 0, slope = 1) + 
      geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0("R: ", cor(truth, imp_avg, method = "pearson") %>% round(2), 
                                                                             "\nRMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2), 
                                                                             "\nMAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2))))
  }
  
  
}





######################################################################



### Clock (linear glmnet age) ###




######################################################################









########

### GOT

########


imputed_c_got5 <- filter_and_impute_multi(wide_data, c('CO', 'CY', 'CM'))
got_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_got5, name = "GOT", color = NULL, imp_num = .x))

got_pred_df <- imputation_df(got_age_analysis)

# final plot
(got_age_predtruth <- predtruth_plot(got_pred_df, name = "GOT"))
ggsave(filename = "age_clock_GOT_pred.png")


got_in_all <- got_age_analysis %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")


########

### Lipids

########

imputed_c_lipids5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM'))
lipids_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_lipids5, name = "Lipids", color = NULL, imp_num = .x))

lipids_pred_df <- imputation_df(lipids_age_analysis)
  
(lipids_age_predtruth <- predtruth_plot(lipids_pred_df, name = "Lipids")) 
ggsave(filename = "age_clock_lipids_pred.png") #14.9 x 8.21


lipids_in_all <- lipids_age_analysis %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")




### Correlation between lipids and GOT residuals

got_avg_resid <- got_pred_df$imp_avg - got_pred_df$truth
lipids_avg_resid <- lipids_pred_df$imp_avg - lipids_pred_df$truth

cor.test(got_avg_resid, lipids_avg_resid, method = "spearman")




########

### Lipids (mice)

########

imputed_c_lipids_mice5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM'), method = "mice")
lipids_age_analysis_mice <- purrr::map(1:5, ~age_control_analysis(imputed_c_lipids_mice5, name = "Lipids", color = NULL, imp_num = .x))

lipids_pred_mice_df <- imputation_df(lipids_age_analysis_mice)

(lipids_age_mice_predtruth <- predtruth_plot(lipids_pred_mice_df, name = "Lipids (mice)")) 

lipids_in_all_mice <- lipids_age_analysis_mice %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")



########

### Lipids: Mice, loo

########







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



imputed_c_combined_amelia5 <- purrr::map2(imputed_c_got5, imputed_c_lipids5, ~combine_got_lipids(.x, .y))
combined_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_combined5, name = "Combined GOT + Lipids", color = NULL, imp_num = .x))

combined_pred_df <- imputation_df(combined_age_analysis)

(combined_age_predtruth <- predtruth_plot(combined_pred_df, name = "Combined GOT + Lipids"))
ggsave(filename = "age_clock_combined_pred.png")


#plot colored by different things
predtruth_plot(combined_pred_df, color = "gender", errorbar = FALSE, name = "Combined GOT + Lipids")
ggsave(filename = "age_clock_combined_pred_gender.png")
predtruth_plot(combined_pred_df, color = "apoe", errorbar = FALSE, name = "Combined GOT + Lipids")
ggsave(filename = "age_clock_combined_pred_apoe.png")
predtruth_plot(combined_pred_df, color = "apoe4", errorbar = FALSE, name = "Combined GOT + Lipids")
ggsave(filename = "age_clock_combined_pred_apoe4.png")




combined_in_all <- combined_age_analysis %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")





########

### Targeted

########

imputed_c_targeted5 <- filter_and_impute_multi(wide_data_targeted, c('CO', 'CY', 'CM'))
targeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_targeted5, name = "Targeted", color = NULL, imp_num = .x))

targeted_pred_df <- imputation_df(targeted_age_analysis)

(targeted_age_predtruth <- predtruth_plot(targeted_pred_df, name = "Targeted")) 
ggsave(filename = "age_clock_targeted_pred.png") #14.9 x 8.21












##################################################################



### Other tests ###



###################################################################





#########################

### Predict on AD/PD ###

########################


# First, we need to create a model fit on the entire dataset
  # no playing around with that loo monkey business anymore
  # We use combined got/lipids

#' Function to prepare the data dn create fits
#' @param imputations is the one element of the output of filter_and_impute_multi used to train
#' @param new_data is one element of the output of filter_and_impute_multi, used to predict
full_model_new_data <- function(imputation, new_data){
  
  true_control_age <- imputation[[1]][, "Age"]
  features <- imputation[[1]][,!colnames(imputation[[1]]) %in% 'Age']
  true_pred_age <- new_data[[1]][,"Age"]
  
  # sometimes the imputation leads to different numbers of columns (since it drops totally na)
  # the fit won't work if they don't have the same sizes, so subset the features to make sure it's same size
  shared_cols <- intersect(colnames(features), colnames(new_data[[1]]))
  features_shared <- features[, colnames(features) %in% shared_cols]
  new_data_shared <- new_data[[1]][, colnames(new_data[[1]]) %in% shared_cols]
  
  #combine the datasets to avoid errors because the factors are different
  combined_data <- rbind(features_shared, new_data_shared)
  
  fit <- fit_glmnet(combined_data[1:nrow(features_shared),], true_age, alpha = 0.5, penalize_age_gender = FALSE, family = "gaussian")
  oos_pred <- predict(fit, newx = combined_data[-(1:nrow(features_shared)),], s = "lambda.min")
  
  
  data_df <- tibble(
    truth = new_data[[1]][,"Age"],
    pred = oos_pred,
    gender = new_data[[1]][,"GenderM"] %>%
      as.factor %>%
      fct_recode(M = '1', F = '0'),
    type = new_data[[2]],
    apoe = new_data[[3]],
    apoe4 = apoe %>% fct_collapse('1' = c('24','34','44'), '0' = c('22', '23', '33'))
  )
  
  
  
  
  list("data_df" = data_df, "fit" = fit)
}


imputed_adpd_got_amelia5 <- filter_and_impute_multi(wide_data, c("AD", "PD"), method = "amelia")
imputed_adpd_lipids_amelia5 <- filter_and_impute_multi(wide_data_lipids, c("AD", "PD"), method = "amelia")

imputed_adpd_combined_amelia5 <- purrr::map2(imputed_adpd_got_amelia5, imputed_adpd_lipids_amelia5, ~combine_got_lipids(.x, .y))


adpd_age_analysis <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_combined_amelia5[[.x]], new_data = imputed_adpd_combined_amelia5[[.x]]))

adpd_pred_df <- imputation_df(adpd_age_analysis)

(adpd_age_predtruth <- predtruth_plot(adpd_pred_df, name = "Combined GOT + Lipids", data_name = "AD/PD")) 
ggsave(filename = "age_c_adpd_pred.png") #14.9 x 8.21



#########################

### univariate regression on Gender ~ Metabolite for each metabolite/lipid ###

########################


c_combined_gender_df <- imputed_c_combined5[[1]][[1]] %>%
  as_tibble() %>%
  rename_all(function(x) str_replace_all(x, "`", "")) %>%
  select(-c('(Intercept)', Age))


c_combined_names <- names(c_combined_gender_df) %>% setdiff('GenderM')
c_combined_gender_p_values <- purrr::map(c_combined_names, ~gender_metabolite_p(c_combined_gender_df, .x)) %>%
  unlist

#bh corrected controls for false discovery rate.
c_combined_gender_p_table <- bind_cols('name' = c_combined_names, 
                            'og_p_value' = c_combined_gender_p_values,
                            'bh_q_value' = p.adjust(c_combined_gender_p_values, method = 'BH')) 

c_combined_gender_p_table$name <- if_else(!str_detect(c_combined_gender_p_table$name, 'Result'), #take advantage of fact that all metabolites have "results" in name
                               c_combined_gender_p_table$name, 
                               str_replace_all(c_combined_gender_p_table$name, 'Result.*', "") %>%
                                 str_replace_all('\\.$', '') %>%
                                 str_replace_all('^\\.+', '') %>%
                                 str_trim() %>%
                                 sapply(function(x) all_matches[match(x, all_matches$Name), 'Metabolite'] %>% deframe))



c_combined_gender_p_table %>% 
  arrange(bh_q_value) 



#########################

### univariate regression on Age ~ Metabolite for each metabolite ###

########################

#' Create table with metabolite/lipid, p value, t value, name, bh corrected p value
#' (uses only the first imputation)
#' @return table with all of the variables with bh p values < 0.01
#' @param data is a imputation list (eg imputed_c_combined5)
#' 
bh_univariate_age <- function(data) {
  df <- data[[1]][[1]] %>%
    as_tibble() %>%
    dplyr::select(-c('(Intercept)', GenderM))
  
  p_values <- df %>%
    names %>%
    setdiff('Age') %>%
    purrr::map(~age_metabolite_p(df, .x)) %>%
    purrr::reduce(rbind) %>%
    as_tibble() %>%
    dplyr::rename('og_p_value' =  `Pr(>|t|)`)
  
  
  p_table <- cbind(p_values,'bh_p_value' = p.adjust(p_values$og_p_value, method = 'BH'))
  
  p_table %>% 
    dplyr::mutate(name = str_replace_all(name, '`', '')) %>% 
    filter(bh_p_value < 0.01)
  
}






## GOT

c_got_univariate_table <- bh_univariate_age(imputed_c_got5)



## Lipids

c_lipids_univariate_table <- bh_univariate_age(imputed_c_lipids5)




## Combined GOT + Lipids

c_combined_univariate_table <- bh_univariate_age(imputed_c_combined5)



## Targeted

c_targeted_univariate_table <- bh_univariate_age(imputed_c_targeted5)

















#############################################

### Gender Logistic Regresssion on combined

############################################

c_gender_logistic <- purrr::map(1:5, ~gender_control_analysis(imputed_c_combined5, imp_num = .x))

c_gender_logistic %>% 
  purrr::map(~.x[[2]]) %>% 
  cowplot::plot_grid(plotlist = .)


c_gender_logistic[[1]][[2]]
ggsave('gender_logistic_c.png')



############################

### Gender vs C ###
## {Lipids, GOT} ##

############################



