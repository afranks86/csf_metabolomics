
#########################################################


### Analysis Functions ####


############################################################

logistic_control_analysis <- function(imputations, varname = "GenderM", imp_num = 1, nlambda = 100){
  
  imputed_c <- imputations[[imp_num]]
  imputed_c_Y <- imputed_c[[1]]
  imputed_c_type <- imputed_c[[2]]
  imputed_c_apoe <- imputed_c[[3]]
  imputed_c_gender <- imputed_c_Y[,varname]
  
  
  
  imputed_c_age <- imputed_c_Y[,'Age']
  imputed_c_features_gender_tmp <- imputed_c_Y %>% 
    as_tibble %>%
    #mutate(Type = imputed_c_type) %>%
    select(-!!sym(varname))
  
  #turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
  imputed_c_features_gender <- model.matrix(~., imputed_c_features_gender_tmp)
  
  
  fitpred_c_loo_gender <- lapply(1:nrow(imputed_c_features_gender), function(x) loo_cvfit_glmnet(x, imputed_c_features_gender, imputed_c_gender, 
                                                                                           alpha = 0.5, family = 'binomial', penalize_age_gender = FALSE, nlambda = nlambda))
  
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
    labs(title = paste0("ROC: ", varname,  "vs C"),
         subtitle = TeX('Untargeted ,$\\alpha = 0.5$, loo'),
         x = 'False Positive Rate',
         y = 'True Positive Rate') + 
    geom_label(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0('AUC:', round(roc_c_gender_loo$auc[1], 3)))

  return(list(importance_c_loo_median_gender, roc_gender_plot, fit_c_loo_gender)) 
}




#' does full glmnet age analysis
#' @param imputations is the output of filter_and_impute_multi()
#' @param name is a string describing the dataset (eg GOT, Lipids, Combined)
#' @param color is a string, the variable what we want the plots to be colored by. options are gender, type, apoe, apoe4 
#' @param imp_num imputation number. 1-5
#' 
#' @return median importance, shapiro test, results table, predtruth plot
age_control_analysis <- function(imputations, name, color = NULL, imp_num = 1, nlambda = 100, ad_indicator = FALSE, pd_indicator = FALSE){
  
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
  
  if(ad_indicator == TRUE){
    imputed_c_features_age_tmp <- imputed_c_features_age_tmp %>%
      mutate(type = imputed_c_type,
             ad = ifelse(imputed_c_type == "AD", 1, 0)) %>%
      select(-type)
  }
  if(pd_indicator == TRUE){
    imputed_c_features_age_tmp <- imputed_c_features_age_tmp %>%
      mutate(type = imputed_c_type,
             pd = ifelse(imputed_c_type == "PD", 1, 0)) %>%
      select(-type)
  }
  
  
  #turn factors into dummy variables
  imputed_c_features_age <- model.matrix(~., imputed_c_features_age_tmp)
  
  
  fitpred_c_loo_age <- lapply(1:nrow(imputed_c_features_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_age, imputed_c_age, 
                                                                                           alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE, nlambda = nlambda))
  
  fit_c_loo_age <- lapply(fitpred_c_loo_age, function(x) x[[1]])
  pred_c_loo_age <- lapply(fitpred_c_loo_age, function(x) x[[2]]) %>%
    unlist
  
  #fit used to get lambda
  cv_fits <- lapply(fitpred_c_loo_age, function(x) x[[3]])
  
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
  
  
  
  return(list(c_loo_age_table, "importance" = importance_c_loo_median_age, shapiro, pred_truth_c, "loo_fits" = fit_c_loo_age, "full_fits" = cv_fits)) 
}


#' create df to combine/average imputations
#' @param analysis is a list created by multiple calls to age_control_analysis (eg. got_age_analysis)
#' @param num_imps is the number of imputations made. either 3 or 5
#' TODO: FIX THIS FUNCTION! code can be much cleaner to allow all num_imps
imputation_df <- function(analysis, num_imps = 5){
  
  if(num_imps == 5){
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
  } else if (num_imps == 3){
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
      mutate(imp_avg = mean(c(imp1, imp2, imp3)),
             imp_min = min(c(imp1, imp2, imp3)),
             imp_max = max(c(imp1, imp2, imp3)))
  } else {
    error("num_imps must be 3 or 5")
  }
  
  
}


#' create ggplot for pred vs truth in these grouped imputation tables
#' @param df is a dataframe with columns truth, imp_avg, imp_min, imp_max, apoe, gender, type
#' @param name is what we want to call the plots in the title (a string)
#' @param color is a string, representing the variable to color points by, one of apoe, gender, type
#' @param errorbar is boolean for whether to plot an errorbar. defaults to true
predtruth_plot <- function(df, pred_name = "imp_avg", name, color = NULL, errorbar = TRUE, data_name = "Control"){
  pred_name <- sym(pred_name)
  if(!is.null(color)){
    color <- sym(color)
  }
  
  if(errorbar){
    ggplot(df) + 
      geom_point(aes(truth, !!pred_name, color = !!color), size = 2.5) + 
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
      geom_point(aes(truth, !!pred_name, color = !!color), size = 2.5) + 
      scale_color_brewer(type = 'qual', palette = 'Set1') +
      labs(title = paste0(data_name, ': True vs Predicted Age'),
           subtitle = paste0(name, " averaged over 5 imputations, alpha = 0.5, loo"),
           x = 'True Age',
           y = 'Predicted Age') + 
      geom_abline(intercept = 0, slope = 1) + 
      geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 0, label = paste0("R: ", cor(truth, !!pred_name, method = "pearson") %>% round(2), 
                                                                             "\nRMSE: ", (truth - !!pred_name)^2 %>% mean %>% sqrt %>% round(2), 
                                                                             "\nMAE: ", (truth - !!pred_name) %>% abs %>% mean %>% round(2))))
  }
  
  
}





######################################################################



### Clock (linear glmnet age) ###
## Note: For the Amelia imputations, we keep empri = ~ 10% of the number of observations
##        This is a loose upper bound sugggested by the amelia paper




######################################################################








########

### GOT

########


imputed_c_got5 <- filter_and_impute_multi(wide_data, c('CO', 'CY', 'CM'), empri = 8)
got_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_got5, name = "GOT", color = NULL, imp_num = .x))

got_pred_df <- imputation_df(got_age_analysis)

# final plot
(got_age_predtruth <- predtruth_plot(got_pred_df, name = "GOT"))
ggsave(filename = "age_clock_GOT_pred.png")


got_in_all <- got_age_analysis %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")




## Create null model as baseline (it's dataset agnositic, since the observations are same)
null_df <- got_pred_df %>%
  ungroup() %>%
  select(truth, gender, apoe, apoe4, type) %>%
  mutate(pred = mean(truth))

(null_age_predtruth <- predtruth_plot(null_df, pred_name = "pred", name = "Mean model", errorbar = FALSE))



########

### Lipids

########

imputed_c_lipids5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM'), empri = 8)
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


# now, this correlation might include some correlation between the methods (as a result of the regularization to the horizontal line)
# to try and capture correlation better, we'll fit linear models for pred ~ truth, and then compare the residuals of those models
  # note: i know it's kinda backwards to do pred ~ truth, but we're just matching the x/y of the plot
got_truthpred_lm <- lm(imp_avg ~ truth, data = got_pred_df)
lipids_truthpred_lm <- lm(imp_avg ~ truth, data = lipids_pred_df)
cor.test(resid(got_truthpred_lm), resid(lipids_truthpred_lm), method = "spearman")


got_pred_df %>%
  ungroup %>%
  mutate(model_line = predict(got_truthpred_lm)) %>%
  predtruth_plot(name = "GOT") + 
  geom_line(aes(truth, model_line), color = "blue")

lipids_pred_df %>%
  ungroup %>%
  mutate(model_line = predict(lipids_truthpred_lm)) %>%
  predtruth_plot(name = "Lipids") + 
  geom_line(aes(truth, model_line), color = "blue")



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


########

### Unargeted

########

imputed_c_untargeted5 <- filter_and_impute_multi(wide_data_untargeted, c('CO', 'CY', 'CM'), empri = 8)
# lambdas were reaching the end of their sequence due to hugeness of data, so we increase nlambda from 100 to 200
untargeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted5, name = "untargeted", color = NULL, imp_num = .x, nlambda = 200))

untargeted_pred_df <- imputation_df(untargeted_age_analysis)

(untargeted_age_predtruth <- predtruth_plot(untargeted_pred_df, name = "Untargeted")) 
ggsave(filename = "age_clock_untargeted_pred.png") #14.9 x 8.21











##################################################################



### Other tests ###



###################################################################





#########################

### Predict on AD/PD ###

########################



################

## Using combined GOT + lipids

################


#' Helper function to bind two datasets together, keeping only shared columns
#' @param  include_metadata A boolean to determine whether to return a list (in same format at imputation1 or imputation2), or just the dataframe
#' @param include_age A boolean to determine whether to include age. This is useful if we are trying to make a feature set/full merge
#' @param add_AD/PD_ind A boollean to determine whether to add indicator columns for AD/PD. only works if include_metadata = T
merge_datasets <- function(imputation1, imputation2, include_metadata = F, include_age = F, add_AD_ind = F, add_PD_ind = FALSE){

  if(include_age == FALSE){
    features1 <- imputation1[[1]][,!colnames(imputation1[[1]]) %in% 'Age']
    features2 <- imputation2[[1]][,!colnames(imputation2[[1]]) %in% 'Age']
  } else if (include_age == TRUE){
    features1 <- imputation1[[1]]
    features2 <- imputation2[[1]]
  }  
  
  
  
  # sometimes the imputation leads to different numbers of columns (since it drops totally na)
  # the fit won't work if they don't have the same sizes, so subset the features to make sure it's same size
  shared_cols <- intersect(colnames(features1), colnames(features2))
  features1_shared <- features1[, colnames(features1) %in% shared_cols]
  features2_shared <- features2[, colnames(features2) %in% shared_cols]
  
  #combine the datasets to avoid errors because the factors are different
  combined_data <- rbind(features1_shared, features2_shared)
  
  if(include_metadata == TRUE){
    # 2:length(imputation1) is number of metadata elements in the list
    # we use unlist(list()) instead of c() to preserve factors
    metadata <- purrr::map(2:length(imputation1), ~unlist(list(imputation1[[.x]], imputation2[[.x]])))
    
    if(add_AD_ind == TRUE){
      ad_indicator <- ifelse(metadata[[1]] == "AD", 1, 0)
      combined_data <- cbind(combined_data, "AD_ind" = ad_indicator)
    }
    if(add_PD_ind == TRUE){
      pd_indicator <- ifelse(metadata[[1]] == "PD", 1, 0)
      combined_data <- cbind(combined_data, "PD_ind" = pd_indicator)
    }
    
    # add return to stop the function from going past the loop
    # NOTE: better practice would be let this function deal with more than 3 metadata args
    return(list(combined_data, metadata[[1]], metadata[[2]], metadata[[3]]))
  }
  
  
  
  combined_data
}

# First, we need to create a model fit on the entire dataset
  # no playing around with that loo monkey business anymore
  # We use combined got/lipids

#' Function to prepare the data dn create fits
#' @param imputations is the one element of the output of filter_and_impute_multi used to train
#' @param new_data is one element of the output of filter_and_impute_multi, used to predict
full_model_new_data <- function(imputation, new_data, nlambda = 100){
  
  true_control_age <- imputation[[1]][, "Age"]
  true_pred_age <- new_data[[1]][,"Age"]
  
  #to avoid any structure shenanigans, merge the datasets so that the columns match
    # (remember that we drop columns that have a certain amount of missingness)
  combined_data <- merge_datasets(imputation, new_data)

  # Make sure to fit only on the training set, and predict only on the test set by indexing combined_data
  fit <- fit_glmnet(combined_data[1:nrow(imputation[[1]]),], true_control_age, alpha = 0.5, penalize_age_gender = FALSE, family = "gaussian", nlambda = nlambda)
  oos_pred <- predict(fit, newx = combined_data[-(1:nrow(imputation[[1]])),], s = "lambda.min")
  
  
  data_df <- tibble(
    truth = true_pred_age,
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


adpd_age_analysis <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_combined_amelia5[[.x]], new_data = imputed_adpd_combined_amelia5[[.x]], nlambda = 200))

adpd_pred_df <- imputation_df(adpd_age_analysis)

(adpd_age_predtruth <- predtruth_plot(adpd_pred_df, name = "Combined GOT + Lipids", data_name = "AD/PD")) 
ggsave(filename = "age_c_adpd_pred.png") #14.9 x 8.21




################

## Using untargeted

################

#empri = .1(num PD + num_AD)
untargeted_adpd_combined_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('AD', 'PD'), empri = 10)


untargeted_adpd_age_analysis <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_untargeted5[[.x]], new_data = untargeted_adpd_combined_amelia5[[.x]], nlambda = 200))

untargeted_adpd_pred_df <- imputation_df(untargeted_adpd_age_analysis)

(untargeted_adpd_age_predtruth <- predtruth_plot(untargeted_adpd_pred_df, name = "Untargeted", data_name = "AD/PD")) 
ggsave(filename = "age_c_adpd_pred_untargeted.png") #14.9 x 8.21




################

## Using untargeted, separate imputation for AD/PD

################

#empri = .1(num PD/AD)
untargeted_ad_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('AD'), empri = 5)
untargeted_pd_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('PD'), empri = 5)

untargeted_adpd_separate_amelia5 <- purrr::map2(untargeted_ad_amelia5, untargeted_pd_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE))

untargeted_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_untargeted5[[.x]], new_data = untargeted_adpd_separate_amelia5[[.x]], nlambda = 200))

untargeted_adpd_pred_df_separate <- imputation_df(untargeted_adpd_age_analysis_separate)

(untargeted_adpd_age_predtruth_separate <- predtruth_plot(untargeted_adpd_pred_df_separate, name = "Untargeted", data_name = "AD/PD, separate", color = "type")) 
ggsave(filename = "age_c_adpd_pred_untargeted_separate.png") #14.9 x 8.21












########

### Untargeted, including ADPD
### We're interested in the coefficient for AD/PD

########

# combine controls with AD/PD imputation. include indicators for AD/PD as predictors
untargeted_all_amelia5 <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, add_AD_ind = TRUE, add_PD_ind = TRUE))


untargeted_all_age_analysis <- purrr::map(1:5, ~age_control_analysis(untargeted_all_amelia5, name = "untargeted (all)", color = NULL, imp_num = .x, nlambda = 200))

untargeted_all_pred_df <- imputation_df(untargeted_all_age_analysis)

(untargeted_all_age_predtruth <- predtruth_plot(untargeted_all_pred_df, name = "Untargeted (all)")) 
ggsave(filename = "age_clock_all_untargeted_pred.png") #14.9 x 8.21




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
    dplyr::mutate(name = str_replace_all(name, '`', ''))
  
}






## GOT

c_got_univariate_table <- bh_univariate_age(imputed_c_got5) %>%
  filter(bh_p_value < 0.01)



## Lipids

c_lipids_univariate_table <- bh_univariate_age(imputed_c_lipids5) %>%
  filter(bh_p_value < 0.01)




## Combined GOT + Lipids

c_combined_univariate_table <- bh_univariate_age(imputed_c_combined5) %>%
  filter(bh_p_value < 0.01)



## Targeted

c_targeted_univariate_table <- bh_univariate_age(imputed_c_targeted5)

c_targeted_univar_sig <- c_targeted_univariate_table %>%
  filter(bh_p_value < 0.01)


## Untargeted

c_untargeted_univariate_table <- bh_univariate_age(imputed_c_untargeted5)

c_untargeted_univar_table_sig <- c_untargeted_univariate_table %>%
  filter(bh_p_value < 0.01)















#############################################

### Gender Logistic Regresssion on untargeted

############################################

c_gender_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_c_untargeted5, varname = "GenderM", imp_num = .x, nlambda = 200))

c_gender_logistic %>% 
  purrr::map(~.x[[2]]) %>% 
  cowplot::plot_grid(plotlist = .)


c_gender_logistic[[1]][[2]]
ggsave('gender_logistic_c.png')




############################

### AD/PD Logisitic Regression on untargeted ###

############################

#### AD First -----------------------------
# we want an AD indicator, but not a PD one
untargeted_all_amelia5_ad_ind <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
                                                                                                                     add_AD_ind = TRUE, add_PD_ind = FALSE))
untargted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(untargeted_all_amelia5_ad_ind, varname ="AD_ind", imp_num = .x, nlambda = 200))

untargeted_ad_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

untargeted_ad_logistic[[1]][[1]]
ggsave("ad_logistic_untargeted.png")


#### PD -------------------------------------
# we want a PD indicator, but not an AD one
untargeted_all_amelia5_pd_ind <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
                                                                                                                      add_AD_ind = FALSE, add_PD_ind = TRUE))
untargted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(untargeted_all_amelia5_pd_ind, varname ="PD_ind", imp_num = .x, nlambda = 200))

untargeted_pd_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

untargeted_pd_logistic[[1]][[1]]
ggsave("pd_logistic_untargeted.png")










############################

### Trying MSEA from MetaboAnalystR on Targeted significant age

############################

## As input, it takes in the names/ids of significant variables
univar_targeted_names <- c_targeted_univariate_table %>%
  filter(bh_p_value < 0.01) %>% select(name) %>% deframe


## We also want to add in a reference list (ie include names on insignificant variables)
## The reference list must be IDs, so we need to map the names to KEGG/HMDB names.
# Post mortem, I'm coming back and renaming ones that I can identify (but that atempt 2 couldn't)
univar_targeted_names_all <- c_targeted_univariate_table %>%
  mutate(name = str_replace_all(name, "_neg|_pos", "") %>% str_to_lower() %>% str_trim(),
         mapped_name = case_when(
          #https://pubchem.ncbi.nlm.nih.gov/compound/12-Ketolithocholic-acid
          name == "3\\?-hydroxy-12 ketolithocholic acid" ~ "12-Ketodeoxycholic acid",
          #https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:27726 / https://www.genome.jp/dbget-bin/www_bget?cpd:C06560
          name == "2-chloro-4,6-diamino-1,3,5-triazine" ~ "Deisopropyldeethylatrazine",
          #http://www.hmdb.ca/metabolites/HMDB0006029
          name == "acetyl-l-glutamine" ~ "N-Acetylglutamine",
          #https://pubchem.ncbi.nlm.nih.gov/compound/Acetylornithine
          name == "acetylornithine" ~  "N-Acetylornithine",
          #https://pubchem.ncbi.nlm.nih.gov/compound/5-Aminovaleric-acid
          name == "amino valerate" ~ "5-Aminopentanoic acid",
          #https://pubchem.ncbi.nlm.nih.gov/compound/N_N-dimethylarginine
          name == "dimethylarginine" ~ "Asymmetric dimethylarginine",
          #https://pubchem.ncbi.nlm.nih.gov/compound/1826
          name == "hiaa" ~ 	"5-Hydroxyindoleacetic acid",
          TRUE ~ name)) %>%
  select(mapped_name)



# Attempt 1: Use our kegg_map.csv
  # .. this one didn't work so well. we only matched ~50% of the the obs
kegg_lookup <- kegg_map %>%
  mutate(name_formatted = str_to_lower(METABOLITE) %>% str_replace_all("\\(.+\\)", "") %>% str_trim)

univar_targeted_names_kegg <- univar_targeted_names_all %>% left_join(kegg_lookup, by = "name_formatted")


# Attempt 2: Use metaboanalyst to help with the mapping
hmdbMap <- InitDataObjects("conc", "pathora", F)
hmdbMap <- Setup.MapData(hmdbMap, deframe(univar_targeted_names_all))
hmdbMap <- CrossReferencing(hmdbMap, "name")
hmdbMap <- CreateMappingResultTable(hmdbMap)

#trying to find any of them manually. Thes ones I could correct are corrected above.
hmdbMap<-PerformDetailMatch(hmdbMap, "Deisopropyldeethylatrazine")
hmdbMap <- GetCandidateList(hmdbMap)

# --
hmdb_mapping <- hmdbMap$dataSet$map.table[,c("Query","Match")] %>%
  as_tibble()

write_delim(select(hmdb_mapping, Match), path = "targeted_names_hmdb_list.txt", delim = "\n", col_names = FALSE)

  




# Using targeted since it's the only analysis with defined mapping
# this one is from the univariate analysis
targeted_sig_names <- univar_targeted_names %>%
  str_replace_all("_neg", "") %>%
  str_replace_all("_pos", "") %>%
  # an alternate name for 12 ketolithocholic acid according to https://pubchem.ncbi.nlm.nih.gov/compound/12-Ketolithocholic-acid
  str_replace("3\\?-Hydroxy-12 Ketolithocholic Acid", "12-Ketodeoxycholic acid")


## Following the package vignette

mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Setup.MapData(mSet, targeted_sig_names)

# upload our reference list (the metabolites we targeted)
mSet<-Setup.HMDBReferenceMetabolome(mSet, "targeted_names_hmdb_list.txt");
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)



## try to match the NAs manually... 

mSet<-PerformDetailMatch(mSet, "HIAA")
mSet <- GetCandidateList(mSet)
# Match found! set to the first candidate
mSet<-SetCandidate(mSet, "HIAA", mSet$name.map$hits.candidate.list[1])
##


#what does this do?? True is the only option in the web-tool version, so i'm keeping it as default
# but it makes a pretty big difference
mSet<-SetMetabolomeFilter(mSet, F)

# Select metabolite set library

#pathway associated library?
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway",2)

# csf associated library?
#mSet<-SetCurrentMsetLib(mSet, "csf", 2)

# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "bar", "png", 72, width=NA)


