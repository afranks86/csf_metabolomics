error_log <- file("analysis_log_102719.Rout", open="wt")
sink(error_log, type = "message")
source(here::here("analysis", "starter.R"))
#########################################################


### Analysis Functions ####


############################################################

#' function to do logistic glmnet with "varname" as predictorr
#' @parm varname is the target variable, if it is in the predictor matrix (so pretty much just gender... I should probably fix this LOL)
#' NOTE: if genderM is not varname, then it will be included as a predictor in the model
#' AD/PD/GBA flags take precedent over varname when determining the target var
#' ie if AD_ind is TRUE, then the only purpose of varname is to set plot names
logistic_control_analysis <- function(imputations, varname = "GenderM", imp_num = 1, nlambda = 100, AD_ind = FALSE, PD_ind = FALSE, GBA_ind = FALSE){
  
  imputed_c <- imputations[[imp_num]]
  imputed_c_Y <- imputed_c[[1]]
  imputed_c_type <- imputed_c[[2]]
  imputed_c_apoe <- imputed_c[[3]]
  imputed_c_gender <- imputed_c_Y[,"GenderM"]
  
  
  
  if(AD_ind & PD_ind){
    stop("Only one of AD_ind and PD_ind must be selected")
  }
  # AD/PD flags take precedent over varname
  if(AD_ind){
    imputed_c_target <- ifelse(imputed_c_type == "AD", 1,0)
  } else if(PD_ind){
    imputed_c_target <- ifelse(imputed_c_type == "PD", 1,0)
  } else if(GBA_ind){
    imputed_c_gba <- imputed_c[[5]]
    imputed_c_target <- ifelse(imputed_c_gba %in% c('E326K Carrier', 'Pathogenic Carrier', 'CT'), 1, 0)
  } else{
    imputed_c_target <- imputed_c_Y[,varname]
  }
  
  
  imputed_c_age <- imputed_c_Y[,'Age']
  
  # if the target variable (varname) is in the predictor matrix, we gotta get rid of it
  # we get rid of intercept here because it gets added back in model.matrix call below.
  if(varname %in% colnames(imputed_c_Y)){
    imputed_c_features_target_tmp <- imputed_c_Y %>% 
      as_tibble %>%
      #mutate(Type = imputed_c_type) %>%
      select(-c(!!sym(varname), "(Intercept)"))
  } else{
    imputed_c_features_target_tmp <- imputed_c_Y %>% 
      as_tibble %>%
      #mutate(Type = imputed_c_type) %>%
      select(-c("(Intercept)"))
  }
  
  
  #turn type into a dummy var (multiple columns. AD is the redundant column (chosen))
  imputed_c_features_target <- model.matrix(~., imputed_c_features_target_tmp)
  
  full_model <- get_full_model(features= imputed_c_features_target, imputed_c_target, alpha = 0.5, family = "binomial", penalize_AD_PD = FALSE, penalize_age_gender = FALSE, nlambda = nlambda)
  fitpred_c_loo_target <- lapply(1:nrow(imputed_c_features_target), function(x) loo_cvfit_glmnet(x, imputed_c_features_target, imputed_c_target, lambda = full_model[[2]], full_fit = full_model[[1]],
                                                                                           alpha = 0.5, family = 'binomial', penalize_age_gender = FALSE, nlambda = nlambda))
  
  fit_c_loo_target <- lapply(fitpred_c_loo_target, function(x) x[[1]])
  pred_c_loo_target <- lapply(fitpred_c_loo_target, function(x) x[[2]]) %>%
    unlist
  
  #some measure of variable importance
  importance_c_loo_target <- lapply(fit_c_loo_target, function(x) importance(x))
  
  importance_c_loo_median_target <- importance_c_loo_target %>% 
    purrr::map(~importance_consolidated_loo(.x)) %>%
    bind_rows() %>%
    #select only the variables that were present in >95% of fits
    select_if(function(x) sum(is.na(x))/length(x) < .05) %>%
    map_dbl(~median(.x, na.rm = T))


  roc_c_target_loo <- fpr_tpr(pred_c_loo_target, imputed_c_target)
  roc_target_plot <- ggplot(roc_c_target_loo) + 
    geom_line(mapping = aes(fpr, tpr)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    theme_minimal() + 
    labs(title = paste0("ROC: ", varname,  "vs C"),
         subtitle = TeX('Untargeted ,$\\alpha = 0.5$, loo'),
         x = 'False Positive Rate',
         y = 'True Positive Rate') + 
    geom_label(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0('AUC:', round(roc_c_target_loo$auc[1], 3)),
               size =6)

  return(list(nonzero = importance_c_loo_median_target, roc = roc_target_plot, fit = fit_c_loo_target, truth = imputed_c_target, pred =pred_c_loo_target)) 
}




#' does full glmnet age analysis
#' @param imputations is the output of filter_and_impute_multi()\
#' @param target is the string with the name of the target variable name. (default Age)
#' @param name is a string describing the dataset (eg GOT, Lipids, Combined)
#' @param color is a string, the variable what we want the plots to be colored by. options are gender, type, apoe, apoe4 
#' @param imp_num imputation number. 1-5
#' 
#' @return median importance, shapiro test, results table, predtruth plot
age_control_analysis <- function(imputations, target = "Age", name, color = NULL, imp_num = 1, nlambda = 100, ad_indicator = FALSE, pd_indicator = FALSE){
  
  if(!is.null(color)){
    color <- sym(color)
  }
  
  imputed_c <- imputations[[imp_num]]
  imputed_c_Y <- imputed_c[[1]]
  imputed_c_type <- imputed_c[[2]]
  imputed_c_apoe <- imputed_c[[3]]
  imputed_c_id <- imputed_c[[4]]
  imputed_c_gender <- imputed_c_Y[,'GenderM']
  
  
  #1/8imputed_c_age <- imputed_c_Y[,'Age']
  imputed_c_age <- imputed_c_Y[,target]
  # readd type as a feature in this analysis
  imputed_c_features_age_tmp <- imputed_c_Y %>% 
    as_tibble %>%
    #mutate(Type = imputed_c_type) %>%
    #1/8select(-Age)
    select(-!!sym(target))
  
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
  
  
  full_model <- get_full_model(features= imputed_c_features_age, imputed_c_age, alpha = 0.5, family = "gaussian", penalize_AD_PD = FALSE, penalize_age_gender = FALSE, nlambda = nlambda)
  # Note: If AD_ind, PD_ind are missing from the dataset, the flag penalize_AD_Pd doesn't do anything
  fitpred_c_loo_age <- lapply(1:nrow(imputed_c_features_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_age, imputed_c_age, lambda = full_model[[2]], full_fit = full_model[[1]],
                                                                                           alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE, penalize_AD_PD = FALSE, nlambda = nlambda))
  
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
                            id = imputed_c_id,
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
    geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, pred, method = "pearson")^2 %>% round(2), 
                                                                           "<br>RMSE: ", (truth - pred)^2 %>% mean %>% sqrt %>% round(2), 
                                                                           "<br>MAE: ", (truth - pred) %>% abs %>% mean %>% round(2))),
               size = 12
               )
  
  
  
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
             type = analysis[[1]][[1]]$type,
             id = analysis[[1]][[1]]$id
      ) %>%
      rowwise() %>%
      mutate(imp_avg = mean(c(imp1, imp2, imp3, imp4, imp5)),
             imp_min = min(c(imp1, imp2, imp3, imp4, imp5)),
             imp_max = max(c(imp1, imp2, imp3, imp4, imp5))) %>%
      ungroup()
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
             type = analysis[[1]][[1]]$type,
             id = analysis[[1]][[1]]$id
      ) %>%
      rowwise() %>%
      mutate(imp_avg = mean(c(imp1, imp2, imp3)),
             imp_min = min(c(imp1, imp2, imp3)),
             imp_max = max(c(imp1, imp2, imp3))) %>%
      ungroup()
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
  
  # Null model is to predict mean
  null_pred <- mean(df$truth) %>% rep(times = length(df$truth))
  rmse_null <- (df$truth - null_pred)^2 %>% mean %>% sqrt %>% round(2)
  mae_null <- (df$truth - null_pred) %>% abs %>% mean %>% round(2)
  
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
      geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, !!pred_name, method = "pearson")^2 %>% round(2), 
                                                                                "<br>RMSE: ", (truth - !!pred_name)^2 %>% mean %>% sqrt %>% round(2), " (",rmse_null,")", 
                                                                                "<br>MAE: ", (truth - !!pred_name) %>% abs %>% mean %>% round(2), " (",mae_null,")")),
                    size = 12)
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
      geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, !!pred_name, method = "pearson")^2 %>% round(2), 
                                                                             "<br>RMSE: ", (truth - !!pred_name)^2 %>% mean %>% sqrt %>% round(2), " (",rmse_null,")", 
                                                                             "<br>MAE: ", (truth - !!pred_name) %>% abs %>% mean %>% round(2), " (",mae_null,")")),
                    size = 12)
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
message("GOT -------------------------------------------")

imputed_c_got5 <- filter_and_impute_multi(wide_data, c('CO', 'CY', 'CM'), empri = 8)
got_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_got5, name = "GOT", color = NULL, imp_num = .x))

got_pred_df <- imputation_df(got_age_analysis)

# final plot
(got_age_predtruth <- predtruth_plot(got_pred_df, name = "GOT"))
ggsave(filename = "age_clock_GOT_pred.png")


got_avg_retained <- got_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

got_in_all <- got_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
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
message("Lipids (amelia) -------------------------------------------")

imputed_c_lipids5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM'), empri = 8)
lipids_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_lipids5, name = "Lipids", color = NULL, imp_num = .x))

lipids_pred_df <- imputation_df(lipids_age_analysis)
  
(lipids_age_predtruth <- predtruth_plot(lipids_pred_df, name = "Lipids")) 
ggsave(filename = "age_clock_lipids_pred.png") #14.9 x 8.21

lipids_avg_retained <- lipids_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_in_all <- lipids_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
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



# ########
# 
# ### Lipids (mice)
# 
# ########
# 
# message("Lipids (mice) -------------------------------------------")
# 
# imputed_c_lipids_mice5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM'), method = "mice")
# lipids_age_analysis_mice <- purrr::map(1:5, ~age_control_analysis(imputed_c_lipids_mice5, name = "Lipids", color = NULL, imp_num = .x))
# 
# lipids_pred_mice_df <- imputation_df(lipids_age_analysis_mice)
# 
# (lipids_age_mice_predtruth <- predtruth_plot(lipids_pred_mice_df, name = "Lipids (mice)")) 
# 
# lipids_in_all_mice <- lipids_age_analysis_mice %>% 
#   purrr::map(~.x[[1]] %>% names) %>% 
#   reduce(intersect) %>%
#   setdiff("(Intercept)")
# 
# 

########

### Lipids: Mice, loo

########






########

### Combined GOT + lipids

########

message("Combined GOT/lipids -------------------------------------------")

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
combined_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_combined_amelia5, name = "Combined GOT + Lipids", color = NULL, imp_num = .x))

combined_pred_df <- imputation_df(combined_age_analysis)

(combined_age_predtruth <- predtruth_plot(combined_pred_df, name = "Combined GOT + Lipids"))
Sys.sleep(5)
ggsave(filename = "age_clock_combined_pred.png")





combined_avg_retained <- combined_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

combined_in_all <- combined_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")





########

### Targeted

########

message("Targeted -------------------------------------------")

imputed_c_targeted5 <- filter_and_impute_multi(wide_data_targeted, c('CO', 'CY', 'CM'))
targeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_targeted5, name = "Targeted", color = NULL, imp_num = .x))
#save(targeted_pred_df, file = "targeted_age_analysis_02032020.RData")
targeted_pred_df <- imputation_df(targeted_age_analysis)


(targeted_age_predtruth <- predtruth_plot(targeted_pred_df, name = "Targeted")) 
ggsave(filename = "age_clock_targeted_pred.png") #14.9 x 8.21

targeted_avg_retained <- targeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

targeted_in_all <- targeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

 
#### wild test
m_targeted <- wide_data_targeted %>% 
  select_if(is.numeric) %>%
  select(-Age)


m_age <- wide_data_untargeted$Age

m_targeted %>%
  select_if(function(x) abs(cor(x, m_age)) > 0.3)

cor(m_age, m_targeted[,1])
for(i in 1:ncol(m_targeted)){
  if(abs(cor(m_targeted[,i], m_age)) > 0.3){
    m_
  }
}

mt_untargeted <- wide_data_untargeted %>%
  select_if(is.numeric) %>%
  select(-Age) %>%
  transmute_all(function(x) scale(x, center = T, scale = T) %>% Hmisc::impute()) 

# mt_imp <- amelia(t(mt_untargeted), empri = 10000000)
mt_imp <- mice(t(mt_untargeted))
# 1:5 %>%
#   purrr::map(~ complete(mt_imp, .x) %>%
#                t %>%
#                as_tibble() %>%
mt_untargeted %>%
  transmute_all(function(x) abs(cor(x, m_age))) %>% 
  slice(1) %>% 
  gather() %>% 
  ggplot() + 
  geom_histogram(aes(value)) #%>%
  #cowplot::plot_grid(plotlist = .)





untar_tiny <- wide_data_untargeted %>% select(setdiff(untargeted_in_at_least_one, "GenderM"), Age, Gender, Type, APOE, Id, GBA_T369M, GBAStatus)
imp_untar_tiny <- filter_and_impute_multi(untar_tiny, c("CO", "CY", "CM"), transpose = F, empri = 200)
untar_tiny_age_analysis <- purrr::map(1:5, ~age_control_analysis(imp_untar_tiny, name = "untargeted tiny", color = NULL, imp_num = .x, nlambda = 200))
untar_tiny_pred_df <- imputation_df(untar_tiny_age_analysis)
(untar_tiny_age_predtruth <- predtruth_plot(untar_tiny_pred_df, name = "untar tiny"))
untar_tiny_avg_retained <- untar_tiny_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untar_tiny_in_all <- untar_tiny_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  reduce(c) %>%
  setdiff("(Intercept)")



untar_tiny_ad_amelia5 <- filter_and_impute_multi(untar_tiny, c('AD'), empri = 100)
untar_tiny_pd_amelia5 <- filter_and_impute_multi(untar_tiny, c('PD'), empri = 100)
untar_tiny_adpd_separate_amelia5 <- purrr::map2(untar_tiny_ad_amelia5, untar_tiny_pd_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE))
untar_tiny_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(imp_untar_tiny[[.x]], new_data = untar_tiny_adpd_separate_amelia5[[.x]], nlambda = 200))
untar_tiny_adpd_pred_df_separate <- imputation_df(untar_tiny_adpd_age_analysis_separate)
(untar_tiny_adpd_age_predtruth <- predtruth_plot(untar_tiny_adpd_pred_df_separate, name = "untar tiny adpd"))

# match
untar_tiny_adpd_matched_controls <- purrr::map(1:nrow(untar_tiny), ~find_control(.x, data = filter(untar_tiny, Type %in% c("AD", "PD")), 
                                                                                           data_control = filter(untar_tiny, Type %in% c("CO", "CY", "CM")))) %>%
  bind_rows(filter(untar_tiny, Type %in% c("AD", "PD"))) %>%
  distinct(.keep_all = T)

# imputation_df() expects a list of lists, but only uses the first element of the inner list. so we add a random "1" to the output. doesn't do anything.
untar_tiny_c_matched_age_analysis_table <- purrr::map(untar_tiny_age_analysis, ~list(.x[[1]] %>% filter(id %in% untar_tiny_adpd_matched_controls$Id), 1))
untar_tiny_c_matched_age_pred_df <- imputation_df(untar_tiny_c_matched_age_analysis_table)

(untar_tiny_c_matched_age_predtruth <- predtruth_plot(untar_tiny_c_matched_age_pred_df, name = "untar_tiny (matched)", color = NULL))

##### random n features

random_n_features_analysis <- function(n){
  untar_random_100_features <- wide_data_untargeted %>%
    names %>%
    setdiff(c("Age", "Gender", "Type", "APOE", "Id", "GBA_T369M", "GBAStatus")) %>%
    base::sample(size = n, replace = F)
  
  untar_tiny_random <- wide_data_untargeted %>% 
    select(untar_random_100_features, Age, Gender, Type, APOE, Id, GBA_T369M, GBAStatus)
  
  imp_untar_tiny_random <- filter_and_impute_multi(untar_tiny_random, c("CO", "CY", "CM"), transpose = T, empri = 200)
  untar_tiny_random_age_analysis <- purrr::map(1:5, ~age_control_analysis(imp_untar_tiny_random, name = "untargeted tiny random 100", color = NULL, imp_num = .x, nlambda = 200))
  untar_tiny_random_pred_df <- imputation_df(untar_tiny_random_age_analysis)
  untar_tiny_random_age_predtruth <- predtruth_plot(untar_tiny_random_pred_df, name = paste("untar tiny", n))
  
  untar_tiny_random_age_predtruth
  
}


random_100_features <- purrr::map(1:3, ~random_n_features_analysis(100))
random_200_features <- purrr::map(1:3, ~random_n_features_analysis(200))
random_500_features <- purrr::map(1:3, ~random_n_features_analysis(500))
random_1000_features <- purrr::map(1:3, ~random_n_features_analysis(1000))
random_2000_features <- purrr::map(1:3, ~random_n_features_analysis(2000))
random_3000_features <- purrr::map(1:3, ~random_n_features_analysis(3000))


cowplot::plot_grid(plotlist = random_100_features)
cowplot::plot_grid(plotlist = random_200_features)
cowplot::plot_grid(plotlist = random_500_features)
cowplot::plot_grid(plotlist = random_1000_features)
cowplot::plot_grid(plotlist = random_2000_features)
cowplot::plot_grid(plotlist = random_3000_features)


# match
untar_tiny_random_adpd_matched_controls <- purrr::map(1:nrow(untar_tiny_random), ~find_control(.x, data = filter(untar_tiny_random, Type %in% c("AD", "PD")), 
                                                                                 data_control = filter(untar_tiny_random, Type %in% c("CO", "CY", "CM")))) %>%
  bind_rows(filter(untar_tiny_random, Type %in% c("AD", "PD"))) %>%
  distinct(.keep_all = T)

# imputation_df() expects a list of lists, but only uses the first element of the inner list. so we add a random "1" to the output. doesn't do anything.
untar_tiny_random_c_matched_age_analysis_table <- purrr::map(untar_tiny_random_age_analysis, ~list(.x[[1]] %>% filter(id %in% untar_tiny_random_adpd_matched_controls$Id), 1))
untar_tiny_random_c_matched_age_pred_df <- imputation_df(untar_tiny_random_c_matched_age_analysis_table)

(untar_tiny_random_c_matched_age_predtruth <- predtruth_plot(untar_tiny_random_c_matched_age_pred_df, name = "untar_tiny_random (matched)", color = NULL))




#######

### Unargeted

########

message("untargeted -------------------------------------------")

imputed_c_untargeted5 <- filter_and_impute_multi(wide_data_untargeted, c('CO', 'CY', 'CM'), empri = 8)
#save(imputed_c_untargeted5, file = "imputed_c_untargeted5_03022020")
# lambdas were reaching the end of their sequence due to hugeness of data, so we increase nlambda from 100 to 200
untargeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted5, name = "untargeted", color = NULL, imp_num = .x, nlambda = 200))
#save(untargeted_age_analysis, file = "untargeted_age_analysis_02032020.RData")

untargeted_pred_df <- imputation_df(untargeted_age_analysis)

(untargeted_age_predtruth <- predtruth_plot(untargeted_pred_df, name = "Untargeted")) 
ggsave(filename = "age_clock_untargeted_pred.png") #14.9 x 8.21

untargeted_avg_retained <- untargeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_in_all <- untargeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

untargeted_in_at_least_one <- untargeted_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")

#plot colored by different things
predtruth_plot(untargeted_pred_df, color = "apoe", errorbar = FALSE, name = "Untargeted")
predtruth_plot(untargeted_pred_df, color = "apoe4", errorbar = FALSE, name = "Untargeted")


########

### Untargeted, select columns with <90% missingness

########


message("Untargeted  subset-------------------------------------------")



# What about doing untargeted analysis on <90% missing?
untargeted_less_90perc_missing <-wide_data_untargeted %>%
  map_dbl(~ sum(is.na(.x))/length(.x)) %>% 
  enframe(value = "perct_missing") %>% filter(perct_missing < .9)

wide_data_untargeted_less90perc_missing <- wide_data_untargeted %>% select_at(vars(untargeted_less_90perc_missing$name))    


# analysis
imputed_c_untargeted_90subset5 <- filter_and_impute_multi(wide_data_untargeted_less90perc_missing , c('CO', 'CY', 'CM'))
untargeted_90subset_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted_90subset5, name = "untargeted_90subset", color = NULL, imp_num = .x))

untargeted_90subset_pred_df <- imputation_df(untargeted_90subset_age_analysis)

(untargeted_90subset_age_predtruth <- predtruth_plot(untargeted_90subset_pred_df, name = "untargeted_90subset")) 
ggsave(filename = "age_clock_untargeted_90subset_pred.png") #14.9 x 8.21

untargeted_90subset_avg_retained <- untargeted_90subset_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_90subset_in_all <- untargeted_90subset_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")










########

### Untargeted, select columns with <10% missingness

########


message("Untargeted  subset-------------------------------------------")



# What about doing untargeted analysis on <10% missing?
untargeted_less_10perc_missing <-wide_data_untargeted %>%
  map_dbl(~ sum(is.na(.x))/length(.x)) %>% 
  enframe(value = "perct_missing") %>% filter(perct_missing < .1)

wide_data_untargeted_less10perc_missing <- wide_data_untargeted %>% select_at(vars(untargeted_less_10perc_missing$name))    


# analysis
imputed_c_untargeted_10subset5 <- filter_and_impute_multi(wide_data_untargeted_less10perc_missing , c('CO', 'CY', 'CM'))
untargeted_10subset_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted_10subset5, name = "untargeted_10subset", color = NULL, imp_num = .x))

untargeted_10subset_pred_df <- imputation_df(untargeted_10subset_age_analysis)

(untargeted_10subset_age_predtruth <- predtruth_plot(untargeted_10subset_pred_df, name = "Untargeted (<10% missing)")) 
ggsave(filename = "age_clock_untargeted_10subset_pred.png") #14.9 x 8.21

untargeted_10subset_avg_retained <- untargeted_10subset_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_10subset_in_all <- untargeted_10subset_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")


########

### Untargeted, select columns with <10% missingness, randomly shuffling the ages to test

########


message("Untargeted  subset, permuted age-------------------------------------------")

wide_data_untargeted_permuted <- wide_data_untargeted %>%
  mutate(Age = sample(Age, replace = F))

# What about doing untargeted analysis on <10% missing?
untargeted_permuted_less_10perc_missing <-wide_data_untargeted_permuted %>%
  map_dbl(~ sum(is.na(.x))/length(.x)) %>% 
  enframe(value = "perct_missing") %>% filter(perct_missing < .1)

wide_data_untargeted_permuted_less10perc_missing <- wide_data_untargeted_permuted %>% select_at(vars(untargeted_less_10perc_missing$name))    


# analysis
imputed_c_untargeted_permuted_10subset5 <- filter_and_impute_multi(wide_data_untargeted_permuted_less10perc_missing , c('CO', 'CY', 'CM'))
untargeted_permuted_10subset_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted_permuted_10subset5, name = "untargeted_permuted_10subset", color = NULL, imp_num = .x))

untargeted_permuted_10subset_pred_df <- imputation_df(untargeted_permuted_10subset_age_analysis)

(untargeted_permuted_10subset_age_predtruth <- predtruth_plot(untargeted_permuted_10subset_pred_df, name = "untargeted_permuted (<10% missing)")) 
ggsave(filename = "age_clock_untargeted_permuted_10subset_pred.png") #14.9 x 8.21

untargeted_permuted_10subset_avg_retained <- untargeted_permuted_10subset_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_permuted_10subset_in_all <- untargeted_permuted_10subset_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")





##################################################################



### Other tests ###



###################################################################





#########################

### Predict on AD/PD ###

########################



################

## Using combined GOT + lipids

################

message("combined got+lipids, ADPD  age-------------------------------------------")

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
  fit <- fit_glmnet(combined_data[1:nrow(imputation[[1]]),], true_control_age, alpha = 0.5, penalize_age_gender = FALSE, penalize_AD_PD = FALSE, family = "gaussian", nlambda = nlambda)
  oos_pred <- predict(fit, newx = combined_data[-(1:nrow(imputation[[1]])),], s = "lambda.min")
  
  
  data_df <- tibble(
    truth = true_pred_age,
    pred = oos_pred,
    gender = new_data[[1]][,"GenderM"] %>%
      as.factor %>%
      fct_recode(M = '1', F = '0'),
    type = new_data[[2]],
    apoe = new_data[[3]],
    id = new_data[[4]],
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




# ################
# 
# ## Using untargeted
# 
# ################
# 
# #empri = .1(num PD + num_AD)
# untargeted_adpd_combined_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('AD', 'PD'), empri = 10)
# 
# 
# untargeted_adpd_age_analysis <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_untargeted5[[.x]], new_data = untargeted_adpd_combined_amelia5[[.x]], nlambda = 200))
# 
# untargeted_adpd_pred_df <- imputation_df(untargeted_adpd_age_analysis)
# 
# (untargeted_adpd_age_predtruth <- predtruth_plot(untargeted_adpd_pred_df, name = "Untargeted", data_name = "AD/PD")) 
# ggsave(filename = "age_c_adpd_pred_untargeted.png") #14.9 x 8.21
# 



################

## Using untargeted, separate imputation for AD/PD

################

message("Untargeted, ADPD age, separate imp -------------------------------------------")

#empri = .1(num PD/AD)
untargeted_ad_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('AD'), empri = 5)
untargeted_pd_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('PD'), empri = 5)

untargeted_adpd_separate_amelia5 <- purrr::map2(untargeted_ad_amelia5, untargeted_pd_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE))

untargeted_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_untargeted5[[.x]], new_data = untargeted_adpd_separate_amelia5[[.x]], nlambda = 200))

untargeted_adpd_pred_df_separate <- imputation_df(untargeted_adpd_age_analysis_separate)

# will use for matched
adpd_null_pred <- mean(untargeted_adpd_pred_df_separate$truth) %>% rep(times = length(untargeted_adpd_pred_df_separate$truth))
adpd_rmse_null <- (untargeted_adpd_pred_df_separate$truth - adpd_null_pred)^2 %>% mean %>% sqrt %>% round(2)
adpd_mae_null <- (untargeted_adpd_pred_df_separate$truth - adpd_null_pred) %>% abs %>% mean %>% round(2)


(untargeted_adpd_age_predtruth_separate <- predtruth_plot(untargeted_adpd_pred_df_separate, name = "Untargeted", data_name = "AD/PD, separate", color = "type")) 
ggsave(filename = "age_c_adpd_pred_untargeted_separate.png") #14.9 x 8.21

untargeted_adpd_pred_df_separate %>% mutate(pred_dist =abs(imp_avg - truth)) %>%
  arrange(desc(pred_dist)) %>% select(-apoe4)





################

## Using untargeted, Comparing metrics for matched controls/ ADPD

################

message("untargeted ADPD matched age -------------------------------------------")

# match controls for AD/PD, join with the ADPD data, and remove duplicates
untargeted_adpd_matched_controls <- purrr::map(1:nrow(wide_data_untargeted), ~find_control(.x, data = filter(wide_data_untargeted, Type %in% c("AD", "PD")), 
                                                                                           data_control = filter(wide_data_untargeted, Type %in% c("CO", "CY", "CM")))) %>%
  bind_rows(filter(wide_data_untargeted, Type %in% c("AD", "PD"))) %>%
  distinct(.keep_all = T)

# imputation_df() expects a list of lists, but only uses the first element of the inner list. so we add a random "1" to the output. doesn't do anything.
untargeted_c_matched_age_analysis_table <- purrr::map(untargeted_age_analysis, ~list(.x[[1]] %>% filter(id %in% untargeted_adpd_matched_controls$Id), 1))
untargeted_c_matched_age_pred_df <- imputation_df(untargeted_c_matched_age_analysis_table)

(untargeted_c_matched_age_predtruth <- predtruth_plot(untargeted_c_matched_age_pred_df, name = "Untargeted (matched)", color = NULL))
ggsave(filename = "age_clock_untargeted_pred_c_matched.png") #14.9 x 8.21


adpd_matched_rsq <- cor(untargeted_c_matched_age_pred_df$truth, untargeted_c_matched_age_pred_df$imp_avg)^2 %>% round(2)
adpd_matched_rmse <- (untargeted_c_matched_age_pred_df$truth - untargeted_c_matched_age_pred_df$imp_avg)^2 %>% mean %>% sqrt %>% round(2)
adpd_matched_mae <- abs(untargeted_c_matched_age_pred_df$truth - untargeted_c_matched_age_pred_df$imp_avg) %>% mean %>% round(2)

# add indicator for points we want to highlight
untargeted_adpd_pred_df_separate <- untargeted_adpd_pred_df_separate %>%
  mutate(mse = (truth - imp_avg)^2) %>%
  group_by(type) %>%
  mutate(biggest_error = ifelse(mse == max(mse), 1,0)) %>%
  ungroup

## Update regular ADPD plot with the stats for this one:
(untargeted_adpd_age_predtruth_separate <- predtruth_plot(untargeted_adpd_pred_df_separate, name = "Untargeted", data_name = "AD/PD, separate", color = "type") +
    geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, imp_avg, method = "pearson")^2 %>% round(2), " (",adpd_matched_rsq,")", 
                                                                              "<br>RMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2), " (",adpd_matched_rmse,")", " (",adpd_rmse_null,")", 
                                                                              "<br>MAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2), " (",adpd_matched_mae,")", " (",adpd_mae_null,")")),
                  size =12) + 
    ggforce::geom_mark_circle(aes(truth, imp_avg, fill = type, filter = biggest_error ==1))) 
ggsave(filename = "age_c_adpd_pred_untargeted_separate.png") #14.9 x 8.21




### residuals plot, faceted by type

# first, we join the control preds with the ad/pd
untargeted_pred_df_with_adpd_oos <- untargeted_c_matched_age_pred_df %>% 
  bind_rows(untargeted_adpd_pred_df_separate) %>%
  mutate(type = fct_relevel(type, "CY", "CM", "CO", "AD", "PD"),
         type_c = fct_collapse(type, `Matched Control` = c("CY", "CM", "CO")) %>%
           fct_relevel("Matched Control", after = 0))

ggplot(untargeted_pred_df_with_adpd_oos) + 
  geom_density_ridges(aes(truth - imp_avg, type_c, fill = type_c), quantile_lines = T, quantiles = 2,
                      jittered_points = T,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7) + 
  labs(title = "Residuals of AD/PD predictions (fit on controls)",
       x = "True Age - Predicted Age",
       y = "Type",
       fill = "Type") + 
  geom_vline(xintercept =0, linetype = "dashed", color = "black")

# include all types instead of collapsed controls
ggplot(untargeted_pred_df_with_adpd_oos) + 
  geom_density_ridges(aes(truth - imp_avg, type, fill = type)) + 
  labs(title = "Residuals of AD/PD predictions (fit on controls)",
       x = "True Age",
       y = "True Age - Predicted Age")


#line graph?
ggplot(untargeted_pred_df_with_adpd_oos) + 
  geom_smooth(aes(truth, truth - imp_avg, color = type_c)) + 
  geom_rug(aes(truth)) +
  labs(title = "Residuals of AD/PD predictions (fit on controls)",
       x = "True Age - Predicted Age",
       y = "Type",
       fill = "Type")



#### asefkjasdklfjhalsdkfjh

################

## Using targeted, separate imputation for AD/PD

################

message("targeted, ADPD age, separate imp -------------------------------------------")

#empri = .1(num PD/AD)
targeted_ad_amelia5 <- filter_and_impute_multi(wide_data_targeted, c('AD'), empri = 10)
targeted_pd_amelia5 <- filter_and_impute_multi(wide_data_targeted, c('PD'), empri = 10)

targeted_adpd_separate_amelia5 <- purrr::map2(targeted_ad_amelia5, targeted_pd_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE))

targeted_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_targeted5[[.x]], new_data = targeted_adpd_separate_amelia5[[.x]], nlambda = 200))

targeted_adpd_pred_df_separate <- imputation_df(targeted_adpd_age_analysis_separate)


targeted_adpd_pred_df_separate <- targeted_adpd_pred_df_separate %>%
  mutate(mse = (truth - imp_avg)^2) %>%
  group_by(type) %>%
  mutate(biggest_error = ifelse(mse == max(mse), 1,0)) %>%
  ungroup

## Update regular ADPD plot with the stats for this one:
(targeted_adpd_age_predtruth_separate <- predtruth_plot(targeted_adpd_pred_df_separate, name = "targeted", data_name = "AD/PD, separate", color = "type") +
    geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, imp_avg, method = "pearson")^2 %>% round(2),  
                                                                              "<br>RMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2),
                                                                              "<br>MAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2))),
                  size =12) + 
    ggforce::geom_mark_circle(aes(truth, imp_avg, fill = type, filter = biggest_error ==1))) 





### asl;khfaklsdhf

########

### Untargeted, including ADPD (ie full model)
### We're interested in the coefficient for AD/PD

########

message("Untargeted ADPD full model -------------------------------------------")

# combine controls with AD/PD imputation. include indicators for AD/PD as predictors
untargeted_all_amelia5 <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, add_AD_ind = TRUE, add_PD_ind = TRUE))
#
#save(untargeted_all_amelia5, file = "untargeted_all_amelia5.RData")
untargeted_all_age_analysis <- purrr::map(1:5, ~age_control_analysis(untargeted_all_amelia5, name = "untargeted (all)", color = NULL, imp_num = .x, nlambda = 200))

untargeted_all_pred_df <- imputation_df(untargeted_all_age_analysis)

(untargeted_all_age_predtruth <- predtruth_plot(untargeted_all_pred_df, name = "Untargeted (all)")) 
ggsave(filename = "age_clock_all_untargeted_pred.png") #14.9 x 8.21



########

### Untargeted MATCHED, including ADPD (ie full model)
### We're interested in the coefficient for AD/PD
### The full full model has huge ADPD coeffs, so we're trying to narrow that down.

########

message("Untargeted ADPD full age model, matched subjects -------------------------------------------")

wide_data_untargeted_matched_c <- purrr::map_df(1:nrow(wide_data_untargeted), ~find_control(.x, data = filter(wide_data_untargeted, Type %in% c("AD", "PD")), 
                                                                                           data_control = filter(wide_data_untargeted, Type %in% c("CO", "CY", "CM")))) %>%
  distinct(.keep_all = TRUE)

imputed_matched_c_untargeted5 <- filter_and_impute_multi(wide_data_untargeted_matched_c, types = c("CO","CY", "CM"), empri = 5)
# combine controls with AD/PD imputation. include indicators for AD/PD as predictors
untargeted_all_matched_amelia5 <- purrr::map2(imputed_matched_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, add_AD_ind = TRUE, add_PD_ind = TRUE))


untargeted_all_matched_age_analysis <- purrr::map(1:5, ~age_control_analysis(untargeted_all_matched_amelia5, name = "untargeted (matched,all types) ", color = NULL, imp_num = .x, nlambda = 200))

untargeted_all_matched_pred_df <- imputation_df(untargeted_all_matched_age_analysis)

(untargeted_all_matched_age_predtruth <- predtruth_plot(untargeted_all_matched_pred_df, name = "Untargeted (all matched)")) 
ggsave(filename = "age_clock_all_matched_untargeted_pred.png") #14.9 x 8.21

untargeted_all_matched_age_analysis[[1]]$importance



################

## Using untargeted, separate imputation for AD/PD
## FIT ON ADPD data, pred on matched control

################

message("Untargeted Age model, fit on ADPD, predict on matched control -------------------------------------------")

untargeted_reverse_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(untargeted_adpd_separate_amelia5[[.x]], new_data = imputed_matched_c_untargeted5[[.x]], nlambda = 200))

untargeted_reverse_adpd_pred_df_separate <- imputation_df(untargeted_reverse_adpd_age_analysis_separate)

(untargeted_reverse_adpd_age_predtruth_separate <- predtruth_plot(untargeted_reverse_adpd_pred_df_separate, name = "Untargeted", data_name = "Control, separate (fit on AD/PD)", color = NULL)) 
ggsave(filename = "age_c_reverse_adpd_pred_untargeted_separate.png") #14.9 x 8.21



#########################

### univariate regression on Gender ~ Metabolite for each metabolite/lipid ###

########################

message("Combined GOT+Lipids univariate logistic regression on Gender ~ Metabolite -------------------------------------------")

c_combined_gender_df <- imputed_c_combined_amelia5[[1]][[1]] %>%
  as_tibble() %>%
  rename_all(function(x) str_replace_all(x, "`", "")) %>%
  select(-c('(Intercept)', Age))


c_combined_names <- names(c_combined_gender_df) %>% setdiff('GenderM')
c_combined_gender_p_values <- purrr::map(c_combined_names, ~age_metabolite_p(c_combined_gender_df, .x, var = "GenderM", family = "binomial")) %>%
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

message("Univariate Age ~ Metabolite (lm) -------------------------------------------")

#' Function to help process unviariate results
#' Create table with metabolite/lipid, p value, t value, name, bh corrected p value (for conc = FALSE)
#' Creates table with metabolite/lipid, variation explained (for conc = TRUE)
#' (using only the first imputation)
#' @return table with all of the variables with bh p values < 0.01
#' @param data is a imputation list (eg imputed_c_combined5)
#' @param imp_num is integer 1-5 (which imputation touse)
#' @param var/family/conc are as in age_metabolite_p. 
#' 
bh_univariate_age <- function(data, var = "Age", family = "gaussian", conc = FALSE) {
  df <- data[[1]][[1]] %>%
    as_tibble() %>%
    dplyr::select(-c('(Intercept)', GenderM)) %>%
    mutate_at(.vars = vars(-one_of("Age", "AD_ind", "PD_ind")), .funs = ~scale(.x, center = T, scale = T))
  
  p_table <- df %>%
    names %>%
    setdiff(c("Age",var)) %>%
    purrr::map(~age_metabolite_p(df, metabolite = .x, var = var, family = family, conc = conc)) %>%
    purrr::reduce(rbind) %>%
    as_tibble()
  
  if(conc == FALSE){
    p_values <- p_table %>%
      dplyr::rename('og_p_value' =  1)
    p_table <- cbind(p_values,'bh_p_value' = p.adjust(p_values$og_p_value, method = 'BH'))
  }
   
  p_table %>% 
    dplyr::mutate(name = str_replace_all(name, '`', ''))
  
}




### Using onlyt he first imputation -----------
## GOT
c_got_univariate_table <- bh_univariate_age(imputed_c_got5) %>%
  filter(bh_p_value < 0.01)


## Lipids
c_lipids_univariate_table <- bh_univariate_age(imputed_c_lipids5) %>%
  filter(bh_p_value < 0.01)


## Combined GOT + Lipids
c_combined_univariate_table <- bh_univariate_age(imputed_c_combined_amelia5) %>%
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

message("Untargeted Gender Logistic -------------------------------------------")

c_gender_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_c_untargeted5, varname = "GenderM", imp_num = .x, nlambda = 200))

c_gender_logistic %>% 
  purrr::map(~.x[[2]]) %>% 
  cowplot::plot_grid(plotlist = .)


c_gender_logistic[[1]][[2]] + 
  labs(title = "GenderM vs C")
ggsave('gender_logistic_c.png')





############################

### Trying MSEA from MetaboAnalystR on Targeted significant age

############################

message("MSEA for age significant (targeted) -------------------------------------------")

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

#univar_targeted_names_kegg <- univar_targeted_names_all %>% left_join(kegg_lookup, by = "name_formatted")


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
mSet<-SetMetabolomeFilter(mSet, T)

# Select metabolite set library

#pathway associated library?
#mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway",2)

# csf associated library?
mSet<-SetCurrentMsetLib(mSet, "csf",2)

# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "bar", "png", 72, width=NA)











############################

### univariate RF on metabolite concentration ~ Age ###
### To see if there's a natural plateu in metabolite concentration as a function of age

############################

message("Univariate Metabolite concentration ~ age (RF untargeted) -------------------------------------------")


# in case we want to use a full dataset
untargeted_all_amelia5_no_ind <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE,
                                                                                                                      add_AD_ind = FALSE, add_PD_ind = FALSE))
#
# using full dataset
untargeted_univariate_full_conc_age_scale_rf <- bh_univariate_age(untargeted_all_amelia5_no_ind, var = "Age", family = "rf", conc = TRUE)


# using just controls
untargeted_univariate_c_conc_age_scale_rf <- bh_univariate_age(imputed_c_untargeted5, var = "Age", family = "rf", conc = TRUE)


#' function to plot a line graph with concentration by age
#' @data is the result of bh_univariate_age()
#' @n is the number of metabolites you want to plot.
#' @names is an optional vector of metabolite names to plot (if null, uses top n var explained)
concentration_rf_plot <- function(data, n, names = NULL){
  
  if(is.null(names)){
    names <- data %>%
      arrange(desc(var_explained)) %>%
      select(name) %>%
      unique %>%
      slice(1:n) %>%
      deframe
  } 
  
  # Plot these 10 metabolites
  conc_vartop10 <- data %>%
    filter(name %in% names) %>%
    # create dummy column for flipped prediction. will fill in later
    mutate(flipped_pred = pred)
  
  ## We are interested in the pattern of these curves, so let's do some transforms
  
  # flip curves that are sloping down, so that everything is sloping up.
  for(var in names) {
    linear_mod <- data %>%
      filter(name == var) %>%
      lm(pred ~ Age, data =.)
    
    if(linear_mod$coefficients["Age"] < 0){
      print(var)
      conc_vartop10 <- conc_vartop10 %>%
        mutate(flipped_pred = ifelse(name == var, -pred, flipped_pred))
    }
  }
  
  # create dummy column to fill in shifted_pred
  conc_vartop10 <- conc_vartop10 %>%
    mutate(shifted_pred = flipped_pred)
  
  # next, shift all curves to end at 0
  for(var in names) {
    # get concentration at endpoint (latest age)
    max_age_conc <- conc_vartop10 %>%
      filter(name == var) %>%
      filter(Age == max(Age)) %>%
      pull(flipped_pred)
    
    # shift everything to get the endpoint at 0.
    if(max_age_conc > 0){
      conc_vartop10 <- conc_vartop10 %>%
        mutate(shifted_pred = ifelse(name == var, flipped_pred - max_age_conc, shifted_pred))
    } else if (max_age_conc <= 0){
      conc_vartop10 <- conc_vartop10 %>%
        mutate(shifted_pred = ifelse(name == var, flipped_pred + max_age_conc, shifted_pred))
    }
  }
  
  #add column with avg
  conc_vartop10 <- conc_vartop10 %>%
    group_by(Age) %>%
    mutate(avg_by_age_raw = mean(pred),
           avg_by_age_shifted = mean(shifted_pred),
           avg_by_age_flipped = mean(flipped_pred)) %>%
    ungroup
  
  # Look at no change
  raw <- conc_vartop10 %>%
    ggplot() + 
    geom_line(aes(group = name, x= Age, y = pred), color = 'gray',size = 1, show.legend = FALSE)  +
    geom_rug(aes(Age, pred),sides = "b")+
    #geom_line(aes(x = Age, y = avg_by_age_raw), size = 1, color = 'red') +
    #gghighlight(avg_by_age, max_highlight = 1, use_group_by = FALSE) +
    labs(title = "Fitted Concentration as a function of Age",
         subtitle = "Untargeted, Random Forest: top 10 most significant, scaled",
         y = "Predicted Concentration")
  
  
  
  # Look at Flipped only
  flipped <- conc_vartop10 %>%
    ggplot() + 
    geom_line(aes(group = name, x= Age, y = flipped_pred), color = 'gray',size = 1, show.legend = FALSE)  +
    geom_smooth(aes(x = Age, y = avg_by_age_flipped), size = 1, color = 'red') +
    geom_rug(aes(Age, flipped_pred),sides = "b")+
    #gghighlight(avg_by_age, max_highlight = 1, use_group_by = FALSE) +
    labs(title = "Fitted Concentration as a function of Age",
         subtitle = "Untargeted, Random Forest: top 10 most significant, scaled, sometimes flipped",
         y = "Predicted Concentration") +
   ylim(c(-2.5, 2.5))
  
  
  # Look at flipped and shifted
  flipped_shifted <- conc_vartop10 %>%
    ggplot() + 
    geom_line(aes(group = name, x= Age, y = shifted_pred), color = 'gray',size = 1, show.legend = FALSE)  +
    geom_smooth(aes(x = Age, y = avg_by_age_shifted), size = 1, color = 'red') +
    geom_rug(aes(Age, shifted_pred),sides = "b")+
    #gghighlight(avg_by_age, max_highlight = 1, use_group_by = FALSE) +
    labs(title = "Fitted Concentration as a function of Age",
         subtitle = "Untargeted, Random Forest: top 10 most significant, scaled, sometimes flipped, endpoint shifted",
         y = "Predicted Concentration") +
    ylim(c(-2,5,2.5))
  
  list("data" = conc_vartop10, "raw" = raw, "flipped" = flipped, "flip_shift" = flipped_shifted)
}



# 
# conc_vartop10 %>%
#   ggplot() + 
#   geom_line(aes(color = name, x= Age, y = truth), size = 1,show.legend = FALSE)  +
#   #geom_line(aes(x = Age, y = avg_by_age), size = 1, color = 'red') +
#   gghighlight(mean(truth), max_highlight = 1) + 
#   labs(title = "True concentration as a function of Age",
#        subtitle = "Untargeted, Random Forest: top 10 most significant",
#        y = "True Concentration")
# ggsave("untargeted_rf_true_conc_age_scale.png")





# Now try using the significant metabolites from the univariate Age ~ Concentration
# We would expect these to defy the trend


untargeted_age_sig_top10_names <- c_untargeted_univar_table_sig %>%
  arrange(bh_p_value) %>%
  slice(1:10) %>%
  pull(name) %>%
  unique


age_vartop10 <- untargeted_univariate_c_conc_age_scale_rf %>%
  filter(name %in% untargeted_age_sig_top10_names) %>%
  group_by(name, Age) %>%
  summarise(pred = median(pred)) %>% 
  ungroup() %>%
  #dummy column to be filled in later
  mutate(flipped_pred = pred)




## REPEAT transforms

# flip curves that are sloping down, so that everything is sloping up.
for(var in untargeted_age_sig_top10_names) {
  linear_mod <- untargeted_univariate_c_conc_age_scale_rf %>%
    filter(name == var) %>%
    lm(pred ~ Age, data =.)
  
  if(linear_mod$coefficients["Age"] < 0){
    print(var)
    age_vartop10 <- age_vartop10 %>%
      mutate(flipped_pred = ifelse(name == var, -pred, flipped_pred))
  }
}

# create dummy column to fill in shifted_pred
age_vartop10 <- age_vartop10 %>%
  mutate(shifted_pred = flipped_pred)

# next, shift all curves to end at 0
for(var in untargeted_age_sig_top10_names) {
  # get concentration at endpoint (latest age)
  max_age_conc <- age_vartop10 %>%
    filter(name == var) %>%
    filter(Age == max(Age)) %>%
    pull(flipped_pred)
  
  # shift everything to get the endpoint at 0.
  if(max_age_conc > 0){
    age_vartop10 <- age_vartop10 %>%
      mutate(shifted_pred = ifelse(name == var, flipped_pred - max_age_conc, shifted_pred))
  } else if (max_age_conc <= 0){
    age_vartop10 <- age_vartop10 %>%
      mutate(shifted_pred = ifelse(name == var, flipped_pred + max_age_conc, shifted_pred))
  }
}

#add column with avg
age_vartop10 <- age_vartop10 %>%
  group_by(Age) %>%
  mutate(avg_by_age_shifted = mean(shifted_pred),
         avg_by_age_flipped = mean(flipped_pred)) %>%
  ungroup






# there's a lot of overlap between this top 10 and the conc ~ age top 10. look at the diff
conc_age_sig_diff <- untargeted_age_sig_top10_names %>%
  setdiff(conc_age_scale_rf_top10)

  
age_vartop10 %>% ggplot() + 
  geom_line(aes(group = name, x= Age, y = flipped_pred), size = 1, color = 'gray',show.legend = FALSE)  +
  geom_smooth(aes(x = Age, y = avg_by_age_flipped)) +
  #gghighlight(name %in% conc_age_sig_diff, max_highlight = 2, use_group_by = FALSE) +
  labs(title = "Fitted Concentration as a function of Age",
       subtitle = "Untargeted, Random Forest: top 10 most significant in age ~ conc, flipped",
       y = "Predicted Concentration")


# 
full_conc_rf_top10 <- concentration_rf_plot(untargeted_univariate_full_conc_age_scale_rf, n = 50) 
c_conc_rf_top10 <- concentration_rf_plot(untargeted_univariate_c_conc_age_scale_rf, n = 10)

# using same metabolites as full in c
full_conc_names <- full_conc_rf_top10$data %>% pull(name) %>% unique
c_conc_rf_fullnames <- concentration_rf_plot(untargeted_univariate_c_conc_age_scale_rf, names = full_conc_names)



######

## Missingness

#####


# Looking by missingness by type and Dataset
percent_missing_by_type <- function(data, name){
  missing <- data %>% 
    group_by(Type) %>%
    group_map(~sum(is.na(.x)) / prod(dim(.x)))
  
  tibble(
    "source" = name,
    "CY" = missing[[4]], 
    "CM" = missing[[2]],
    "CO" = missing[[3]],
    "AD" = missing[[1]],
    "PD" = missing[[5]]
  )
}

data_sources <- list(wide_data, wide_data_lipids, wide_data_targeted, wide_data_untargeted)
data_names <- list("GOT", "Lipids", "Targeted", "Untargeted")
pct_missing_data <- purrr::map2(data_sources, data_names, ~percent_missing_by_type(.x, .y)) %>%
  bind_rows %>%
  gather(key = "Type", value = "pct_missing", -source) %>%
  mutate(Type = factor(Type, levels = c("CY", "CM", "CO", "AD", "PD")))


ggplot(pct_missing_data) + 
  geom_col(aes(x = source, y= pct_missing, fill = Type), position = 'dodge') +
  labs(title = "Percent Missingness by dataset by type",
       x = "Dataset",
       y = "Percent Missing")



# Run chi squared test on missingness between the five types

#get sample size by type
n_by_type_untargeted <- wide_data_untargeted %>%
  group_by(Type) %>% 
  tally() %>%
  spread(key = 'Type', value = 'n') %>%
  rename_all(function(x) paste('n', x, sep = '_'))

missingness_by_type_untargeted_counts <- wide_data_untargeted %>%
  group_by(Type) %>%
  group_map(~ map_int(.x, function(y) sum(is.na(y)))) %>%
  set_names(wide_data_untargeted$Type %>% droplevels %>% levels) %>%
  lapply(function(x) enframe(x, name = 'name', value = 'num_missing')) 

missingness_by_type_untargeted_pct <- wide_data_untargeted %>% 
  group_by(Type) %>%
  group_map(~ map_dbl(.x, function(y) round(sum(is.na(y))/length(y), digits = 3))) %>%
  set_names(wide_data_untargeted$Type %>% droplevels %>% levels) %>%
  lapply(function(x) enframe(x, name = 'name', value = 'pct_missing'))



missingness_by_type_all <- purrr::reduce(missingness_by_type_untargeted_counts, inner_join, by = 'name') %>%
  #set the names to show type
  set_names(c('name', paste('num_missing',levels(droplevels(wide_data_untargeted$Type)), sep = '_'))) %>%
  dplyr::filter(!(name %in% c('GBAStatus', 'cognitive_status', 'GBA_T369M'))) %>%
  cbind(n_by_type_untargeted) %>%
  dplyr::mutate(
    pct_missing_CY = num_missing_CY/n_CY,
    pct_missing_CM = num_missing_CM/n_CM,
    pct_missing_CO = num_missing_CO/n_CO,
    pct_missing_AD = num_missing_AD/n_AD,
    pct_missing_PD = num_missing_PD/n_PD
  ) %>%
  #filter(reduce(list(pct_missing_AD, pct_missing_CM,pct_missing_CO,pct_missing_CY,pct_missing_PD), `==`)) %>%
  rowwise() %>%
  #mutate(p_value = (prop.test(x =  str_subset(names(.), 'num_missing'), n = str_subset(names(.), 'n_')))$p.value)
  dplyr::mutate(p_value = (prop.test(x = c(num_missing_AD, num_missing_CM, num_missing_CO, num_missing_CY, num_missing_PD), 
                                     n = c(n_AD,n_CM,n_CO,n_CY,n_PD)))$p.value) %>%
  cbind('bh_q_value' = p.adjust(.$p_value, method = 'BH')) %>%
  dplyr::filter(bh_q_value < 0.01) %>%
  gather('Type', 'pct_missing', contains('pct_missing')) %>%
  mutate(Type = factor(Type, levels = c("pct_missing_CY","pct_missing_CM","pct_missing_CO","pct_missing_AD","pct_missing_PD")))


ggplot(missingness_by_type_all, aes(pct_missing, name)) +
  geom_point(aes(color = Type), size = 3, position = position_jitter(width = 0.01, height = 0,seed = 1)) + 
  scale_color_manual(labels = c("CY", "CM", "CO", "AD", "PD"), values = c("lightskyblue", "dodgerblue", "blue", "darkgreen", "purple")) +
  theme(axis.text.x = element_text(angle= 90, hjust =1)) + 
  labs(title = 'Percent Missingness by Type',
       subtitle = 'GOT + Lipids, filtered using 2-tailed Pearson Chi squared test, BH q < 0.05',
       x = "Percent Missing",
       y = "Name")


# correlation between missingness and the significant metabolites?
untargeted_10subset_in_all



#### Fitting a series of univariate missingness ~ metabolite + Gender tests ------------------------------------

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

na_univar_table <- function(data){
  missingness_some_na_names_untargeted <- data %>% 
    #leaving only metabolites/lipids
    dplyr::select(-one_of("Type", "Gender", "Age", "APOE", "GBAStatus", "Id", 
                          "GBA_T369M")) %>%
    #need at least one NA, but not all NA
    dplyr::select_if(function(x) any(is.na(x)) & !all(is.na(x))) %>%
    names
  
  form_missing_genderAge <- formula(missing~Gender + Age)
  missingness_overall_p_values_untargeted <- purrr::map(missingness_some_na_names_untargeted, ~missingness_metabolite_p(data, form_missing_genderAge, .x)[[1]]) %>%
    unlist
  missingness_age_p_values_untargeted <- purrr::map(missingness_some_na_names_untargeted, ~missingness_metabolite_p(data, form_missing_genderAge, .x)[[2]]) %>%
    unlist
  missingness_gender_p_values_untargeted <- purrr::map(missingness_some_na_names_untargeted, ~missingness_metabolite_p(data, form_missing_genderAge, .x)[[3]]) %>%
    unlist
  
  
  #bh corrected controls for false discovery rate.
  missingness_p_table_untargeted <- bind_cols('name' = missingness_some_na_names_untargeted, 
                                              'og_overall_p_value' = missingness_overall_p_values_untargeted,
                                              'og_age_p_value' = missingness_age_p_values_untargeted,
                                              'og_gender_value' = missingness_gender_p_values_untargeted,
                                              'bh_overall_q_value' = p.adjust(missingness_overall_p_values_untargeted, method = 'BH'),
                                              'bh_age_q_value' = p.adjust(missingness_age_p_values_untargeted, method = 'BH'),
                                              'bh_gender_q_value' = p.adjust(missingness_gender_p_values_untargeted, method = 'BH'))
  
}




#### Fitting a missingness model with a 1/0 matrix ---------------------------------------






sink(type="message")
close(error_log)
