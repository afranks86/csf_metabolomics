## test of reduced feature modeling
source(here::here("analysis", "starter.R"))

# attempt 1:
# for each obs
# get missingness cols
# take them out of the rest of the data
# impute rest of the data
# fit model

# attempt 2: 
# for each obs
# impute
# get missingness cols
# take them out of the imputed data
# fit model


#' Fit an elastic net model
#' @param imp_num is an integer, one of 1:length(imps)
#' @param imps is the output of filter_and_imput_multi()
#' @param target is a string,  the name of a variable in imps[[imp_num]]$Y, used as response
#' @param family is a string, either gaussian or binomial
#' @param alpha is a numeric between 0 and 1, the elastic net param
#' @param nlambda is a positive number, how far glmnet searches for values of lambda
#' @param *_indicator... deprecated?
#' @param penalize_* are logicals, indicating whether we should penalize them in the fitting process.
train_reduced_model <- function(imp_num, imps, target = "Age", 
                                family = "gaussian", 
                                alpha = 0.5, nlambda = 200, 
                                AD_ind= F, PD_ind = F,
                                GBA_ind = F, APOE4_ind = F,
                                penalize_age_gender = F,
                                penalize_AD_PD = T,
                                types
                                ){
  message(str_glue('imp_num is {imp_num}'))
  imputed_c <- imps[[imp_num]]
  imputed_c_Y <- imputed_c[[1]]
  imputed_c_type <- imputed_c[[2]]
  imputed_c_apoe <- imputed_c[[3]]
  imputed_c_id <- imputed_c[[4]]
  imputed_c_gender <- imputed_c_Y[,'GenderM']
  
  # if types is null, we include all types
  if(!is.null(types)){
    types_index <- which(imputed_c_type %in% types)
    imputed_c_Y <- imputed_c_Y[types_index,]
    imputed_c_type <- imputed_c_type[types_index]
    imputed_c_apoe <- imputed_c_apoe[types_index]
    imputed_c_gender <- imputed_c_gender[types_index]
  }
  
  
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
  } else if(APOE4_ind){
    imputed_c_target <- imputed_c_apoe %>%
      as.character() %>%
      str_detect("4")
  } else{
    imputed_c_target <- imputed_c_Y[,target]
  }
  
  if(family == "binomial"){
    # calculate what fraction of the total each class has
    freq_frac <- table(imputed_c_target)/length(imputed_c_target)
    # assign 1 - that value to a "weights" vector
    weights <- 1 - freq_frac[as.character(imputed_c_target)]    
  } else{
    weights <- rep(1, length(imputed_c_target))
  }
  
  # readd type as a feature in this analysis
  imputed_c_features_tmp <- imputed_c_Y %>% 
    as_tibble %>%
    #mutate(Type = imputed_c_type) %>%
    #1/8select(-Age)
    # we'll add back intercept in model.matrix below.
    select(-any_of(c(target, "(Intercept)")))
  
  
  #turn factors into dummy variables
  imputed_c_features <- model.matrix(~., imputed_c_features_tmp)
  
  #set penalty factors. (1 for each is default)
  p_factors <- rep(1, ncol(imputed_c_features))
  #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
  if(penalize_age_gender == FALSE){
    age_gender_index <- which(colnames(imputed_c_features) %in% c('Age', 'GenderM'))
    p_factors[age_gender_index] <- 0
  }
  if(penalize_AD_PD == FALSE){
    ad_pd_index <- which(colnames(imputed_c_features) %in% c("AD_ind", "PD_ind"))
    p_factors[ad_pd_index] <- 0
  }
  
  fit <- cv.glmnet(imputed_c_features, imputed_c_target, 
                  family = family, alpha = alpha, 
                  standardize = TRUE, penalty.factor = p_factors, 
                  nlambda = nlambda, weights = weights
                  )
  
  list(fit = fit,
       feature_order = colnames(imputed_c_features)
       )
}



#' Elastic net regression on an n-1 sized dataset
#' where index is the row left out of data. Default settings fits the age prediction model
#' @param index is an integer in 1:nrow(data)
#' @param data is a dataframe, one of wide_data_* variants
#' @param target is a string, the name of the reponse variable in data
#' @param types is a string vector, some subset of CO, CY, CM, PD, AD. these will be used for analysis
#' @param empri0 is the starting empirical prior for imputation. will be incremented by 5 until code runs
#' @param all other params are exclusively for train_reduced_model
reduced_fit <- function(index, data, target = 'Age', 
                        types = c("CO", "CY", "CM"), 
                        empri0 = 50, family = "gaussian", 
                        AD_ind= F, PD_ind = F,
                        GBA_ind = F, APOE4_ind = F,
                        penalize_age_gender = F,
                        penalize_AD_PD = T){
  message(str_glue("index = {index}"))
  train <- data[-index,]
  test <- data[index,]
  
  test_meta <- test %>%
    select(any_of(metadata_cols)) %>%
    mutate(GBAStatus = as.factor(case_when(GBA_T369M == 'CT' ~ 'CT', 
                                           TRUE ~ GBAStatus))) %>%
    mutate(AD_ind = ifelse(Type == "AD", 1, 0),
           PD_ind = ifelse(Type == "PD", 1, 0),
           GBA_ind = ifelse(GBAStatus %in% c('E326K Carrier', 'Pathogenic Carrier', 'CT'), 1, 0),
           APOE4_ind = ifelse(str_detect(as.character(APOE), "4"), 1, 0)
           )
  
  train_used <- train %>%
    filter(Type %in% types)
    
  test_x <- test %>% 
    select(-any_of(metadata_cols)) 
  
  
  missingness_cols <- names(test_x[,which(is.na(test_x))])

  train_reduced <- train %>%
    select(-all_of(missingness_cols))
  test_reduced <- test %>%
    select(-all_of(missingness_cols))
  
  # if Age is not the reponse, include it in the model
  if(target != "Age"){
    unused_cols <- setdiff(metadata_cols, c('Gender', 'Age'))
  } else{
    unused_cols <- setdiff(metadata_cols, c('Gender'))
  }
  
  test_reduced_x <- test_reduced %>%
    select(-any_of(unused_cols)) %>%
    #select(-any_of(metadata_cols)) %>%
    # scale each column using the same scaling done on the training data
    future_map2_dfc(names(.), ~scale_data(.x, train_used, .y)) %>%
    model.matrix(~., .)
  
  train_imps <- try(filter_and_impute_multi(train_reduced, types, scale = T, empri = empri0), silent = T)
  while(class(train_imps) == "try-error"){
    empri0 <- empri0 + 5
    train_imps <- try(
      filter_and_impute_multi(train_reduced, types, scale = T, empri = empri0), 
      silent = T
      )
  }
  train_fits <- furrr::future_map(1:5, ~train_reduced_model(.x, train_imps, 
                                                            target = target,
                                                            family = family,
                                                            AD_ind = AD_ind,
                                                            PD_ind = PD_ind,
                                                            GBA_ind = GBA_ind,
                                                            APOE4_ind = APOE4_ind,
                                                            types = types,
                                                            penalize_age_gender = penalize_age_gender,
                                                            penalize_AD_PD = penalize_AD_PD
                                                            ))
  
  # making sure the column ordering is consistent between train and test
  # and converting result back into a matrix of appropriate size.
  test_reduced_x <- matrix(test_reduced_x[,train_fits[[1]]$feature_order], 
                           nrow = 1, 
                           dimnames = list("", train_fits[[1]]$feature_order) 
                           )
  
  test_preds <- purrr::map_dbl(1:5, ~predict(train_fits[[.x]]$fit, newx = test_reduced_x, s = "lambda.min", type = "response"))
  
  pred_df <- test_preds %>%
    enframe() %>%
    pivot_wider(names_from = "name", values_from = "value", names_prefix = "imp") %>%
    mutate(imp_avg = mean(c_across(starts_with('imp')))) %>%
    cbind(test_meta)
    
    
  list(
    fit = train_fits, 
    pred = pred_df,
    empri_final = empri0
    )
}



#' Fit all n leave one out gaussian elastic net aging models
#' @param data is one of the wide_data_* variant
#' @returns fit_age has each individual fit, a list where each output is an element of reduced_fit_age()
#' @returns pred_df is a dataframe with a summary of predictions
#' @returns figure is a truth v predicted plot
reduced_age_analysis <- function(data){
  set.seed(1)
  c_data <- data %>%
    filter(Type %in% c("CO", "CM", "CY"))
  data_reduced <- furrr::future_map(1:nrow(c_data), ~reduced_fit(.x, c_data))
  
  
  reduced_pred_df <- data_reduced %>%
    purrr::map_df(~.x$pred) %>%
    rowwise() %>%
    mutate(imp_avg = mean(c_across(starts_with('imp'))),
           imp_min = min(c_across(starts_with('imp'))),
           imp_max = max(c_across(starts_with('imp')))
    ) %>%
    ungroup() %>%
    rename_with(~str_to_lower(.x), everything()) %>%
    rename("truth" = age)
  
  reduced_plot <- predtruth_plot(reduced_pred_df, name = "Reduced")
  
  list(fit_age = data_reduced,
       pred_df = reduced_pred_df,
       figure = reduced_plot
  )
  
}





reduced_logistic_analysis <- function(data, target = 'GenderM', AD_ind = F, PD_ind = F, GBA_ind = F, APOE4_ind = F, types){
  set.seed(1)
  
  data_reduced <- furrr::future_map(1:nrow(data), ~reduced_fit(.x, data, 
                                                               family = 'binomial', target = target,
                                                               AD_ind = AD_ind, PD_ind = PD_ind,
                                                               GBA_ind = GBA_ind, APOE4_ind = APOE4_ind,
                                                               types = types
                                                               ))
  
  
  reduced_pred_df <- data_reduced %>%
    purrr::map_df(~.x$pred) %>%
    rowwise() %>%
    mutate(imp_avg = mean(c_across(starts_with('imp'))),
           imp_min = min(c_across(starts_with('imp'))),
           imp_max = max(c_across(starts_with('imp')))
    ) %>%
    ungroup() %>%
    rename("truth" = target) %>%
    rename_with(~str_to_lower(.x), -contains('_ind'))
    
  
  roc_mean <- fpr_tpr(reduced_pred_df$imp_avg, reduced_pred_df$truth)
  pr_mean <- fpr_tpr(reduced_pred_df$imp_avg, reduced_pred_df$truth, x = "rec", y = "prec") %>%
    mutate(pos_class_rate = mean(reduced_pred_df$truth))
  
  roc_min <- fpr_tpr(reduced_pred_df$imp_min, reduced_pred_df$truth)
  pr_min <- fpr_tpr(reduced_pred_df$imp_min, reduced_pred_df$truth, x = "rec", y = "prec") %>%
    mutate(pos_class_rate = mean(reduced_pred_df$truth))
  
  roc_max <- fpr_tpr(reduced_pred_df$imp_max, reduced_pred_df$truth)
  pr_max <- fpr_tpr(reduced_pred_df$imp_max, reduced_pred_df$truth, x = "rec", y = "prec") %>%
    mutate(pos_class_rate = mean(reduced_pred_df$truth))
  
  roc_df <- bind_rows("mean" = roc_mean,
                      "max" = roc_max,
                      "min" = roc_min,
                      .id = "imp")
  
  pr_df <- bind_rows("mean" = pr_mean,
                      "max" = pr_max,
                      "min" = pr_min,
                      .id = "imp")
  
  list(fit_binom = data_reduced,
       predtruth_df = reduced_pred_df,
       roc_df = roc_df,
       pr_df = pr_df
       )
  
}


targeted_reduced <- reduced_age_analysis(wide_data_targeted)
untargeted_reduced <- reduced_age_analysis(wide_data_untargeted)
got_reduced <- reduced_age_analysis(wide_data)
lipids_reduced <- reduced_age_analysis(wide_data_lipids)

targeted_reduced <- readRDS(here('aging_output_files', 'targeted_reduced.Rds'))
untargeted_reduced <- readRDS(here('aging_output_files', 'untargeted_reduced.Rds'))
got_reduced <- readRDS(here('aging_output_files', 'got_reduced.Rds'))
lipids_reduced <- readRDS(here('aging_output_files', 'lipids_reduced.Rds'))

saveRDS(targeted_reduced, file = here('aging_output_files', 'targeted_reduced.Rds'))
saveRDS(untargeted_reduced, file = here('aging_output_files', 'untargeted_reduced.Rds'))
saveRDS(lipids_reduced, file = here('aging_output_files', 'lipids_reduced.Rds'))
saveRDS(got_reduced, file = here('aging_output_files', 'got_reduced.Rds'))


#' Function to get the median coefficient across five imputations for one of the leave one out models
#' @param index is a natural number between 1 and the number of observations used in reduced_output
#' @param reduced_output is the output of reduced_age_analysis()
#' @return a vector of median retained coefficients across the five imputations for models fit on that loo model
retained_one_obs <- function(index, reduced_output){
  purrr::map(1:5, ~importance(reduced_output$fit_age[[index]]$fit[[.x]]$fit, metabolites = F) %>%
               enframe(value = str_glue("imp{.x}"))) %>%
    reduce(full_join, by = "name") %>%
    rowwise() %>%
    transmute(name, median_coef = median(c_across(starts_with('imp')), na.rm = T)) %>%
    ungroup() %>%
    deframe()
}


#' look at the retained coefs for fits across each of n leave one out models
#' select the ones that show up >95% of the time.
#' @param reduced_output is the output of reduced_age_analysis
retained_over_fits <- function(reduced_output){
  purrr::map(1:85, ~retained_one_obs(.x, reduced_output)) %>%
    bind_rows() %>%
    #select only the variables that were present in >95% of fits
    select_if(function(x) sum(is.na(x))/length(x) < .05) %>%
    map_dbl(~median(.x, na.rm = T))
}

targeted_retained <- retained_over_fits(targeted_reduced)  
untargeted_retained <- retained_over_fits(untargeted_reduced)  
got_retained <- retained_over_fits(got_reduced)
lipids_retained <- retained_over_fits(lipids_reduced)







####

targeted_ad_reduced <- reduced_logistic_analysis(data = wide_data_targeted, 
                                                 target ="AD_ind", AD_ind = T, 
                                                 types = c("AD", "CO", "CM", "CY")
                                                 )
got_ad_reduced <- reduced_logistic_analysis(data = wide_data, 
                                            target ="AD_ind", AD_ind = T, 
                                            types = c("AD", "CO", "CM", "CY")
)
untargeted_ad_reduced <- reduced_logistic_analysis(data = wide_data_untargeted, 
                                                 target ="AD_ind", AD_ind = T, 
                                                 types = c("AD", "CO", "CM", "CY")
                                                 )


lipids_ad_reduced <- reduced_logistic_analysis(data = wide_data_lipids, 
                                            target ="AD_ind", AD_ind = T, 
                                            types = c("AD", "CO", "CM", "CY")
                                            )


####

targeted_pd_reduced <- reduced_logistic_analysis(data = wide_data_targeted, 
                                                 target ="PD_ind", PD_ind = T, 
                                                 types = c("PD", "CO", "CM", "CY")
                                                 )
got_pd_reduced <- reduced_logistic_analysis(data = wide_data, 
                                            target ="PD_ind", PD_ind = T, 
                                            types = c("PD", "CO", "CM", "CY")
)


untargeted_pd_reduced <- reduced_logistic_analysis(data = wide_data_untargeted, 
                                                 target ="PD_ind", PD_ind = T, 
                                                 types = c("PD", "CO", "CM", "CY")
                                                 )
lipids_pd_reduced <- reduced_logistic_analysis(data = wide_data_lipids, 
                                                 target ="PD_ind", PD_ind = T, 
                                                 types = c("PD", "CO", "CM", "CY")
                                               )



  
