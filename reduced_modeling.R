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
                                types, post_select = F,
                                take_log = F, include_age_gender = T,
                                adpd_detrend = F
                                ){
  message(str_glue('imp_num is {imp_num}'))
  imputed_c <- imps[[imp_num]]
  imputed_c_Y <- imputed_c[[1]]
  imputed_c_type <- imputed_c[[2]]
  imputed_c_apoe <- imputed_c[[3]]
  imputed_c_id <- imputed_c[[4]]
  imputed_c_gender <- imputed_c_Y[,'GenderM']
  adpd_coef <- NULL
  
  # if types is null, we include all types
  if(!is.null(types)){
    # adpd_detrend required us to impute all types in reduced_fit()
    # save AD/PD data separately if this is the case
    if(adpd_detrend){
      adpd_ind <- which(imputed_c_type %in% c("AD", "PD"))
      imputed_adpd_Y <- imputed_c_Y[adpd_ind,]
      adpd_adind <- imputed_c_type[adpd_ind] == "AD"
    }
    
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
    
    if(all(target == 'Age', take_log)){
      imputed_c_target <- log(imputed_c_target)
    }
  }
  
  ## FIT AD vs PD model on n - 1 to remove 
  # condition on 
  if(adpd_detrend){
    if(include_age_gender == FALSE){
      age_gender_index_adpd <- which(colnames(imputed_adpd_Y) %in% c('Age', 'GenderM'))
      imputed_adpd_Y <- imputed_adpd_Y[,-age_gender_index_adpd]
    }
    adpd_mod <- cv.glmnet(imputed_adpd_Y, adpd_adind,
                       family = "binomial", alpha = alpha, 
                       standardize = TRUE, nlambda = nlambda
                       )
    
    adpd_coef <- coef(adpd_mod, s = "lambda.min") %>%
      as.matrix() %>%
      as_tibble(rownames= "id") %>%
      set_names(c("id", "beta")) %>%
      # what to do with intercepts?. this removes duplicate intercept
      filter(beta != 0)
    
    # orthogonalize dataset by ad/pd beta
    imputed_c_Y <- purrr::map_dfr(1:nrow(imputed_c_Y), ~sub_coef(.x, imputed_c_Y, adpd_coef, meta = T)) %>%
      as.matrix()
    
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
  
  
  if(include_age_gender == FALSE){
    age_gender_index <- which(colnames(imputed_c_features) %in% c('Age', 'GenderM'))
    imputed_c_features <- imputed_c_features[,-age_gender_index]
  }
  
  #set penalty factors. (1 for each is default)
  p_factors <- rep(1, ncol(imputed_c_features))
  
  #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
  if(penalize_age_gender == FALSE & include_age_gender == TRUE){
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
  
  if(post_select){
    fit_coefs <- coef(fit, s = 'lambda.min') %>%
      as.matrix() 
    nonzero_coefs <- names(fit_coefs[which(fit_coefs != 0),])
    
    features_post <- imputed_c_features[,which(colnames(imputed_c_features) %in% nonzero_coefs)]
    
    if(penalize_age_gender == FALSE){
      age_gender_index <- which(colnames(features_post) %in% c('Age', 'GenderM'))
      p_factors[age_gender_index] <- 0
    }
    if(penalize_AD_PD == FALSE){
      ad_pd_index <- which(colnames(features_post) %in% c("AD_ind", "PD_ind"))
      p_factors[ad_pd_index] <- 0
    }
    
    fit_post <- cv.glmnet(features_post, imputed_c_target, 
                     family = family, alpha = alpha, 
                     standardize = TRUE, penalty.factor = p_factors, 
                     nlambda = nlambda, weights = weights
    )
    
    # adpd coef is null unless adpd_detrend is TRUE
    return(list('fit' = fit,
                feature_order = colnames(features_post),
                'fit_post' = fit_post,
                adpd_coef = adpd_coef
                ))
  }
  # adpd coef is null unless adpd_detrend is TRUE
  list(fit = fit,
       feature_order = colnames(imputed_c_features),
       adpd_coef = adpd_coef
       )
}



#' Elastic net regression on an n-1 sized dataset
#' where index is the row left out of data. Default settings fits the age prediction model
#' @param index is an integer in 1:nrow(data)
#' @param data is a dataframe, one of wide_data_* variants
#' @param target is a string, the name of the reponse variable in data
#' @param types is a string vector, some subset of CO, CY, CM, PD, AD. these will be used for analysis
#' @param empri0 is the starting empirical prior for imputation. will be incremented by 5 until code runs
#' @param remove_cor_age removes features with cor < .3 with age, following horvath dog paper
#' @param all other params are exclusively for train_reduced_model
reduced_fit <- function(index, data, target = 'Age', 
                        types = c("CO", "CY", "CM"), 
                        empri0 = 50, num = 5, family = "gaussian", 
                        AD_ind= F, PD_ind = F,
                        GBA_ind = F, APOE4_ind = F,
                        penalize_age_gender = F,
                        penalize_AD_PD = T,
                        post_select = F,
                        remove_cor_age = F,
                        include_led = F,
                        take_log = F, include_age_gender = T,
                        adpd_detrend = F){
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
    
  test_x <- test %>% 
    select(-any_of(metadata_cols)) 
  
  
  missingness_cols <- names(test_x[,which(is.na(test_x))])

  train_reduced <- train %>%
    select(-all_of(missingness_cols))
    # mutate(across(all_of(missingness_cols), ~if_else(is.na(.x), 0, 1)))
  
  #could just set these all to 0, but keeping same format as sanity check
  test_reduced <- test %>%
    select(-all_of(missingness_cols))
    # mutate(across(all_of(missingness_cols), ~if_else(is.na(.x), 0, 1)))
  
  
  # if Age is not the reponse, include it in the model
  if(target != "Age"){
    unused_cols <- setdiff(metadata_cols, c('Gender', 'Age'))
  } else{
    unused_cols <- setdiff(metadata_cols, c('Gender'))
  }
  
  if(include_age_gender == F){
    unused_cols <- metadata_cols
  }
  
  
  # for scaling, filter out types (filter_and_impute_multi will do this for everything else)
  train_used <- train_reduced %>%
    filter(Type %in% types)
  
  # remove cols with <.3 correlation (following horvath dog paper)
  if(remove_cor_age){
    # make sure gender isn't there and age is there
    train_used_numeric <- train_used %>%
      select(-any_of(c(unused_cols, 'Gender') %>% 
                       setdiff('Age'))
             )
    
    # completely NA cols and cols with < .3 correlation with age
    low_cor_cols <- train_used_numeric %>%
      select(where(~is.na(cor(.x, train_used_numeric$Age, use = 'na.or.complete', method = 'spearman')) || abs(cor(.x, train_used_numeric$Age, use = 'complete.obs', method = 'spearman')) < .3)) %>%
      names()
    
    train_reduced <- train_reduced %>%
      select(-all_of(low_cor_cols))
    test_reduced <- test_reduced %>%
      select(-all_of(low_cor_cols))
  }
  
  test_reduced_x <- test_reduced %>%
    select(-any_of(unused_cols)) %>%
    #select(-any_of(metadata_cols)) %>%
    # scale each column using the same scaling done on the training data
    furrr::future_map2_dfc(names(.), ~scale_data(.x, train_used, .y)) %>%
    model.matrix(~., .)
  
  
  # check for missing values
  
  train_reduced_metabolites <- train_reduced %>%
    select(-any_of(metadata_cols))
  
  if(sum(is.na(train_reduced_metabolites)) > 0){
    # don't transpose before amelia if matrix is long
    # subtract 2 for intercept and gender
    if(ncol(test_reduced_x) - 2 < nrow(train_reduced)){
      transpose <- F
    } else{
      transpose <- T
    }
    
    ## if we are orthogonalizing by ad/pd regression, then we need to impute the ENTIRE dataset
    # which is different from what we were doing before (i.e. if we're doing PD v C, we didn't impute AD)
    # ....so create a new "types" vector to include both AD and PD
    if(adpd_detrend){
      if(!("AD" %in% types)){
        types_imp <- c(types, "AD")
      }
      if(!("PD" %in% types)){
        types_imp <- c(types, "PD")
      }
    } else{
      types_imp <- types
    }
    
    train_imps <- try(filter_and_impute_multi(train_reduced, types_imp, scale = T, empri = empri0, num = num, transpose = transpose, include_led = include_led), silent = T)
    while(class(train_imps) == "try-error"){
      empri0 <- empri0 + 10
      message(str_glue('empri0 = {empri0}'))
      train_imps <- try(
        filter_and_impute_multi(train_reduced, types_imp, scale = T, empri = empri0, num = num), 
        silent = T
      )
    }
  } else{
    # if there are no missing values, no imputation is needed, so set impute = F
    train_imps <- filter_and_impute_multi(train_reduced, types_imp,  scale = T, empri = empri0, impute = F)
    
  }
  
  
  # with_progress({
  #   p <- progressor(steps = length(train_imps))
  #   train_fits <- furrr::future_map(1:length(train_imps), ~{
  #     p()
  #     Sys.sleep(.2)
  #     train_reduced_model(.x, train_imps,
  #                         target = target,
  #                         family = family,
  #                         AD_ind = AD_ind,
  #                         PD_ind = PD_ind,
  #                         GBA_ind = GBA_ind,
  #                         APOE4_ind = APOE4_ind,
  #                         types = types,
  #                         penalize_age_gender = penalize_age_gender,
  #                         penalize_AD_PD = penalize_AD_PD,
  #                         post_select = post_select,
  #                         take_log = take_log, include_age_gender = include_age_gender,
  #                         adpd_detrend = adpd_detrend
  #     )
  #   }, .options = furrr_options(seed = T))
  # })
  
  train_fits <- purrr::map(1:length(train_imps), ~train_reduced_model(.x, train_imps,
                                                            target = target,
                                                            family = family,
                                                            AD_ind = AD_ind,
                                                            PD_ind = PD_ind,
                                                            GBA_ind = GBA_ind,
                                                            APOE4_ind = APOE4_ind,
                                                            types = types,
                                                            penalize_age_gender = penalize_age_gender,
                                                            penalize_AD_PD = penalize_AD_PD,
                                                            post_select = post_select,
                                                            take_log = take_log, include_age_gender = include_age_gender,
                                                            adpd_detrend = adpd_detrend
  ))
  
  
  # making sure the column ordering is consistent between train and test, and cut down to post select cols
  # and converting result back into a matrix of appropriate size.
  
  test_reduced_x <- purrr::map(train_fits, 
                               ~{
                                 # if adpd_detrend, then need to orthogonalize the test row by beta
                                 if(adpd_detrend){
                                   new_test_reduced_x_tmp <- sub_coef(1, test_reduced_x, .x$adpd_coef, meta = T)
                                   new_test_reduced_x <- matrix(new_test_reduced_x_tmp, nrow = 1, dimnames = list("", names(new_test_reduced_x_tmp)))
                                 } else{
                                   new_test_reduced_x <- test_reduced_x
                                 }
                                 matrix(new_test_reduced_x[,.x$feature_order], 
                                        nrow = 1, 
                                        dimnames = list("", .x$feature_order)
                                 )
                               }
  )
  
  if(post_select){
    test_preds <- purrr::map_dbl(1:length(train_fits), ~predict(train_fits[[.x]]$fit_post, newx = test_reduced_x[[.x]], s = "lambda.min", type = "response"))
  } else{
    # # making sure the column ordering is consistent between train and test
    # # and converting result back into a matrix of appropriate size.
    # test_reduced_x <- matrix(test_reduced_x[,train_fits[[1]]$feature_order], 
    #                          nrow = 1, 
    #                          dimnames = list("", train_fits[[1]]$feature_order) 
    # )
    test_preds <- purrr::map_dbl(1:length(train_fits), ~predict(train_fits[[.x]]$fit, newx = test_reduced_x[[.x]], s = "lambda.min", type = "response"))
  }
  
  if(take_log){
    test_preds <- exp(test_preds)
  }
  
  pred_df <- test_preds %>%
    enframe() %>%
    pivot_wider(names_from = "name", values_from = "value", names_prefix = "imp") %>%
    mutate(imp_avg = mean(c_across(starts_with('imp'))),
           na_cols = n_distinct(missingness_cols)) %>%
    cbind(test_meta)
    
    
  list(
    fit = train_fits, 
    pred = pred_df,
    empri_final = empri0
    )
}



#' Fit all n leave one out gaussian elastic net aging models
#' @param data is one of the wide_data_* variant
#' @param remove_cor_age and post_select are as in reduced_fit(.)
#' @returns fit_age has each individual fit, a list where each output is an element of reduced_fit_age()
#' @returns pred_df is a dataframe with a summary of predictions
#' @returns figure is a truth v predicted plot
reduced_age_analysis <- function(data, target = 'Age', types = c('CO', 'CM', 'CY'), post_select = F, remove_cor_age = F, include_led = F, take_log = F){
  set.seed(1)
  c_data <- data %>%
    filter(Type %in% types)
  data_reduced <- furrr::future_map(1:nrow(c_data), ~reduced_fit(.x, c_data, types = types, target = target, post_select = post_select, remove_cor_age = remove_cor_age, include_led = include_led, take_log = take_log),
                                    .options = furrr_options(seed = TRUE))
  
  
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





reduced_logistic_analysis <- function(data, target = 'GenderM', AD_ind = F, PD_ind = F, GBA_ind = F, APOE4_ind = F, types, empri0 = 50, num = 5, penalize_age_gender = F, include_age_gender = T, adpd_detrend = F){
  set.seed(1)
  
  message(str_glue('target is {target} and ncol is {ncol(data)}'))
  data_to_use <- data #%>%
    #filter(Type %in% types)
  
  types_index <- which(data$Type %in% types)
  
  with_progress({
    p <- progressor(steps = length(types_index))
    data_reduced <- furrr::future_map(types_index, ~{
      p()
      Sys.sleep(.2)
      reduced_fit(.x, data_to_use,
        family = 'binomial', target = target,
        AD_ind = AD_ind, PD_ind = PD_ind,
        GBA_ind = GBA_ind, APOE4_ind = APOE4_ind,
        types = types, empri0 = empri0, num = num,
        penalize_age_gender = penalize_age_gender,
        include_age_gender = include_age_gender,
        adpd_detrend = adpd_detrend
      )
      }, .options = furrr_options(seed = T))
  })
  
  # data_reduced <- purrr::map(types_index, ~reduced_fit(.x, data_to_use,
  #                                                              family = 'binomial', target = target,
  #                                                              AD_ind = AD_ind, PD_ind = PD_ind,
  #                                                              GBA_ind = GBA_ind, APOE4_ind = APOE4_ind,
  #                                                              types = types, empri0 = empri0, num = num,
  #                                                              penalize_age_gender = penalize_age_gender,
  #                                                              include_age_gender = include_age_gender,
  #                                                              adpd_detrend = adpd_detrend
  # ))
  
  
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

roc_fig <- function(roc_df, predtruth_df){
  #brier_score <- with(predtruth_df, brier(imp_avg, truth))
  
  
  ggplot(roc_df) + 
    geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    labs(#title = "ROC: AD vs C",
      #subtitle = TeX('Untargeted'),
      x = 'False Positive Rate',
      y = 'True Positive Rate') + 
    geom_richtext(x = Inf, y = -Inf, 
                  hjust = 1, vjust = 0, 
                  size = rel(7),# label.padding = unit(1, "lines"),
                  label = paste0('AUC:', round(mean(roc_df$auc), 3) 
                                 #'<br>Brier:', round(brier_score, 3)
                  )
    ) +
    theme(
      axis.text = element_text(size = rel(1.25)),
      axis.title = element_text(size = rel(2)),
      plot.title = element_text(size = rel(2.5))
    )
}


#' Function to get the median coefficient across five imputations for one of the leave one out models
#' @param index is a natural number between 1 and the number of observations used in reduced_output
#' @param reduced_output is the output of reduced_age_analysis()
#' @return a vector of median retained coefficients across the five imputations for models fit on that loo model
retained_one_obs <- function(index, reduced_output){
  # quick check in case there is only 1 imputation
  number_of_imputations <- length(reduced_output[[1]][[1]]$fit)
  purrr::map(1:number_of_imputations, ~importance(reduced_output[[1]][[index]]$fit[[.x]]$fit, metabolites = F) %>%
               enframe(value = str_glue("imp{.x}"))) %>%
    reduce(full_join, by = "name") %>%
    rowwise() %>%
    transmute(name, median_coef = median(c_across(starts_with('imp')), na.rm = T)) %>%
    ungroup() %>%
    deframe()
}


#' look at the retained coefs for fits across each of n leave one out models
#' select the ones that show up > 1 - pct_mis% of the time.
#' @param reduced_output is the output of reduced_age_analysis
#' @param pct_mis -- only select variables that are present in 1 - pct_mis% of the fits
#' @param only_return_pct ignores pct_mis and just shows the percentage of missingness for each feature
retained_over_fits <- function(reduced_output, pct_mis = 0.05, only_return_pct = FALSE){
  
  if(only_return_pct){
    return(
      purrr::map(1:length(reduced_output[[1]]), ~retained_one_obs(.x, reduced_output)) %>%
        bind_rows() %>%
        map_dbl(~sum(is.na(.x) / length(.x))) %>%
        enframe(value = "pct_mis")
    )
  }
  
  purrr::map(1:length(reduced_output[[1]]), ~retained_one_obs(.x, reduced_output)) %>%
    bind_rows() %>%
    select_if(function(x) sum(is.na(x))/length(x) <= pct_mis) %>%
    map_dbl(~median(exp(.x), na.rm = T)) %>%
    enframe(value = "median_OR")
}


#' get lambda value for 5 imputations for a single observation
#' @param reduced_obj is an output of reduced_age_analysis()
#' @param index is 1:85 in the case of controls, (i.e. the row number of the left-out obs)
#' @return vector of 5 lambda values, corresponding to the 5 imputations for the given index
get_lambda_median <- function(reduced_obj, index){
  reduced_obj[[1]][[index]]$fit %>%
    purrr::map_dbl(~.x$fit$lambda.min)
}



targeted_reduced <- reduced_age_analysis(wide_data_targeted)
untargeted_reduced <- reduced_age_analysis(wide_data_untargeted)
got_reduced <- reduced_age_analysis(wide_data)
lipids_reduced <- reduced_age_analysis(wide_data_lipids)

### combined profile
wide_data_combined <- wide_data %>%
  left_join(select(wide_data_targeted, -any_of(setdiff(metadata_cols, 'Id'))), by = 'Id') %>%
  left_join(select(wide_data_untargeted, -any_of(setdiff(metadata_cols, 'Id'))), by = 'Id') %>%
  left_join(select(wide_data_lipids, -any_of(setdiff(metadata_cols, 'Id'))), by = 'Id')

combined_reduced <- reduced_age_analysis(wide_data_combined)

targeted_reduced <- readRDS(here('aging_output_files', 'targeted_reduced.Rds'))
untargeted_reduced <- readRDS(here('aging_output_files', 'untargeted_reduced.Rds'))
got_reduced <- readRDS(here('aging_output_files', 'got_reduced.Rds'))
lipids_reduced <- readRDS(here('aging_output_files', 'lipids_reduced.Rds'))
combined_reduced <- readRDS(here('aging_output_files', 'combined_reduced.Rds'))

saveRDS(targeted_reduced, file = here('aging_output_files', 'targeted_reduced.Rds'))
saveRDS(untargeted_reduced, file = here('aging_output_files', 'untargeted_reduced.Rds'))
saveRDS(lipids_reduced, file = here('aging_output_files', 'lipids_reduced.Rds'))
saveRDS(got_reduced, file = here('aging_output_files', 'got_reduced.Rds'))

saveRDS(combined_reduced, file = here('aging_output_files', 'combined_reduced.Rds'))
ggsave('age_clock_targeted_reduced_pred.png',
       plot = predtruth_plot(targeted_reduced$pred_df, name = "Targeted"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('age_clock_untargeted_reduced_pred.png',
       plot = predtruth_plot(untargeted_reduced$pred_df, name = "Untargeted"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('age_clock_GOT_reduced_pred.png',
       plot = predtruth_plot(got_reduced$pred_df, name = "GOT"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('age_clock_lipids_reduced_pred.png',
       plot = predtruth_plot(lipids_reduced$pred_df, name = "Lipids"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('age_clock_combined_reduced_pred.png',
       plot = predtruth_plot(combined_reduced$pred_df, name = "Combined"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)

figure_3 <- predtruth_plot(targeted_reduced$pred_df, name = "Targeted") /
  predtruth_plot(untargeted_reduced$pred_df, name = "Untargeted") /
  predtruth_plot(got_reduced$pred_df, name = "GOT") /
  predtruth_plot(lipids_reduced$pred_df, name = "Lipids")

ggsave('figure_3.tiff',
       device = 'tiff',
       type = 'cairo',
       plot = figure_3,
       path = here('plots', 'aging_figs', 'main_figs'),
       width = 83, height = 200, dpi = 300, units = 'mm')



targeted_retained <- retained_over_fits(targeted_reduced)  
untargeted_retained <- retained_over_fits(untargeted_reduced)  
got_retained <- retained_over_fits(got_reduced)
lipids_retained <- retained_over_fits(lipids_reduced)




targeted_lambdas <- map(1:85, ~get_lambda_median(targeted_reduced, .x)) %>%
  reduce(`c`)
lipids_lambdas <- map(1:85, ~get_lambda_median(lipids_reduced, .x)) %>%
  reduce(`c`)
untargeted_lambdas <- map(1:85, ~get_lambda_median(untargeted_reduced, .x)) %>%
  reduce(`c`)
got_lambdas <- map(1:85, ~get_lambda_median(got_reduced, .x)) %>%
  reduce(`c`)



### plot only matched controls ####

untargeted_adpd_matched_controls <- readRDS(here("aging_output_files", "untargeted_adpd_matched_controls.Rds"))
untargeted_c_matched_age_reduced <- untargeted_reduced$pred_df %>%
  filter(id %in% untargeted_adpd_matched_controls$Id)

ggsave(filename = "age_clock_untargeted_reduced_pred_matched.png",
       plot = predtruth_plot(untargeted_c_matched_age_reduced, name = "Untargeted (matched)"),
       width = 83,
       height = 100,
       path = here("plots", "aging_figs", 'main_figs'),
       units = 'mm',
       dpi = 300
       )


####### Untargeted Variations #################

####### rawscaled
untargeted_reduced_rawscaled <- reduced_age_analysis(wide_data_untargeted_rawscaled)
saveRDS(untargeted_reduced_rawscaled, file = here('aging_output_files', 'untargeted_reduced_rawscaled.Rds'))

####### models for each sex
wide_data_untargeted_male <- wide_data_untargeted %>%
  filter(Gender == 'M')
wide_data_untargeted_female <- wide_data_untargeted %>%
  filter(Gender == 'F')

untargeted_reduced_male <- reduced_age_analysis(wide_data_untargeted_male)
untargeted_reduced_female <- reduced_age_analysis(wide_data_untargeted_female)
saveRDS(untargeted_reduced_male, file = here('aging_output_files', 'untargeted_reduced_male.Rds'))
saveRDS(untargeted_reduced_female, file = here('aging_output_files', 'untargeted_reduced_female.Rds'))

# get random samples of same size. check 3 samples to get average
untargeted_reduced_n_male <- purrr::map(1:3, 
                                          ~reduced_age_analysis(
                                            wide_data_untargeted %>%
                                            filter(Type %in% c('CY', 'CM', 'CO')) %>%
                                            sample_n(44, replace = F)
                                            )
                                          )

untargeted_reduced_n_female <- purrr::map(1:3, 
                                          ~reduced_age_analysis(
                                            wide_data_untargeted %>%
                                              filter(Type %in% c('CY', 'CM', 'CO')) %>%
                                              sample_n(41, replace = F)
                                          )
)
saveRDS(untargeted_reduced_n_male, file = here('aging_output_files', 'untargeted_reduced_n_male.Rds'))
saveRDS(untargeted_reduced_n_female, file = here('aging_output_files', 'untargeted_reduced_n_female.Rds'))

cowplot::plot_grid(plotlist = purrr::map(wide_data_untargeted_n_male, ~.x$figure))
cowplot::plot_grid(plotlist = purrr::map(wide_data_untargeted_n_female, ~.x$figure))


####### post selection

targeted_reduced_post <- reduced_age_analysis(wide_data_targeted, post_select = T)
untargeted_reduced_post <- reduced_age_analysis(wide_data_untargeted, post_select = T)
got_reduced_post <- reduced_age_analysis(wide_data, post_select = T)
lipids_reduced_post <- reduced_age_analysis(wide_data_lipids, post_select = T)

targeted_reduced_post <- readRDS(file = here('aging_output_files', 'targeted_reduced_post.Rds'))
untargeted_reduced_post <- readRDS(file = here('aging_output_files', 'untargeted_reduced_post.Rds'))
lipids_reduced_post <- readRDS(file = here('aging_output_files', 'lipids_reduced_post.Rds'))
got_reduced_post <- readRDS(file = here('aging_output_files', 'got_reduced_post.Rds'))

saveRDS(targeted_reduced_post, file = here('aging_output_files', 'targeted_reduced_post.Rds'))
saveRDS(untargeted_reduced_post, file = here('aging_output_files', 'untargeted_reduced_post.Rds'))
saveRDS(lipids_reduced_post, file = here('aging_output_files', 'lipids_reduced_post.Rds'))
saveRDS(got_reduced_post, file = here('aging_output_files', 'got_reduced_post.Rds'))


ggsave('post_age_clock_targeted_reduced_pred.png',
       plot = predtruth_plot(targeted_reduced_post$pred_df, name = "Targeted Post"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('post_age_clock_untargeted_reduced_pred.png',
       plot = predtruth_plot(untargeted_reduced_post$pred_df, name = "Untargeted Post"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('post_age_clock_GOT_reduced_pred.png',
       plot = predtruth_plot(got_reduced_post$pred_df, name = "GOT Post"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('post_age_clock_lipids_reduced_pred.png',
       plot = predtruth_plot(lipids_reduced_post$pred_df, name = "Lipids Post"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)


##### cor > .3 (dog paper)
targeted_reduced_cor <- reduced_age_analysis(wide_data_targeted, remove_cor_age = T)
untargeted_reduced_cor <- reduced_age_analysis(wide_data_untargeted, remove_cor_age = T)
got_reduced_cor <- reduced_age_analysis(wide_data, remove_cor_age = T)
lipids_reduced_cor <- reduced_age_analysis(wide_data_lipids, remove_cor_age = T)

targeted_reduced_cor <- readRDS(file = here('aging_output_files', 'targeted_reduced_cor.Rds'))
untargeted_reduced_cor <- readRDS(file = here('aging_output_files', 'untargeted_reduced_cor.Rds'))
lipids_reduced_cor <- readRDS(file = here('aging_output_files', 'lipids_reduced_cor.Rds'))
got_reduced_cor <- readRDS(file = here('aging_output_files', 'got_reduced_cor.Rds'))

saveRDS(targeted_reduced_cor, file = here('aging_output_files', 'targeted_reduced_cor.Rds'))
saveRDS(untargeted_reduced_cor, file = here('aging_output_files', 'untargeted_reduced_cor.Rds'))
saveRDS(lipids_reduced_cor, file = here('aging_output_files', 'lipids_reduced_cor.Rds'))
saveRDS(got_reduced_cor, file = here('aging_output_files', 'got_reduced_cor.Rds'))


ggsave('cor_age_clock_targeted_reduced_pred.png',
       plot = predtruth_plot(targeted_reduced_cor$pred_df, name = "Targeted"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('cor_age_clock_untargeted_reduced_pred.png',
       plot = predtruth_plot(untargeted_reduced_cor$pred_df, name = "Untargeted"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('cor_age_clock_GOT_reduced_pred.png',
       plot = predtruth_plot(got_reduced_cor$pred_df, name = "GOT"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)
ggsave('cor_age_clock_lipids_reduced_pred.png',
       plot = predtruth_plot(lipids_reduced_cor$pred_df, name = "Lipids"),
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)


### taking the log of age first
targeted_reduced_takelog <- reduced_age_analysis(wide_data_targeted, remove_cor_age = T, take_log = T)
untargeted_reduced_takelog <- reduced_age_analysis(wide_data_untargeted, remove_cor_age = T, take_log = T)
got_reduced_takelog <- reduced_age_analysis(wide_data, remove_cor_age = T, take_log = T)
lipids_reduced_takelog <- reduced_age_analysis(wide_data_lipids, remove_cor_age = T, take_log = T)

saveRDS(targeted_reduced_takelog, file = here('aging_output_files', 'targeted_reduced_takelog.Rds'))
saveRDS(untargeted_reduced_takelog, file = here('aging_output_files', 'untargeted_reduced_takelog.Rds'))
saveRDS(lipids_reduced_takelog, file = here('aging_output_files', 'lipids_reduced_takelog.Rds'))
saveRDS(got_reduced_takelog, file = here('aging_output_files', 'got_reduced_takelog.Rds'))







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

targeted_ad_reduced <- readRDS(file = here('aging_output_files', 'targeted_ad_reduced.Rds'))
untargeted_ad_reduced <- readRDS(file = here('aging_output_files', 'untargeted_ad_reduced.Rds'))
lipids_ad_reduced <- readRDS(file = here('aging_output_files', 'lipids_ad_reduced.Rds'))
got_ad_reduced <- readRDS(file = here('aging_output_files', 'got_ad_reduced.Rds'))

# saveRDS(targeted_ad_reduced, file = here('aging_output_files', 'targeted_ad_reduced.Rds'))
# saveRDS(untargeted_ad_reduced, file = here('aging_output_files', 'untargeted_ad_reduced.Rds'))
# saveRDS(lipids_ad_reduced, file = here('aging_output_files', 'lipids_ad_reduced.Rds'))
# saveRDS(got_ad_reduced, file = here('aging_output_files', 'got_ad_reduced.Rds'))



untargeted_ad_c_roc <- roc_fig(untargeted_ad_reduced$roc_df, untargeted_ad_reduced$predtruth_df) + 
  labs(title = 'Untargeted')
targeted_ad_c_roc <- roc_fig(targeted_ad_reduced$roc_df, targeted_ad_reduced$predtruth_df) + 
  labs(title = 'Targeted')
lipids_ad_c_roc <- roc_fig(lipids_ad_reduced$roc_df, lipids_ad_reduced$predtruth_df) + 
  labs(title = 'Lipids')
got_ad_c_roc <- roc_fig(got_ad_reduced$roc_df, got_ad_reduced$predtruth_df) + 
  labs(title = 'GOT')


# get null model from adpd_glmnet.R
meta_ad_c_auc <- round(base_ad_c_roc$auc[1], 3)

roc_ad_reduced_grid <- untargeted_ad_c_roc + targeted_ad_c_roc + lipids_ad_c_roc + 
  plot_annotation(title = "AD vs C",theme = theme(title = element_text(size = rel(1.5)))#,
                  #subtitle = str_glue("(Age, sex only AUC: {meta_ad_c_auc})")) +
  ) +
 plot_layout(ncol = 3, nrow = 1)

ggsave('roc_ad_reduced_grid.png',
       plot = roc_ad_reduced_grid,
       path = here("plots", "adpd_figs"), 
       width = 16, height = 4)


ggsave('roc_ad_reduced_untargeted.png',
       plot = roc_fig(untargeted_ad_reduced$roc_df, untargeted_ad_reduced$predtruth_df) + 
         labs(title = 'Untargeted: AD vs C'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)
ggsave('roc_ad_reduced_targeted.png',
       plot = roc_fig(targeted_ad_reduced$roc_df, targeted_ad_reduced$predtruth_df) + 
         labs(title = 'Targeted: AD vs C'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)
ggsave('roc_ad_reduced_lipids.png',
       plot = roc_fig(lipids_ad_reduced$roc_df, lipids_ad_reduced$predtruth_df) + 
         labs(title = 'Lipids: AD vs C'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)
ggsave('roc_ad_reduced_got.png',
       plot = roc_fig(got_ad_reduced$roc_df, got_ad_reduced$predtruth_df) + 
         labs(title = 'GOT: AD vs C'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)

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
# saveRDS(targeted_pd_reduced, file = here('aging_output_files', 'targeted_pd_reduced.Rds'))
# saveRDS(untargeted_pd_reduced, file = here('aging_output_files', 'untargeted_pd_reduced.Rds'))
# saveRDS(lipids_pd_reduced, file = here('aging_output_files', 'lipids_pd_reduced.Rds'))
# saveRDS(got_pd_reduced, file = here('aging_output_files', 'got_pd_reduced.Rds'))


targeted_pd_reduced <- readRDS(file = here('aging_output_files', 'targeted_pd_reduced.Rds'))
untargeted_pd_reduced <- readRDS(file = here('aging_output_files', 'untargeted_pd_reduced.Rds'))
lipids_pd_reduced <- readRDS(file = here('aging_output_files', 'lipids_pd_reduced.Rds'))
got_pd_reduced <- readRDS(file = here('aging_output_files', 'got_pd_reduced.Rds'))



untargeted_pd_c_roc <- roc_fig(untargeted_pd_reduced$roc_df, untargeted_pd_reduced$predtruth_df) + 
  labs(title = 'Untargeted')
targeted_pd_c_roc <- roc_fig(targeted_pd_reduced$roc_df, targeted_pd_reduced$predtruth_df) + 
  labs(title = 'Targeted')
lipids_pd_c_roc <- roc_fig(lipids_pd_reduced$roc_df, lipids_pd_reduced$predtruth_df) + 
  labs(title = 'Lipids')
got_pd_c_roc <- roc_fig(got_pd_reduced$roc_df, got_pd_reduced$predtruth_df) + 
  labs(title = 'GOT')

# get null model from adpd_glmnet.R (base_ad_c...)
meta_pd_c_auc <- round(base_pd_c_roc$auc[1], 3)

roc_pd_reduced_grid <- untargeted_pd_c_roc + targeted_pd_c_roc + lipids_pd_c_roc + 
  plot_annotation(title = "PD vs C", theme = theme(title = element_text(size = rel(1.5)))#,
                  #subtitle = str_glue("(Age, sex only: {meta_pd_c_auc})")) +
  )+
  plot_layout(ncol = 3, nrow = 1)

ggsave('roc_pd_reduced_grid.png',
       plot = roc_pd_reduced_grid,
       path = here("plots", "adpd_figs"), 
       width = 16, height = 4)


#### Untargeted PD has almost perfect prediction... remove all metabolites retained in the full model and see what happens

untargeted_pd_full <- readRDS(file = here('ad_pd', 'untargeted_pd_full.Rds'))


untar_pd_full_retained <- get_importance_tables(untargeted_pd_full) %>%
  filter(!(Name %in% c('Age', 'GenderM', '(Intercept)'))) %>%
  pull(Name)

wide_data_untargeted_pd_rd2 <- wide_data_untargeted %>%
  select(-all_of(contains(untar_pd_full_retained)))

untargeted_pd_reduced_rd2 <- readRDS(file = here('ad_pd', 'untargeted_pd_reduced_rd2.Rds'))

untargeted_pd_reduced_rd2 <- reduced_logistic_analysis(data = wide_data_untargeted_pd_rd2, 
                                                   target ="PD_ind", PD_ind = T, 
                                                   types = c("PD", "CO", "CM", "CY")
)
saveRDS(untargeted_pd_reduced_rd2, file = here('ad_pd', 'untargeted_pd_reduced_rd2.Rds'))

ggsave('roc_pd_untargeted_rd2.png',
       plot = roc_fig(untargeted_pd_reduced_rd2$roc_df) + 
         labs(title = 'Untargeted: PD vs Controls',
              subtitle = 'Round 2'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)

##### PD vs AD
targeted_adpd_reduced <- readRDS(file = here('aging_output_files', 'targeted_adpd_reduced.Rds'))
untargeted_adpd_reduced <- readRDS(file = here('aging_output_files', 'untargeted_adpd_reduced.Rds'))
lipids_adpd_reduced <- readRDS(file = here('aging_output_files', 'lipids_adpd_reduced.Rds'))
got_adpd_reduced <- readRDS(file = here('aging_output_files', 'got_adpd_reduced.Rds'))

targeted_adpd_reduced <- reduced_logistic_analysis(data = wide_data_targeted, 
                                                 target ="AD_ind", AD_ind = T, 
                                                 types = c("AD", "PD")
)
got_adpd_reduced <- reduced_logistic_analysis(data = wide_data, 
                                            target ="AD_ind", AD_ind = T, 
                                            types = c("AD", "PD")
)
untargeted_adpd_reduced <- reduced_logistic_analysis(data = wide_data_untargeted, 
                                                   target ="AD_ind", AD_ind = T, 
                                                   types = c("AD", "PD")
)

lipids_adpd_reduced <- reduced_logistic_analysis(data = wide_data_lipids, 
                                               target ="AD_ind", AD_ind = T, 
                                               types = c("AD", "PD")
)

# saveRDS(targeted_adpd_reduced, file = here('aging_output_files', 'targeted_adpd_reduced.Rds'))
# saveRDS(untargeted_adpd_reduced, file = here('aging_output_files', 'untargeted_adpd_reduced.Rds'))
# saveRDS(lipids_adpd_reduced, file = here('aging_output_files', 'lipids_adpd_reduced.Rds'))
# saveRDS(got_adpd_reduced, file = here('aging_output_files', 'got_adpd_reduced.Rds'))
 

untargeted_adpd_roc <- roc_fig(untargeted_adpd_reduced$roc_df, untargeted_adpd_reduced$predtruth_df) + 
  labs(title = 'Untargeted')
targeted_adpd_roc <- roc_fig(targeted_adpd_reduced$roc_df, targeted_adpd_reduced$predtruth_df) + 
  labs(title = 'Targeted')
lipids_adpd_roc <- roc_fig(lipids_adpd_reduced$roc_df, lipids_adpd_reduced$predtruth_df) + 
  labs(title = 'Lipids')
got_adpd_roc <- roc_fig(got_adpd_reduced$roc_df, got_adpd_reduced$predtruth_df) + 
  labs(title = 'GOT')


# ggsave('roc_adpd_untargeted.png',
#        plot = roc_fig(untargeted_adpd_reduced$roc_df) + 
#          labs(title = 'Untargeted: AD vs PD'),
#        path = here("plots", "adpd_figs"), 
#        width = 14, height = 10)
# ggsave('roc_adpd_targeted.png',
#        plot = roc_fig(targeted_adpd_reduced$roc_df) + 
#          labs(title = 'Targeted: AD vs PD'),
#        path = here("plots", "adpd_figs"), 
#        width = 14, height = 10)
# ggsave('roc_adpd_lipids.png',
#        plot = roc_fig(lipids_adpd_reduced$roc_df) + 
#          labs(title = 'Lipids: AD vs PD'),
#        path = here("plots", "adpd_figs"), 
#        width = 14, height = 10)
# ggsave('roc_adpd_got.png',
#        plot = roc_fig(got_adpd_reduced$roc_df) + 
#          labs(title = 'GOT: AD vs PD'),
#        path = here("plots", "adpd_figs"), 
#        width = 14, height = 10)


# get null model
meta_adpd_auc <- round(base_ad_pd_roc$auc[1], 3)
roc_adpd_reduced_grid <- untargeted_adpd_roc + targeted_adpd_roc + lipids_adpd_roc + 
  plot_annotation(title = "AD vs PD", theme = theme(title = element_text(size = rel(1.5)))#,
                  #subtitle = str_glue("(Age, sex only: {meta_adpd_auc})")) +
  )+
  plot_layout(ncol = 3, nrow = 1)

ggsave('roc_adpd_reduced_grid.png',
       plot = roc_adpd_reduced_grid,
       path = here("plots", "adpd_figs"), 
       width = 16, height = 4)




roc_biggrid_w_agesex <-  
  untargeted_ad_c_roc + targeted_ad_c_roc + lipids_ad_c_roc + 
  untargeted_pd_c_roc + targeted_pd_c_roc + lipids_pd_c_roc +
  untargeted_adpd_roc + targeted_adpd_roc + lipids_adpd_roc +
  plot_layout(ncol = 4, nrow= 3)


#### APOE

targeted_apoe_reduced <- reduced_logistic_analysis(data = wide_data_targeted, 
                                                 target ="APOE4_ind", APOE4_ind = T, 
                                                 types = c("AD", 'PD', "CO", "CM", "CY")
)
got_apoe_reduced <- reduced_logistic_analysis(data = wide_data, 
                                            target ="APOE4_ind", APOE4_ind = T, 
                                            types = c("AD", "PD", "CO", "CM", "CY")
)
# note high empri0
untargeted_apoe_reduced <- reduced_logistic_analysis(data = wide_data_untargeted, 
                                                   target ="APOE4_ind", APOE4_ind = T, 
                                                   types = c("AD", "PD", "CO", "CM", "CY"),
                                                   empri0 = 200
)

lipids_apoe_reduced <- reduced_logistic_analysis(data = wide_data_lipids, 
                                               target ="APOE4_ind", APOE4_ind = T, 
                                               types = c("AD", "PD", "CO", "CM", "CY")
)

saveRDS(targeted_apoe_reduced, file = here('aging_output_files', 'targeted_apoe_reduced.Rds'))
saveRDS(untargeted_apoe_reduced, file = here('aging_output_files', 'untargeted_apoe_reduced.Rds'))
saveRDS(lipids_apoe_reduced, file = here('aging_output_files', 'lipids_apoe_reduced.Rds'))
saveRDS(got_apoe_reduced, file = here('aging_output_files', 'got_apoe_reduced.Rds'))


### GBA

targeted_GBA_reduced <- reduced_logistic_analysis(data = wide_data_targeted, 
                                                   target ="GBA_ind", GBA_ind = T, 
                                                   types = c('PD')
)
got_GBA_reduced <- reduced_logistic_analysis(data = wide_data, 
                                              target ="GBA_ind", GBA_ind = T, 
                                              types = c("PD")
)
untargeted_GBA_reduced <- reduced_logistic_analysis(data = wide_data_untargeted, 
                                                     target ="GBA_ind", GBA_ind = T, 
                                                     types = c("PD")
)

lipids_GBA_reduced <- reduced_logistic_analysis(data = wide_data_lipids, 
                                                 target ="GBA_ind", GBA_ind = T, 
                                                 types = c("PD")
)

saveRDS(targeted_GBA_reduced, file = here('aging_output_files', 'targeted_GBA_reduced.Rds'))
saveRDS(untargeted_GBA_reduced, file = here('aging_output_files', 'untargeted_GBA_reduced.Rds'))
saveRDS(lipids_GBA_reduced, file = here('aging_output_files', 'lipids_GBA_reduced.Rds'))
saveRDS(got_GBA_reduced, file = here('aging_output_files', 'got_GBA_reduced.Rds'))


### LED

targeted_led_reduced <- reduced_age_analysis(wide_data_targeted_panuc, 
                                             target = 'led',
                                             types = c('PD'),
                                             include_led = T)
untargeted_led_reduced <- reduced_age_analysis(wide_data_untargeted_panuc, 
                                               target = 'led',
                                               types = c('PD'),
                                               include_led = T)
got_led_reduced <- reduced_age_analysis(wide_data_panuc, 
                                        target = 'led',
                                        types = c('PD'),
                                        include_led = T)
lipids_led_reduced <- reduced_age_analysis(wide_data_lipids_panuc, 
                                           target = 'led',
                                           types = c('PD'),
                                           include_led = T)
saveRDS(targeted_led_reduced, file = here('aging_output_files', 'targeted_led_reduced.Rds'))
saveRDS(untargeted_led_reduced, file = here('aging_output_files', 'untargeted_led_reduced.Rds'))
saveRDS(lipids_led_reduced, file = here('aging_output_files', 'lipids_led_reduced.Rds'))
saveRDS(got_led_reduced, file = here('aging_output_files', 'got_led_reduced.Rds'))


####### change everything to missingnes indicators



targeted_ad_reduced_naind <- readRDS(file = here('aging_output_files', 'targeted_ad_reduced_naind.Rds'))
untargeted_ad_reduced_naind <- readRDS(file = here('aging_output_files', 'untargeted_ad_reduced_naind.Rds'))
lipids_ad_reduced_naind <- readRDS(file = here('aging_output_files', 'lipids_ad_reduced_naind.Rds'))
got_ad_reduced_naind <- readRDS(file = here('aging_output_files', 'got_ad_reduced_naind.Rds'))


targeted_ad_reduced_naind <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
                                                 target ="AD_ind", AD_ind = T, 
                                                 types = c("AD", "CO", "CM", "CY")
)
got_ad_reduced_naind <- reduced_logistic_analysis(data = wide_data_got_naind, 
                                            target ="AD_ind", AD_ind = T, 
                                            types = c("AD", "CO", "CM", "CY")
)
untargeted_ad_reduced_naind <- reduced_logistic_analysis(data = wide_data_untargeted_naind, 
                                                   target ="AD_ind", AD_ind = T, 
                                                   types = c("AD", "CO", "CM", "CY")
)

lipids_ad_reduced_naind <- reduced_logistic_analysis(data = wide_data_lipids_naind, 
                                               target ="AD_ind", AD_ind = T, 
                                               types = c("AD", "CO", "CM", "CY")
)

saveRDS(targeted_ad_reduced_naind, file = here('aging_output_files', 'targeted_ad_reduced_naind.Rds'))
saveRDS(untargeted_ad_reduced_naind, file = here('aging_output_files', 'untargeted_ad_reduced_naind.Rds'))
saveRDS(lipids_ad_reduced_naind, file = here('aging_output_files', 'lipids_ad_reduced_naind.Rds'))
saveRDS(got_ad_reduced_naind, file = here('aging_output_files', 'got_ad_reduced_naind.Rds'))
####

targeted_pd_reduced_naind <- readRDS(file = here('aging_output_files', 'targeted_pd_reduced_naind.Rds'))
untargeted_pd_reduced_naind <- readRDS(file = here('aging_output_files', 'untargeted_pd_reduced_naind.Rds'))
lipids_pd_reduced_naind <- readRDS(file = here('aging_output_files', 'lipids_pd_reduced_naind.Rds'))
got_pd_reduced_naind <- readRDS(file = here('aging_output_files', 'got_pd_reduced_naind.Rds'))


targeted_pd_reduced_naind <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
                                                 target ="PD_ind", PD_ind = T, 
                                                 types = c("PD", "CO", "CM", "CY")
)
got_pd_reduced_naind <- reduced_logistic_analysis(data = wide_data_got_naind, 
                                            target ="PD_ind", PD_ind = T, 
                                            types = c("PD", "CO", "CM", "CY")
)



untargeted_pd_reduced_naind <- reduced_logistic_analysis(data = wide_data_untargeted_naind, 
                                                   target ="PD_ind", PD_ind = T, 
                                                   types = c("PD", "CO", "CM", "CY")
)

retained_over_fits(untargeted_pd_reduced_naind)

lipids_pd_reduced_naind <- reduced_logistic_analysis(data = wide_data_lipids_naind, 
                                               target ="PD_ind", PD_ind = T, 
                                               types = c("PD", "CO", "CM", "CY")
)
saveRDS(targeted_pd_reduced_naind, file = here('aging_output_files', 'targeted_pd_reduced_naind.Rds'))
saveRDS(untargeted_pd_reduced_naind, file = here('aging_output_files', 'untargeted_pd_reduced_naind.Rds'))
saveRDS(lipids_pd_reduced_naind, file = here('aging_output_files', 'lipids_pd_reduced_naind.Rds'))
saveRDS(got_pd_reduced_naind, file = here('aging_output_files', 'got_pd_reduced_naind.Rds'))


# near perfect prediction with untargeted pd naindicators, so let's try again with round 2

untargeted_pd_full_naind <- readRDS(file = here('ad_pd', 'untargeted_pd_full_naind.Rds'))


untar_pd_full_retained_naind <- get_importance_tables(untargeted_pd_full_naind) %>%
  filter(!(Name %in% c('Age', 'GenderM', '(Intercept)'))) %>%
  pull(Name)

wide_data_untargeted_pd_rd2_naind <- wide_data_untargeted_naind %>%
  select(-all_of(contains(untar_pd_full_retained_naind)))

untargeted_pd_reduced_rd2_naind <- readRDS(file = here('ad_pd', 'untargeted_pd_reduced_rd2_naind.Rds'))

untargeted_pd_reduced_rd2_naind <- reduced_logistic_analysis(data = wide_data_untargeted_pd_rd2_naind, 
                                                       target ="PD_ind", PD_ind = T, 
                                                       types = c("PD", "CO", "CM", "CY")
)
saveRDS(untargeted_pd_reduced_rd2_naind, file = here('ad_pd', 'untargeted_pd_reduced_rd2_naind.Rds'))

ggsave('roc_pd_untargeted_rd2_naind.png',
       plot = roc_fig(untargeted_pd_reduced_rd2_naind$roc_df) + 
         labs(title = 'Untargeted: PD vs Controls',
              subtitle = 'Missingness Indicators, Round 2'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)


#################################

#  possible replacement of above -- penalize age/sex in na indicators
# idea: the reason we don't penalize is because we don't want the model to pick up on age/sex effects in the metabolome
# but if we dont' assume that age/sex are associated with missingness, then removing them will give a better sense of how models perform using metabolites only
targeted_ad_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'targeted_ad_reduced_naind_nopenalize.Rds'))
untargeted_ad_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'untargeted_ad_reduced_naind_nopenalize.Rds'))
lipids_ad_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'lipids_ad_reduced_naind_nopenalize.Rds'))
got_ad_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'got_ad_reduced_naind_nopenalize.Rds'))


targeted_ad_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
                                                       target ="AD_ind", AD_ind = T, 
                                                       types = c("AD", "CO", "CM", "CY"),
                                                       penalize_age_gender = T
)
got_ad_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_got_naind, 
                                                  target ="AD_ind", AD_ind = T, 
                                                  types = c("AD", "CO", "CM", "CY"),
                                                  penalize_age_gender = T
)
untargeted_ad_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_untargeted_naind, 
                                                         target ="AD_ind", AD_ind = T, 
                                                         types = c("AD", "CO", "CM", "CY"),
                                                         penalize_age_gender = T
)

lipids_ad_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_lipids_naind, 
                                                     target ="AD_ind", AD_ind = T, 
                                                     types = c("AD", "CO", "CM", "CY"),
                                                     penalize_age_gender = T
)

saveRDS(targeted_ad_reduced_naind_nopenalize, file = here('aging_output_files', 'targeted_ad_reduced_naind_nopenalize.Rds'))
saveRDS(untargeted_ad_reduced_naind_nopenalize, file = here('aging_output_files', 'untargeted_ad_reduced_naind_nopenalize.Rds'))
saveRDS(lipids_ad_reduced_naind_nopenalize, file = here('aging_output_files', 'lipids_ad_reduced_naind_nopenalize.Rds'))
saveRDS(got_ad_reduced_naind_nopenalize, file = here('aging_output_files', 'got_ad_reduced_naind_nopenalize.Rds'))
####

targeted_pd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'targeted_pd_reduced_naind_nopenalize.Rds'))
untargeted_pd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'untargeted_pd_reduced_naind_nopenalize.Rds'))
lipids_pd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'lipids_pd_reduced_naind_nopenalize.Rds'))
got_pd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'got_pd_reduced_naind_nopenalize.Rds'))


targeted_pd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
                                                       target ="PD_ind", PD_ind = T, 
                                                       types = c("PD", "CO", "CM", "CY"),
                                                       penalize_age_gender = T
)
got_pd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_got_naind, 
                                                  target ="PD_ind", PD_ind = T, 
                                                  types = c("PD", "CO", "CM", "CY"),
                                                  penalize_age_gender = T
)



untargeted_pd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_untargeted_naind, 
                                                         target ="PD_ind", PD_ind = T, 
                                                         types = c("PD", "CO", "CM", "CY"),
                                                         penalize_age_gender = T
)

#retained_over_fits(untargeted_pd_reduced_naind_nopenalize)

lipids_pd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_lipids_naind, 
                                                     target ="PD_ind", PD_ind = T, 
                                                     types = c("PD", "CO", "CM", "CY"),
                                                     penalize_age_gender = T
)
saveRDS(targeted_pd_reduced_naind_nopenalize, file = here('aging_output_files', 'targeted_pd_reduced_naind_nopenalize.Rds'))
saveRDS(untargeted_pd_reduced_naind_nopenalize, file = here('aging_output_files', 'untargeted_pd_reduced_naind_nopenalize.Rds'))
saveRDS(lipids_pd_reduced_naind_nopenalize, file = here('aging_output_files', 'lipids_pd_reduced_naind_nopenalize.Rds'))
saveRDS(got_pd_reduced_naind_nopenalize, file = here('aging_output_files', 'got_pd_reduced_naind_nopenalize.Rds'))

#####

targeted_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'targeted_adpd_reduced_naind_nopenalize.Rds'))
untargeted_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'untargeted_adpd_reduced_naind_nopenalize.Rds'))
lipids_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'lipids_adpd_reduced_naind_nopenalize.Rds'))
got_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'got_adpd_reduced_naind_nopenalize.Rds'))


targeted_adpd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
                                                                  target ="AD_ind", AD_ind = T, 
                                                                  types = c("PD", "AD"),
                                                                  penalize_age_gender = T
)
got_adpd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_got_naind, 
                                                             target ="AD_ind", AD_ind = T, 
                                                             types = c("PD", "AD"),
                                                             penalize_age_gender = T
)



untargeted_adpd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_untargeted_naind, 
                                                                    target ="AD_ind", AD_ind = T, 
                                                                    types = c("PD", "AD"),
                                                                    penalize_age_gender = T
)

#retained_over_fits(untargeted_adpd_reduced_naind_nopenalize)

lipids_adpd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_lipids_naind, 
                                                                target ="AD_ind", AD_ind = T, 
                                                                types = c("PD", "AD"),
                                                                penalize_age_gender = T
)
saveRDS(targeted_adpd_reduced_naind_nopenalize, file = here('aging_output_files', 'targeted_adpd_reduced_naind_nopenalize.Rds'))
saveRDS(untargeted_adpd_reduced_naind_nopenalize, file = here('aging_output_files', 'untargeted_adpd_reduced_naind_nopenalize.Rds'))
saveRDS(lipids_adpd_reduced_naind_nopenalize, file = here('aging_output_files', 'lipids_adpd_reduced_naind_nopenalize.Rds'))
saveRDS(got_adpd_reduced_naind_nopenalize, file = here('aging_output_files', 'got_adpd_reduced_naind_nopenalize.Rds'))





########

# Detrend sex and age to get AUC without effect of age / ex

########




### AD vs c------------------
targeted_ad_reduced_detrend <- reduced_logistic_analysis(data = wide_data_targeted_detrend, 
                                                 target ="AD_ind", AD_ind = T, 
                                                 types = c("AD", "CO", "CM", "CY"),
                                                 include_age_gender = F
)
got_ad_reduced_detrend <- reduced_logistic_analysis(data = wide_data_detrend, 
                                            target ="AD_ind", AD_ind = T, 
                                            types = c("AD", "CO", "CM", "CY"),
                                            include_age_gender = F
)
untargeted_ad_reduced_detrend <- reduced_logistic_analysis(data = wide_data_untargeted_detrend, 
                                                   target ="AD_ind", AD_ind = T, 
                                                   types = c("AD", "CO", "CM", "CY"),
                                                   include_age_gender = F
)

lipids_ad_reduced_detrend <- reduced_logistic_analysis(data = wide_data_lipids_detrend, 
                                               target ="AD_ind", AD_ind = T, 
                                               types = c("AD", "CO", "CM", "CY"),
                                               include_age_gender = F
)


targeted_ad_c_roc_detrend <- roc_fig(targeted_ad_reduced_detrend$roc_df, targeted_ad_reduced_detrend$predtruth_df) + 
  labs(title = 'Targeted')
lipids_ad_c_roc_detrend <- roc_fig(lipids_ad_reduced_detrend$roc_df, lipids_ad_reduced_detrend$predtruth_df) + 
  labs(title = 'Lipids')
untargeted_ad_c_roc_detrend <- roc_fig(untargeted_ad_reduced_detrend$roc_df, untargeted_ad_reduced_detrend$predtruth_df) + 
  labs(title = 'AD vs Controls')

targeted_ad_c_roc_detrend + lipids_ad_c_roc_detrend + untargeted_ad_c_roc_detrend +
  plot_annotation(title = "AD vs Controls")

saveRDS(targeted_ad_reduced_detrend, file = here('aging_output_files', 'targeted_ad_reduced_detrend.Rds'))
saveRDS(untargeted_ad_reduced_detrend, file = here('aging_output_files', 'untargeted_ad_reduced_detrend.Rds'))
saveRDS(lipids_ad_reduced_detrend, file = here('aging_output_files', 'lipids_ad_reduced_detrend.Rds'))
saveRDS(got_ad_reduced_detrend, file = here('aging_output_files', 'got_ad_reduced_detrend.Rds'))


targeted_ad_reduced_detrend <- readRDS(file = here('aging_output_files', 'targeted_ad_reduced_detrend.Rds'))
untargeted_ad_reduced_detrend <- readRDS(file = here('aging_output_files', 'untargeted_ad_reduced_detrend.Rds'))
lipids_ad_reduced_detrend <- readRDS(file = here('aging_output_files', 'lipids_ad_reduced_detrend.Rds'))


### PD vs c ------------------

targeted_pd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_targeted_detrend, 
                                                 target ="PD_ind", PD_ind = T, 
                                                 types = c("PD", "CO", "CM", "CY"),
                                                 include_age_gender = F
)
got_pd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_detrend, 
                                            target ="PD_ind", PD_ind = T, 
                                            types = c("PD", "CO", "CM", "CY"),
                                            include_age_gender = F
)


untargeted_pd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_untargeted_detrend, 
                                                   target ="PD_ind", PD_ind = T, 
                                                   types = c("PD", "CO", "CM", "CY"),
                                                   include_age_gender = F
)
lipids_pd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_lipids_detrend, 
                                               target ="PD_ind", PD_ind = T, 
                                               types = c("PD", "CO", "CM", "CY"),
                                               include_age_gender = F
)

targeted_pd_c_roc_detrend <- roc_fig(targeted_pd_reduced_detrend$roc_df, targeted_pd_reduced_detrend$predtruth_df) + 
  labs(title = 'Targeted')
lipids_pd_c_roc_detrend <- roc_fig(lipids_pd_reduced_detrend$roc_df, lipids_pd_reduced_detrend$predtruth_df) + 
  labs(title = 'Lipids')
untargeted_pd_c_roc_detrend <- roc_fig(untargeted_pd_reduced_detrend$roc_df, untargeted_pd_reduced_detrend$predtruth_df) + 
  labs(title = 'PD vs Controls')

targeted_pd_c_roc_detrend + lipids_pd_c_roc_detrend + untargeted_pd_c_roc_detrend

# 
# saveRDS(targeted_pd_reduced_detrend, file = here('aging_output_files', 'targeted_pd_reduced_detrend.Rds'))
# saveRDS(untargeted_pd_reduced_detrend, file = here('aging_output_files', 'untargeted_pd_reduced_detrend.Rds'))
# saveRDS(lipids_pd_reduced_detrend, file = here('aging_output_files', 'lipids_pd_reduced_detrend.Rds'))
# saveRDS(got_pd_reduced_detrend, file = here('aging_output_files', 'got_pd_reduced_detrend.Rds'))

targeted_pd_reduced_detrend <- readRDS(file = here('aging_output_files', 'targeted_pd_reduced_detrend.Rds'))
untargeted_pd_reduced_detrend <- readRDS(file = here('aging_output_files', 'untargeted_pd_reduced_detrend.Rds'))
lipids_pd_reduced_detrend <- readRDS(file = here('aging_output_files', 'lipids_pd_reduced_detrend.Rds'))

lipids_pd_reduced_detrend$predtruth_df %>% 
  inner_join(panuc_data, by = c("id" = "subject_id")) %$% 
  cor(disease_duration_onset, imp_avg, method = "spearman")

### AD vs PD ---------------------

targeted_adpd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_targeted_detrend, 
                                                   target ="AD_ind", AD_ind = T, 
                                                   types = c("AD", "PD"),
                                                   include_age_gender = F
)
got_adpd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_detrend, 
                                              target ="AD_ind", AD_ind = T, 
                                              types = c("AD", "PD"),
                                              include_age_gender = F
)
untargeted_adpd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_untargeted_detrend, 
                                                     target ="AD_ind", AD_ind = T, 
                                                     types = c("AD", "PD"),
                                                     include_age_gender = F
)

lipids_adpd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_lipids_detrend, 
                                                 target ="AD_ind", AD_ind = T, 
                                                 types = c("AD", "PD"),
                                                 include_age_gender = F
)

targeted_adpd_c_roc_detrend <- roc_fig(targeted_adpd_reduced_detrend$roc_df, targeted_adpd_reduced_detrend$predtruth_df) + 
  labs(title = 'Targeted')
lipids_adpd_c_roc_detrend <- roc_fig(lipids_adpd_reduced_detrend$roc_df, lipids_adpd_reduced_detrend$predtruth_df) + 
  labs(title = 'Lipids')
untargeted_adpd_c_roc_detrend <- roc_fig(untargeted_adpd_reduced_detrend$roc_df, untargeted_adpd_reduced_detrend$predtruth_df) + 
  labs(title = 'AD vs PD')

targeted_adpd_c_roc_detrend + lipids_adpd_c_roc_detrend + untargeted_adpd_c_roc_detrend

# saveRDS(targeted_adpd_reduced_detrend, file = here('aging_output_files', 'targeted_adpd_reduced_detrend.Rds'))
# saveRDS(untargeted_adpd_reduced_detrend, file = here('aging_output_files', 'untargeted_adpd_reduced_detrend.Rds'))
# saveRDS(lipids_adpd_reduced_detrend, file = here('aging_output_files', 'lipids_adpd_reduced_detrend.Rds'))
# saveRDS(got_adpd_reduced_detrend, file = here('aging_output_files', 'got_adpd_reduced_detrend.Rds'))

targeted_adpd_reduced_detrend <- readRDS(file = here('aging_output_files', 'targeted_adpd_reduced_detrend.Rds'))
untargeted_adpd_reduced_detrend <- readRDS(file = here('aging_output_files', 'untargeted_adpd_reduced_detrend.Rds'))
lipids_adpd_reduced_detrend <- readRDS(file = here('aging_output_files', 'lipids_adpd_reduced_detrend.Rds'))


figure_3 <- untargeted_ad_c_roc_detrend + untargeted_pd_c_roc_detrend + untargeted_adpd_c_roc_detrend +
  plot_annotation(title = "ROC curves using the untargeted profile", theme = theme(title = element_text(size = rel(2))))

ggsave("figure_3.tiff",
       type = 'cairo',
       device = 'tiff',
       plot = figure_3,
       path = here("plots", "adpd_figs", "main_figs"),
       width = 83*3,
       height = 100,
       units = 'mm', dpi = 300)





### which metabolites are retained in EVERY SINGLE model (n models * 5 imputations)

# PD
tar_pd_detrend_retained <- retained_over_fits(targeted_pd_reduced_detrend, pct_mis= 0)
untar_pd_detrend_retained <- retained_over_fits(untargeted_pd_reduced_detrend, pct_mis= 0)
lipids_pd_detrend_retained <- retained_over_fits(lipids_pd_reduced_detrend, pct_mis= 0)

# untar didn't have any retained.. check the distn of pct_retained
untar_pd_detrend_retained <- retained_over_fits(untargeted_pd_reduced_detrend, only_return_pct = T)

(pd_detrend_retained_summary <- bind_rows(
  Targeted = tar_pd_detrend_retained,
  Lipids = lipids_pd_detrend_retained,
  .id = "profile"
) %>%
  filter(name != "(Intercept)") %>%
  group_by(profile) %>%
  arrange(desc(abs(log(median_OR)))) %>%
  ungroup() %>%
  knitr::kable()
)

# AD
tar_ad_detrend_retained <- retained_over_fits(targeted_ad_reduced_detrend, pct_mis= 0)
untar_ad_detrend_retained <- retained_over_fits(untargeted_ad_reduced_detrend, pct_mis= 0)
lipids_ad_detrend_retained <- retained_over_fits(lipids_ad_reduced_detrend, pct_mis= 0)


(ad_detrend_retained_summary <- bind_rows(
  Targeted = tar_ad_detrend_retained,
  Lipids = lipids_ad_detrend_retained,
  .id = "profile"
) %>%
    filter(name != "(Intercept)") %>%
    group_by(profile) %>%
    arrange(desc(abs(log(median_OR)))) %>%
    ungroup() %>%
    knitr::kable()
)

# AD v PD
tar_adpd_detrend_retained <- retained_over_fits(targeted_adpd_reduced_detrend, pct_mis= 0)
untar_adpd_detrend_retained <- retained_over_fits(untargeted_adpd_reduced_detrend, pct_mis= 0)
lipids_adpd_detrend_retained <- retained_over_fits(lipids_adpd_reduced_detrend, pct_mis= 0)

(adpd_detrend_retained_summary <- bind_rows(
  Targeted = tar_adpd_detrend_retained,
  Lipids = lipids_adpd_detrend_retained,
  .id = "profile"
) %>%
    filter(name != "(Intercept)") %>%
    group_by(profile) %>%
    arrange(desc(abs(log(median_OR)))) %>%
    ungroup() %>%  
    knitr::kable()
)








#### adpd_detrend, e.g. removing signal between AD and PD before fitting on controls

### AD
targeted_ad_reduced_detrend_orthog <- reduced_logistic_analysis(data = wide_data_targeted_detrend, 
                                                         target ="AD_ind", AD_ind = T, 
                                                         types = c("AD", "CO", "CM", "CY"),
                                                         include_age_gender = F, adpd_detrend = T
)


# saveRDS(targeted_ad_reduced_detrend_orthog, file = here('aging_output_files', 'targeted_ad_reduced_detrend_orthog.Rds'))

### PD
targeted_pd_reduced_detrend_orthog <- reduced_logistic_analysis(data = wide_data_targeted_detrend, 
                                                         target ="PD_ind", PD_ind = T, 
                                                         types = c("PD", "CO", "CM", "CY"),
                                                         include_age_gender = F, adpd_detrend = T, num = 3
)


untargeted_pd_reduced_detrend_orthog <- reduced_logistic_analysis(data = wide_data_untargeted_detrend, 
                                                           target ="PD_ind", PD_ind = T, 
                                                           types = c("PD", "CO", "CM", "CY"),
                                                           include_age_gender = F, adpd_detrend = T, empri0 = 1000, num = 3
)
lipids_pd_reduced_detrend_orthog <- reduced_logistic_analysis(data = wide_data_lipids_detrend, 
                                                       target ="PD_ind", PD_ind = T, 
                                                       types = c("PD", "CO", "CM", "CY"),
                                                       include_age_gender = F, adpd_detrend = T
)


saveRDS(targeted_pd_reduced_detrend_orthog, file = here('aging_output_files', 'targeted_pd_reduced_detrend_orthog.Rds'))
saveRDS(untargeted_pd_reduced_detrend_orthog, file = here('aging_output_files', 'untargeted_pd_reduced_detrend_orthog.Rds'))
saveRDS(lipids_pd_reduced_detrend_orthog, file = here('aging_output_files', 'lipids_pd_reduced_detrend_orthog.Rds'))



targeted_ad_c_roc_detrend_orthog <- roc_fig(targeted_ad_reduced_detrend_orthog$roc_df, targeted_ad_reduced_detrend_orthog$predtruth_df) + 
  labs(title = 'Targeted')

targeted_pd_c_roc_detrend_orthog <- roc_fig(targeted_pd_reduced_detrend_orthog$roc_df, targeted_pd_reduced_detrend_orthog$predtruth_df) + 
  labs(title = 'Targeted')


a <- which(targeted_pd_reduced_detrend_orthog$predtruth_df$type != "AD")
aa <- targeted_pd_reduced_detrend_orthog$predtruth_df %>%
  filter(type != "AD")
fpr_tpr(aa$imp1, aa$type == "PD")



