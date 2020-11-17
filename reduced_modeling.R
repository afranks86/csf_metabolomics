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
                                take_log = F
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
    
    if(all(target == 'Age', take_log)){
      imputed_c_target <- log(imputed_c_target)
    }
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
    
    return(list('fit' = fit,
                feature_order = colnames(features_post),
                'fit_post' = fit_post
                ))
  }
  
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
#' @param remove_cor_age removes features with cor < .3 with age, following horvath dog paper
#' @param all other params are exclusively for train_reduced_model
reduced_fit <- function(index, data, target = 'Age', 
                        types = c("CO", "CY", "CM"), 
                        empri0 = 50, family = "gaussian", 
                        AD_ind= F, PD_ind = F,
                        GBA_ind = F, APOE4_ind = F,
                        penalize_age_gender = F,
                        penalize_AD_PD = T,
                        post_select = F,
                        remove_cor_age = F,
                        include_led = F,
                        take_log = F){
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
    map2_dfc(names(.), ~scale_data(.x, train_used, .y)) %>%
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
    
    
    train_imps <- try(filter_and_impute_multi(train_reduced, types, scale = T, empri = empri0, transpose = transpose, include_led = include_led), silent = T)
    while(class(train_imps) == "try-error"){
      empri0 <- empri0 + 5
      message(str_glue('empri0 = {empri0}'))
      train_imps <- try(
        filter_and_impute_multi(train_reduced, types, scale = T, empri = empri0), 
        silent = T
      )
    }
  } else{
    # if there are no missing values, no imputation is needed, so set impute = F
    train_imps <- filter_and_impute_multi(train_reduced, types, scale = T, empri = empri0, impute = F)
    
  }
  
  
  with_progress({
    p <- progressor(steps = length(train_imps))
    train_fits <- furrr::future_map(1:length(train_imps), ~{
      p()
      Sys.sleep(.2)
      train_reduced_model(.x, train_imps,
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
                          take_log = take_log
      )
    })
  })
  
  # train_fits <- purrr::map(1:5, ~train_reduced_model(.x, train_imps, 
  #                                                           target = target,
  #                                                           family = family,
  #                                                           AD_ind = AD_ind,
  #                                                           PD_ind = PD_ind,
  #                                                           GBA_ind = GBA_ind,
  #                                                           APOE4_ind = APOE4_ind,
  #                                                           types = types,
  #                                                           penalize_age_gender = penalize_age_gender,
  #                                                           penalize_AD_PD = penalize_AD_PD
  # ))
  
  
  
  if(post_select){
    # making sure the column ordering is consistent between train and test, and cut down to post select cols
    # and converting result back into a matrix of appropriate size.
    test_reduced_x <- purrr::map(train_fits, 
      ~matrix(test_reduced_x[,.x$feature_order], 
                             nrow = 1, 
                             dimnames = list("", .x$feature_order)
              )
    )
    test_preds <- purrr::map_dbl(1:length(train_fits), ~predict(train_fits[[.x]]$fit_post, newx = test_reduced_x[[.x]], s = "lambda.min", type = "response"))
  } else{
    # making sure the column ordering is consistent between train and test
    # and converting result back into a matrix of appropriate size.
    test_reduced_x <- matrix(test_reduced_x[,train_fits[[1]]$feature_order], 
                             nrow = 1, 
                             dimnames = list("", train_fits[[1]]$feature_order) 
    )
    test_preds <- purrr::map_dbl(1:length(train_fits), ~predict(train_fits[[.x]]$fit, newx = test_reduced_x, s = "lambda.min", type = "response"))
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
  data_reduced <- furrr::future_map(1:nrow(c_data), ~reduced_fit(.x, c_data, types = types, target = target, post_select = post_select, remove_cor_age = remove_cor_age, include_led = include_led, take_log = take_log))
  
  
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





reduced_logistic_analysis <- function(data, target = 'GenderM', AD_ind = F, PD_ind = F, GBA_ind = F, APOE4_ind = F, types, empri0 = 50, penalize_age_gender = F){
  set.seed(1)
  
  message(str_glue('target is {target} and ncol is {ncol(data)}'))
  data_to_use <- data %>%
    filter(Type %in% types)
  
  # with_progress({
  #   p <- progressor(steps = nrow(data))
  #   data_reduced <- furrr::future_map(1:nrow(data), ~{
  #     p()
  #     Sys.sleep(.2)
  #     reduced_fit(.x, data_to_use,
  #       family = 'binomial', target = target,
  #       AD_ind = AD_ind, PD_ind = PD_ind,
  #       GBA_ind = GBA_ind, APOE4_ind = APOE4_ind,
  #       types = types
  #       )
  #     })
  # })
  
  data_reduced <- purrr::map(1:nrow(data_to_use), ~reduced_fit(.x, data_to_use,
                                                               family = 'binomial', target = target,
                                                               AD_ind = AD_ind, PD_ind = PD_ind,
                                                               GBA_ind = GBA_ind, APOE4_ind = APOE4_ind,
                                                               types = types, empri0 = empri0, 
                                                               penalize_age_gender = penalize_age_gender
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
#' select the ones that show up >95% of the time.
#' @param reduced_output is the output of reduced_age_analysis
retained_over_fits <- function(reduced_output){
  purrr::map(1:length(reduced_output[[1]]), ~retained_one_obs(.x, reduced_output)) %>%
    bind_rows() %>%
    #select only the variables that were present in >95% of fits
    select_if(function(x) sum(is.na(x))/length(x) < .05) %>%
    map_dbl(~median(.x, na.rm = T))
}

targeted_retained <- retained_over_fits(targeted_reduced)  
untargeted_retained <- retained_over_fits(untargeted_reduced)  
got_retained <- retained_over_fits(got_reduced)
lipids_retained <- retained_over_fits(lipids_reduced)



### plot only matched controls ####

untargeted_adpd_matched_controls <- readRDS(here("aging_output_files", "untargeted_adpd_matched_controls.Rds"))
untargeted_c_matched_age_reduced <- untargeted_reduced$pred_df %>%
  filter(id %in% untargeted_adpd_matched_controls$Id)

ggsave(filename = "age_clock_untargeted_reduced_pred_matched.png",
       plot = predtruth_plot(untargeted_c_matched_age_reduced, name = "Untargeted (matched)"),
       width = 14,
       height = 10,
       path = here("plots", "aging_figs")) #14.9 x 8.21


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



roc_fig <- function(roc_df){
  ggplot(roc_df) + 
    geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    labs(#title = "ROC: AD vs C",
         #subtitle = TeX('Untargeted'),
         x = 'False Positive Rate',
         y = 'True Positive Rate') + 
    geom_text(x = Inf, y = -Inf, 
              hjust = 1, vjust = -0.5, 
              size = 20,# label.padding = unit(1, "lines"),
              label = paste0('AUC:', round(mean(roc_df$auc), 3)))
}



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

saveRDS(targeted_ad_reduced, file = here('aging_output_files', 'targeted_ad_reduced.Rds'))
saveRDS(untargeted_ad_reduced, file = here('aging_output_files', 'untargeted_ad_reduced.Rds'))
saveRDS(lipids_ad_reduced, file = here('aging_output_files', 'lipids_ad_reduced.Rds'))
saveRDS(got_ad_reduced, file = here('aging_output_files', 'got_ad_reduced.Rds'))
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
saveRDS(targeted_pd_reduced, file = here('aging_output_files', 'targeted_pd_reduced.Rds'))
saveRDS(untargeted_pd_reduced, file = here('aging_output_files', 'untargeted_pd_reduced.Rds'))
saveRDS(lipids_pd_reduced, file = here('aging_output_files', 'lipids_pd_reduced.Rds'))
saveRDS(got_pd_reduced, file = here('aging_output_files', 'got_pd_reduced.Rds'))


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

saveRDS(targeted_adpd_reduced, file = here('aging_output_files', 'targeted_adpd_reduced.Rds'))
saveRDS(untargeted_adpd_reduced, file = here('aging_output_files', 'untargeted_adpd_reduced.Rds'))
saveRDS(lipids_adpd_reduced, file = here('aging_output_files', 'lipids_adpd_reduced.Rds'))
saveRDS(got_adpd_reduced, file = here('aging_output_files', 'got_adpd_reduced.Rds'))
 

ggsave('roc_adpd_untargeted.png',
       plot = roc_fig(untargeted_adpd_reduced$roc_df) + 
         labs(title = 'Untargeted: AD vs PD'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)
ggsave('roc_adpd_targeted.png',
       plot = roc_fig(targeted_adpd_reduced$roc_df) + 
         labs(title = 'Targeted: AD vs PD'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)
ggsave('roc_adpd_lipids.png',
       plot = roc_fig(lipids_adpd_reduced$roc_df) + 
         labs(title = 'Lipids: AD vs PD'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)
ggsave('roc_adpd_got.png',
       plot = roc_fig(got_adpd_reduced$roc_df) + 
         labs(title = 'GOT: AD vs PD'),
       path = here("plots", "adpd_figs"), 
       width = 14, height = 10)



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
