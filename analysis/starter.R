######################################

### Loading libraries, reading in data
### Also defining many helper functions used in both aging and AD/PD analysis

######################################

library(ggExtra)
library(magrittr)
library(ggridges)
#library(microbenchmark)
#library(mvtnorm)
#library(mvnfast)
#library(rstiefel)
#library(mgCov)
library(Amelia) # imputation
library(modelr)
#library(robust)
library(Matrix)
#library(car)
#library(patchwork)
library(glmnet)
library(ROCR)
library(latex2exp)
library(rsample)
library(yardstick)
library(here)
library(gridExtra)
library(gbm3)
library(MetaboAnalystR)
library(gghighlight)
library(ggtext)
library(ggforce)
library(denoiseR)
library(ropls)
library(vctrs)
library(scales)
library(furrr)
library(progressr)
library(snow)
library(tidyverse)
library(patchwork)
plan(multisession, workers = 4)
#cluster_amelia <- makeSOCKcluster(c("localhost", "localhost", "localhost", "localhost"))

source(here("analysis/utility.R"))


theme_set(theme_bw(base_size = 12) #+
              #theme(axis.title.x = element_text(size = 24),
              #      axis.title.y = element_text(size = 24))
          )


############################

### Helper Functions ###

############################

# Looking by missingness by age and Dataset
percent_missing_by_age <- function(data, name){
    missing <- data %>% 
        mutate(age_group = cut_interval(Age, n = 3, dig.lab = 2)) %>%
        group_by(age_group) %>%
        group_map(~list("age_group" = as.character(deframe(.y)),
                        "pct_missing" = sum(is.na(.x)) / prod(dim(.x)))
                  )
    
    tibble(
        "source" = name,
        !!(missing[[1]]$age_group) := missing[[1]]$pct_missing,
        !!(missing[[2]]$age_group) := missing[[2]]$pct_missing,
        !!(missing[[3]]$age_group) := missing[[3]]$pct_missing
    )
}

#' quick function to collapse the different control categories into one (ie cy, cm, co)
#' mostly useful for ad/pd analysis where we don't care as much about control age 
#' @param data is a dataframe, with `Type` as a column, with factors "AD", "PD", "CO", "CY', "CM" in Type
collapse_controls <- function(data){
    data %>%
        mutate(Type = fct_collapse(Type, "Controls" = c("CY", "CM", "CO")) %>%
                   fct_relevel("AD", after = 1))
}

#' function to convert observations that are 3 MAD away from the median
#' @param data is a dataframe, expected as "wide_data_*"
#' @param metadata is a dataframe in long format with the metadata info (ie subject_data)
modify_outliers <- function(data, metadata){
    data %>% mutate_at(vars(-one_of(names(metadata))), 
                  function(x) ifelse(abs(x) < 3*mad(x, na.rm = TRUE), x, NA)
                  )
}


#' filter type and impute using amelia
#' type is a vector of strings, one of the levels of type
filter_and_impute <- function(data, types, method = "amelia"){
    filtered <- data %>%
        filter(Type %in% types)
    
    #drop unused levels
    type <- filtered$Type %>%
        droplevels %>%
        #relevel factors alphabetically so it's at least consistent
        fct_relevel(sort(levels(.)))
    
    
    apoe <- filtered$APOE %>%
        droplevels
    
    id <- filtered$Id
    
    filtered_features <- filtered %>%
        dplyr::select(-one_of("Type2", "Type", "Gender", "Age", "APOE", "Batch",
                              #"Data File",  (not found in dataset, so removed)
                              "Index", "GBAStatus",  "Id",
                              "GBA_T369M", "cognitive_status"))
        

    #0 doesn't make any sense, so we replace all zeroes with NA, and remove totally NA columns and rows
    Y <- filtered_features %>% 
        mutate_all(~replace(., .==0, NA)) %>%
        #remove columns that are ALL NA
        select_if(function(x) any(!is.na(x))) %>%
        janitor::remove_empty('rows') %>%
        as.matrix()
    
    
    # clean up column names
    Y_colnames <- colnames(Y) %>%
        str_replace_all('`', '') %>%
        str_replace_all('\\\\', '')
    
    #do imputation
    if(method == "amelia"){
        Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
    } else if(method == "mice"){
        Yt <- mice(t(Y), m = 1, seed = 1) %>% 
            mice::complete(1)
    } else if(method == "ada"){
        # tranpose should be unnecessary here
        Y <- imputeada(t(Y))$completeObs
    }
    
    
    #our functions take matrices, and we add age/gender (1 = F, 2 = M)
    Y_tmp <- t(Yt) %>% 
        as_tibble %>%
        set_colnames(Y_colnames) %>%
        mutate(Age = filtered$Age,
               Gender = filtered$Gender)
    #convert gender to a dummy var (1 if male, 0 if female)
    Y <- model.matrix(~., Y_tmp)
    
    
    
    #keep gba for potential gba analysis. GBA is only non-NA for PD
    if('PD' %in% types){
        gba <- filtered %>%
            mutate(GBAStatus = as.factor(case_when(GBA_T369M == 'CT' ~ 'CT', 
                                                   TRUE ~ GBAStatus))) %>%
            select(GBAStatus) %>%
            deframe
        message("list order is Y, type, apoe, gba, id")
        return(list(Y, type, apoe,gba, id))
    }
    message("list order is Y, type, apoe, id")
    return(list(Y, type, apoe, id))
    
    
}



#' function to clean each of the imputations
#' @param Yt is a transposed features matrix (transposed to impute since we need long)
imputed_data_cleaning <- function(Yt, Y_colnames, age, gender, transpose = T, led = NULL){
    

    if(transpose){
        #our functions take matrices, and we add age/gender (1 = F, 2 = M)
        Y_tmp <- t(Yt) %>% 
            as_tibble %>%
            set_colnames(Y_colnames)
    }else{
        Y_tmp <- Yt %>% 
            as_tibble %>%
            set_colnames(Y_colnames)
    }
    
    
    # if the method is amelia, age/gender are thrown out before imptuation. 
        # if this is the case, add it back in.
    if(!("Age" %in% Y_colnames) | !( ("GenderM" %in% Y_colnames) | ("Gender" %in% Y_colnames) )){
        Y_tmp <- Y_tmp %>%
            mutate(Age = age,
                   Gender = gender)
        }
    
    
    # if we're interested in the led column, add it in now
    if(!is.null(led)){
        Y_tmp <- Y_tmp %>%
            mutate(led = led)
    }

    
    # this only does anything if the method for imputation is mice. amelia already has this done.
    # it keep only features and target
    Y <- Y_tmp %>% 
        dplyr::select(-any_of(c("Type2", "Type", "APOE23", "APOE24", "APOE33", "Batch",
                          #"Data File",  (not found in dataset, so removed)
                          "Index", "GBAStatus",  "Id",
                          "GBA_T369M", "cognitive_status"))) %>% 
        #convert factors to dummmies (1 if male, 0 if female)
        model.matrix.lm(~., ., na.action = "na.pass")
    

    colnames(Y) <- str_replace_all(colnames(Y), '`', '') %>%
        str_replace_all('\\\\', '')
    
    Y

}

#' check if a vector is numeric binary (specifically, if it only contains 1 or 0)
#' 
is_binary <- function(vec){
    all(sort(unique(vec)) %in% c(0,1))
}

#' scale column (with name) of df_to_scale by subtracting the mean and dividng sd of same col in df_for_ref 
scale_data <- function(col_to_scale, df_for_ref, name){
    
    metadata_cols <- c("Type2", "Type", "Gender", "Age", "APOE", "Batch", "Index", 
                       "GBAStatus", "Id", "GBA_T369M", "cognitive_status", "subject_id", 
                       "no_meds_reported", "led", "GenderM")
    # if it's metadata, return the original column without any changes.
    if(name %in% metadata_cols){
        return(col_to_scale)
    }
    
    col_for_ref <- df_for_ref %>% pull(name)
    
    # we need at least 2 observations in the ref col to get a standard deviation
    # so if the whole column has < 2 observations in df_for_ref, then return the original column without any changes
    if(sum(!is.na(col_for_ref)) < 2){
        return(col_to_scale)
    }
    
    # if the columns is binary (contains only 1 or 0), don't scale
    if(is_binary(col_for_ref)){
        return(col_to_scale)
    }
        
    
    (col_to_scale - mean(col_for_ref, na.rm = T)) / sd(col_for_ref, na.rm = T)
}


#' filter type and impute using amelia/mice. Removes columns with > 90% missingness
#' type is a vector of strings, one of the levels of type
#' empri is the argument to amelia to set ridge prior. makes it easier to converge
#' @param AD_ind/PD_ind optional flag to add ad/pd indicators to the features.
#' @param na_threshold is a number between 0 and 1. columns with more than na_threshold missingness are removed
filter_and_impute_multi <- function(data, types, method = "amelia", num = 5, empri = 100, AD_ind = FALSE, PD_ind = FALSE, transpose = T, impute = T, replace_zeroes = T, include_led = F, scale = T, na_threshold = 0.5){
    set.seed(1)
    filtered <- data %>%
        filter(Type %in% types)
    
    #drop unused levels
    type <- filtered$Type %>%
        droplevels %>%
        #relevel factors alphabetically so it's at least consistent
        fct_relevel(sort(levels(.)))
    
    
    apoe <- filtered$APOE %>%
        droplevels
    
    id <- filtered$Id
    age <- filtered$Age
    gender <- filtered$Gender
    
    if(include_led == TRUE){
        led <- filtered$led
    } else {
        led <- NULL
    }
    
    
    
    # removing columns we don't want in the imputation
    # amelia requires all data to be roughly normal (so no metadata)
    # but mice can take anything
    ####NH EDIT 11/19. maybe it's best to keep everything consistent..
    
    if(method == "amelia"){
        filtered_features <- filtered %>%
            dplyr::select(-any_of(c("Type2", "Type", "Gender", "Age", "APOE", "Batch",
                                  #"Data File",  (not found in dataset, so removed)
                                  "Index", "GBAStatus",  "Id",
                                  "GBA_T369M", "cognitive_status",
                                  #found in panuc
                                  "subject_id", "no_meds_reported", "led")))
    } else if(method == "mice"){
        filtered_features <- filtered %>%
            # don't need type since we have age.
            #dplyr::select(-one_of("Type2", "Type",  "Index",  "Id"))
            dplyr::select(-any_of(c("Type2", "Type", "Gender", "Age", "APOE", "Batch",
                                  #"Data File",  (not found in dataset, so removed)
                                  "Index", "GBAStatus",  "Id",
                                  "GBA_T369M", "cognitive_status",
                                  #found in panuc
                                  "subject_id", "no_meds_reported", "led")))
    }
    
    
    # if 0s are still in the data and we want to replace them with NA (if it's not binary var)
    # this is because 0s are data artifacts-- they wouldn't have shown up in our data if they were 0
    if(replace_zeroes){
        filtered_features <- filtered_features %>%
            mutate_if(~!is_binary(.x), function(x) replace(x, x==0, NA))
    }
    
    
    message(str_glue('scaling...'))
    # we want to make sure scaling the data happens on controls, and then is applied to ad/pd
    # so this if statement checks if we're imputing AD/PD, while making sure that at least some controls are in the data (ndistinct > 2)
    if(all("PD" %in% types | "AD" %in% types & n_distinct(data$Type) > 2)){
        
        c_data <- data %>%
            filter(Type %in% c("CO", "CY", "CM"))
        
        
        Y_filtered <- filtered_features %>% 
            #remove columns that are ALL NA
            # select_if(function(x) any(!is.na(x))) %>%
            # change columns with >50% NA to missingness indicators
            # mutate_if(function(x) sum(is.na(x))/nrow(filtered_features) > .5,
            #           ~ifelse(is.na(.x), 1, 0)) %>%
            # remove columns with > 50% NA
            select_if(function(x) sum(is.na(x))/nrow(filtered_features) <= na_threshold) %>%
            # remove rows that are totally NA
            janitor::remove_empty('rows')
        
        
        Y <- Y_filtered %>%
            # scale each column using the same scaling done on the controls
            furrr::future_map2_dfc(names(.), ~scale_data(.x, c_data, .y)) %>%
            # purrr::map2_dfc(names(.), ~scale_data(.x, c_data, .y)) %>%
            # need to use model.matrix to get the factors to turn into dummies (otherwise they're treated as numeric)
            # it would be nice we if could keep the factors for mice, but it won't work because we need to transpose
            # need to use model.matrix.lm to get the na.action argument to ignore na
            model.matrix.lm(~., ., na.action = "na.pass") %>%
            # remove the intercept column for imputation (don't worry, it'll come back)
            .[,-1]
    } else{
        
        Y <- filtered_features
        if(scale){
            # scale the data before imputation so that imputed values don't impact scale
            Y <- Y %>% 
                mutate_if(~!is_binary(.x), function(x) scale(x, center = T, scale = T))
        }
        
        Y <- Y %>% 
            #remove columns that are ALL NA
            #select_if(function(x) any(!is.na(x))) %>%
            # mutate_if(function(x) sum(is.na(x))/nrow(filtered_features) > .5,
            #           ~ifelse(is.na(.x), 1, 0)) %>%
            # remove columns with > 50% NA
            select_if(function(x) sum(is.na(x))/nrow(filtered_features) <= na_threshold) %>%
            # remove rows that are totally NA
            janitor::remove_empty('rows') %>%
            # need to use model.matrix to get the factors to turn into dummies (otherwise they're treated as numeric)
            # it would be nice we if could keep the factors for mice, but it won't work because we need to transpose
            # need to use model.matrix.lm to get the na.action argument to ignore na
            model.matrix.lm(~., ., na.action = "na.pass") %>%
            # remove the intercept column for imputation (don't worry, it'll come back)
            .[,-1]
        
        
    }
    
        
    
    
    # clean up column names
    Y_colnames <- colnames(Y) %>%
        str_replace_all('`', '') %>%
        str_replace_all('\\\\', '')
    
    if(impute){
        message(str_glue('starting to impute'))

        #do imputation
        if(method == "amelia"){
            if(transpose){
                #Yt_list <- amelia(t(Y), m = num, empri = empri, p2s = 0, parallel = 'snow', ncpus = 4, cl = cluster_amelia)$imputations
                Yt_list <- amelia(t(Y), m = num, empri = empri, p2s = 2)$imputations
                imputed_Y <- 1:num %>% paste0('imp', .) %>%
                    purrr::map(~imputed_data_cleaning(Yt_list[[.x]], Y_colnames = Y_colnames, age = age, gender = gender, led = led))
            } else{
                Yt_list <- amelia(Y, m = num, empri = empri, p2s = 2)$imputations
                imputed_Y <- 1:num %>% paste0('imp', .) %>%
                    purrr::map(~imputed_data_cleaning(Yt_list[[.x]], Y_colnames = Y_colnames, age = age, gender = gender, transpose = F, led = led))
            }
            
            message(str_glue('done imputing!'))
        } else if (method == "mice"){
            if(transpose){
                #need to change names on the tranposed dataset so that they aren't just numbers. v for variable
                Yt <- t(Y) %>% set_colnames(paste0("V",colnames(.)))
                # Note: mice by default removes collinear values for imputation.
                # by setting remove.collinar = FALSE, it's ignoring their relationship
                Yt_list <- mice(Yt, m = num, seed = 1, remove.collinear = FALSE)
                imputed_Y <- 1:num %>% furrr::future_map(~imputed_data_cleaning(Yt_list %>% mice::complete(.x), Y_colnames = Y_colnames, age = age, gender = gender, led = led))
            } else {
                Yt <- Y %>% set_colnames(paste0("V",colnames(.)))
                # Note: mice by default removes collinear values for imputation.
                # by setting remove.collinar = FALSE, it's ignoring their relationship
                Yt_list <- mice(Yt, m = num, seed = 1, remove.collinear = FALSE)
                imputed_Y <- 1:num %>% furrr::future_map(~imputed_data_cleaning(Yt_list %>% mice::complete(.x), Y_colnames = Y_colnames, age = age, gender = gender, transpose = F, led = led))
            }
        
                          
        }
    } else {
        # this is for the odd case when we're messing with missingness and don't have any. but we still want data in the same form.
        num <- 1
        imputed_Y <- imputed_data_cleaning(Y, Y_colnames = Y_colnames, age = age, gender = gender, transpose = F, led = led)
    }
    
    #keep gba for potential gba analysis. GBA is only non-NA for PD
    if('PD' %in% types){
        gba <- filtered %>%
            mutate(GBAStatus = as.factor(case_when(GBA_T369M == 'CT' ~ 'CT', 
                                                   TRUE ~ GBAStatus))) %>%
            dplyr::select(GBAStatus) %>%
            deframe
        message("list order is Y, type, apoe, id, gba")
        if(impute){
            return(1:num %>% purrr::map(~list(Y = imputed_Y[[.x]], Type = type, APOE = apoe, ID = id, GBA = gba)))
        } else {
            return(1:num %>% purrr::map(~list(Y = imputed_Y, Type = type, APOE = apoe, ID = id, GBA = gba)))
        }
        
    }
    message("list order is Y, type, apoe, id")
    if (impute){
        return(1:num %>% purrr::map(~list(Y = imputed_Y[[.x]], Type = type, APOE = apoe, ID = id)))
    } else {
        return(1:num %>% purrr::map(~list(Y = imputed_Y, Type = type, APOE = apoe, ID = id)))
    }
    
}







#' Quick helper to pull out held out row of a filter_and_impute matrix, and the corresponding metadata
#' 
#' @param imputation is the output of filter_and_impute (list with 4 elements: matrix, type, apoe, id)
#' @param age is an int, the single held out age
held_out_imputed <- function(imputation, index, age){
    imputation[[1]][index,"Age"] <- age
    # first pull out the indexed row of the matrix
    Y <- imputation[[1]][index,]
    metadata <- 2:4 %>% purrr::map(~imputation[[.x]][index])
    
    list("Y" = Y, "type" = metadata[[1]], "apoe" = metadata[[2]], "id" = metadata[[3]])
}






#' Attempt 1 at the loo mice
#' 
#'  Step 1: remove age from 1 observation, 
#'  Step 2: impute
#'  Step 3: retain the imputed row from the 1 observation
#'  Step 4: (outside this function) combine all these rows into a matrix
loo_filter_and_impute <- function(index, data, method = 'mice'){
    #make sure the features only have controls
    c_rows <- data %>% 
        filter(Type %in% c("CO", "CM", "CY"))
    
    held_out_row <- c_rows[index,]
    
    # save the age for the held out row and remove that obs
    loo_pred <- held_out_row$Age
    c_rows[index, "Age"] <- NA
    
    # impute
    loo_imputed <- filter_and_impute_multi(c_rows, types = c("CO", "CM", "CY"), method = method, num = 3)
    
    
    # pull out only the held out row of each imputation to save, with all the metadata
    loo_imputed %>%
        purrr::map(~held_out_imputed(.x, index, age = loo_pred))
    
        
}








#' Attempt 2 at loo mice
#' Impute with leave one out, then do loo_cvfit_glmnet
#' Step 1: pick observation (by index paramter)
#' Step 2: remove age from this observation
#' Step 3: impute using mice (with all metadata)
#' Step 4: fit model, and predict on 1

loo_filter_impute_fitpred <- function(index, data, method = "mice"){
    
}






















#does loocv for logistic regression, using deviance loss by default (check?), eval at different elastic net alphas
#if penalize_age_gender is false, set penalty coeff to 0
fit_glmnet <- function(features, labels, alpha, penalize_age_gender = TRUE, penalize_AD_PD = TRUE, family= 'binomial', nlambda = 200, weights = rep(1, length(labels))){
    set.seed(1)
    #set foldid so that same folds every time (it's loocv so it's just an ordering)
    foldid <- sample(nrow(features))
    if(family == "binomial"){
        # calculate what fraction of the total each class has
        freq_frac <- table(labels)/length(labels)
        # assign 1 - that value to a "weights" vector
        weights <- 1 - freq_frac[as.character(labels)]    
    }
    
    #set penalty factors. (1 for each is default)
    p_factors <- rep(1, ncol(features))
    #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
    if(penalize_age_gender == FALSE){
        age_gender_index <- which(colnames(features) %in% c('Age', 'GenderM'))
        p_factors[age_gender_index] <- 0
    }
    if(penalize_AD_PD == FALSE){
        ad_pd_index <- which(colnames(features) %in% c("AD_ind", "PD_ind"))
        p_factors[ad_pd_index] <- 0
    }
    
    
    fit <- cv.glmnet(features, labels, family = family, alpha = alpha, standardize = TRUE, penalty.factor = p_factors, 
                     nlambda = nlambda, weights = weights)

    
    
    return(fit)
}

# # function to fit glmnet on n-1 observations and predict on 1, and do n times.
# # to use when you have a single lambda value you want to use (ie min lambda on full dataset)
# # index is the one observation to remove
# # only pass in binomial or gaussian
# # DEPRECATED!! I leave this around for the times we do EDA with it in gotms_age, lipids_glmnet, ect..
# #https://stats.stackexchange.com/questions/304440/building-final-model-in-glmnet-after-cross-validation
# loo_pred_glmnet <- function(lambda, index, features, labels, alpha, penalize_age_gender = TRUE, family = 'binomial'){
#     #features and label, leaving out one observation for training
#     loo_features <- features[-index,]
#     loo_labels <- labels[-index]
#     
#     #features and labels on the held out observation
#     
#     new_features <- features[index,] %>% matrix(nrow = 1, dimnames = list('', colnames(loo_features) ))
#     new_label <- labels[index]
#     
#     #set penalty factors. (1 for each is default)
#     p_factors <- rep(1, ncol(loo_features))
#     #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
#     if(penalize_age_gender == FALSE){
#         age_gender_index <- which(colnames(loo_features) %in% c('Age', 'GenderM'))
#         p_factors[age_gender_index] <- 0
#     }
#     
#     #documentation says we should avoid passing in single value of lambda. how else to do?
#     #fit on n-1
#     fit <- glmnet(loo_features, loo_labels, family = family, lambda = lambda, alpha = alpha, standardize = TRUE, penalty.factor = p_factors, gruoped = FALSE)
#     
#     #predict on 1
#     pred <- predict(fit, newx = new_features, type = 'response', s = lambda)
#     return(pred)
# }

# # alternative to loo_pred_glmnet
# # does loocv cross validation on each fit (With n-1 obs) and returns the fit. 
# #leads to n different lambda.mins. idea is to compare fits to see how similar they are
# loo_cvfit_glmnet <- function(index, features, labels, alpha, penalize_age_gender = TRUE, family = 'binomial'){
#     #features and label, leaving out one observation for training
#     loo_features <- features[-index,]
#     loo_labels <- labels[-index]
#     
#     #features and labels on the held out observation (converting it to the right type and giving names)
#         #ie when doing loo, length(labels) = 1
#     new_features <- features[index,] %>% matrix(nrow = length(index), dimnames = list(1:length(index), colnames(loo_features) ))
#     new_label <- labels[index]
#     
#     
#     #fit on n-1
#     #note that this does cv, so we will have n lambda.mins
#     fit <- fit_glmnet(loo_features, loo_labels, family = family, alpha = alpha, penalize_age_gender = penalize_age_gender)
#     
#     #predict on 1
#     pred <- predict(fit, newx = new_features, type = 'response', s = 'lambda.min')
#     return(list(fit, pred))
# }

#' quick function to adapt the filter_and_impute_multi() output format to work with fit_glmnet()
#' @param imps is the output of filter_and_impute_multi()
#' @param imp_num is an integer, one of 1:length(imps)
full_age_model_from_imputed <- function(imp_num, imps){
    full_matrix <- imps[[imp_num]][[1]]
    age_index <- which(colnames(full_matrix) == 'Age')
    features <- full_matrix[,-age_index]
    target <- full_matrix[,age_index]
    
    fit_glmnet(features, target, 
               alpha = 0.5, penalize_age_gender = FALSE, family = 'gaussian')
    
}

# get_full_model <- function(features, labels, alpha, penalize_age_gender = TRUE, penalize_AD_PD = TRUE, family = 'binomial', nlambda = 100){
#     full_fit <- fit_glmnet(features = features, labels = labels, alpha = alpha, family = family, penalize_age_gender = penalize_age_gender, penalize_AD_PD = penalize_AD_PD, nlambda = nlambda)
#     lambda <- full_fit$lambda.1se
#     
#     list(full_fit, lambda)
# }
    
# Fit a model on the full data to find lambda, and use that lambda for leave one out.
#' This deprecates loo_pred_glmnet
# we set penalize_age_gender/penalize_ADPD to be TRUE by default just so it's guaranteed to work with data that doesn't have these columns
    # but whenever these columns are present, the arguments should be false.
    # but even if we set the penalize tags as FALSE when the columns aren't there, everything will work (the penalize tag will be ignored)
loo_cvfit_glmnet <- function(index, features, labels, lambda = NULL, full_fit = NULL, alpha, penalize_age_gender = TRUE, penalize_AD_PD = TRUE, family = 'binomial', nlambda = 100, weights = rep(1, length(labels) - 1), post_select = FALSE){
    
    set.seed(1)

    #features and label, leaving out one observation for training
    loo_features <- features[-index,]
    loo_labels <- labels[-index]
    
    if(family == "binomial"){
        # calculate what fraction of the total each class has
        freq_frac <- table(loo_labels)/length(loo_labels)
        # assign 1 - that value to a "weights" vector
        weights <- 1 - freq_frac[as.character(loo_labels)]    
    }
    
    
    #features and labels on the held out observation
    
    new_features <- features[index,] %>% matrix(nrow = 1, dimnames = list('', colnames(loo_features) ))
    new_label <- labels[index]
    
    #set penalty factors. (1 for each is default)
    p_factors <- rep(1, ncol(loo_features))
    #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
    if(penalize_age_gender == FALSE){
        age_gender_index <- which(colnames(loo_features) %in% c('Age', 'GenderM'))
        p_factors[age_gender_index] <- 0
    }
    if(penalize_AD_PD == FALSE){
        ad_pd_index <- which(colnames(loo_features) %in% c("AD_ind", "PD_ind"))
        p_factors[ad_pd_index] <- 0
    }
    
    #documentation says we should avoid passing in single value of lambda. how else to do?
    #fit on n-1
    # fit <- glmnet(loo_features, loo_labels, family = family, lambda = lambda, alpha = alpha, standardize = TRUE, penalty.factor = p_factors)
    fit <- cv.glmnet(loo_features, loo_labels, family = family, alpha = alpha, standardize = TRUE, penalty.factor = p_factors, 
                     nlambda = nlambda, weights = weights)
    
    
    
    #predict on 1
    pred <- predict(fit, newx = new_features, type = 'response', s = "lambda.min")
    return(list(fit,pred, "cvfit" = full_fit))
}











#function to return fpr, tpr given prediction and true label
#label ordering goes c(negative class, positive class)
#' @param pred are predictions
#' @param label is the true values
#' @param y is the metric on y axis of desired plot (fpr for roc)
#' @param x is the metric on x axis of desired plot (fpr for roc)
#' Can also make precision recall curve by letting x = rec (or tpr), y = prec (or ppv)
#' summary measure options from ROCR: measure = "f"; "sens" ; "spec", "ppv", "npv"
fpr_tpr <- function(pred, label, y = 'tpr', x = 'fpr', label_ordering = NULL){
    rocpred <- ROCR::prediction(pred, label, label.ordering = label_ordering)
    rocfpr_tpr <- ROCR::performance(rocpred, measure = y, x.measure = x)
    
    # if plotting roc, get auc. if plotting precision-recall, get roprc
    if(y == 'tpr'){
        rocauc <- ROCR::performance(rocpred, measure = 'auc')
    } else{
        rocauc <- ROCR::performance(rocpred, measure = 'aucpr')
    }
    
    return(tibble(x = deframe(rocfpr_tpr@x.values), 
                  y = deframe(rocfpr_tpr@y.values),
                  auc = as.numeric(rocauc@y.values)))
}

## get idea of variable importance. note that this relies on our variables being standardized before 
# also glmnet doesn't give us a good idea of standard error because it fits using coord descent (algorithmic, not statistical)
# to fix, maybe refit model in glm?
# numbers are not on any kind of scale, but give some idea of relative importance\
#metabolites is a flag for determining whether we need to map names to their metabolite name.
importance <- function(fit, metabolites = TRUE){
    coefficients <- coef(fit, s = 'lambda.min') %>% 
        as.matrix
    #order the coefficients for weak measure of importance, and remove coefs with 0 coeff
    #coefficients_sorted <- coefficients[order(abs(coefficients), decreasing = TRUE) & abs(coefficients) > 0,]
    coefficients_sorted_with_zeroes <- coefficients[order(abs(coefficients), decreasing = TRUE),]
    coefficients_sorted <- coefficients_sorted_with_zeroes[coefficients_sorted_with_zeroes != 0]
    #doing basic cleaning that we want done, regardless of whether it's metabolite, lipid, or anything else
    names(coefficients_sorted) <- names(coefficients_sorted) %>%
        str_replace_all('`', '') %>%
        str_replace_all('\\\\', '')
    
    if(metabolites == TRUE){
        #map metabolites to their names. keep the names of gender and age, since they aren't metabolites
        names(coefficients_sorted) <- if_else(#names(coefficients_sorted) %in% c('Age', 'GenderM', 'TypeCM', 'TypeCO', 'TypeCY', 'TypePD', "APOE23", "APOE24", "APOE33", "APOE34", "APOE44", '(Intercept)',names(wide_data_lipids)),
                                              !str_detect(names(coefficients_sorted), 'Result'), #take advantage of fact that all metabolites have "results" in name
                                              names(coefficients_sorted), 
                                              str_replace_all(names(coefficients_sorted), 'Result.*', "") %>%
                                                  str_replace_all('\\.$', '') %>%
                                                  str_replace_all('^\\.+', '') %>%
                                                  str_trim() %>%
                                                  sapply(function(x) metabolite_lookup_table[match(x, metabolite_lookup_table$Name), 'Metabolite'] %>% deframe))
    }
    
    return(coefficients_sorted)
    
}

#' helper function to pull out significant variables
#' 
#' After doing leave one out prediction, we're left with a ton of almost identical fits,
#' each fit having its own set of retained variables. We pull out all variables with at least 95% presence,
#' and take the median of each of their coefficients
#' 
#' @param importance is a list of named numeric vectors, ie the output of importance()
importance_consolidated_loo <- function(retained){
    retained %>% 
        enframe %>% 
        dplyr::group_by(name) %>% 
        #this summary is for when we map 2 metabolites to the same name (within a list)
        dplyr::summarise(value = median(value)) %>% 
        ungroup() %>% 
        spread(name, value)
}

#' Quick function to get a table of mean coefficients of a full model
#' @param full_model is output of map(1:5, logistic_control_analysis, ...) with flag full_model = T
#' @param drop_missing if TRUE, will only look at coefs that appear in all 5 imputations
get_importance_tables <- function(full_model, drop_missing = F){
    # if the imputation number is not 5, then full_model is a single instance of logistic_control_analysis().
    if(length(full_model) != 5){
        importance_table <- importance(full_model) %>% 
            enframe(value = "Coef") %>%
            arrange(desc(abs(Coef))) %>%
            filter(Coef != 0) %>%
            select('Name' = name, Coef)
    } else{
        if(drop_missing){
            importance_table <- full_model %>% 
                purrr::map(~importance(.x) %>% enframe(value = "coef")) %>%
                reduce(full_join, by = "name") %>%
                rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
                drop_na() %>%
                #mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
                rowwise() %>%
                mutate(name = str_replace_all(name, "_pos|_neg", ""),
                       mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
                       median_coef = median(c(imp1, imp2, imp3, imp4, imp5)),
                       sd_coef = sd(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2)
                       ) %>%
                arrange(desc(abs(mean_coef))) %>%
                select("Name" = name, "Avg Coef" = mean_coef, "sd" = sd_coef) %>%
                ungroup() %>%
                filter(`Avg Coef` != 0)
        } else{
            importance_table <- full_model %>% 
                purrr::map(~importance(.x) %>% enframe(value = "coef")) %>%
                reduce(full_join, by = "name") %>%
                rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
                mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
                rowwise() %>%
                mutate(name = str_replace_all(name, "_pos|_neg", ""),
                       mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
                       median_coef = median(c(imp1, imp2, imp3, imp4, imp5)),
                       sd_coef = sd(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2)
                       ) %>%
                arrange(desc(abs(mean_coef))) %>%
                select("Name" = name, "Avg Coef" = mean_coef, "sd" = sd_coef) %>%
                ungroup() %>%
                filter(`Avg Coef` != 0)
        }
        }
        
    return(importance_table)
}


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


#########################################################


### Analysis Functions ####


############################################################

#' function to do logistic glmnet with "varname" as predictorr
#' @parm varname is the target variable, if it is in the predictor matrix (so pretty much just gender... I should probably fix this LOL)
#' @types a string vector of which types to keep (some subset of AD, PD, CO, CY, CM)
#' NOTE: if genderM is not varname, then it will be included as a predictor in the model
#' AD/PD/GBA flags take precedent over varname when determining the target var
#' ie if AD_ind is TRUE, then the only purpose of varname is to set plot names
logistic_control_analysis <- function(imputations, varname = "GenderM", imp_num = 1, nlambda = 100, AD_ind = FALSE, PD_ind = FALSE, GBA_ind = FALSE, APOE4_ind = FALSE, types = NULL, full_model = FALSE, penalize_age_gender = F, include_age_gender = T, alpha = 0.5){
    
    imputed_c <- imputations[[imp_num]]
    imputed_c_Y <- imputed_c[[1]]
    imputed_c_type <- imputed_c[[2]]
    imputed_c_apoe <- imputed_c[[3]]
    imputed_c_gender <- imputed_c_Y[,"GenderM"]
    
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
    
    if(include_age_gender == FALSE){
        age_gender_index <- which(colnames(imputed_c_features_target) %in% c('Age', 'GenderM'))
        imputed_c_features_target <- imputed_c_features_target[,-age_gender_index]
    }
    
    
    if(full_model){
        full_fit <- fit_glmnet(features = imputed_c_features_target, labels = imputed_c_target, 
                               alpha = alpha, family = 'binomial', 
                               penalize_age_gender = penalize_age_gender, 
                               nlambda = nlambda)
        return(full_fit)
    }
    
    #full_model <- get_full_model(features= imputed_c_features_target, imputed_c_target, alpha = 0.5, family = "binomial", penalize_AD_PD = FALSE, penalize_age_gender = FALSE, nlambda = nlambda)
    fitpred_c_loo_target <- lapply(1:nrow(imputed_c_features_target), function(x) loo_cvfit_glmnet(x, imputed_c_features_target, imputed_c_target, 
                                                                                                   alpha = alpha, family = 'binomial', penalize_age_gender = penalize_age_gender, 
                                                                                                   nlambda = nlambda))
    
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
    
    
    # get ROC curve
    roc_c_target_loo <- fpr_tpr(pred_c_loo_target, imputed_c_target)
    # roc_target_plot <- ggplot(roc_c_target_loo) + 
    #     geom_line(mapping = aes(x, y)) + 
    #     geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    #     theme_minimal() + 
    #     labs(title = paste0("ROC: ", varname,  "vs C"),
    #          subtitle = TeX('Untargeted ,$\\alpha = 0.5$, loo'),
    #          x = 'False Positive Rate',
    #          y = 'True Positive Rate') + 
    #     geom_label(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0('AUC:', round(roc_c_target_loo$auc[1], 3)),
    #                size =12)
    
    # get precision recall curve
    precision_recall <- fpr_tpr(pred_c_loo_target, imputed_c_target, x = "rec", y = "prec") %>%
        mutate(pos_class_rate = mean(imputed_c_target))
    # precision_recall_plot <- ggplot(precision_recall) + 
    #     geom_line(mapping = aes(x, y)) + 
    #     #geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    #     theme_minimal() + 
    #     labs(title = paste0("Precision Recall: ", varname,  "vs C"),
    #          subtitle = TeX('Untargeted ,$\\alpha = 0.5$, loo'),
    #          x = 'Recall',
    #          y = 'Precision') 
    #     #geom_label(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0('AUC:', round(roc_c_target_loo$auc[1], 3)),
    #     #           size =12)
    # 
    return(list(nonzero = importance_c_loo_median_target, 
                #roc_plot = roc_target_plot, 
                fit = fit_c_loo_target, 
                truth = imputed_c_target, 
                pred =pred_c_loo_target, 
                roc_df = roc_c_target_loo,
                pr_df = precision_recall)) 
}



#' does full glmnet age analysis
#' @param imputations is the output of filter_and_impute_multi()\
#' @param target is the string with the name of the target variable name. (default Age)
#' @param name is a string describing the dataset (eg GOT, Lipids, Combined)
#' @param color is a string, the variable what we want the plots to be colored by. options are gender, type, apoe, apoe4 
#' @param imp_num imputation number. 1-5
#' 
#' @return median importance, shapiro test, results table, predtruth plot
age_control_analysis <- function(imputations, target = "Age", name, color = NULL, imp_num = 1, nlambda = 100, ad_indicator = FALSE, pd_indicator = FALSE, map_got = TRUE, post_select = FALSE){
    
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
        # we'll add back in intercept in model.matrix below
        select(-all_of(c(target, "(Intercept)")))
    
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
    
    
    #full_model <- get_full_model(features= imputed_c_features_age, imputed_c_age, alpha = 0.5, family = "gaussian", penalize_AD_PD = FALSE, penalize_age_gender = FALSE, nlambda = nlambda)
    # Note: If AD_ind, PD_ind are missing from the dataset, the flag penalize_AD_Pd doesn't do anything
    # fitpred_c_loo_age <- lapply(1:nrow(imputed_c_features_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_age, imputed_c_age, lambda = full_model[[2]], full_fit = full_model[[1]],
    #                                                                                          alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE, penalize_AD_PD = FALSE, nlambda = nlambda))
    fitpred_c_loo_age <- lapply(1:nrow(imputed_c_features_age), function(x) loo_cvfit_glmnet(x, imputed_c_features_age, imputed_c_age, alpha = 0.5, 
                                                                                             family = 'gaussian', penalize_age_gender = FALSE, penalize_AD_PD = FALSE, nlambda = nlambda, post_select = post_select))
    
    
    fit_c_loo_age <- lapply(fitpred_c_loo_age, function(x) x[[1]])
    pred_c_loo_age <- lapply(fitpred_c_loo_age, function(x) x[[2]]) %>%
        unlist
    
    #fit used to get lambda
    cv_fits <- lapply(fitpred_c_loo_age, function(x) x[[3]])
    
    #some measure of variable importance. metabolites = TRUE tries to map GOT to names
    importance_c_loo_age <- lapply(fit_c_loo_age, function(x) importance(x, metabolites = map_got))
    
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

#' from https://stackoverflow.com/questions/12643391/how-to-remove-leading-0-in-a-numeric-r-variable
#' removes leading 0s for JoG formatting
numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) }

#' create ggplot for pred vs truth in these grouped imputation tables
#' @param df is a dataframe with columns truth, imp_avg, imp_min, imp_max, apoe, gender, type, id
#' @param name is what we want to call the plots in the title (a string)
#' @param color is a string, representing the variable to color points by, one of apoe, gender, type
#' @param errorbar is boolean for whether to plot an errorbar. defaults to true
#' @param position is the position of the textbox with error metrics. defaults to topleft, can also be topright, bottom right
predtruth_plot <- function(df, pred_name = "imp_avg", name, color = NULL, errorbar = TRUE, data_name = "Control", position = "topleft"){
    pred_name <- sym(pred_name)
    if(!is.null(color)){
        color <- sym(color)
    }
    
    if(position == "topleft"){
        box_pos <- tibble(x = -Inf, y = Inf, hjust = 0, vjust = 1)
    } else if(position == "topright"){
        box_pos <- tibble(x = Inf, y = Inf, hjust = 0, vjust = 1)
    } else if(position == "botright"){
        box_pos <- tibble(x = Inf, y = -Inf, hjust= 0, vjust = -1)
    } else{
        stop("position must be topleft, topright, or botright")
    }
    
    
    
    # Null model is to predict mean
    null_pred <- mean(df$truth) %>% rep(times = length(df$truth))
    rmse_null <- (df$truth - null_pred)^2 %>% mean %>% sqrt %>% round(2)
    mae_null <- (df$truth - null_pred) %>% abs %>% mean %>% round(2)
    
    if(errorbar){
        # need to manually keep width constant. ow error bar length is determined by number of points at that x
        df <- df %>%
            group_by(truth) %>%
            mutate(width = 2 * n()) %>%
            ungroup()
        
        ggplot(df) + 
            geom_errorbar(aes(x = truth, ymin = imp_min, ymax = imp_max, group = id, width = width), position = position_dodge(.3),  alpha = 0.7, size = 0.3) + 
            geom_point(aes(truth, !!pred_name, color = !!color, group = id), size = 0.75, position = position_dodge(.3)) +
            scale_color_viridis_d() +
            labs(title = name,
                 x = 'Chronological Age',
                 y = 'Predicted Age') + 
            geom_abline(intercept = 0, slope = 1) + 
            expand_limits(x = 15, y = c(15, 100)) + 
            geom_richtext(aes(x = box_pos$x, y = box_pos$y, hjust = box_pos$hjust, vjust = box_pos$vjust, label = paste0("R^2: ", cor(truth, !!pred_name, method = "pearson")^2 %>% round(2) %>% numformat(),
                                                                                                                         "<br>RMSE: ", (truth - !!pred_name)^2 %>% mean %>% sqrt %>% round(2), " [",rmse_null,"]",
                                                                                                                         "<br>MAE: ", (truth - !!pred_name) %>% abs %>% mean %>% round(2), " [",mae_null,"]")),
                          size = rel(2))
    } else {
        # same as above, but just without errorbar
        ggplot(df) + 
            geom_point(aes(truth, !!pred_name, color = !!color), size = 1, position = position_dodge(width = .2)) + 
            scale_color_viridis_d() +
            labs(title = name,
                 #subtitle = paste0(name, " averaged over 5 imputations, alpha = 0.5, loo"),
                 x = 'Chronological Age',
                 y = 'Predicted Age') + 
            geom_abline(intercept = 0, slope = 1) + 
            expand_limits(x = 15, y = c(15, 100)) + 
            geom_richtext(aes(x = box_pos$x, y = box_pos$y, hjust = box_pos$hjust, vjust = box_pos$vjust, label = paste0("R^2: ", cor(truth, !!pred_name, method = "pearson")^2 %>% round(2) %>% numformat(), 
                                                                                                                         "<br>RMSE: ", (truth - !!pred_name)^2 %>% mean %>% sqrt %>% round(2), " [",rmse_null,"]", 
                                                                                                                         "<br>MAE: ", (truth - !!pred_name) %>% abs %>% mean %>% round(2), " [",mae_null,"]")),
                          size = rel(2))
    }
    
    
}


#' Function to do post selection inference
#' ie a second round of elastic net regression using only the features that are nonzero in any of the 5 imputations
#' @param analysis is the output of age_control_analysis()
#' @param data is the output of filter_and_impute_multi()
#' @param name is just the dataset we're using (eg lipids, untargeted, targeted, GOT)
post_select <- function(data, analysis, name){
    # first, pull out the features that are in any of the imputations
    feature_subset <- analysis %>%
        purrr::map(~.x[[2]] %>% names) %>% 
        reduce(c) %>%
        setdiff("(Intercept)") %>%
        unique()
    
    #make sure to include the target variable: Age
    feature_subset <- c(feature_subset, "Age")
    
    # change the y component of each of the feature matrices to only include our feature subset
    data_subsetted <- data %>%
        purrr::map(~.x %>% 
                       purrr::list_modify(Y = .x[[1]][,colnames(.x[[1]]) %in% feature_subset])
        )
    
    
    
    #re-run the analysis on this new subsetted dataframe
    purrr::map(1:5, ~age_control_analysis(data_subsetted, name = name, color = NULL, imp_num = .x))
    
}



#' Function to do univariate glm's.
#' @param metabolite is a string name in data
#' @param var is the other variable (either predicting/a predictor of metabolite)
#' @param conc is whether we're looking at concentration or not. this puts metabolite on the lhs
#' @param get_coefs is whether to return full coefficient estimates, or just the t/p values
age_metabolite_p <- function(data, metabolite, var = "Age", family = "gaussian", conc = FALSE,
                             get_coefs = FALSE){
    
    if(conc == FALSE){
        form <- paste0(var, "~ `", metabolite, "`") %>%
            as.formula
    } else if (conc == TRUE){
        form <- paste0("`", metabolite, "` ~ ", var) %>% 
            as.formula()
    } else {
        stop("conc must be true or false")
    }
    
    # which summary items to keep (if get_coefs = TRUE, then include beta/sd)
    summary_cols <- if(get_coefs) 1:4 else 3:4 
    
    metabolite <- sym(metabolite)
    var = sym(var)
    df <- data %>%
        as_tibble() %>%
        select(!!var, !!metabolite)
    
    if(family %in% c("gaussian", "binomial")){
        fit <- glm(form, data = df, family = family)
        
        #get the t-value and p score (excluding intercept)
        tryCatch({
            enframe(summary(fit)$coefficients[-1,summary_cols]) %>% spread(key = name, value = value) %>%
                cbind('name' = rlang::as_string(metabolite))
        },
        error = function(w) cat('metabolite is ', metabolite, "\n")
        )
    } else if (family == "rf") {
        fit <- randomForest::randomForest(form, data = df)
        
        #predictions. this method guarantees that the results are aggregated by age.
        preds <- predict(fit, tibble(Age = pull(df, !!var)))
        
        # variation explained. this is equiv to what randomForest spits out as "% variation explained"
        pseudo_r2 <- 1 - sum((fit$y - preds)^2) / sum((fit$y - mean(fit$y))^2)
        tibble("name" = rlang::as_string(metabolite),
               "truth" = fit$y,
               "pred" = preds,
               !!var := df %>% select(!!var) %>% deframe,
               "var_explained" = pseudo_r2)
        
    }
    
    
    
    

    
}

#' Function to help process unviariate results
#' Create table with metabolite/lipid, p value, t value, name, bh corrected p value (for conc = FALSE)
#' Creates table with metabolite/lipid, variation explained (for conc = TRUE)
#' (using only the first imputation)
#' @return table with all of the variables with bh p values < 0.01
#' @param data is a imputation list (eg imputed_c_combined5)
#' @param imp_num is integer 1-5 (which imputation touse)
#' @param var/family/conc/get_coefs are as in age_metabolite_p. 
#' @param types_index is a numeric vector of the row indices associated with the types of interest (eg if we're trying to classify pd against controls, types_index would the index of all PD and controls)

bh_univariate_age <- function(data, var = "Age", family = "gaussian", conc = FALSE, imp_num = 1, scale = F, types_index,
                              get_coefs = F) {
    message(str_glue("imp number is {imp_num}"))
    if(scale){
        df <- data[[imp_num]][[1]] %>%
            as_tibble() %>%
            dplyr::slice(types_index) %>%
            dplyr::select(-any_of(c('(Intercept)', 'GenderM'))) %>%
            mutate_at(.vars = vars(-any_of(c("Age", "AD_ind", "PD_ind", "GBA_ind"))), .funs = ~scale(.x, center = T, scale = T))
    } else{
        df <- data[[imp_num]][[1]] %>%
            as_tibble() %>%
            slice(types_index) %>%
            dplyr::select(-any_of(c('(Intercept)', 'GenderM')))
    }
    
    
    p_table <- df %>%
        names %>%
        setdiff(c("Age",var)) %>%
        purrr::map(~age_metabolite_p(df, metabolite = .x, var = var, family = family, conc = conc, get_coefs = get_coefs)) %>%
        purrr::reduce(rbind) %>%
        as_tibble()
    
    if(conc == FALSE){
        p_values <- p_table %>%
            dplyr::rename('og_p_value' =  `Pr(>|z|)`)
        p_table <- cbind(p_values,'bh_p_value' = p.adjust(p_values$og_p_value, method = 'BH'))
    }
    
    p_table %>% 
        dplyr::mutate(name = str_replace_all(name, '`', ''))
    
}





#' Do ADPD logistic regression on subset of predictors
#' 
#' We're concerned that our analysis might just be picking up on drugs. 
#' We can test for this by using only a subset of predictors.
#' The subset can be chosen explicitly (eg excluding the most significant) or randomly
#' Note: We impute after randomizing, so we don't get leakage from other features.
#' Note: We remove all columns that have >10% missingness before subsetting
#' @param data is one of the wide_data_* variants from starter.R
#' @param varname A string, as in logistic_control_analysis()
#' @param features a string, the feature to use
#' @param nfeatures An int. if not NULL, ignore features and select n random features.
#' @param AD/PD/GBA_ind Chooses target var. Only one should be true.
adpd_subset_analysis <- function(data, varname, features, nfeatures = NULL, AD_ind = FALSE, PD_ind = FALSE, GBA_ind = FALSE){
    # get column names with <10% missing
    colnames_less10perct <- data %>% 
        map_dbl(~ sum(is.na(.x))/nrow(data) < .1) %>%
        names() %>%
        append(c("GBAStatus", "GBA_T369M"))
    
    if(is.null(nfeatures)){
        # our feature set will be this intersected with features param
        final_features <- intersect(colnames_less10perct, features)
    } else{
        # otherwise take a random sampling
        final_features <- sample(colnames_less10perct, size = nfeatures)
    }
    
    
    wide_df <- data %>% select(final_features)
    imputation <- filter_and_impute_multi(wide_df, c('CO', 'CY', 'CM', "AD", "PD"))
    
    
    purrr::map(1:5, ~logistic_control_analysis(imputation, varname =varname, imp_num = .x, nlambda = 200, AD_ind = AD_ind, PD_ind = PD_ind, GBA_ind = GBA_ind))
    
    
}

compute_deviance_resid <- function(obs, pred){
    sqrt(-2*(obs*log(pred) + (1-obs)*log(1-pred)))* ifelse(obs > pred,1,-1)
}

#' Quick function to recreate mummichog pathway output
#' path is the path to the mummichog output file
mummichog_plot <- function(path){
    if(!file.exists(path)){
        stop("path does not exist")
    }
    pathway_table <- read_tsv(path) %>%
        mutate(neg_log10_p = -log10(`p-value`),
               pathway = fct_reorder(pathway, neg_log10_p)) %>%
        filter(overlap_size > 0) %>%
        top_n(10, `neg_log10_p`)
    
    ggplot(pathway_table) +
        geom_col(aes(pathway, neg_log10_p)) + 
        geom_hline(yintercept = -log10(0.05), lty = "dashed") +
        labs(y = TeX("-log_{10} $p$-values (BH corrected)"), 
             x = "") + 
        theme(axis.text.y = element_text(size = rel(1.75))) + 
        coord_flip() 
}

#' Function to compute brier score
#' @param pred is a vector of predicted probabilities (1 = positive case)
#' @param label is the true label (1 = positive case)
brier <- function(pred, label){
    mean((pred - label)^2)
}


### Functions for iterative, reduced elastic net models

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
        imap_dfr(~tibble(name = .y,
                         mean_OR = mean(exp(.x), na.rm = T), 
                         sd_coef = sd(.x, na.rm = T)
        )
        )
    # map_dbl(~mean(exp(.x), na.rm = T)) %>%
    # enframe(value = "mean_OR")
}


#' get lambda value for 5 imputations for a single observation
#' @param reduced_obj is an output of reduced_age_analysis()
#' @param index is 1:85 in the case of controls, (i.e. the row number of the left-out obs)
#' @return vector of 5 lambda values, corresponding to the 5 imputations for the given index
get_lambda_median <- function(reduced_obj, index){
    reduced_obj[[1]][[index]]$fit %>%
        purrr::map_dbl(~.x$fit$lambda.min)
}



############################

### Reading in Data ###

############################
set.seed(1)
# data_path <- '~/course/ND_Metabolomics/'
# processed_files <- dir(path = data_path, pattern="^preprocessed_gotms_data*")
# ## Most recent file
# load(max(file.path(data_path, processed_files[grep("-20+", processed_files)])))
# load(file.path(data_path, 'data', 'got-ms',"identification_map.RData"))

data_path <- file.path('E:', 'Projects', 'metabolomics', 'ND_Metabolomics')
processed_files <- dir(path = file.path(data_path, 'analysis'), pattern="^preprocessed_gotms_data*")
## Most recent file
load(max(file.path(data_path, 'analysis', processed_files[grep("-20+", processed_files)])), verbose = T)
load(file.path(data_path, 'data', 'got-ms',"identification_map.RData"), verbose = T)


qc_long_got <- QC_long

wide_data_got_og <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend",
                          "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=Abundance)
dim(wide_data_got_og)


wide_data_rawscaled <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "Abundance", "Trend",
                          "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=RawScaled)


wide_data_raw <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("RawScaled", "Abundance", "Trend",
                          "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=Raw)

subject_data_detrend <- readRDS(here::here("output_files", "subject_data_got_detrend.Rds"))
wide_data_got_detrend <- subject_data_detrend %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Abundance",
                          "RunIndex", "Name","Data File", "Code")) %>%
    spread(key=Metabolite, value=Abundance_as)


### code snippit from Alex for better metabolite_matching
# subject_data %>% filter(Metabolite == "176 Results") %>%
#     select(Mode) %>% table

#subject_data$Metabolite %>% unique

feature_names <- all_matches %>% ungroup %>%
    dplyr::rename(MetId=Name) %>% 
    mutate(pre_diff = abs(`Precursor Ion` - PreOrig)) %>%
    mutate(BestMatch = ifelse(pre_diff < 0.02, Metabolite, "Unknown")) %>%
    mutate(MetaboliteName = paste(BestMatch, PreOrig, ProdOrig, MetId, Mode, sep="_")) %>%
    dplyr::select(-Metabolite)

level_key <- feature_names$MetaboliteName
names(level_key) <- paste(feature_names$MetId, feature_names$Mode, sep="_")
level_key <- level_key[!duplicated(names(level_key))]

### End code snippit

metabolite_lookup_table <- feature_names %>% 
    dplyr::select("Name" = MetId, "Metabolite" =  MetaboliteName)







### NOTE: For lipids, need to load in new files. the GOT subject_data df is overwritten with the below lines!!
processed_files_lipids <- dir(path = file.path(data_path, 'analysis'), pattern="^preprocessed_lipid_data*")
## Most recent file
load(max(file.path(data_path, 'analysis', processed_files_lipids[grep("-20+", processed_files_lipids)])))


# processed_files_lipids <- dir(path = data_path, pattern="^preprocessed_lipid_data*")
# load(max(file.path(data_path, processed_files_lipids[grep("-20+", processed_files_lipids)])))

qc_long_lipids <- QC_long

wide_data_lipids_og <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Name", "Mode", "RunIndex")) %>%
    spread(key=Lipid, value=Abundance) %>%
    as_tibble(.name_repair = 'minimal')

subject_data_detrend <- readRDS(here::here("output_files", "subject_data_lipids_detrend.Rds"))
wide_data_lipids_detrend <- subject_data_detrend %>%     
    filter(!(Type %in% c("Other"))) %>%
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Abundance",
                          "RunIndex", "Name","Data File", "Code", "Mode")) %>%
    spread(key=Lipid, value=Abundance_as) %>%
    as_tibble(.name_repair = 'minimal')

### untargeted data
processed_files_untargeted <- dir(path = file.path(data_path, 'analysis'), pattern="^preprocessed_untargeted_data*")
load(max(file.path(data_path, 'analysis', processed_files_untargeted[grep("-20+", processed_files_untargeted)])), verbose = T)

raw_data_untargeted <- subject_data

qc_long_untargeted <- QC_long

wide_data_untargeted_og <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    select(Age, Gender, APOE, Type, Metabolite, Abundance, GBAStatus, GBA_T369M, Id) %>%
    #mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "RawScaled", "Trend",
    #                       "RunIndex", "Name", 'Id')) %>%
    spread(key=Metabolite, value=Abundance)

wide_data_untargeted_raw <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    dplyr::select(Age, Gender, APOE, Type, Metabolite, Raw, GBAStatus, GBA_T369M, Id) %>%
    #mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "Abundance", "Trend",
    #                       "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=Raw)



wide_data_untargeted_rawscaled <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    dplyr::select(Age, Gender, APOE, Type, Metabolite, RawScaled, GBAStatus, GBA_T369M, Id) %>%
    #mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "Abundance", "Trend",
    #                       "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=RawScaled)

subject_data_detrend <- readRDS(here::here("output_files", "subject_data_untargeted_detrend.Rds"))

wide_data_untargeted_detrend <- subject_data_detrend %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    dplyr::select(Age, Gender, APOE, Type, Metabolite, Abundance_as, GBAStatus, GBA_T369M, Id) %>%
    # mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Abundance",
    #                       "RunIndex", "Name","Data File", "Code")) %>%
    spread(key=Metabolite, value=Abundance_as)
    

wide_data_untargeted_detrend_nolev <- subject_data_detrend %>%
    filter(!between(`m/z`, 196, 199)) %>%
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    dplyr::select(Age, Gender, APOE, Type, Metabolite, Abundance_as, GBAStatus, GBA_T369M, Id) %>%
    # mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Abundance",
    #                       "RunIndex", "Name","Data File", "Code")) %>%
    spread(key=Metabolite, value=Abundance_as)


wide_data_untargeted_detrend_nolev_noenta <- subject_data_detrend %>%
    filter(!between(`m/z`, 192, 202)) %>% #l-dopa mass is 197.19
    filter(!between(`m/z`, 300, 310)) %>% #entacapone mass is 305.29
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    dplyr::select(Age, Gender, APOE, Type, Metabolite, Abundance_as, GBAStatus, GBA_T369M, Id) %>%
    # mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Abundance",
    #                       "RunIndex", "Name","Data File", "Code")) %>%
    spread(key=Metabolite, value=Abundance_as)



# Targeted data
load(file.path(data_path, 'data', 'preprocessed_csf_data.RData'), verbose =  T)
raw_data_targeted <- subject_data
qc_long_targeted <- QC_long
wide_data_targeted_og <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    select(Age, Gender, APOE, Type, Metabolite, Abundance, GBAStatus, GBA_T369M, Id) %>%
    #mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "RawScaled", "Trend",
    #                       "RunIndex", "Name", 'Id')) %>%
    spread(key=Metabolite, value=Abundance)


#create foldid so we can test different alphas on same sets
#set.seed(1)
#foldid <- sample(nrow(imputed_pd_co_y))
subject_data_detrend <- readRDS(here::here("output_files", "subject_data_targeted_detrend.Rds"))

wide_data_targeted_detrend <- subject_data_detrend %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Abundance",
                          "RunIndex", "Name","Data File", "Code")) %>%
    spread(key=Metabolite, value=Abundance_as)


#######################################

### Coering outliers to NAs ###

#######################################

wide_data_targeted <- wide_data_targeted_og %>% modify_outliers(subject_data)
wide_data <- wide_data_got_og %>% modify_outliers(subject_data)
wide_data_lipids <- wide_data_lipids_og %>% modify_outliers(subject_data)
wide_data_untargeted <- wide_data_untargeted_og %>% modify_outliers(subject_data)

wide_data_untargeted_detrend <- wide_data_untargeted_detrend %>% modify_outliers(subject_data)
wide_data_targeted_detrend <- wide_data_targeted_detrend %>% modify_outliers(subject_data)
wide_data_detrend <- wide_data_got_detrend %>% modify_outliers(subject_data)
wide_data_lipids_detrend <- wide_data_lipids_detrend %>% modify_outliers(subject_data)



#### Misc operations ####

# get names of all the possible metadata columns in case we want to do some subsetting
metadata_cols <- c("Type2", "Type", "Gender", "Age", "APOE", "Batch",
                      #"Data File",  (not found in dataset, so removed)
                      "Index", "GBAStatus",  "Id",
                      "GBA_T369M", "cognitive_status",
                      #found in panuc
                      "subject_id", "no_meds_reported", "led", "Code")


#get the types in our dataset (ie AD, PD, CO, ..)
all_types <- wide_data$Type %>% 
    unique %>%
    as.character

#smaller wide_data_untargeted with columns with < .1 NA
untargeted_less_10perct_na_columns <- wide_data_untargeted %>% 
    map_dbl(~ sum(is.na(.x))/nrow(wide_data_untargeted) < .1) %>%
    names



# Get a clean, processed version of the data -- ie 
    # 1. remove features with >10% missingness
    # 2. Center, scale the numeric features
untargeted_c_processed <- wide_data_untargeted %>%
    filter(Type %in% c("CO", "CY", 'CM')) %>%
    select_if(~sum(is.na(.x))/nrow(wide_data_untargeted) < .1) %>%
    mutate_at(vars(-any_of(metadata_cols)), ~as.vector(scale(.x, center = TRUE, scale = TRUE)))

targeted_c_processed <- wide_data_targeted %>%
    filter(Type %in% c("CO", "CY", 'CM')) %>%
    select_if(~sum(is.na(.x))/nrow(wide_data_targeted) < .1) %>%
    mutate_at(vars(-any_of(metadata_cols)), ~as.vector(scale(.x, center = TRUE, scale = TRUE)))


untargeted_all_processed <- wide_data_untargeted %>%
    select_if(~sum(is.na(.x))/nrow(wide_data_untargeted) < .1) %>%
    mutate_at(vars(-any_of(metadata_cols)), ~as.vector(scale(.x, center = TRUE, scale = TRUE)))

untargeted_all_processed_detrend <- wide_data_untargeted_detrend %>%
    select_if(~sum(is.na(.x))/nrow(wide_data_untargeted_detrend) < .1) %>%
    mutate_at(vars(-any_of(metadata_cols)), ~as.vector(scale(.x, center = TRUE, scale = TRUE)))


wide_data_combined <- wide_data %>%
    left_join(select(wide_data_targeted, -any_of(setdiff(metadata_cols, 'Id'))), by = 'Id') %>%
    left_join(select(wide_data_untargeted, -any_of(setdiff(metadata_cols, 'Id'))), by = 'Id') %>%
    left_join(select(wide_data_lipids, -any_of(setdiff(metadata_cols, 'Id'))), by = 'Id')

combined_c_processed <- wide_data_combined %>%
    filter(Type %in% c("CO", "CY", 'CM')) %>%
    select_if(~sum(is.na(.x))/nrow(wide_data_combined) < .1) %>%
    mutate_at(vars(-any_of(metadata_cols)), ~as.vector(scale(.x, center = TRUE, scale = TRUE)))



#### PANUC data
panuc_data <- readxl::read_xlsx(file.path(data_path, 'data', 'PANUC', 'panuc-0133-2019_07_09.xlsx'))

wide_data_targeted_panuc <- panuc_data %>% 
    select(Id = subject_id, led, no_meds_reported) %>%
    inner_join(wide_data_targeted, by = "Id") %>%
    #drop na for now  (8 obs)
    filter(!is.na(led))

wide_data_untargeted_panuc <- panuc_data %>% 
    select(Id = subject_id, led, no_meds_reported) %>%
    inner_join(wide_data_untargeted, by = "Id") %>%
    #drop na for now  (8 obs)
    filter(!is.na(led))

wide_data_panuc <- panuc_data %>% 
    select(Id = subject_id, led, no_meds_reported) %>%
    inner_join(wide_data, by = "Id") %>%
    #drop na for now  (8 obs)
    filter(!is.na(led))

wide_data_lipids_panuc <- panuc_data %>% 
    select(Id = subject_id, led, no_meds_reported) %>%
    inner_join(wide_data_lipids, by = "Id") %>%
    #drop na for now  (8 obs)
    filter(!is.na(led))



#######################################

### Creating matched PD/AD C dataset ###
## Useful for controlling for outliers (eg a problem with doing PD vs CO is that there are some very young people with PD, messing with the importance of age)

#######################################


wide_data_pd <- wide_data %>% 
    filter(Type == 'PD')

wide_data_ad <- wide_data %>%
    filter(Type == 'AD')

wide_data_control <- wide_data %>%
    filter(Type %in% c('CO', 'CY', 'CM'))

wide_data_pd_lipids <- wide_data_lipids %>%
    filter(Type == 'PD')

wide_data_ad_lipids <- wide_data_lipids %>%
    filter(Type == 'AD')

wide_data_control_lipids <- wide_data_lipids %>%
    filter(Type %in% c('CO', 'CY', 'CM'))


###
#test on a single row (useful for checking to make sure it's working)
#test_row <- wide_data_pd[5,]
#want same gender, closest age, closest batch
# wide_data_control %>%
#     filter(Gender == test_row$Gender) %>%
#     filter(abs(Age - test_row$Age) == min(abs(Age - test_row$Age))) %>%
#     filter(abs(Batch - test_row$Batch) == min(abs(Batch - test_row$Batch)))
###

#write above into a function
#' @param data is the smaller dataset you want to find a match for
#' @param data_control is the larger dataset you find to find a math from
find_control <- function(row_num, data, data_control){
    #select a row to match
    subject <- data[row_num,]
    # find a control that matches
    match <- data_control %>%
        filter(Gender == subject$Gender) %>%
        filter(abs(Age - subject$Age) == min(abs(Age - subject$Age))) 
    
    if("Batch" %in% names(data_control)){
        match <- match %>%
            filter(abs(Batch - subject$Batch) == min(abs(Batch - subject$Batch)))
    }
     
    match %>% slice(1) #if there are multiple matches for all three of these, then pick the first one (this should be random)
}

# apply function to get pd matched controls
    #the lapply returns a list, where each element is a row, a control match
    # bind_rows will join the 1 row list together into a df, and combine it with wide_data_pd
wide_data_matched_pd_c <- lapply(1:nrow(wide_data_pd), function(x) find_control(x, data = wide_data_pd, data_control = wide_data_control)) %>% 
    bind_rows(wide_data_pd)


## do the same for ad
wide_data_matched_ad_c <- lapply(1:nrow(wide_data_ad), function(x) find_control(x, data = wide_data_ad, data_control = wide_data_control)) %>%
    bind_rows(wide_data_ad)

### do the same for lipids
wide_data_lipids_matched_pd_c <- lapply(1:nrow(wide_data_pd_lipids), function(x) find_control(x, wide_data_pd_lipids, wide_data_lipids)) %>%
    bind_rows(wide_data_pd_lipids)
wide_data_lipids_matched_ad_c <- lapply(1:nrow(wide_data_ad_lipids), function(x) find_control(x, wide_data_ad_lipids, wide_data_lipids)) %>%
    bind_rows(wide_data_ad_lipids)
 
# versions of the data where abundances are replaced with missingness indicators (0 if missing, 1 if present)

wide_data_targeted_naind <- wide_data_targeted %>%
    mutate(across(-any_of(metadata_cols), ~ifelse(is.na(.x), 0, 1)))
wide_data_got_naind <- wide_data %>%
    mutate(across(-any_of(metadata_cols), ~ifelse(is.na(.x), 0, 1)))
wide_data_lipids_naind <- wide_data_lipids %>%
    mutate(across(-any_of(metadata_cols), ~ifelse(is.na(.x), 0, 1)))
wide_data_untargeted_naind <- wide_data_untargeted %>%
    mutate(across(-any_of(metadata_cols), ~ifelse(is.na(.x), 0, 1)))


                           


