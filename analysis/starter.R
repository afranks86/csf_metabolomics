###########################

### Starter code from Prof Franks ###

###########################

# library(tidyverse)
# library(ggridges)
# library(gbm3)
# library(patchwork)
# library(scales)
# library(quantreg)
# library(glmnet)
# library(mice)
# 
# source("utility.R")
# 
# processed_files <- dir(pattern="^preprocessed_lipid_data*")
# 
# 
# processed_files <- dir(pattern="^preprocessed_gotms_data*")
# 
# 
# processed_files <- dir(pattern="^preprocessed_untargeted_data*")
# 
# ## load("preprocessed_gotms_data.RData")
# 
# ## Most recent file
# load(max(processed_files[grep("-20+", processed_files)]))
# 
# 
# wide_data <- subject_data %>%     
#     filter(!(Type %in% c("Other"))) %>%
#     mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
#     dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Batch", "Name",
#                           "Id", "Index", "GBAStatus", "GBA_T369M",
#                           "cognitive_status")) %>%
#     spread(key=Lipid, value=Abundance)
# 
# 
# Yage <- wide_data %>% dplyr::select(Age) %>% as.matrix
# 
# Ytype <- wide_data %>% dplyr::select(Type) %>% as.matrix
# 
# X <- wide_data %>% dplyr::select(-one_of(c("Age", "Gender", "Type", "APOE",
#                                            "Mode", "RunIndex", "Type2"))) %>%
#     as.matrix
# 
# apply(X, 2, function(x) mean(is.na(x))) %>% summary
# X[, 100]
# 
# 
# X <- X[, colMeans(is.na(X)) < 0.2]
# colnames(X) <- make.names(colnames(X))
# mice_X <- mice(X)
# 
# Xcomp <- complete(mice_X) %>% as.matrix
# 
# res <- cv.glmnet(x=Xcomp, y=Yage, family="gaussian", alpha=0.5)
# coef(res, res$lambda.min)
# plot(res)





###########################

### nathan's universal script for all his glmnet analyses ###

###########################


library(tidyverse)
library(ggExtra)
library(magrittr)
library(microbenchmark)
#library(mvtnorm)
#library(mvnfast)
#library(rstiefel)
#library(mgCov)
library(Amelia)
library(modelr)
#library(robust)
library(ggridges)
library(Matrix)
#library(car)
#library(patchwork)
library(glmnet)
library(ROCR)
library(latex2exp)
library(rsample)
library(yardstick)
library(here)
library(CCA)
library(gridExtra)
library(gbm3)
library(MetaboAnalystR)
source(here("analysis/utility.R"))


theme_set(theme_bw(base_size = 20))

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
load(max(file.path(data_path, 'analysis', processed_files[grep("-20+", processed_files)])))
load(file.path(data_path, 'data', 'got-ms',"identification_map.RData"), verbose = T)

wide_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend",
                          "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=Abundance)
dim(wide_data)


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



### code snippit from Alex for better metabolite_matching
subject_data %>% filter(Metabolite == "176 Results") %>%
    select(Mode) %>% table

subject_data$Metabolite %>% unique

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
    select("Name" = MetId, "Metabolite" =  MetaboliteName)






### NOTE: For lipids, need to load in new files. the GOT subject_data df is overwritten with the below lines!!

processed_files_lipids <- dir(path = file.path(data_path, 'analysis'), pattern="^preprocessed_lipid_data*")
## Most recent file
load(max(file.path(data_path, 'analysis', processed_files_lipids[grep("-20+", processed_files_lipids)])))


# processed_files_lipids <- dir(path = data_path, pattern="^preprocessed_lipid_data*")
# load(max(file.path(data_path, processed_files_lipids[grep("-20+", processed_files_lipids)])))



wide_data_lipids <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Name", "Mode", "RunIndex")) %>%
    spread(key=Lipid, value=Abundance) %>%
    as_tibble(.name_repair = 'minimal')



wide_data_combined <- wide_data %>%
    #remove duplicate columns
    select(-c(Age, Type, Gender, Batch, Index, GBAStatus, GBA_T369M, cognitive_status, APOE, Type2)) %>%
    inner_join(wide_data_lipids, by = 'Id') 


### untargeted data
processed_files_untargeted <- dir(path = file.path(data_path, 'analysis'), pattern="^preprocessed_untargeted_data*")
load(max(file.path(data_path, 'analysis', processed_files_untargeted[grep("-20+", processed_files_untargeted)])), verbose = T)

raw_data_untargeted <- subject_data

wide_data_untargeted <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    select(Age, Gender, APOE, Type, Metabolite, Abundance, GBAStatus, GBA_T369M, Id) %>%
    #mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "RawScaled", "Trend",
    #                       "RunIndex", "Name", 'Id')) %>%
    spread(key=Metabolite, value=Abundance)


untargeted_columns <- wide_data_untargeted %>% 
    map_dbl(~ sum(is.na(.x))/nrow(wide_data_untargeted) < .9) %>%
    names

#smaller wide_data_untargeted with columns with < .9 NA
wide_data_untargeted_dropped <- wide_data_untargeted %>%
    select(untargeted_columns)



# Targeted data
load(file.path(data_path, 'data', 'preprocessed_csf_data.RData'), verbose =  T)
raw_data_targeted <- subject_data
wide_data_targeted <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    select(Age, Gender, APOE, Type, Metabolite, Abundance, GBAStatus, GBA_T369M, Id) %>%
    #mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    # dplyr::select(-one_of("Raw", "RawScaled", "Trend",
    #                       "RunIndex", "Name", 'Id')) %>%
    spread(key=Metabolite, value=Abundance)


#create foldid so we can test different alphas on same sets
set.seed(1)
#foldid <- sample(nrow(imputed_pd_co_y))


############################

### Helper Functions ###

############################


#' filter type and impute using amelia
#' col is a vector of strings, one of the levels of type
filter_and_impute <- function(data, types){
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
    #Y[Y==0] <- NA
    
    #do imputation
    Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
    #our functions take matrices, and we add age/gender (1 = F, 2 = M)
    Y_tmp <- t(Yt) %>% 
        as_tibble %>%
        mutate(Age = filtered$Age,
               Gender = filtered$Gender)
    #convert gender to a dummy var (1 if male, 0 if female)
    Y <- model.matrix(~., Y_tmp)
    
    colnames(Y) <- str_replace_all(colnames(Y), '`', '') %>%
        str_replace_all('\\\\', '')
    
    
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

#does loocv for logistic regression, using deviance loss by default (check?), eval at different elastic net alphas
#if penalize_age_gender is false, set penalty coeff to 0
fit_glmnet <- function(features, labels, alpha, penalize_age_gender = TRUE, family= 'binomial'){
    set.seed(1)
    #set foldid so that same folds every time (it's loocv so it's just an ordering)
    foldid <- sample(nrow(features))
    
    #set penalty factors. (1 for each is default)
    p_factors <- rep(1, ncol(features))
    #set age/gender pfactors to 0 if penalize_age_gender flag is true. (assumes age/gender is very important)
    if(penalize_age_gender == FALSE){
        age_gender_index <- which(colnames(features) %in% c('Age', 'GenderM'))
        p_factors[age_gender_index] <- 0
    }
    #grouped = false is already enforced for small folds. I'm just writing it explicitly to avoid the warning.
    if(family == 'binomial'){
        fit <- cv.glmnet(features, labels, family = 'binomial', 
                         type.measure = 'deviance', nfolds = nrow(features),
                         foldid = foldid, alpha = alpha, standardize = TRUE, penalty.factor = p_factors, grouped = FALSE)
    }
    else if(family == 'gaussian'){
        fit <- cv.glmnet(features, labels, family = 'gaussian', 
                         type.measure = 'mse', nfolds = nrow(features),
                         foldid = foldid, alpha = alpha, standardize = TRUE, penalty.factor = p_factors, grouped = FALSE)
        
    }
    
    
    return(fit)
}

# function to fit glmnet on n-1 observations and predict on 1, and do n times.
# to use when you have a single lambda value you want to use (ie min lambda on full dataset)
# index is the one observation to remove
# only pass in binomial or gaussian
#https://stats.stackexchange.com/questions/304440/building-final-model-in-glmnet-after-cross-validation
loo_pred_glmnet <- function(lambda, index, features, labels, alpha, penalize_age_gender = TRUE, family = 'binomial'){
    #features and label, leaving out one observation for training
    loo_features <- features[-index,]
    loo_labels <- labels[-index]
    
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
    
    #documentation says we should avoid passing in single value of lambda. how else to do?
    #fit on n-1
    fit <- glmnet(loo_features, loo_labels, family = family, lambda = lambda, alpha = alpha, standardize = TRUE, penalty.factor = p_factors, gruoped = FALSE)
    
    #predict on 1
    pred <- predict(fit, newx = new_features, type = 'response', s = lambda)
    return(pred)
}

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

# Fit a model on the full data to find lambda, and use that lambda for leave one out.
loo_cvfit_glmnet <- function(index, features, labels, alpha, penalize_age_gender = TRUE, family = 'binomial'){
    full_fit <- fit_glmnet(features = features, labels = labels, alpha = alpha, family = family, penalize_age_gender = penalize_age_gender)
    lambda <- full_fit$lambda.min
    
    
    #features and label, leaving out one observation for training
    loo_features <- features[-index,]
    loo_labels <- labels[-index]
    
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
    
    #documentation says we should avoid passing in single value of lambda. how else to do?
    #fit on n-1
    fit <- glmnet(loo_features, loo_labels, family = family, lambda = lambda, alpha = alpha, standardize = TRUE, penalty.factor = p_factors)
    
    #predict on 1
    pred <- predict(fit, newx = new_features, type = 'response', s = lambda)
    return(list(fit,pred))
}




#' Impute with leave one out, then do loo_cvfit_glmnet
#' 
impute_c_loo_cvfit_glmnet <- function(index, og_data, alpha, penalize_age_gender = TRUE, family = 'binomial', imp_num = 1){
    #make sure the features only have controls
    c_rows <- og_data %>% 
        filter(Type %in% c("CO", "CM", "CY"))
    
    imputation_num <- paste0("imp", imp_num)
    
    loo_rows <- c_rows[-index,]
    new_rows <- c_rows[index,] 
    
    new_apoe <- new_rows$APOE
    new_type <- new_rows$Type
    new_id <- new_rows$Id
    new_gender <- new_rows$Gender
    new_age <- new_rows$Age
    
    new_feature <- new_rows %>% 
        dplyr::mutate(GenderM = ifelse(Gender == "M", 1, 0)) %>%
        dplyr::select(-one_of("Type2", "Type",  "APOE", "Gender", "Batch",
                              #"Data File",  (not found in dataset, so removed)
                              "Index", "GBAStatus",  "Id",
                              "GBA_T369M", "cognitive_status")) %>%
        as.matrix(nrow = 1) 
    
    # split the first imputation by type
    loo_imputed_separate_list <- purrr::map(c("CO", "CM", "CY"), ~filter_and_impute(loo_rows, .x))
    
    
    #creates a list of 3 matrices  (one for CO, CM, and CY)
    # combine the results. rbind.fill.matrix coerces NAs in the columns that don't have a match,
    # so we throw these columns out. this is different from the non-separated version, where the NAs are filled
    loo_imputed_separate_Y <- loo_imputed_separate_list %>% 
        purrr::map(~.x[[1]]) %>% 
        plyr::rbind.fill.matrix() %>%
        .[,!apply(is.na(.), 2, any)] %>%
        .[,colnames(.) != "(Intercept)"]
    

    # now we impute again (this time, there are only NA's in the held out row)
    full_separate_Y <- plyr::rbind.fill.matrix(loo_imputed_separate_Y, new_feature) 
    full_imputed_separate_Y <- amelia(t(full_separate_Y), m = 3, empri = 100)$imputations[[imputation_num]] %>% t
    
    # see how different imputations are from one another
    test <- amelia(t(full_separate_Y), m = 5, empri = 100)$imputations %>%
        purrr::map(~t(.x) %>%
                as_tibble(.name_repair = 'minimal') %>%
                select_if(function(x) any(!is.na(x))) %>%
                janitor::remove_empty('rows'))

    # # see how different the imputations are from the old version
    
    
    
    # list/unlist is a cool trick to preserve factors
    full_imputed_separate_type <- loo_imputed_separate_list %>% purrr::map(~.x[[2]]) %>% list(., new_type) %>% unlist 
    full_imputed_separate_apoe <- loo_imputed_separate_list %>% purrr::map(~.x[[3]]) %>% list(., new_apoe) %>% unlist
    full_imputed_separate_id <- loo_imputed_separate_list %>% purrr::map(~.x[[4]]) %>% list(., new_id) %>% unlist
    full_imputed_separate_age <- full_imputed_separate_Y[,'Age']
    full_imputed_separate_gender <- full_imputed_separate_Y[,'GenderM'] %>% as.factor %>%
        fct_recode(M = '1', F = '0')
    
    
    #zeroes don't make sense, so we replace them with NA
    # then we remove all columns that are completely na
    # then we remove any completely empty rows
    full_imputed_features_separate_age_tmp <- full_imputed_separate_Y %>% 
        as_tibble(.name_repair = 'minimal') %>%
        select(-c(Age)) %>%
        #mutate_all(~replace(., .==0, NA)) %>%
        #remove columns that are ALL NA
        select_if(function(x) any(!is.na(x))) %>%
        janitor::remove_empty('rows')
    
    # turn all factors into dummies, add intercept
    full_imputed_features_separate_age <- model.matrix(~., full_imputed_features_separate_age_tmp)
    
    #run the model (note that since we binded the new row to the end, the "index" should be the last row)
    fitpred <- loo_cvfit_glmnet(nrow(full_imputed_features_separate_age), full_imputed_features_separate_age, full_imputed_separate_age, 
                     alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE)
    
    message("list order is fitpred, truth, type, apoe, gender")
    list(fitpred[[1]], fitpred[[2]], new_age, new_type, new_apoe, new_gender)
    
}



#function to return fpr, tpr given prediction and true label
#label ordering goes c(negative class, positive class)
fpr_tpr <- function(pred, label, label_ordering = NULL){
    rocpred <- ROCR::prediction(pred, label, label.ordering = label_ordering)
    rocfpr_tpr <- ROCR::performance(rocpred, measure = 'tpr', x.measure = 'fpr')
    rocauc <- ROCR::performance(rocpred, measure = 'auc')
    return(tibble(fpr = deframe(rocfpr_tpr@x.values), 
                  tpr = deframe(rocfpr_tpr@y.values),
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


#get the types in our dataset (ie AD, PD, CO, ..)
all_types <- wide_data$Type %>% 
    unique %>%
    as.character



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
test_row <- wide_data_pd[5,]
#want same gender, closest age, closest batch
wide_data_control %>%
    filter(Gender == test_row$Gender) %>%
    filter(abs(Age - test_row$Age) == min(abs(Age - test_row$Age))) %>%
    filter(abs(Batch - test_row$Batch) == min(abs(Batch - test_row$Batch)))
###

#write above into a function
find_control <- function(row_num, data){
    #select a row to match
    subject <- data[row_num,]
    # find a control that matches
    wide_data_control %>%
        filter(Gender == subject$Gender) %>%
        filter(abs(Age - subject$Age) == min(abs(Age - subject$Age))) %>%
        filter(abs(Batch - subject$Batch) == min(abs(Batch - subject$Batch))) %>%
        slice(1) #if there are multiple matches for all three of these, then pick the first one (this should be random)
}

# apply function to get pd matched controls
    #the lapply returns a list, where each element is a row, a control match
    # bind_rows will join the 1 row list together into a df, and combine it with wide_data_pd
wide_data_matched_pd_c <- lapply(1:nrow(wide_data_pd), function(x) find_control(x, data = wide_data_pd)) %>% 
    bind_rows(wide_data_pd)


## do the same for ad
wide_data_matched_ad_c <- lapply(1:nrow(wide_data_ad), function(x) find_control(x, data = wide_data_ad)) %>%
    bind_rows(wide_data_ad)

### do the same for lipids
wide_data_lipids_matched_pd_c <- lapply(1:nrow(wide_data_pd_lipids), function(x) find_control(x, wide_data_pd_lipids)) %>%
    bind_rows(wide_data_pd_lipids)
wide_data_lipids_matched_ad_c <- lapply(1:nrow(wide_data_ad_lipids), function(x) find_control(x, wide_data_ad_lipids)) %>%
    bind_rows(wide_data_ad_lipids)
 



                           


