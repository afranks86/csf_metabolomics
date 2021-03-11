source(here::here("analysis", "starter.R"))



## descriptive summary table
panuc_data <- read_csv('E:/Projects/metabolomics/ND_Metabolomics/data/PANUC/panuc.csv')
ad_c_tracking <- read_csv('E:/Projects/metabolomics/ND_Metabolomics/data/ADTracking.csv') %>%
  select(LPAge, Sex = Gender, Type = PrimaryDx, OnsetAge) %>%
  mutate(Type = fct_collapse(Type,'Control' = c('CO', 'CY', 'CM')))
pd_tracking <- read_csv('E:/Projects/metabolomics/ND_Metabolomics/data/PDTracking.csv') %>%
  left_join(panuc_data, by = c('PaNUCID' = 'subject_id')) %>%
  transmute(LPAge = AgeAtDraw, Sex, race, Type = `Disease Diagnosis`, 
         cognitive_status, OnsetAge = ageatonset,
         moca_score, 
         disease_duration_onset) %>%
  filter(Type == 'PD') %>%
  mutate(Sex = case_when(
    Sex == 'Female' ~ "F",
    Sex == 'Male' ~ 'M',
    TRUE ~ Sex
    )
  )

summary_subject_info <- ad_c_tracking %>%
  bind_rows(pd_tracking) %>%
  mutate(Type = fct_relevel(Type, 'Control', after = 0L)) %>%
  group_by(Type) %>%
  summarize(
    num_subjects = n(),
    num_female = sum(Sex == "F"),
    num_white = sum(race == 'White'),
    min_LPage = min(LPAge), mean_LPage = mean(LPAge), max_LPage = max(LPAge),
    min_OnsetAge = min(OnsetAge), mean_OnsetAge = mean(OnsetAge), max_OnsetAge = max(OnsetAge),
    has_dementia = sum(cognitive_status == 'Dementia'),
    no_cognitive_impairment = sum(cognitive_status == 'No cognitive impairment')
  ) %>%
  mutate(across(-c('Type', 'num_subjects'), ~round(.x, digits = 2))) %>%
  transmute(Type, num_subjects,
            num_female,
            num_white,
            LPAge = str_glue('{mean_LPage} [{min_LPage}, {max_LPage}]'),
            OnsetAge = str_glue('{mean_OnsetAge} [{min_OnsetAge}, {max_OnsetAge}]'),
            has_dementia,
            no_cognitive_impairment)

# write_csv(summary_subject_info, path = here('ad_pd', 'summary_subject_info.csv'))


# transpose
tar_c_index <- which(imputed_all_targeted5[[1]][[2]] %in% c("CO", "CM", "CY")) %>%
  sample(size = 56)
tar_c <- imputed_all_targeted5[[1]][[1]][tar_c_index,]

tar_ad_index <- which(imputed_all_targeted5[[1]][[2]] == "AD") %>%
  sample(size = 56)
tar_ad <- imputed_all_targeted5[[1]][[1]][tar_ad_index,]

tar_pd_index <- which(imputed_all_targeted5[[1]][[2]] == "PD")
tar_pd <- imputed_all_targeted5[[1]][[1]][tar_pd_index,]

cca_lambda_c_ad <- estim.regul(tar_c, tar_ad, plt = T)
cca_ad_c <- rcc(tar_c, tar_ad, cca_lambda_c_ad$lambda1, cca_lambda_c_ad$lambda2)
plt.cc(cca_ad_c)
ggplot(tibble(cor = cca_ad_c$cor) %>% rownames_to_column() %>% slice(1:10)) + 
  geom_col(aes(rowname, cor))

cca_lambda_c_pd <- estim.regul(tar_c, tar_pd, plt = T)
cca_c_pd <- rcc(tar_c, tar_pd, cca_lambda_c_pd$lambda1, cca_lambda_c_pd$lambda2)
plt.cc(cca_c_pd)
ggplot(tibble(cor = cca_c_pd$cor) %>% rownames_to_column() %>% slice(1:10)) + 
  geom_col(aes(rowname, cor))


cca_lambda_ad_pd <- estim.regul(tar_ad, tar_pd, plt = T)
cca_ad_pd <- rcc(tar_ad, tar_pd, cca_lambda_ad_pd$lambda1, cca_lambda_ad_pd$lambda2)
plt.cc(cca_ad_pd)
ggplot(tibble(cor = cca_ad_c$cor) %>% rownames_to_column() %>% slice(1:10)) + 
  geom_col(aes(rowname, cor))



got_c_index <- which(imputed_all_got5[[1]][[2]] %in% c("CO", "CM", "CY")) %>%
  sample(size = 56)
got_c <- imputed_all_got5[[1]][[1]][got_c_index,]

got_ad_index <- which(imputed_all_got5[[1]][[2]] == "AD") %>%
  sample(size = 56)
got_ad <- imputed_all_got5[[1]][[1]][got_ad_index,]

got_pd_index <- which(imputed_all_got5[[1]][[2]] == "PD")
got_pd <- imputed_all_got5[[1]][[1]][got_pd_index,]

got_cca_lambda_c_ad <- estim.regul(got_c, got_ad, plt = T)
got_cca_ad_c <- rcc(got_c, got_ad, got_cca_lambda_c_ad$lambda1, got_cca_lambda_c_ad$lambda2)

got_cca_lambda_c_pd <- estim.regul(got_c, got_pd, plt = T)
got_cca_c_pd <- rcc(got_c, got_pd, got_cca_lambda_c_pd$lambda1, got_cca_lambda_c_pd$lambda2)

got_cca_lambda_ad_pd <- estim.regul(got_ad, got_pd, plt = T)
got_cca_ad_pd <- rcc(got_ad, got_pd, got_cca_lambda_ad_pd$lambda1, got_cca_lambda_ad_pd$lambda2)
plt.cc(cca_ad_pd)


# lda multiclass? 
# cca on ad/control, project control data onto first few components.

## Before we get into it, here's a quick distribution of age. could've picked any dataset -- same subjects
age_dist_ad_pd <- ggplot(wide_data_targeted %>% collapse_controls(), aes(Age)) +
  geom_histogram(bins = 10) +
  facet_wrap(~Type) + 
  labs(x = "Age",
       y = "Count",
       title = "Distribution of Age")

# ggsave("age_dist.png",
#        plot = age_dist_ad_pd,
#        path = here("plots", "adpd_figs"),
#        width= 14,
#        height = 10)




## Missing data overview plot

data_sources <- list(wide_data, wide_data_lipids, wide_data_targeted, wide_data_untargeted)
data_names <- list("GOT", "Lipids", "Targeted", "Untargeted")
pct_missing_data_adpd <- purrr::map2(data_sources, data_names, ~percent_missing_by_type(.x, .y)) %>%
  bind_rows %>%
  gather(key = "Type", value = "pct_missing", -source) %>%
  collapse_controls()


missingness_overview_adpd_plot <- ggplot(pct_missing_data_adpd) + 
  geom_col(aes(x = source, y= pct_missing, fill = Type), position = 'dodge') +
  labs(title = "Percent Missingness By Dataset",
       x = "Dataset",
       y = "Percent Missing")
ggsave("missingness_overview_adpd_bar.png",
       plot = missingness_overview_adpd_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)



#### PCA ---------

# note: untargeted_processed has only features with < 10% missingness, and scaled/centered features
# note: nipals handles missing values by using null weights, so no need to impute first.

pca_untargeted_x <- untargeted_all_processed %>% 
  dplyr::select(-any_of(metadata_cols)) 

pca_all_untargeted <- pca_untargeted_x %>%
  opls(algoC = "nipals")


pca_scores_all_untargeted <- pca_all_untargeted@scoreMN %>%
  as_tibble() %>%
  bind_cols(untargeted_all_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id)) %>%
  collapse_controls()

pc1_all_var_explained <- pca_all_untargeted@modelDF["p1","R2X"] %>%
  scales::percent()
pc2_all_var_explained <- pca_all_untargeted@modelDF["p2","R2X"] %>%
  scales::percent()

# # look at median age for pc1 > 0 and pc1 <0
# pca_scores_c_untargeted %>%
#   transmute(Age, p1, 
#             p1_greater_than_zero = p1 > 0) %>%
#   group_by(p1_greater_than_zero) %>%
#   summarize(p1_median_age = median(Age))

pca_type_untargeted <- ggplot(pca_scores_all_untargeted, aes(p1, p2, color = Type)) + 
  geom_point(size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "Principal Component Analysis",
       subtitle = "Colored by Phenotype",
       x = str_glue("PC1 ({pc1_all_var_explained})"),
       y = str_glue("PC2 ({pc2_all_var_explained})"))
ggsave("pca_adpd_untar.png",
       plot = pca_type_untargeted,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


### PLS-DA on phenotype

plsda_all_untargeted <- opls(pca_untargeted_x, algoC = "nipals",
       y = as.character(collapse_controls(untargeted_all_processed)$Type), predI = 2)

plsda_all_untargeted <- readRDS(here("aging_output_files", "plsda_all_untargeted.Rds"))
# saveRDS(plsda_all_untargeted, here("aging_output_files", "plsda_all_untargeted.Rds"))

plsda_pc1_varx_explained <- plsda_all_untargeted@modelDF["p1","R2X"] %>%
  scales::percent()
plsda_pc2_varx_explained <- plsda_c_untargeted@modelDF["p2","R2X"] %>%
  scales::percent()


plsda_scores_all_untargeted <- plsda_all_untargeted@scoreMN %>%
  as_tibble() %>%
  bind_cols(untargeted_all_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id)) %>%
  collapse_controls()

plsda_type_untargeted <- ggplot(plsda_scores_all_untargeted, aes(p1, p2, color = Type)) + 
  geom_point(size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "Partial Least Squares Discriminant Analysis",
       subtitle = "by phenotype",
       x = str_glue("PC1 ({plsda_pc1_varx_explained})"),
       y = str_glue("PC2 ({plsda_pc2_varx_explained})")) + 
  scale_color_viridis_d()

ggsave("plsda_type_untar.png",
       plot = plsda_type_untargeted,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


#### not stratifying the imputation at all
## most conservative way to impute

imputed_all_untargeted5 <- filter_and_impute_multi(wide_data_untargeted, c('CO', 'CY', 'CM', "AD", "PD"))

imputed_all_untargeted5 <- readRDS(here("ad_pd", "imputed_all_untargeted5.Rds"))
# saveRDS(imputed_all_untargeted5, here("ad_pd", "imputed_all_untargeted5.Rds"))

# # removing all columns with >10% missingness (to help imputation)
# colnames_less10perct_untargeted <- wide_data_untargeted %>% 
#   map_dbl(~ sum(is.na(.x))/nrow(wide_data_untargeted) < .1) %>%
#   names()
# 
# imputed_less10perct_untargeted <-wide_data_untargeted %>%
#   select_at(vars(colnames_less10perct_untargeted)) %>%
#   filter_and_impute_multi(types = c("CO", "CY", "CM", "AD", "PD"))
# 
# 
# # removing all columns with NA. still run filter and impute to get data in right format?
# subset_nona_untargeted <-wide_data_untargeted %>%
#   select_if(~ !any(is.na(.))) %>%
#   filter_and_impute_multi(types = c("CO", "CY", "CM", "AD", "PD"))
# 



######

### PCA On lipids

# Idea: do plsda and see if the weights correspond to lipid classes

#####

# get a cleaned version of the lipids data for pca
lipids_all_processed <- wide_data_lipids %>%
  select_if(~sum(is.na(.x))/nrow(wide_data_lipids) < .5) %>%
  mutate_at(vars(-any_of(metadata_cols)), ~as.vector(scale(.x, center = TRUE, scale = TRUE)))

plsda_all_lipids <- opls(lipids_all_processed %>% select(-any_of(metadata_cols)), algoC = "nipals",
                             y = as.character(collapse_controls(lipids_all_processed)$Type), predI = 2)

plsda_pc1_varx_explained_l <- plsda_all_lipids@modelDF["p1","R2X"] %>%
  scales::percent()
plsda_pc2_varx_explained_l <- plsda_c_lipids@modelDF["p2","R2X"] %>%
  scales::percent()


plsda_scores_all_lipids <- plsda_all_lipids@scoreMN %>%
  as_tibble() %>%
  bind_cols(lipids_all_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id)) %>%
  collapse_controls()

plsda_type_lipids <- ggplot(plsda_scores_all_lipids, aes(p1, p2, color = Type)) + 
  geom_point(size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "Partial Least Squares Discriminant Analysis",
       subtitle = "by phenotype",
       x = str_glue("PC1 ({plsda_pc1_varx_explained})"),
       y = str_glue("PC2 ({plsda_pc2_varx_explained})")) + 
  scale_color_viridis_d()




############################

### AD/PD Logisitic Regression, baseline model using only sex and age ###

############################
message("Untargeted ADPD logistic -------------------------------------------")

# AD vs C ---------------

base_ad_c_df <- wide_data_targeted %>%
  filter(Type %in% c('AD', 'CO', 'CM', 'CY')) %>%
  transmute(AD_ind = Type == 'AD',
            male_ind = Gender == 'M',
            Age)

base_ad_c_out <- purrr::map_df(1:nrow(base_ad_c_df),
                               ~{
                                 loo_obs <- slice(base_ad_c_df, .x)
                                 train_obs <- slice(base_ad_c_df, -.x)
                                 # calculate what fraction of the total each class has
                                 freq_frac <- table(train_obs$AD_ind)/length(train_obs$AD_ind)
                                 # assign 1 - that value to a "weights" vector
                                 weights_train <- 1 - freq_frac[as.character(train_obs$AD_ind)]    
                                 base_ad_c_mod <- glm(AD_ind ~ male_ind + Age, 
                                                      data = train_obs, 
                                                      family = 'quasibinomial',
                                                      weights = weights_train)
                                 tibble(pred = predict(base_ad_c_mod, select(loo_obs, -AD_ind), type = 'response'),
                                        label = loo_obs$AD_ind)
                               })
  

base_ad_c_roc <- fpr_tpr(pred = base_ad_c_out$pred, 
                          label = base_ad_c_out$label)

base_ad_c_brier <- brier(pred = base_ad_c_out$pred, 
                         label = base_ad_c_out$label)


# PD vs C ---------------

base_pd_c_df <- wide_data_targeted %>%
  filter(Type %in% c('PD', 'CO', 'CM', 'CY')) %>%
  transmute(PD_ind = Type == 'PD',
            male_ind = Gender == 'M',
            Age)

base_pd_c_out <- purrr::map_df(1:nrow(base_pd_c_df),
                               ~{
                                 loo_obs <- slice(base_pd_c_df, .x)
                                 train_obs <- slice(base_pd_c_df, -.x)
                                 # calculate what fraction of the total each class has
                                 freq_frac <- table(train_obs$PD_ind)/length(train_obs$PD_ind)
                                 # assign 1 - that value to a "weights" vector
                                 weights_train <- 1 - freq_frac[as.character(train_obs$PD_ind)]    
                                 base_pd_c_mod <- glm(PD_ind ~ male_ind + Age, 
                                                      data = train_obs, 
                                                      family = 'quasibinomial',
                                                      weights = weights_train)
                                 tibble(pred = predict(base_pd_c_mod, select(loo_obs, -PD_ind), type = 'response'),
                                        label = loo_obs$PD_ind)
                               })


base_pd_c_roc <- fpr_tpr(pred = base_pd_c_out$pred, 
                         label = base_pd_c_out$label)

base_pd_c_brier <- brier(pred = base_pd_c_out$pred, 
                         label = base_pd_c_out$label)
# AD vs PD --------------
base_ad_pd_df <- wide_data_targeted %>%
  filter(Type %in% c('AD', 'PD')) %>%
  transmute(AD_ind = Type == 'AD',
            male_ind = Gender == 'M',
            Age)

base_ad_pd_out <- purrr::map_df(1:nrow(base_ad_pd_df),
                               ~{
                                 loo_obs <- slice(base_ad_pd_df, .x)
                                 train_obs <- slice(base_ad_pd_df, -.x)
                                 # calculate what fraction of the total each class has
                                 freq_frac <- table(train_obs$AD_ind)/length(train_obs$AD_ind)
                                 # assign 1 - that value to a "weights" vector
                                 weights_train <- 1 - freq_frac[as.character(train_obs$AD_ind)]    
                                 base_ad_pd_mod <- glm(AD_ind ~ male_ind + Age, 
                                                      data = train_obs, 
                                                      family = 'quasibinomial',
                                                      weights = weights_train)
                                 tibble(pred = predict(base_ad_pd_mod, select(loo_obs, -AD_ind), type = 'response'),
                                        label = loo_obs$AD_ind)
                               })


base_ad_pd_roc <- fpr_tpr(pred = base_ad_pd_out$pred, 
                         label = base_ad_pd_out$label)

base_ad_pd_brier <- brier(pred = base_ad_pd_out$pred, 
                         label = base_ad_pd_out$label)
############################

### AD/PD Logisitic Regression on untargeted ###

############################




#### AD First -----------------------------



# we want an AD indicator, but not a PD one
# untargeted_all_amelia5_ad_ind <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = TRUE, add_PD_ind = FALSE))
# untargeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(untargeted_all_amelia5_ad_ind, varname ="AD_ind", imp_num = .x, nlambda = 200))
untargeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("AD", "CO", "CM", "CY")))
untargeted_ad_logistic <- readRDS(here("ad_pd", "untargeted_ad_logistic.Rds"))
# saveRDS(untargeted_ad_logistic, here("ad_pd", "untargeted_ad_logistic.Rds"))

# untargeted_ad_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)

# untargeted_ad_logistic[[1]][[2]] + 
#   labs(title = "ROC: AD vs {Controls, AD}")
# ggsave("ad_logistic_untargeted.png")

untargeted_ad_roc_table <- untargeted_ad_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

untargeted_ad_roc_plot <- ggplot(untargeted_ad_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: AD vs C",
       subtitle = TeX('Untargeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(untargeted_ad_roc_table$auc), 3)))


ggsave("roc_ad_untargeted.png",
       plot = untargeted_ad_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


untargeted_ad_pr_table <- untargeted_ad_logistic %>%
  purrr::map(~.x$pr_df) %>%
  bind_rows(.id = "imp")

untargeted_ad_pr_plot <- ggplot(untargeted_ad_pr_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 0, linetype = 2) +
  labs(title = "ROC: AD vs C",
       subtitle = TeX('Untargeted'),
       x = 'Recall',
       y = 'Precision') #+ 
  # geom_text(x = Inf, y = -Inf, 
  #           hjust = 1, vjust = -0.5, 
  #           size = 20,# label.padding = unit(1, "lines"),
  #           label = paste0('AUC:', round(mean(untargeted_ad_roc_table$auc), 3)))


# diagnostic plots
heatmap.fit(untargeted_ad_logistic[[1]]$truth, pred = untargeted_ad_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(untargeted_ad_logistic[[1]]$truth, pred = untargeted_ad_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = untargeted_ad_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "AD (untargeted) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

untargeted_ad_avg_retained <- untargeted_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_ad_in_all <- untargeted_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

untargeted_ad_in_at_least_one <- untargeted_ad_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


untar_ad_retained_table <- tibble(
  data = "Untargeted",
  response = "AD",
  avg_retained = untargeted_ad_avg_retained,
  num_in_all = length(untargeted_ad_in_all),
  num_in_any = length(untargeted_ad_in_at_least_one)
)

#relationship between retained ad and retained age . None
intersect(untargeted_in_all, untargeted_ad_in_all)

test_untargeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("AD", "PD")))
test_untargeted_ad_roc_table <- test_untargeted_ad_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

test_untargeted_ad_pr_table <- test_untargeted_ad_logistic %>%
  purrr::map(~.x$pr_df) %>%
  bind_rows(.id = "imp")

test_untargeted_ad_pr_plot <- ggplot(test_untargeted_ad_pr_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = unique(), slope = 0, linetype = 2) +
  labs(title = "ROC: AD vs C",
       subtitle = TeX('Untargeted'),
       x = 'Recall',
       y = 'Precision')

ggsave(here("nathan_test_delete_pls.png"), test_untargeted_ad_pr_plot)
#### PD -------------------------------------
# we want a PD indicator, but not an AD one
# untargeted_all_amelia5_pd_ind <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = FALSE, add_PD_ind = TRUE))
# untargeted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(untargeted_all_amelia5_pd_ind, varname ="PD_ind", imp_num = .x, nlambda = 200))

untargeted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T, types = c("CY", "CO", "CM", "PD")))
untargeted_pd_logistic <- readRDS(here("ad_pd", "untargeted_pd_logistic.Rds"))

# saveRDS(untargeted_pd_logistic, here("ad_pd", "untargeted_pd_logistic.Rds"))

# untargeted_pd_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)

# untargeted_pd_logistic[[1]][[2]] + 
#   labs(title = "ROC: PD vs {Controls, AD}")
# ggsave("pd_logistic_untargeted.png")

untargeted_pd_roc_table <- untargeted_pd_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

untargeted_pd_roc_plot <- ggplot(untargeted_pd_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('Untargeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(untargeted_pd_roc_table$auc), 3)))

ggsave("roc_pd_untargeted.png",
       plot = untargeted_pd_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)

untargeted_pd_avg_retained <- untargeted_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_pd_in_all <- untargeted_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

untargeted_pd_in_at_least_one <- untargeted_pd_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


untar_pd_retained_table <- tibble(
  data = "Untargeted",
  response = "PD",
  avg_retained = untargeted_pd_avg_retained,
  num_in_all = length(untargeted_pd_in_all),
  num_in_any = length(untargeted_pd_in_at_least_one)
)


# diagnostic plots
heatmap.fit(untargeted_pd_logistic[[1]]$truth, pred = untargeted_pd_logistic[[1]]$pred)
dev_resid_pd <- compute_deviance_resid(untargeted_pd_logistic[[1]]$truth, pred = untargeted_pd_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = untargeted_pd_logistic[[1]]$pred, y = dev_resid_pd)) +
  labs(title = "PD (untargeted) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")


# relationship between retained pd and retained ad?
intersect(untargeted_ad_in_all, untargeted_pd_in_all)

#relationship between retained pd and retained age . None
intersect(untargeted_in_all, untargeted_pd_in_all)



############################

### AD/PD Logisitic Regression on targeted ###

############################

message("targeted ADPD logistic -------------------------------------------")

# note: not tranposed since n > p
imputed_all_targeted5 <- filter_and_impute_multi(wide_data_targeted, c('CO', 'CY', 'CM', "AD", "PD"), transpose = F, empri = 18)
imputed_all_targeted5 <- readRDS(here("ad_pd", "imputed_all_targeted5.Rds"))
# saveRDS(imputed_all_targeted5, here("ad_pd", "imputed_all_targeted5.Rds"))

# # removing all columns with >10% missingness (to help imputation)
# colnames_less10perct_targeted <- wide_data_targeted %>% 
#   map_dbl(~ sum(is.na(.x))/nrow(wide_data_targeted) < .1) %>%
#   names()
# 
# # probably need to fix to lower empri.. but won't converge ow.
# imputed_less10perct_targeted <-wide_data_targeted %>%
#   select_at(vars(colnames_less10perct_targeted)) %>%
#   filter_and_impute_multi(types = c("CO", "CY", "CM", "AD", "PD"), transpose = F, empri = 20)


#### AD First -----------------------------
# we want an AD indicator, but not a PD one
# targeted_all_amelia5_ad_ind <- purrr::map2(imputed_c_targeted5, targeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = TRUE, add_PD_ind = FALSE))
# targeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(targeted_all_amelia5_ad_ind, varname ="AD_ind", imp_num = .x, nlambda = 200))
targeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("CO", "CM", "CY", "AD")))
targeted_ad_logistic <- readRDS(here("ad_pd", "targeted_ad_logistic.Rds"))
# saveRDS(targeted_ad_logistic, here("ad_pd", "targeted_ad_logistic.Rds"))

# # compare with all
# targeted_ad_logistic_all <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))

targeted_ad_roc_table <- targeted_ad_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

targeted_ad_roc_plot <- ggplot(targeted_ad_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: AD vs C",
       subtitle = TeX('Targeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(targeted_ad_roc_table$auc), 3)))

ggsave("roc_ad_targeted.png",
       plot = targeted_ad_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)

# targeted_ad_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# 
# targeted_ad_logistic[[1]][[2]] + 
#   labs(title = "ROC: AD vs {Controls, AD}",
#        subtitle = TeX('Targeted ,$\\alpha = 0.5$'))
# ggsave("ad_logistic_targeted.png")



# diagnostic plots
heatmap.fit(targeted_ad_logistic[[1]]$truth, pred = targeted_ad_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(targeted_ad_logistic[[1]]$truth, pred = targeted_ad_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = targeted_ad_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "AD (targeted) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 
# get average coefs across imputations. also 
targeted_ad_coefs <- targeted_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(name = str_replace_all(name, "_pos|_neg", ""),
         mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
         median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef))) %>%
  select("Name" = name, "Avg Coef" = mean_coef) %>%
  ungroup() %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM")))

# write separate tables for positive and negative coefs
targeted_ad_coefs %>%
  filter(`Avg Coef` > 0) %>%
  write_csv(here("ad_pd", "pos_coef_table_targeted_ad.csv"))


targeted_ad_coefs %>%
  filter(`Avg Coef` < 0) %>%
  write_csv(here("ad_pd", "neg_coef_table_targeted_ad.csv"))



targeted_ad_avg_retained <- targeted_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

targeted_ad_in_all <- targeted_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

targeted_ad_in_at_least_one <- targeted_ad_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


tar_ad_retained_table <- tibble(
  data = "Targeted",
  response = "AD",
  avg_retained = targeted_ad_avg_retained,
  num_in_all = length(targeted_ad_in_all),
  num_in_any = length(targeted_ad_in_at_least_one)
)






#relationship between retained ad and retained age . None
intersect(targeted_in_all, targeted_ad_in_all)



#### PD -------------------------------------
# we want a PD indicator, but not an AD one
# targeted_all_amelia5_pd_ind <- purrr::map2(imputed_c_targeted5, targeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = FALSE, add_PD_ind = TRUE))
# targeted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(targeted_all_amelia5_pd_ind, varname ="PD_ind", imp_num = .x, nlambda = 200))

targeted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD")))
targeted_pd_logistic <- readRDS(here("ad_pd", "targeted_pd_logistic.Rds"))
# saveRDS(targeted_pd_logistic, here("ad_pd", "targeted_pd_logistic.Rds"))

# targeted_pd_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)

targeted_pd_roc_table <- targeted_pd_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

targeted_pd_roc_plot <- ggplot(targeted_pd_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('Targeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(targeted_pd_roc_table$auc), 3)))

ggsave("roc_pd_targeted.png",
       plot = targeted_pd_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


# targeted_pd_logistic[[1]][[2]] + 
#   labs(title = "ROC: PD vs {Controls, AD}",
#        subtitle = TeX('Targeted ,$\\alpha = 0.5$'))
# ggsave("pd_logistic_targeted.png")

targeted_pd_coefs <- targeted_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(name = str_replace_all(name, "_pos|_neg", ""),
         mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
         median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef))) %>%
  select("Name" = name, "Avg Coef" = mean_coef) %>%
  ungroup() %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM")))

# write separate tables for positive and negative coefs
targeted_pd_coefs %>%
  filter(`Avg Coef` > 0) %>%
  write_csv(here("ad_pd", "pos_coef_table_targeted_pd.csv"))

targeted_pd_coefs %>%
  filter(`Avg Coef` < 0) %>%
  write_csv(here("ad_pd", "neg_coef_table_targeted_pd.csv"))

targeted_pd_avg_retained <- targeted_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

targeted_pd_in_all <- targeted_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

targeted_pd_in_at_least_one <- targeted_pd_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


tar_pd_retained_table <- tibble(
  data = "Targeted",
  response = "PD",
  avg_retained = targeted_pd_avg_retained,
  num_in_all = length(targeted_pd_in_all),
  num_in_any = length(targeted_pd_in_at_least_one)
)

# diagnostic plots
heatmap.fit(targeted_pd_logistic[[1]]$truth, pred = targeted_pd_logistic[[1]]$pred)
dev_resid_pd <- compute_deviance_resid(targeted_pd_logistic[[1]]$truth, pred = targeted_pd_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = targeted_pd_logistic[[1]]$pred, y = dev_resid_pd)) +
  labs(title = "PD (targeted) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")


# relationship between retained pd and retained ad?
intersect(targeted_ad_in_all, targeted_pd_in_all)

#relationship between retained pd and retained age . None
intersect(targeted_in_all, targeted_pd_in_all)


############################

### AD/PD Logisitic Regression on lipids ###

############################

message("lipids ADPD logistic -------------------------------------------")


imputed_all_lipids5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM', "AD", "PD"), empri = 250)
imputed_all_lipids5 <- readRDS(here("ad_pd", "imputed_all_lipids5.Rds"))
# saveRDS(imputed_all_lipids5, here("ad_pd", "imputed_all_lipids5.Rds"))

# 
# imputed_less10perct_lipids_mice <-wide_data_lipids %>%
#   select_at(vars(colnames_less10perct_lipids)) %>%
#   filter_and_impute_multi(types = c("CO", "CY", "CM", "AD", "PD"),  method = "mice")

                                                                                                           


#### AD First -----------------------------
# we want an AD indicator, but not a PD one
# lipids_all_amelia5_ad_ind <- purrr::map2(imputed_c_lipids5, lipids_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = TRUE, add_PD_ind = FALSE))
# lipids_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(lipids_all_amelia5_ad_ind, varname ="AD_ind", imp_num = .x, nlambda = 200))


lipids_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c('CO', "CM", "CY", "AD")))
lipids_ad_logistic <- readRDS(here("ad_pd", "lipids_ad_logistic.Rds"))
# saveRDS(lipids_ad_logistic, here("ad_pd", "lipids_ad_logistic.Rds"))

# lipids_ad_logistic_allmetadata <- lipids_ad_logistic
# lipids_ad_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# lipids_ad_logistic[[1]][[2]] + 
#   labs(title = "ROC: AD vs {Controls, AD}",
#        subtitle = TeX('Lipids ,$\\alpha = 0.5$'))
# ggsave("ad_logistic_lipids.png")


lipids_ad_roc_table <- lipids_ad_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

lipids_ad_roc_plot <- ggplot(lipids_ad_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: AD vs C",
       subtitle = TeX('Lipids'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(lipids_ad_roc_table$auc), 3)))

ggsave("roc_ad_lipids.png",
       plot = lipids_ad_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


# diagnostic plots
heatmap.fit(lipids_ad_logistic[[1]]$truth, pred = lipids_ad_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(lipids_ad_logistic[[1]]$truth, pred = lipids_ad_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = lipids_ad_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "AD (lipids) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

lipids_ad_coefs <- lipids_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(name = str_replace_all(name, "_pos|_neg", ""),
         mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
         median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef))) %>%
  select("Name" = name, "Avg Coef" = mean_coef) %>%
  ungroup() %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM")))

# write separate tables for positive and negative coefs
lipids_ad_coefs %>%
  filter(`Avg Coef` > 0) %>%
  write_csv(here("ad_pd", "pos_coef_table_lipids_ad.csv"))

lipids_ad_coefs %>%
  filter(`Avg Coef` < 0) %>%
  write_csv(here("ad_pd", "neg_coef_table_lipids_ad.csv"))

lipids_ad_avg_retained <- lipids_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_ad_in_all <- lipids_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

lipids_ad_in_at_least_one <- lipids_ad_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


lipids_ad_retained_table <- tibble(
  data = "Lipids",
  response = "AD",
  avg_retained = lipids_ad_avg_retained,
  num_in_all = length(lipids_ad_in_all),
  num_in_any = length(lipids_ad_in_at_least_one)
)

#relationship between retained ad and retained age . None
intersect(lipids_in_all, lipids_ad_in_all)



#### PD -------------------------------------
# we want a PD indicator, but not an AD one
# lipids_all_amelia5_pd_ind <- purrr::map2(imputed_c_lipids5, lipids_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = FALSE, add_PD_ind = TRUE))
# lipids_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(lipids_all_amelia5_pd_ind, varname ="PD_ind", imp_num = .x, nlambda = 200))

lipids_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD")))
lipids_pd_logistic <- readRDS(here("ad_pd", "lipids_pd_logistic.Rds"))
# saveRDS(lipids_pd_logistic, here("ad_pd", "lipids_pd_logistic.Rds"))

# lipids_pd_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# lipids_pd_logistic[[1]][[2]] + 
#   labs(title = "ROC: PD vs {Controls, AD}",
#        subtitle = TeX('Lipids ,$\\alpha = 0.5$'))
# ggsave("pd_logistic_lipids.png")

lipids_pd_roc_table <- lipids_pd_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

lipids_pd_roc_plot <- ggplot(lipids_pd_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('Lipids'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(lipids_pd_roc_table$auc), 3)))

ggsave("roc_pd_lipids.png",
       plot = lipids_pd_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


lipids_pd_coefs <- lipids_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(name = str_replace_all(name, "_pos|_neg", ""),
         mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
         median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef))) %>%
  select("Name" = name, "Avg Coef" = mean_coef) %>%
  ungroup() %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM")))

# write separate tables for positive and negative coefs
lipids_pd_coefs %>%
  filter(`Avg Coef` > 0) %>%
  write_csv(here("ad_pd", "pos_coef_table_lipids_pd.csv"))

lipids_pd_coefs %>%
  filter(`Avg Coef` < 0) %>%
  write_csv(here("ad_pd", "neg_coef_table_lipids_pd.csv"))


lipids_pd_avg_retained <- lipids_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_pd_in_all <- lipids_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

lipids_pd_in_at_least_one <- lipids_pd_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


lipids_pd_retained_table <- tibble(
  data = "Lipids",
  response = "PD",
  avg_retained = lipids_pd_avg_retained,
  num_in_all = length(lipids_pd_in_all),
  num_in_any = length(lipids_pd_in_at_least_one)
)

# diagnostic plots
heatmap.fit(lipids_pd_logistic[[1]]$truth, pred = lipids_pd_logistic[[1]]$pred)
dev_resid_pd <- compute_deviance_resid(lipids_pd_logistic[[1]]$truth, pred = lipids_pd_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = lipids_pd_logistic[[1]]$pred, y = dev_resid_pd)) +
  labs(title = "PD (lipids) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")






# relationship between retained pd and retained ad?
intersect(lipids_ad_in_all, lipids_pd_in_all)

#relationship between retained pd and retained age . None
intersect(lipids_in_all, lipids_pd_in_all)







############################

### AD/PD Logisitic Regression on got ###

############################

message("got ADPD logistic -------------------------------------------")


imputed_all_got5 <- filter_and_impute_multi(wide_data, c('CO', 'CY', 'CM', "AD", "PD"), empri = 250)
imputed_all_got5 <- readRDS(here("ad_pd", "imputed_all_got5.Rds"))
# saveRDS(imputed_all_got5, here("ad_pd", "imputed_all_got5.Rds"))



#### AD First -----------------------------
# we want an AD indicator, but not a PD one
# got_all_amelia5_ad_ind <- purrr::map2(imputed_c_got5, got_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = TRUE, add_PD_ind = FALSE))
# got_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(got_all_amelia5_ad_ind, varname ="AD_ind", imp_num = .x, nlambda = 200))


got_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_got5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("CO", "CY", "CM", "AD")))
got_ad_logistic <- readRDS(here("ad_pd", "got_ad_logistic.Rds"))
# saveRDS(got_ad_logistic, here("ad_pd", "got_ad_logistic.Rds"))


got_ad_roc_table <- got_ad_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

got_ad_roc_plot <- ggplot(got_ad_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: AD vs C",
       subtitle = TeX('GOT-MS'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(got_ad_roc_table$auc), 3)))

ggsave("roc_ad_got.png",
       plot = got_ad_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


# diagnostic plots
heatmap.fit(got_ad_logistic[[1]]$truth, pred = got_ad_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(got_ad_logistic[[1]]$truth, pred = got_ad_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = got_ad_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "AD (got) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

got_ad_coefs <- got_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef)))

got_ad_avg_retained <- got_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

got_ad_in_all <- got_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

got_ad_in_at_least_one <- got_ad_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


got_ad_retained_table <- tibble(
  data = "GOT",
  response = "AD",
  avg_retained = got_ad_avg_retained,
  num_in_all = length(got_ad_in_all),
  num_in_any = length(got_ad_in_at_least_one)
)




#### PD -------------------------------------
# we want a PD indicator, but not an AD one
# got_all_amelia5_pd_ind <- purrr::map2(imputed_c_got5, got_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = FALSE, add_PD_ind = TRUE))
# got_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(got_all_amelia5_pd_ind, varname ="PD_ind", imp_num = .x, nlambda = 200))

got_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_got5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD")))
got_pd_logistic <- readRDS(here("ad_pd", "got_pd_logistic.Rds"))
saveRDS(got_pd_logistic, here("ad_pd", "got_pd_logistic.Rds"))

# got_pd_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# got_pd_logistic[[1]][[2]] + 
#   labs(title = "ROC: PD vs {Controls, AD}",
#        subtitle = TeX('got ,$\\alpha = 0.5$'))
# ggsave("pd_logistic_got.png")

got_pd_roc_table <- got_pd_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

got_pd_roc_plot <- ggplot(got_pd_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('GOT-MS'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(got_pd_roc_table$auc), 3)))

ggsave("roc_pd_got.png",
       plot = got_pd_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


got_pd_coefs <- got_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef)))


got_pd_avg_retained <- got_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

got_pd_in_all <- got_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

got_pd_in_at_least_one <- got_pd_logistic %>%
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")


got_pd_retained_table <- tibble(
  data = "got",
  response = "PD",
  avg_retained = got_pd_avg_retained,
  num_in_all = length(got_pd_in_all),
  num_in_any = length(got_pd_in_at_least_one)
)




###########################

# Get full models for coefs

############################

targeted_pd_full <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T))
targeted_ad_full <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("CO", "CM", "CY", "AD"), full_model = T))
lipids_pd_full <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T))
lipids_ad_full <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("CO", "CM", "CY", "AD"), full_model = T))

targeted_adpd_full <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("PD", "AD"), full_model = T))
lipids_adpd_full <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T, types = c("PD", "AD"), full_model = T))


untargeted_pd_full <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T))


targeted_pd_full <- readRDS(file = here('ad_pd', 'targeted_pd_full.Rds'))
targeted_ad_full <- readRDS(file = here('ad_pd', 'targeted_ad_full.Rds'))
lipids_pd_full <- readRDS(file = here('ad_pd', 'lipids_pd_full.Rds'))
lipids_ad_full <- readRDS(file = here('ad_pd', 'lipids_ad_full.Rds'))
targeted_adpd_full <- readRDS(file = here('ad_pd', 'targeted_adpd_full.Rds'))
lipids_adpd_full <- readRDS(file = here('ad_pd', 'lipids_adpd_full.Rds'))


saveRDS(targeted_pd_full, file = here('ad_pd', 'targeted_pd_full.Rds'))
saveRDS(targeted_ad_full, file = here('ad_pd', 'targeted_ad_full.Rds'))
saveRDS(lipids_pd_full, file = here('ad_pd', 'lipids_pd_full.Rds'))
saveRDS(lipids_ad_full, file = here('ad_pd', 'lipids_ad_full.Rds'))

saveRDS(targeted_adpd_full, file = here('ad_pd', 'targeted_adpd_full.Rds'))
saveRDS(lipids_adpd_full, file = here('ad_pd', 'lipids_adpd_full.Rds'))

saveRDS(untargeted_pd_full, file = here('ad_pd', 'untargeted_pd_full.Rds'))





full_targeted_ad_coefs <- get_importance_tables(targeted_ad_full, drop_missing = T)

# write separate tables for positive and negative coefs
# exponentiate to get odds ratios associated with a unit increase in the predictor.
full_targeted_ad_coefs %>%
  filter(`Avg Coef` > 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "pos_OR_table_targeted_ad_full.csv"))

full_targeted_ad_coefs %>%
  filter(`Avg Coef` < 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "neg_OR_table_targeted_ad_full.csv"))

full_targeted_ad_coefs %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  mutate(name = Name %>% str_to_lower() %>% str_trim(),
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
  select(mapped_name) %>%
  write_tsv(here("ad_pd", "msea_ad_names_multivar_sig_full.txt"))


full_targeted_pd_coefs <- get_importance_tables(targeted_pd_full, drop_missing = TRUE)
# write separate tables for positive and negative coefs
full_targeted_pd_coefs %>%
  filter(`Avg Coef` > 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "pos_OR_table_targeted_pd_full.csv"))

full_targeted_pd_coefs %>%
  filter(`Avg Coef` < 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "neg_OR_table_targeted_pd_full.csv"))

full_targeted_pd_coefs %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  mutate(name = Name %>% str_to_lower() %>% str_trim(),
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
  select(mapped_name) %>%
  write_tsv(here("ad_pd", "msea_pd_names_multivar_sig_full.txt"))

#
full_lipids_ad_coefs <- get_importance_tables(lipids_ad_full, drop_missing = T)

# write separate tables for positive and negative coefs
full_lipids_ad_coefs %>%
  filter(`Avg Coef` > 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "pos_OR_table_lipids_ad_full.csv"))

full_lipids_ad_coefs %>%
  filter(`Avg Coef` < 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "neg_OR_table_lipids_ad_full.csv"))


full_lipids_pd_coefs <- get_importance_tables(lipids_pd_full, drop_missing = T)

# write separate tables for positive and negative coefs
full_lipids_pd_coefs %>%
  filter(`Avg Coef` > 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "pos_OR_table_lipids_pd_full.csv"))

full_lipids_pd_coefs %>%
  filter(`Avg Coef` < 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "neg_OR_table_lipids_pd_full.csv"))


### ad vs pd
full_targeted_adpd_coefs <- get_importance_tables(targeted_adpd_full, drop_missing = T)

# write separate tables for positive and negative coefs
full_targeted_adpd_coefs %>%
  filter(`Avg Coef` > 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "pos_OR_table_targeted_adpd_full.csv"))

full_targeted_adpd_coefs %>%
  filter(`Avg Coef` < 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "neg_OR_table_targeted_adpd_full.csv"))

full_targeted_adpd_coefs %>%
  filter(Name %in% c("(Intercept)", "Age", "GenderM"))

# lipids
full_lipids_adpd_coefs <- get_importance_tables(lipids_adpd_full, drop_missing = T)

full_lipids_adpd_coefs %>%
  filter(`Avg Coef` > 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "pos_OR_table_lipids_adpd_full.csv"))

full_lipids_adpd_coefs %>%
  filter(`Avg Coef` < 0) %>%
  transmute(Name, OR = round(exp(`Avg Coef`), digits = 2)) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(here("ad_pd", "neg_OR_table_lipids_adpd_full.csv"))

full_lipids_adpd_coefs %>%
  filter(Name %in% c("(Intercept)", "Age", "GenderM"))

#### full na indicator models

# get na_indicator data in form suitable for logistic_analysis function. doesn't actually impute anything
all_targeted_naind <- filter_and_impute_multi(wide_data_targeted_naind, c('CO', 'CY', 'CM', "AD", "PD"), transpose = F, impute =F)
all_targeted_naind <- readRDS(here("ad_pd", "all_targeted_naind.Rds"))
# saveRDS(all_targeted_naind, here("ad_pd", "all_targeted_naind.Rds"))

all_lipids_naind <- filter_and_impute_multi(wide_data_lipids_naind, c('CO', 'CY', 'CM', "AD", "PD"), transpose = F, impute =F)
all_lipids_naind <- readRDS(here("ad_pd", "all_lipids_naind.Rds"))
# saveRDS(all_lipids_naind, here("ad_pd", "all_lipids_naind.Rds"))

all_untargeted_naind <- filter_and_impute_multi(wide_data_untargeted_naind, c('CO', 'CY', 'CM', "AD", "PD"), transpose = F, impute =F)
all_untargeted_naind <- readRDS(here("ad_pd", "all_untargeted_naind.Rds"))
# saveRDS(all_untargeted_naind, here("ad_pd", "all_untargeted_naind.Rds"))


targeted_pd_full_naind <- logistic_control_analysis(all_targeted_naind, varname ="PD_ind", imp_num = 1, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T)
targeted_ad_full_naind <- logistic_control_analysis(all_targeted_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("CO", "CM", "CY", "AD"), full_model = T)
lipids_pd_full_naind <- logistic_control_analysis(all_lipids_naind, varname ="PD_ind", imp_num = 1, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T)
lipids_ad_full_naind <- logistic_control_analysis(all_lipids_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("CO", "CM", "CY", "AD"), full_model = T)

targeted_adpd_full_naind <- logistic_control_analysis(all_targeted_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("PD", "AD"), full_model = T)
lipids_adpd_full_naind <- logistic_control_analysis(all_lipids_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("PD", "AD"), full_model = T)

untargeted_pd_full_naind <- logistic_control_analysis(all_untargeted_naind, varname ="PD_ind", imp_num = 1, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T)

targeted_pd_full_naind <- readRDS(file = here('ad_pd', 'targeted_pd_full_naind.Rds'))
targeted_ad_full_naind <- readRDS(file = here('ad_pd', 'targeted_ad_full_naind.Rds'))
lipids_pd_full_naind <- readRDS(file = here('ad_pd', 'lipids_pd_full_naind.Rds'))
lipids_ad_full_naind <- readRDS(file = here('ad_pd', 'lipids_ad_full_naind.Rds'))
targeted_adpd_full_naind <- readRDS(file = here('ad_pd', 'targeted_adpd_full_naind.Rds'))
lipids_adpd_full_naind <- readRDS(file = here('ad_pd', 'lipids_adpd_full_naind.Rds'))

untargeted_pd_full_naind <- readRDS(file = here('ad_pd', 'untargeted_pd_full_naind.Rds'))



saveRDS(targeted_pd_full_naind, file = here('ad_pd', 'targeted_pd_full_naind.Rds'))
saveRDS(targeted_ad_full_naind, file = here('ad_pd', 'targeted_ad_full_naind.Rds'))
saveRDS(lipids_pd_full_naind, file = here('ad_pd', 'lipids_pd_full_naind.Rds'))
saveRDS(lipids_ad_full_naind, file = here('ad_pd', 'lipids_ad_full_naind.Rds'))
saveRDS(targeted_adpd_full_naind, file = here('ad_pd', 'targeted_adpd_full_naind.Rds'))
saveRDS(lipids_adpd_full_naind, file = here('ad_pd', 'lipids_adpd_full_naind.Rds'))

saveRDS(untargeted_pd_full_naind, file = here('ad_pd', 'untargeted_pd_full_naind.Rds'))

get_importance_tables(lipids_pd_full_naind) %>%
  mutate(Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  transmute(Name, OR = round(exp(Coef), digits = 2)) %>%
  write_csv(here("ad_pd", "OR_table_lipids_pd_naind.csv"))
  

get_importance_tables(targeted_ad_full_naind) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  mutate(Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  transmute(Name, OR = round(exp(Coef), digits = 2)) %>%
  write_csv(here("ad_pd", "OR_table_targeted_ad_naind.csv"))


get_importance_tables(lipids_adpd_full_naind) %>%
  mutate(Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  transmute(Name, OR = round(exp(Coef), digits = 2)) %>%
  write_csv(here("ad_pd", "OR_table_lipids_adpd_naind.csv"))




###### same thing without penalty
targeted_pd_full_naind_nopenalize <- logistic_control_analysis(all_targeted_naind, varname ="PD_ind", imp_num = 1, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T, penalize_age_gender = T)
targeted_ad_full_naind_nopenalize <- logistic_control_analysis(all_targeted_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("CO", "CM", "CY", "AD"), full_model = T, penalize_age_gender = T)
lipids_pd_full_naind_nopenalize <- logistic_control_analysis(all_lipids_naind, varname ="PD_ind", imp_num = 1, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T, penalize_age_gender = T)
lipids_ad_full_naind_nopenalize <- logistic_control_analysis(all_lipids_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("CO", "CM", "CY", "AD"), full_model = T, penalize_age_gender = T)
targeted_adpd_full_naind_nopenalize <- logistic_control_analysis(all_targeted_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("PD", "AD"), full_model = T, penalize_age_gender = T)
lipids_adpd_full_naind_nopenalize <- logistic_control_analysis(all_lipids_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("PD", "AD"), full_model = T, penalize_age_gender = T)

targeted_pd_full_naind_nopenalize <- readRDS(file = here('ad_pd', 'targeted_pd_full_naind_nopenalize.Rds'))
targeted_ad_full_naind_nopenalize <- readRDS(file = here('ad_pd', 'targeted_ad_full_naind_nopenalize.Rds'))
lipids_pd_full_naind_nopenalize <- readRDS(file = here('ad_pd', 'lipids_pd_full_naind_nopenalize.Rds'))
lipids_ad_full_naind_nopenalize <- readRDS(file = here('ad_pd', 'lipids_ad_full_naind_nopenalize.Rds'))
targeted_adpd_full_naind_nopenalize <- readRDS(file = here('ad_pd', 'targeted_adpd_full_naind_nopenalize.Rds'))
lipids_adpd_full_naind_nopenalize <- readRDS(file = here('ad_pd', 'lipids_adpd_full_naind_nopenalize.Rds'))


untargeted_pd_full_naind_nopenalize <- logistic_control_analysis(all_untargeted_naind, varname ="PD_ind", imp_num = 1, nlambda = 200, PD_ind = T, types = c("CO", "CM", "CY", "PD"), full_model = T, penalize_age_gender = T)
untargeted_adpd_full_naind_nopenalize <- logistic_control_analysis(all_untargeted_naind, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("AD", "PD"), full_model = T, penalize_age_gender = T)

saveRDS(targeted_pd_full_naind_nopenalize, file = here('ad_pd', 'targeted_pd_full_naind_nopenalize.Rds'))
saveRDS(targeted_ad_full_naind_nopenalize, file = here('ad_pd', 'targeted_ad_full_naind_nopenalize.Rds'))
saveRDS(lipids_pd_full_naind_nopenalize, file = here('ad_pd', 'lipids_pd_full_naind_nopenalize.Rds'))
saveRDS(lipids_ad_full_naind_nopenalize, file = here('ad_pd', 'lipids_ad_full_naind_nopenalize.Rds'))
saveRDS(targeted_adpd_full_naind_nopenalize, file = here('ad_pd', 'targeted_adpd_full_naind_nopenalize.Rds'))
saveRDS(lipids_adpd_full_naind_nopenalize, file = here('ad_pd', 'lipids_adpd_full_naind_nopenalize.Rds'))

saveRDS(untargeted_pd_full_naind_nopenalize, file = here('ad_pd', 'untargeted_pd_full_naind_nopenalize.Rds'))
saveRDS(untargeted_adpd_full_naind_nopenalize, file = here('ad_pd', 'untargeted_adpd_full_naind_nopenalize.Rds'))

get_importance_tables(lipids_pd_full_naind_nopenalize) %>%
  mutate(Coef = round(Coef, digits = 2),
         Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  write_csv(here("ad_pd", "coef_table_lipids_pd_naind_nopenalize.csv"))

get_importance_tables(lipids_ad_full_naind_nopenalize) %>%
  mutate(Coef = round(Coef, digits = 2),
         Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  write_csv(here("ad_pd", "coef_table_lipids_ad_naind_nopenalize.csv"))

get_importance_tables(lipids_adpd_full_naind_nopenalize) %>%
  mutate(Coef = round(Coef, digits = 2),
         Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  write_csv(here("ad_pd", "coef_table_lipids_adpd_naind_nopenalize.csv"))


get_importance_tables(targeted_ad_full_naind_nopenalize) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  mutate(Coef = round(Coef, digits = 2),
         Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  write_csv(here("ad_pd", "coef_table_targeted_ad_naind_nopenalize.csv"))

get_importance_tables(targeted_pd_full_naind_nopenalize) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  mutate(Coef = round(Coef, digits = 2),
         Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  write_csv(here("ad_pd", "coef_table_targeted_pd_naind_nopenalize.csv"))

get_importance_tables(targeted_adpd_full_naind_nopenalize) %>%
  filter(!(Name %in% c("(Intercept)"))) %>%
  mutate(Coef = round(Coef, digits = 2),
         Name = str_replace_all(Name, "_pos|_neg", "")) %>%
  write_csv(here("ad_pd", "coef_table_targeted_adpd_naind_nopenalize.csv"))



###########################

### Univariate AD/PD Logistic Regresion on untargeted

### For use with mummichogcommand: `mummichog -f {input_file} -o {output_name} -c 0.05 -m {neg/pos}`

##########################

mz_retention_untargeted <- raw_data_untargeted %>%
  group_by(Metabolite, Mode) %>%
  slice(1) %>%
  dplyr::select(Metabolite, Mode, `m/z`, `Retention time (min)`)



######## AD ----------------------------------------------------------

# add ad indicator
untargeted_all_amelia5_ad_ind <- imputed_all_untargeted5 %>%
  purrr::map(~cbind(.x[[1]], "AD_ind" = ifelse(imputed_all_untargeted5[[1]][[2]] == "AD", 1, 0)) %>% list())

ad_c_index <- which(imputed_all_untargeted5[[1]][[2]] %in% c('AD', 'CO', 'CY', 'CM'))

untargeted_univar_ad_logistic_list <- purrr::map(1:5, ~bh_univariate_age(untargeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = ad_c_index))

untargeted_univar_ad_logistic <- untargeted_univar_ad_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()

untargeted_univar_ad_logistic <- readRDS(here("ad_pd", "untargeted_univar_ad_logistic.Rds"))
# saveRDS(untargeted_univar_ad_logistic, here("ad_pd", "untargeted_univar_ad_logistic.Rds"))

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_ad_untargeted <- untargeted_univar_ad_logistic %>% 
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2) %>%
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  # dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 'z-score' = `z value`, Metabolite, Mode)

mummichog_ad_untargeted %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("ad_pd", 'mummichog_ad_neg.txt'))

mummichog_ad_untargeted %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("ad_pd", 'mummichog_ad_pos.txt'))


mummichog_ad_pos_plot <- mummichog_plot(here("ad_pd", "mummichog_ad_pos", "tables", "mcg_pathwayanalysis_mummichog_ad_pos.tsv")) +
  labs(title = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(1.75)))
ggsave("mummichog_ad_pos.png",
       plot = mummichog_ad_pos_plot, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 14)



######## PD ----------------------------------------------------------

# add pd indicator
untargeted_all_amelia5_pd_ind <- imputed_all_untargeted5 %>%
  purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(imputed_all_untargeted5[[1]][[2]] == "PD", 1, 0)) %>% list())

pd_c_index <- which(imputed_all_untargeted5[[1]][[2]] %in% c('PD', 'CO', 'CY', 'CM'))

untargeted_univar_pd_logistic_list <- purrr::map(1:5, ~bh_univariate_age(untargeted_all_amelia5_pd_ind, types_index = pd_c_index, var = "PD_ind", family = "binomial", conc = FALSE, imp_num = .x))

untargeted_univar_pd_logistic <- untargeted_univar_pd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()

untargeted_univar_pd_logistic <- readRDS(here("ad_pd", "untargeted_univar_pd_logistic.Rds"))
# saveRDS(untargeted_univar_pd_logistic, here("ad_pd", "untargeted_univar_pd_logistic.Rds"))

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_pd_untargeted <- untargeted_univar_pd_logistic %>% 
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2) %>%
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  #dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 'z-score' = `z value`, Metabolite, Mode)

mummichog_pd_untargeted %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("ad_pd", 'mummichog_pd_neg.txt'))

mummichog_pd_untargeted %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("ad_pd", 'mummichog_pd_pos.txt'))

# pretty much null i think
mummichog_pd_neg_plot <- mummichog_plot(here("ad_pd", "mummichog_pd_neg", "tables", "mcg_pathwayanalysis_mummichog_pd_neg.tsv")) +
  labs(title = "Negative Mode")
ggsave("mummichog_pd_neg.png",
       plot = mummichog_pd_neg_plot, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 18)

mummichog_pd_pos_plot <- mummichog_plot(here("ad_pd", "mummichog_pd_pos", "tables", "mcg_pathwayanalysis_mummichog_pd_pos.tsv")) +
  labs(title = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(1.75)))
ggsave("mummichog_pd_pos.png",
       plot = mummichog_pd_pos_plot, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 14)

######## PD vs AD ----------------------------------------------------------


pd_ad_index <- which(imputed_all_untargeted5[[1]][[2]] %in% c('PD', 'AD'))

untargeted_univar_adpd_logistic_list <- purrr::map(1:5, ~bh_univariate_age(untargeted_all_amelia5_pd_ind, types_index = pd_ad_index, var = "PD_ind", family = "binomial", conc = FALSE, imp_num = .x))

untargeted_univar_adpd_logistic <- untargeted_univar_adpd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()

untargeted_univar_adpd_logistic <- readRDS(here("ad_pd", "untargeted_univar_adpd_logistic.Rds"))
# saveRDS(untargeted_univar_adpd_logistic, here("ad_pd", "untargeted_univar_adpd_logistic.Rds"))

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_adpd_untargeted <- untargeted_univar_adpd_logistic %>% 
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2) %>%
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  #dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 'z-score' = `z value`, Metabolite, Mode)

mummichog_adpd_untargeted %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("ad_pd", 'mummichog_adpd_neg.txt'))

mummichog_adpd_untargeted %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("ad_pd", 'mummichog_adpd_pos.txt'))

# pretty much null i think
mummichog_adpd_neg_plot <- mummichog_plot(here("ad_pd", "mummichog_adpd_neg", "tables", "mcg_pathwayanalysis_mummichog_adpd_neg.tsv")) +
  labs(title = "Negative Mode")
ggsave("mummichog_adpd_neg.png",
       plot = mummichog_adpd_neg_plot, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 18)

mummichog_adpd_pos_plot <- mummichog_plot(here("ad_pd", "mummichog_adpd_pos", "tables", "mcg_pathwayanalysis_mummichog_adpd_pos.tsv")) +
  labs(title = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(1.75)))
ggsave("mummichog_adpd_pos.png",
       plot = mummichog_adpd_pos_plot, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 14)


######## PD missingness indicator ----------------------------------------------------------

# add pd indicator
untargeted_all_amelia5_pd_naind <- all_untargeted_naind %>%
  purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(all_untargeted_naind[[1]][[2]] == "PD", 1, 0)) %>% list())

pd_c_na_index <- which(all_untargeted_naind[[1]][[2]] %in% c('PD', 'CO', 'CY', 'CM'))
untargeted_univar_pd_logistic_naind <- bh_univariate_age(untargeted_all_amelia5_pd_naind, types_index = pd_c_na_index, var = "PD_ind", family = "binomial", conc = FALSE, imp_num = 1, scale =F)

untargeted_univar_pd_logistic_naind <- readRDS(here("ad_pd", "untargeted_univar_pd_logistic_naind.Rds"))
# saveRDS(untargeted_univar_pd_logistic_naind, here("ad_pd", "untargeted_univar_pd_logistic_naind.Rds"))

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_pd_untargeted_naind <- untargeted_univar_pd_logistic_naind %>% 
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2) %>%
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  # dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 'z-score' = `z value`, Metabolite, Mode)

mummichog_pd_untargeted_naind %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("ad_pd", 'mummichog_pd_neg_naind.txt'))

mummichog_pd_untargeted_naind %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("ad_pd", 'mummichog_pd_pos_naind.txt'))

# # pretty much null i think
# mummichog_pd_neg_plot_naind <- mummichog_plot(here("ad_pd", "mummichog_pd_neg", "tables", "mcg_pathwayanalysis_mummichog_pd_neg.tsv")) +
#   labs(title = "Negative Mode")
# ggsave("mummichog_pd_neg.png",
#        plot = mummichog_pd_neg_plot, 
#        path = here("plots", "adpd_figs"),
#        width = 20,
#        height = 18)

mummichog_pd_pos_plot_naind <- mummichog_plot(here("ad_pd", "mummichog_pd_pos_naind", "tables", "mcg_pathwayanalysis_mummichog_pd_pos_naind.tsv")) +
  labs(title = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(1.75)))
ggsave("mummichog_pd_pos_naind.png",
       plot = mummichog_pd_pos_plot_naind, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 14)

######## AD missingness indicator ----------------------------------------------------------

# add ad indicator
untargeted_all_amelia5_ad_naind <- all_untargeted_naind %>%
  purrr::map(~cbind(.x[[1]], "AD_ind" = ifelse(all_untargeted_naind[[1]][[2]] == "AD", 1, 0)) %>% list())

ad_c_na_index <- which(all_untargeted_naind[[1]][[2]] %in% c('AD', 'CO', 'CY', 'CM'))
untargeted_univar_ad_logistic_naind <- bh_univariate_age(untargeted_all_amelia5_ad_naind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = 1, scale =F, types_index <- ad_c_na_index)

untargeted_univar_ad_logistic_naind <- readRDS(here("ad_pd", "untargeted_univar_ad_logistic_naind.Rds"))
# saveRDS(untargeted_univar_ad_logistic_naind, here("ad_pd", "untargeted_univar_ad_logistic_naind.Rds"))

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_ad_untargeted_naind <- untargeted_univar_ad_logistic_naind %>% 
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2) %>%
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  # dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 'z-score' = `z value`, Metabolite, Mode)

# both mummichog results are null.
mummichog_ad_untargeted_naind %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("ad_pd", 'mummichog_ad_neg_naind.txt'))

mummichog_ad_untargeted_naind %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("ad_pd", 'mummichog_ad_pos_naind.txt'))

mummichog_ad_pos_plot_naind <- mummichog_plot(here("ad_pd", "mummichog_ad_pos_naind", "tables", "mcg_pathwayanalysis_mummichog_ad_pos_naind.tsv")) +
  labs(title = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(1.75)))
ggsave("mummichog_adpd_pos_naind.png",
       plot = mummichog_adpd_pos_plot_naind, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 14)
######## AD/PD (classify ad against pd) missingness indicator ----------------------------------------------------------

ad_pd_na_index <- which(all_untargeted_naind[[1]][[2]] %in% c('AD', 'PD'))

untargeted_univar_adpd_logistic_naind <- bh_univariate_age(untargeted_all_amelia5_adpd_naind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = 1, scale =F, types_index <- ad_pd_na_index)

untargeted_univar_adpd_logistic_naind <- readRDS(here("ad_pd", "untargeted_univar_adpd_logistic_naind.Rds"))
# saveRDS(untargeted_univar_adpd_logistic_naind, here("ad_pd", "untargeted_univar_adpd_logistic_naind.Rds"))

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_adpd_untargeted_naind <- untargeted_univar_adpd_logistic_naind %>% 
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2) %>%
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  # dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 'z-score' = `z value`, Metabolite, Mode)

# both mummichog results are null.
mummichog_adpd_untargeted_naind %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("ad_pd", 'mummichog_adpd_neg_naind.txt'))

mummichog_adpd_untargeted_naind %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("ad_pd", 'mummichog_adpd_pos_naind.txt'))

# pretty much null i think
mummichog_adpd_neg_plot_naind <- mummichog_plot(here("ad_pd", "mummichog_adpd_neg", "tables", "mcg_pathwayanalysis_mummichog_adpd_neg.tsv")) +
  labs(title = "Negative Mode")
ggsave("mummichog_adpd_neg.png",
       plot = mummichog_adpd_neg_plot,
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 18)

mummichog_adpd_pos_plot_naind <- mummichog_plot(here("ad_pd", "mummichog_adpd_pos_naind", "tables", "mcg_pathwayanalysis_mummichog_adpd_pos_naind.tsv")) +
  labs(title = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(1.75)))
ggsave("mummichog_adpd_pos_naind.png",
       plot = mummichog_adpd_pos_plot_naind, 
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 14)

############################

### Univariate AD/PD Logisitic Regression on targeted ###
### for use with MSEA

############################

message("Targeted Univariate ADPD Logistic -------------------------------------------")

#### AD First -----------------------------
# we want an AD indicator, but not a PD one


targeted_all_amelia5_ad_ind <- imputed_all_targeted5 %>%
  purrr::map(~cbind(.x[[1]], "AD_ind" = ifelse(imputed_all_targeted5[[1]][[2]] == "AD", 1, 0)) %>% list())

ad_c_index_targeted <- which(imputed_all_targeted5[[1]][[2]] %in% c('AD', "CO", "CM", "CY"))

# Create table with bh-corrected p values
targeted_univar_ad_logistic_list <- purrr::map(1:5, ~bh_univariate_age(targeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = ad_c_index_targeted))
targeted_univar_ad_logistic <- targeted_univar_ad_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()

targeted_univar_ad_logistic <- readRDS(here("ad_pd", "targeted_univar_ad_logistic.Rds"))
# saveRDS(targeted_univar_ad_logistic, here("ad_pd", "targeted_univar_ad_logistic.Rds"))


targeted_univar_ad_logistic %>%
  filter(bh_p_value < 0.05)

msea_ad_table <-targeted_univar_ad_logistic %>%
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
           TRUE ~ name))

# the list of significant ones -- basically none
univar_ad_names_sig <- msea_ad_table %>%
  filter(bh_p_value < 0.05) %>%
  select(mapped_name)
write_tsv(univar_ad_names_sig, here("ad_pd", "msea_ad_names_sig.txt"))




#### PD -------------------------------------
# # we want a PD indicator, but not an AD one
# targeted_all_amelia5_pd_ind <- purrr::map2(imputed_c_targeted5, targeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE,
# add_AD_ind = FALSE, add_PD_ind = TRUE))

# add ad_indicator. note that we have to add the extra "list()" at the end because this removes metadata, and bh_univariate_age expects a list of lists
# targeted_all_amelia5_pd_ind <- imputed_less10perct_targeted %>%
#   purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(imputed_less10perct_targeted[[1]][[2]] == "PD", 1, 0)) %>% list())
targeted_all_amelia5_pd_ind <- imputed_all_targeted5 %>%
  purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(imputed_all_targeted5[[1]][[2]] == "PD", 1, 0)) %>% list())

pd_c_index_targeted <- which(imputed_all_targeted5[[1]][[2]] %in% c('PD', "CO", "CM", "CY"))
# Create table with bh-corrected p values
targeted_univar_pd_logistic_list <- purrr::map(1:5, ~bh_univariate_age(targeted_all_amelia5_pd_ind, var = "PD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = pd_c_targeted_index))

# for each metabolite, use the smallest p value out of the 5 imputations. this is very generous
targeted_univar_pd_logistic <- targeted_univar_pd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  slice(1) %>%
  ungroup()

targeted_univar_pd_logistic <- readRDS(here("ad_pd", "targeted_univar_pd_logistic.Rds"))
# saveRDS(targeted_univar_pd_logistic, here("ad_pd", "targeted_univar_pd_logistic.Rds"))

targeted_univar_pd_logistic %>%
  filter(bh_p_value < 0.05)



msea_pd_table <-targeted_univar_pd_logistic %>%
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
           TRUE ~ name))

# the list of significant ones
univar_pd_names_sig <- msea_pd_table %>%
  filter(bh_p_value < 0.05) %>%
  select(mapped_name)
write_tsv(univar_pd_names_sig, here("ad_pd", "msea_pd_names_sig.txt"))


#### ADPD classification -----------------------------
# classify AD against PD (ad is positive class)

ad_pd_index_targeted <- which(imputed_all_targeted5[[1]][[2]] %in% c('AD', "PD"))

# Create table with bh-corrected p values
targeted_univar_adpd_logistic_list <- purrr::map(1:5, ~bh_univariate_age(targeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = ad_pd_index_targeted))
targeted_univar_adpd_logistic <- targeted_univar_adpd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()

targeted_univar_adpd_logistic <- readRDS(here("ad_pd", "targeted_univar_adpd_logistic.Rds"))
# saveRDS(targeted_univar_adpd_logistic, here("ad_pd", "targeted_univar_adpd_logistic.Rds"))


targeted_univar_adpd_logistic %>%
  filter(bh_p_value < 0.05)

msea_adpd_table <-targeted_univar_adpd_logistic %>%
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
           TRUE ~ name))

# the list of significant ones 
univar_adpd_names_sig <- msea_adpd_table %>%
  filter(bh_p_value < 0.05) %>%
  select(mapped_name)
write_tsv(univar_ad_names_sig, here("ad_pd", "msea_adpd_names_sig.txt"))
##


## We need a reference list (ie include names on insignificant variables). can use either ad or pd
univar_targeted_names_ref <- msea_pd_table %>%
  select(mapped_name)
write_tsv(univar_targeted_names_ref, here("ad_pd", "msea_names_ref.txt"))


## run through MSEA, get an output csv, and read it in here to plot

# CSF disease library ---
msea_pd_csf <- read_csv(here("ad_pd", "msea_ora_csf_pd_result.csv")) %>%
  mutate(bh_p = p.adjust(`Raw p`, method = "BH"),
         set = fct_reorder(X1, hits/total)) %>%
  filter(total > 2) %>%
  top_n(10, hits/total) %>%
  select("Set" = X1, "Total" = total, "Expected" = expected, "Hits" = hits) %>%
  arrange(desc(Hits/Total))

write_csv(msea_pd_csf, here("ad_pd", "msea_csf_pd_result_top10.csv"))

# SMPDB library -----
msea_pd_smpdb <- read_csv(here("ad_pd", "msea_ora_smpdb_pd_result.csv")) %>%
  mutate(bh_p = p.adjust(`Raw p`, method = "BH"),
         set = fct_reorder(X1, hits/total)) %>%
  filter(total > 2) %>%
  top_n(10, hits/total) %>%
  select("Set" = X1, "Total" = total, "Expected" = expected, "Hits" = hits) %>%
  arrange(desc(Hits/Total))

write_csv(msea_pd_smpdb, here("ad_pd", "msea_smpdb_pd_result_top10.csv"))





#####################

### Knockoffs

##################

my_glmnet_stat_binom <- function(X, X_k, y){
  stat.glmnet_coefdiff(X, X_k, y, nlambda = 200, alpha = 0.5, family = "binomial")
}

my_lasso_stat_binom <- function(X, X_k, y){
  stat.glmnet_coefdiff(X, X_k, y, nlambda = 200, family = "binomial")
}


# ad logistic

ad_c_index <- which(imputed_all_targeted5[[1]]$Type != "PD")

targeted_ad_c_x <- imputed_all_targeted5[[1]]$Y[ad_c_index,-which(colnames(imputed_all_targeted5[[1]]$Y) %in% c("(Intercept)"))]
targeted_ad_c_y <- fct_collapse(fct_drop(imputed_all_targeted5[[1]]$Type[ad_c_index]), "C" = c("CO", "CM", "CY"))
# try with glmnet
knockoff_tar <- knockoff.filter(targeted_ad_c_x, targeted_ad_c_y, knockoffs = create.second_order, statistic = my_glmnet_stat_binom, fdr = 0.05)
# try with lasso
#knockoff_tar <- knockoff.filter(targeted_ad_c_x, targeted_ad_c_y,knockoffs = create.second_order, statistic = my_lasso_stat_binom, fdr = 0.05)

print(knockoff_tar$selected)


# pd logistic

pd_c_index <- which(imputed_all_targeted5[[1]]$Type != "AD")

targeted_pd_c_x <- imputed_all_targeted5[[1]]$Y[pd_c_index,-which(colnames(imputed_all_targeted5[[1]]$Y) %in% c("(Intercept)"))]
targeted_pd_c_y <- fct_collapse(fct_drop(imputed_all_targeted5[[1]]$Type[pd_c_index]), "C" = c("CO", "CM", "CY"))
knockoff_tar <- knockoff.filter(targeted_pd_c_x, targeted_pd_c_y, knockoffs = create.second_order, statistic = my_glmnet_stat_binom, fdr = 0.05)
print(knockoff_tar$selected)




##################

### MSEA for AD/PD

#################

message("MSEA for ad/pd significant (targeted) -------------------------------------------")

## AD ----------------------------------
# using bh p < -0.05
ad_sig_names <- targeted_univar_ad_logistic %>%
  filter(bh_p_value < 0.05) %>%
  select(name) %>%
  deframe() %>%
  str_replace_all("_neg", "") %>%
  str_replace_all("_pos", "")



mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Setup.MapData(mSet, ad_sig_names)

# upload our reference list (the metabolites we targeted)
mSet<-Setup.HMDBReferenceMetabolome(mSet, "targeted_names_hmdb_list.txt");
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)



#what does this do?? True is the only option in the web-tool version, so i'm keeping it as default
# but it makes a pretty big difference
mSet<-SetMetabolomeFilter(mSet, T)

# Select metabolite set library

#pathway associated library?
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway",2)

# csf associated library?
#mSet<-SetCurrentMsetLib(mSet, "csf", 2)

# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ad_ora_0_", "bar", "png", 72, width=NA)


## PD ----------------------------------
# using bh p < -0.05
pd_sig_names <- targeted_univar_pd_logistic %>%
  filter(bh_p_value < 0.05) %>%
  arrange(bh_p_value) %>%
  select(name) %>%
  deframe() %>%
  str_replace_all("_neg", "") %>%
  str_replace_all("_pos", "")



mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Setup.MapData(mSet, pd_sig_names)

# upload our reference list (the metabolites we targeted)
mSet<-Setup.HMDBReferenceMetabolome(mSet, "targeted_names_hmdb_list.txt");
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)

## try to match the NAs manually... 

mSet<-PerformDetailMatch(mSet, "6-Methyl-DL-tryptophan")
mSet <- GetCandidateList(mSet)

# try the other one.
mSet<-PerformDetailMatch(mSet, "N-Acetylethanolamine")
mSet <- GetCandidateList(mSet)


# results: matches are found, but neither of them look right.

##---

#what does this do?? True is the only option in the web-tool version, so i'm keeping it as default
# but it makes a pretty big difference
mSet<-SetMetabolomeFilter(mSet, T)

# Select metabolite set library

#pathway associated library?
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway",2)

# csf associated library?
mSet<-SetCurrentMsetLib(mSet, "csf", 2)

# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "pd_ora_0_", "bar", "png", 72, width=NA)






###########

## (within pd) GBA Logistic

##########

# quick exploration of GBA carriers
wide_data_untargeted %>% 
  mutate(GBAStatus = fct_collapse(GBAStatus, 
                                  Carrier = c('E326K Carrier', 'Pathogenic Carrier', 'CT'))) %>%
  filter(Type == "PD") %>%
  count(GBAStatus)#, Gender)

wide_data_untargeted %>% 
  mutate(GBAStatus = fct_collapse(GBAStatus, 
                                  Carrier = c('E326K Carrier', 'Pathogenic Carrier', 'CT'))) %>% 
  filter(Type == "PD") %>%
  ggplot() + 
  geom_histogram(aes(Age)) +
  facet_wrap(~Gender)

wide_data_untargeted %>% 
  mutate(GBAStatus = fct_collapse(GBAStatus, 
                                  Carrier = c('E326K Carrier', 'Pathogenic Carrier', 'CT'))) %>% 
  filter(Type == "PD") %>%
  ggplot() + 
  geom_histogram(aes(Age)) +
  facet_wrap(~GBAStatus)





# isolating just the PD subjects
imputed_all_pd_index <- which(imputed_all_untargeted5[[1]][[2]] == "PD")
imputed_all_pd_features <- purrr::map(imputed_all_untargeted5, ~.x[[1]][imputed_all_pd_index,])
imputed_all_pd_metadata <- purrr::map(1:5, function(x) purrr::map(2:length(imputed_all_untargeted5), function(y) imputed_all_untargeted5[[x]][[y]][imputed_all_pd_index]))

# get data in form of filter_and_impute_multi output
imputed_all_pd <- purrr::map(1:5, ~append(imputed_all_pd_metadata[[.x]], list(imputed_all_pd_features[[.x]]), after = 0))


untargeted_gba_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_pd, varname ="GBA Status", imp_num = .x, nlambda = 200, GBA_ind = T))
untargeted_gba_logistic <- readRDS(here("ad_pd", "untargeted_gba_logistic.Rds"))
# saveRDS(untargeted_gba_logistic, here("ad_pd", "untargeted_gba_logistic.Rds"))


untargeted_gba_roc_table <- untargeted_gba_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

untargeted_gba_roc_plot <- ggplot(untargeted_gba_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: GBA Carriers vs Non-Carriers",
       subtitle = TeX('Untargeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(untargeted_gba_roc_table$auc), 3)))

ggsave("roc_gba_untargeted.png",
       plot = untargeted_gba_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


# untargeted_gba_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# untargeted_gba_logistic[[1]][[2]] + 
#   labs(title = "ROC: GBA carriers vs non-carriers")
# ggsave("gba_logistic_untargeted.png")


untargeted_gba_avg_retained <- untargeted_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_gba_in_all <- untargeted_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")


########

## GBA logistic, targeted
#######
# # isolating just the PD subjects
# targeted_imputed_all_pd_index <- which(imputed_all_targeted5[[1]][[2]] == "PD")
# targeted_imputed_all_pd_features <- purrr::map(imputed_all_targeted5, ~.x[[1]][targeted_imputed_all_pd_index,])
# targeted_imputed_all_pd_metadata <- purrr::map(1:5, function(x) purrr::map(2:length(imputed_all_targeted5), function(y) imputed_all_targeted5[[x]][[y]][targeted_imputed_all_pd_index]))
# 
# # get data in form of filter_and_impute_multi output
# targeted_imputed_all_pd <- purrr::map(1:5, ~append(targeted_imputed_all_pd_metadata[[.x]], list(targeted_imputed_all_pd_features[[.x]]), after = 0))

# do a new imputation on just the PD subjects
targeted_imputed_pd <- filter_and_impute_multi(wide_data_targeted, types = c("PD"), empri = 5)

targeted_gba_logistic <- purrr::map(1:5, ~logistic_control_analysis(targeted_imputed_pd, varname ="GBA Status", imp_num = .x, nlambda = 200, GBA_ind = T))
targeted_gba_logistic <- readRDS(here("ad_pd", "targeted_gba_logistic.Rds"))
# saveRDS(targeted_gba_logistic, here("ad_pd", "targeted_gba_logistic.Rds"))


targeted_gba_roc_table <- targeted_gba_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

targeted_gba_roc_plot <- ggplot(targeted_gba_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: GBA Carriers vs Non-Carriers",
       subtitle = TeX('Targeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(targeted_gba_roc_table$auc), 3)))

ggsave("roc_gba_targeted.png",
       plot = targeted_gba_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)



# targeted_gba_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# targeted_gba_logistic[[1]][[2]] + 
#   labs(title = "ROC: GBA carriers vs non-carriers",
#        subtitle = TeX('Targeted ,$\\alpha = 0.5$'))
# ggsave("gba_logistic_targeted.png")


targeted_gba_coefs <- targeted_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(name = str_replace_all(name, "_pos|_neg", ""),
         mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
         median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef))) %>%
  select("Name" = name, "Avg Coef" = mean_coef) %>%
  ungroup() %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM")))

# write separate tables for positive and negative coefs
targeted_gba_coefs %>%
  filter(`Avg Coef` > 0) %>%
  write_csv(here("ad_pd", "pos_coef_table_targeted_gba.csv"))

targeted_gba_coefs %>%
  filter(`Avg Coef` < 0) %>%
  write_csv(here("ad_pd", "neg_coef_table_targeted_gba.csv"))



targeted_gba_avg_retained <- targeted_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

targeted_gba_in_all <- targeted_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")



##### 
### MSEA, GBA
#####

targeted_all_amelia5_gba_ind <- imputed_all_targeted5 %>%
  purrr::map(~cbind(.x[[1]], "GBA_ind" = ifelse(imputed_all_targeted5[[1]][[5]] %in% c("Pathogenic Carrier", "CT", "E326K Carrier"), 1, 0)) %>% list())

# Create table with bh-corrected p values
targeted_univar_gba_logistic <- bh_univariate_age(targeted_all_amelia5_gba_ind, var = "GBA_ind", family = "binomial", conc = FALSE)
targeted_univar_gba_logistic <- readRDS(here("ad_pd", "targeted_univar_gba_logistic.Rds"))
# saveRDS(targeted_univar_gba_logistic, here("ad_pd", "targeted_univar_gba_logistic.Rds"))


# using bh p < -0.05
gba_sig_names <- targeted_univar_gba_logistic %>%
  filter(bh_p_value < 0.05) %>%
  arrange(bh_p_value) %>%
  select(name) %>%
  deframe() %>%
  str_replace_all("_neg", "") %>%
  str_replace_all("_pos", "")



mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Setup.MapData(mSet, pd_sig_names)

# upload our reference list (the metabolites we targeted)
mSet<-Setup.HMDBReferenceMetabolome(mSet, "targeted_names_hmdb_list.txt");
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)

## try to match the NAs manually... 

mSet<-PerformDetailMatch(mSet, "6-Methyl-DL-tryptophan")
mSet <- GetCandidateList(mSet)

# try the other one.
mSet<-PerformDetailMatch(mSet, "N-Acetylethanolamine")
mSet <- GetCandidateList(mSet)


# results: matches are found, but neither of them look right.

##---

#what does this do?? True is the only option in the web-tool version, so i'm keeping it as default
# but it makes a pretty big difference
mSet<-SetMetabolomeFilter(mSet, T)

# Select metabolite set library

#pathway associated library?
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway",2)

# csf associated library?
mSet<-SetCurrentMsetLib(mSet, "csf", 2)

# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "GBA_ora_0_", "bar", "png", 72, width=NA)


########

## GBA logistic, Lipids
#######
# isolating just the PD subjects
lipids_imputed_pd_index <- which(imputed_all_lipids5[[1]][[2]] == "PD")
lipids_imputed_pd_features <- purrr::map(imputed_all_lipids5, ~.x[[1]][lipids_imputed_pd_index,])
lipids_imputed_pd_metadata <- purrr::map(1:5, function(x) purrr::map(2:length(imputed_all_lipids5), function(y) imputed_all_lipids5[[x]][[y]][lipids_imputed_pd_index]))

# get data in form of filter_and_impute_multi output
lipids_imputed_pd <- purrr::map(1:5, ~append(lipids_imputed_pd_metadata[[.x]], list(lipids_imputed_pd_features[[.x]]), after = 0))


lipids_gba_logistic <- purrr::map(1:5, ~logistic_control_analysis(lipids_imputed_pd, varname ="GBA Status", imp_num = .x, nlambda = 200, GBA_ind = T))
lipids_gba_logistic <- readRDS(here("ad_pd", "lipids_gba_logistic.Rds"))
# saveRDS(lipids_gba_logistic, here("ad_pd", "lipids_gba_logistic.Rds"))


lipids_gba_roc_table <- lipids_gba_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

lipids_gba_roc_plot <- ggplot(lipids_gba_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: GBA Carriers vs Non-Carriers",
       subtitle = TeX('Lipids'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(lipids_gba_roc_table$auc), 3)))

ggsave("roc_gba_lipids.png",
       plot = lipids_gba_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)

# lipids_gba_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# lipids_gba_logistic[[1]][[2]] + 
#   labs(title = "ROC: GBA carriers vs non-carriers",
#        subtitle = TeX('Lipids ,$\\alpha = 0.5$'))
# ggsave("gba_logistic_lipids.png")


lipids_gba_coefs <- lipids_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(name = str_replace_all(name, "_pos|_neg", ""),
         mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5)) %>% round(2),
         median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef))) %>%
  select("Name" = name, "Avg Coef" = mean_coef) %>%
  ungroup() %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM")))

# write separate tables for positive and negative coefs
lipids_gba_coefs %>%
  filter(`Avg Coef` > 0) %>%
  write_csv(here("ad_pd", "pos_coef_table_lipids_gba.csv"))

lipids_gba_coefs %>%
  filter(`Avg Coef` < 0) %>%
  write_csv(here("ad_pd", "neg_coef_table_lipids_gba.csv"))

lipids_gba_avg_retained <- lipids_gba_logistic %>% 
   purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_gba_in_all <- lipids_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")


########

## GBA logistic, got
#######
# isolating just the PD subjects
got_imputed_pd_index <- which(imputed_all_got5[[1]][[2]] == "PD")
got_imputed_pd_features <- purrr::map(imputed_all_got5, ~.x[[1]][got_imputed_pd_index,])
got_imputed_pd_metadata <- purrr::map(1:5, function(x) purrr::map(2:length(imputed_all_got5), function(y) imputed_all_got5[[x]][[y]][got_imputed_pd_index]))

# get data in form of filter_and_impute_multi output
got_imputed_pd <- purrr::map(1:5, ~append(got_imputed_pd_metadata[[.x]], list(got_imputed_pd_features[[.x]]), after = 0))


got_gba_logistic <- purrr::map(1:5, ~logistic_control_analysis(got_imputed_pd, varname ="GBA Status", imp_num = .x, nlambda = 200, GBA_ind = T))
got_gba_logistic <- readRDS(here("ad_pd", "got_gba_logistic.Rds"))
# saveRDS(got_gba_logistic, here("ad_pd", "got_gba_logistic.Rds"))


got_gba_roc_table <- got_gba_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

got_gba_roc_plot <- ggplot(got_gba_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: GBA Carriers vs Non-Carriers",
       subtitle = TeX('GOT-MS'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(got_gba_roc_table$auc), 3)))

ggsave("roc_gba_got.png",
       plot = got_gba_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)

# got_gba_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# got_gba_logistic[[1]][[2]] + 
#   labs(title = "ROC: GBA carriers vs non-carriers",
#        subtitle = TeX('got ,$\\alpha = 0.5$'))
# ggsave("gba_logistic_got.png")


got_gba_coefs <- got_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(median_coef)))


got_gba_avg_retained <- got_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

got_gba_in_all <- got_gba_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

##################

### APOE4 logistic, untargeted

#################

# quick idea of what the counts looks like
wide_data_untargeted %>% 
  mutate(APOE4 = APOE %>% 
           as.character %>% 
           str_detect("4"),
         Type = fct_collapse(Type, "Control" = c("CO", "CM", "CY"))) %>% 
  count(APOE4, Type, sort = T) %>%
  ggplot() +
  geom_col(aes(Type, n, fill = APOE4), position = "dodge")

untargeted_apoe_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted5, varname ="APOE_ind", imp_num = .x, nlambda = 200, APOE4_ind = T))
untargeted_apoe_logistic <- readRDS(here("ad_pd", "untargeted_apoe_logistic.Rds"))
# saveRDS(untargeted_apoe_logistic, here("ad_pd", "untargeted_apoe_logistic.Rds"))


untargeted_apoe_roc_table <- untargeted_apoe_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

untargeted_apoe_roc_plot <- ggplot(untargeted_apoe_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: APOE4",
       subtitle = TeX('Untargeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(untargeted_apoe_roc_table$auc), 3)))

ggsave("roc_apoe_untargeted.png",
       plot = untargeted_apoe_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)

# 
# untargeted_apoe_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# untargeted_apoe_logistic[[1]][[2]] + 
#   labs(title = "ROC: APOE4",
#        subtitle = "Targeted: Controls, AD, PD",
#        size = 12)
# ggsave("apoe_logistic_untargeted.png")


# diagnostic plots
heatmap.fit(untargeted_apoe_logistic[[1]]$truth, pred = untargeted_apoe_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(untargeted_apoe_logistic[[1]]$truth, pred = untargeted_apoe_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = untargeted_apoe_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "apoe (untargeted) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

untargeted_apoe_avg_retained <- untargeted_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_apoe_in_all <- untargeted_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

#relationship between retained apoe and retained age . None
intersect(untargeted_in_all, untargeted_apoe_in_all)



##################

### APOE4 logistic, Targeted

#################

targeted_apoe_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="APOE_ind", imp_num = .x, nlambda = 200, APOE4_ind = T))

targeted_apoe_logistic <- readRDS(here("ad_pd", "targeted_apoe_logistic.Rds"))
# saveRDS(targeted_apoe_logistic, here("ad_pd", "targeted_apoe_logistic.Rds"))


targeted_apoe_roc_table <- targeted_apoe_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

targeted_apoe_roc_plot <- ggplot(targeted_apoe_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: APOE4",
       subtitle = TeX('Targeted'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(targeted_apoe_roc_table$auc), 3)))

ggsave("roc_apoe_targeted.png",
       plot = targeted_apoe_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)


# targeted_apoe_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)
# 
# targeted_apoe_logistic[[1]][[2]] + 
#   labs(title = "ROC: APOE4",
#        subtitle = "Targeted: Controls, AD, PD",
#        size = 12)
# ggsave("apoe_logistic_targeted.png")


# diagnostic plots
heatmap.fit(targeted_apoe_logistic[[1]]$truth, pred = targeted_apoe_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(targeted_apoe_logistic[[1]]$truth, pred = targeted_apoe_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = targeted_apoe_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "apoe (targeted) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

targeted_apoe_avg_retained <- targeted_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

targeted_apoe_in_all <- targeted_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

#relationship between retained apoe and retained age . None
intersect(targeted_in_all, targeted_apoe_in_all)


##################

### APOE4 logistic, LIPIDS

#################

lipids_apoe_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="APOE_ind", imp_num = .x, nlambda = 200, APOE4_ind = T))

lipids_apoe_logistic <- readRDS(here("ad_pd", "lipids_apoe_logistic.Rds"))
# saveRDS(lipids_apoe_logistic, here("ad_pd", "lipids_apoe_logistic.Rds"))


lipids_apoe_roc_table <- lipids_apoe_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

lipids_apoe_roc_plot <- ggplot(lipids_apoe_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: APOE4",
       subtitle = TeX('Lipids'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(lipids_apoe_roc_table$auc), 3)))

ggsave("roc_apoe_lipids.png",
       plot = lipids_apoe_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)



lipids_apoe_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

lipids_apoe_logistic[[1]][[2]] + 
  labs(title = "ROC: APOE4",
       subtitle = "Lipids: Controls, AD, PD",
       size = 12)
ggsave("apoe_logistic_lipids.png")


# diagnostic plots
heatmap.fit(lipids_apoe_logistic[[1]]$truth, pred = lipids_apoe_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(lipids_apoe_logistic[[1]]$truth, pred = lipids_apoe_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = lipids_apoe_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "apoe (lipids) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

lipids_apoe_avg_retained <- lipids_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_apoe_in_all <- lipids_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

#relationship between retained apoe and retained age . None
intersect(lipids_in_all, lipids_apoe_in_all)




##################

### APOE4 logistic, got

#################

got_apoe_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_got5, varname ="APOE_ind", imp_num = .x, nlambda = 200, APOE4_ind = T))

got_apoe_logistic <- readRDS(here("ad_pd", "got_apoe_logistic.Rds"))
# saveRDS(got_apoe_logistic, here("ad_pd", "got_apoe_logistic.Rds"))


got_apoe_roc_table <- got_apoe_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

got_apoe_roc_plot <- ggplot(got_apoe_roc_table) + 
  geom_line(mapping = aes(x, y, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: APOE4",
       subtitle = TeX('GOT-MS'),
       x = 'False Positive Rate',
       y = 'True Positive Rate') + 
  geom_text(x = Inf, y = -Inf, 
            hjust = 1, vjust = -0.5, 
            size = 20,# label.padding = unit(1, "lines"),
            label = paste0('AUC:', round(mean(got_apoe_roc_table$auc), 3)))

ggsave("roc_apoe_got.png",
       plot = got_apoe_roc_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)



got_apoe_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

got_apoe_logistic[[1]][[2]] + 
  labs(title = "ROC: APOE4",
       subtitle = "got: Controls, AD, PD",
       size = 12)
ggsave("apoe_logistic_got.png")


# diagnostic plots
heatmap.fit(got_apoe_logistic[[1]]$truth, pred = got_apoe_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(got_apoe_logistic[[1]]$truth, pred = got_apoe_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = got_apoe_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "apoe (got) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

got_apoe_avg_retained <- got_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

got_apoe_in_all <- got_apoe_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

#relationship between retained apoe and retained age . None
intersect(got_in_all, got_apoe_in_all)






############

### AD/PD logistic 
## Looking only at a subset of features to see how robust the predictions are

#########

### AD ----
nonsig_ad_names <- setdiff(untargeted_ad_in_all, c("Age", "GenderM")) %>% 
  setdiff(names(wide_data_untargeted), .)
untargeted_ad_logistic_nonsig <- adpd_subset_analysis(wide_data_untargeted, "AD", features = nonsig_ad_names, AD_ind = T)

untargeted_ad_logistic_nonsig %>%
  purrr::map(~.x[[2]] + labs(title = "ROC: AD vs C, nonsignificant")) %>%
  cowplot::plot_grid(plotlist = .)


### PD --- (removes all nonzero that showed up in all imputations above. 35 predictors)
nonsig_pd_names <- setdiff(untargeted_pd_in_all, c("Age", "GenderM")) %>% 
  setdiff(names(wide_data_untargeted), .)
untargeted_pd_logistic_nonsig <- adpd_subset_analysis(wide_data_untargeted, "PD nonsig ", features = nonsig_pd_names, PD_ind = T)

untargeted_pd_logistic_nonsig %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)


### PD --- (removes all nonzero that showed up in the nonsig_pd_names above)
untargeted_pd_in_all_round2 <- untargeted_pd_logistic_nonsig %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

nonsig_pd_names_round2 <- setdiff(untargeted_pd_in_all_round2, c("Age", "GenderM")) %>% 
  setdiff(nonsig_pd_names, .)
untargeted_pd_logistic_nonsig_round2 <- adpd_subset_analysis(wide_data_untargeted, "PD nonsig round 2", features = nonsig_pd_names_round2, PD_ind = T)

untargeted_pd_logistic_nonsig_round2 %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)


### PD --- (removes all nonzero that showed up in the nonsig_pd_names above) round 3
untargeted_pd_in_all_round3 <- untargeted_pd_logistic_nonsig_round2 %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

nonsig_pd_names_round3 <- setdiff(untargeted_pd_in_all_round3, c("Age", "GenderM")) %>% 
  setdiff(nonsig_pd_names_round2, .)
untargeted_pd_logistic_nonsig_round3 <- adpd_subset_analysis(wide_data_untargeted, "PD nonsig round 3", features = nonsig_pd_names_round3, PD_ind = T)

untargeted_pd_logistic_nonsig_round3 %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)



#########

## Predict missingness given AD/PD

#########


univar_missing_logistic <- function(metabolite){
  metabolite <- sym(metabolite)
  df <- wide_data_untargeted %>% 
    select(!!metabolite, Type) %>%
    mutate(Type = fct_collapse(Type, control = c("CY", "CO", "CM")) %>%
             droplevels() %>%
             fct_relevel("control", "AD", "PD"),
           missing_ind = ifelse(is.na(!!metabolite), 1, 0))
  
  fit <- glm(missing_ind ~ Type, data = df, family = "binomial")
  pred <- predict(fit, type = "response")
  
  fit_coefs <- fit$coefficients
  p_vals <- summary(fit)$coefficients[,"Pr(>|z|)"]
    
  results_df <- tibble(name = as.character(metabolite),
                       coef_intercept = fit_coefs["(Intercept)"],
                       coef_AD = fit_coefs["TypeAD"],
                       coef_PD = fit_coefs["TypePD"],
                       p_intercept = p_vals["(Intercept)"],
                       p_AD = p_vals["TypeAD"],
                       p_PD = p_vals["TypePD"],
                       truth = df$missing_ind,
                       pred = pred
                       )
  return(list(results_df, fit = fit)) 
}
untargeted_missing_ADPD <- wide_data_untargeted %>%
  names %>%
  setdiff(c("Age", "Gender", "APOE", "Type", "GBAStatus", "GBA_T369M", 
            "Id")) %>%
  purrr::map(~univar_missing_logistic(.x))
    
untargeted_missing_ADPD <- readRDS(here("ad_pd", "untargeted_missing_ADPD.Rds"))
# saveRDS(untargeted_missing_ADPD, here("ad_pd", "untargeted_missing_ADPD.Rds"))


untargeted_missing_adpd_sig_names <- untargeted_missing_ADPD %>% purrr::map(~.x[[1]]) %>%
  bind_rows() %>%
  distinct() %>%
  mutate(#name = purrr::map_chr(untargeted_missing_ADPD, ~.x[[1]] %>% as.character),
         bh_p_AD = p.adjust(p_AD, method = 'BH'),
         bh_p_PD = p.adjust(p_PD, method = "BH")) %>%
  filter(bh_p_AD < 0.05 | bh_p_PD < 0.05) %>%
  pull(name) %>%
  unique()

wide_data_untargeted_matched_c_ids <- purrr::map_df(1:nrow(wide_data_untargeted), ~find_control(.x, data = filter(wide_data_untargeted, Type %in% c("AD", "PD")), 
                                                                                            data_control = filter(wide_data_untargeted, Type %in% c("CO", "CY", "CM")))) %>%
  distinct(.keep_all = TRUE) %>%
  pull(Id)

missingness_by_type_untargeted_sig_pct <- wide_data_untargeted %>% 
  filter(Type %in% c("AD", "PD") | Id %in% wide_data_untargeted_matched_c_ids) %>%
  select(all_of(untargeted_missing_adpd_sig_names), Type) %>%
  mutate(Type = fct_collapse(Type, `Matched Controls` = c("CY", "CO", "CM")) %>%
           droplevels() %>%
           fct_relevel("Matched Controls", "AD", "PD")) %>%
  group_by(Type) %>%
  group_map(~ map_dbl(.x, function(y) round(sum(is.na(y))/length(y), digits = 3))) %>%
  set_names(c("Matched Controls", "AD", "PD")) %>%
  lapply(function(x) enframe(x, name = 'name', value = 'pct_missing')) %>%
  bind_rows(.id = "Type")

  


ggplot(missingness_by_type_untargeted_sig_pct, aes(pct_missing, name)) +
  geom_point(aes(color = Type), size = 3, position = position_jitter(width = 0.01, height = 0,seed = 1)) + 
  #scale_color_manual(labels = c("CY", "CM", "CO", "AD", "PD"), values = c("lightskyblue", "dodgerblue", "blue", "darkgreen", "purple")) +
  theme(axis.text.x = element_text(angle= 90, hjust =1)) + 
  labs(title = 'Percent Missingness by Type',
       subtitle = 'Untargeted, filtered using univar logistic regression, BH q < 0.05',
       x = "Percent Missing",
       y = "Name")  


untar_missingness_adpd_plot <- ggplot(missingness_by_type_untargeted_sig_pct, aes(fct_reorder(name, pct_missing), pct_missing)) +
  geom_col(aes(fill = fct_relevel(Type, "AD", after = 1)), position = "dodge", width = 0.75) + 
  #scale_color_manual(labels = c("CY", "CM", "CO", "AD", "PD"), values = c("lightskyblue", "dodgerblue", "blue", "darkgreen", "purple")) +
  #theme(axis.text.x = element_text(angle= 90, hjust =1)) + 
  theme(axis.text.x = element_blank()) + 
  labs(title = 'Percent Missingness by Type',
       subtitle = 'Untargeted, FDR < 0.05',
       y = "Percent Missing",
       x = "Name",
       fill = "Type")  

ggsave("untar_missingness_adpd.png",
       plot = untar_missingness_adpd_plot,
       path = here("plots", "adpd_figs"),
       width = 20,
       height = 10)

#### outlier analysis
# looking at outlier seen in the full control trained model fit on ADPD


#' For a given row of data, see which covariates in the row are more than 3 MAD away from their column
check_if_outlier <- function(row, data){
  # create emtpy dataframe to fill in
  outlier_df <- tibble(name = character(0), 
                       outlier_val = numeric(0),
                       col_abs_median = numeric(0))
  
  # for each feature
  for (i in 1:ncol(row)){
    # get the value for the outlier row (single value)
    # and the value for the whole dataframe (vector)
    row_value <- pull(row, i)
    data_col_value <- pull(data, i)
    # if the row is numeric and non-NA
    if (is.numeric(row_value) & !is.na(row_value)){
      # compute median +- 3 MAD
      upper_threshold <- median(data_col_value, na.rm = T) + 3*mad(data_col_value, na.rm =T)
      lower_threshold <- median(data_col_value, na.rm = T) - 3*mad(data_col_value, na.rm = T)
      # save the information about the covariate if the value is beyond 3 MAD
      if (abs(row_value) >= upper_threshold | abs(row_value) <= lower_threshold){
        outlier_df <- outlier_df %>% add_row(name = names(row)[i],
                                              outlier_val = row_value,
                                              col_abs_median = median(abs(data_col_value), na.rm = T))
      }
    }
  }
  outlier_df
}

# Looking at the outlier from the untargeted ADPD analysis (model train on full controls)
pathogenic_outlier_untargeted <- check_if_outlier(filter(wide_data_untargeted, Id == "PWA11-0318"),
                                                  data = wide_data_untargeted)

pathogenic_outlier_untargeted_PD <- check_if_outlier(filter(wide_data_untargeted, Id == "PWA11-0318"),
                                       data = wide_data_untargeted %>% filter(Type == "PD"))

pathogenic_outlier_untargeted_carrier <- check_if_outlier(filter(wide_data_untargeted, Id == "PWA11-0318"),
                                                  data = wide_data_untargeted %>% filter(GBAStatus %in% c("Pathogenic Carrier",)))



pathogenic_outlier_targeted <- check_if_outlier(row = filter(wide_data_targeted, Id == "PWA11-0318"),
                                       data = wide_data_targeted)

pathogenic_outlier_lipids <- check_if_outlier(row = filter(wide_data_lipids, Id == "PWA11-0318"),
                                                data = wide_data_lipids)



##########

# changing >50% NA columns to just missingness indicators

######

#' function to change column with lots of missingness to an indicator for missingness
#' @param col is a vector, a column in the data
#' @param na_thresh is the threshold at which we decide to change a column
#' ie default 0.5 says change a column to indicator if it has >.5 missingness
change_to_ind <- function(col, na_thresh = 0.5){
  # make sure to only do the operation if >50% NA
  if(sum(is.na(col)) / length(col) > na_thresh){
    ifelse(is.na(col), 0, 1)
  } else {
    col
  }
}

wide_data_untargeted_na_ind <- wide_data_untargeted %>%
  # convert 0 to NA since we aggreed 0s were a kind of missingness
  mutate_if(is.numeric, ~ifelse(.x == 0, NA, .x)) %>%
  # rename and change the columns that we're changing
  rename_if(function(x) is.numeric(x) & (sum(is.na(x)) / length(x)) > 0.5, ~paste0(.x, "_ind")) %>%
  mutate_if(is.numeric, ~change_to_ind(.x))


wide_data_lipids_na_ind <- wide_data_lipids %>%
  # convert 0 to NA since we aggreed 0s were a kind of missingness
  mutate_if(is.numeric, ~ifelse(.x == 0, NA, .x)) %>%
  # rename and change the columns that we're changing
  rename_if(function(x) is.numeric(x) & (sum(is.na(x)) / length(x)) > 0.5, ~paste0(.x, "_ind")) %>%
  mutate_if(is.numeric, ~change_to_ind(.x))



# now impute the rest of the missingness. note the extra replace_zeroes = F argument since we changed zeroes to be informative in the missing columns
imputed_all_untargeted_naind5 <- filter_and_impute_multi(wide_data_untargeted_na_ind, types = c('CO', 'CY', 'CM', "AD", "PD"), empri = 20, replace_zeroes = FALSE)
imputed_all_lipids_naind5 <- filter_and_impute_multi(wide_data_lipids_na_ind, types = c('CO', 'CY', 'CM', "AD", "PD"), empri = 200, replace_zeroes = FALSE)


### Untargeted, AD ----------------
untargeted_naind_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted_naind5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))
untargeted_naind_ad_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

untargeted_naind_ad_logistic[[1]][[2]] + 
  labs(title = "ROC: AD vs {Controls, AD}",
       subtitle = TeX('Untargeted with NA indicators (>.5 NA) ,$\\alpha = 0.5$'))
ggsave("ad_logistic_untargeted_naind.png")


# diagnostic plots
heatmap.fit(untargeted_naind_ad_logistic[[1]]$truth, pred = untargeted_naind_ad_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(untargeted_naind_ad_logistic[[1]]$truth, pred = untargeted_naind_ad_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = untargeted_naind_ad_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "AD (untargeted_naind) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

untargeted_naind_ad_coefs <- untargeted_naind_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(median_coef)))

untargeted_naind_ad_avg_retained <- untargeted_naind_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_naind_ad_in_all <- untargeted_naind_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")



### Untargeted, PD ----------------

untargeted_naind_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted_naind5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T))
untargeted_naind_pd_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

untargeted_naind_pd_logistic[[1]][[2]] + 
  labs(title = "ROC: PD vs {Controls, AD}",
       subtitle = TeX('Untargeted with NA indicators (>.5 NA) ,$\\alpha = 0.5$'))
ggsave("pd_logistic_untargeted_naind.png")


# diagnostic plots
heatmap.fit(untargeted_naind_pd_logistic[[1]]$truth, pred = untargeted_naind_pd_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(untargeted_naind_pd_logistic[[1]]$truth, pred = untargeted_naind_pd_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = untargeted_naind_pd_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "PD (untargeted_naind) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

untargeted_naind_pd_coefs <- untargeted_naind_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(median_coef)))

untargeted_naind_pd_avg_retained <- untargeted_naind_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_naind_pd_in_all <- untargeted_naind_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")



### lipids, AD ----------------
lipids_naind_ad_logistic <- purrr::map(1:5,~ logistic_control_analysis(imputed_all_lipids_naind5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))
lipids_naind_ad_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

lipids_naind_ad_logistic[[1]][[2]] + 
  labs(title = "ROC: AD vs {Controls, AD}",
       subtitle = TeX('lipids with NA indicators (>.5 NA) ,$\\alpha = 0.5$'))
ggsave("ad_logistic_lipids_naind.png")


# diagnostic plots
heatmap.fit(lipids_naind_ad_logistic[[1]]$truth, pred = lipids_naind_ad_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(lipids_naind_ad_logistic[[1]]$truth, pred = lipids_naind_ad_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = lipids_naind_ad_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "AD (lipids_naind) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

lipids_naind_ad_coefs <- lipids_naind_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(median_coef)))

lipids_naind_ad_avg_retained <- lipids_naind_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_naind_ad_in_all <- lipids_naind_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")



### lipids, PD ----------------
lipids_naind_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids_naind5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T))
lipids_naind_pd_logistic %>%
  purrr::map(~.x[[2]]) %>%
  cowplot::plot_grid(plotlist = .)

lipids_naind_pd_logistic[[1]][[2]] + 
  labs(title = "ROC: PD vs {Controls, AD}",
       subtitle = TeX('lipids with NA indicators (>.5 NA) ,$\\alpha = 0.5$'))
ggsave("pd_logistic_lipids_naind.png")


# diagnostic plots
heatmap.fit(lipids_naind_pd_logistic[[1]]$truth, pred = lipids_naind_pd_logistic[[1]]$pred)
dev_resid <- compute_deviance_resid(lipids_naind_pd_logistic[[1]]$truth, pred = lipids_naind_pd_logistic[[1]]$pred)
ggplot() + geom_point(aes(x = lipids_naind_pd_logistic[[1]]$pred, y = dev_resid)) +
  labs(title = "PD (lipids_naind) Residuals plot",
       x = "Fitted Value",
       y = "Deviance")

#-- 

lipids_naind_pd_coefs <- lipids_naind_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(median_coef)))

lipids_naind_pd_avg_retained <- lipids_naind_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_naind_pd_in_all <- lipids_naind_pd_logistic %>% 
  purrr::map(~.x[[1]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")






