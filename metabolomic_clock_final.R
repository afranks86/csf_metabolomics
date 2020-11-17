#error_log <- file("analysis_log_102719.Rout", open="wt")
#sink(error_log, type = "message")
source(here::here("analysis", "starter.R"))
  


#########################################

# Data Exploration (PCA and PLS-DA)

##########################################


## Before we get into it, here's a quick distribution of age. could've picked any dataset -- same subjects
age_dist_c <- qplot(wide_data_targeted %>% 
        filter(Type %in% c("CO", "CM", "CY")) %>%
        pull(Age),
      binwidth = 5) + 
  labs(x = "Age",
       y = "Count",
       title = "Controls: Distribution of Age") + 
  scale_fill_viridis_c()

ggsave("age_dist_c.png",
       plot = age_dist_c,
       path = here("plots", "aging_figs"),
       width= 14,
       height = 10)

### Batch Balance --------------
subject_data_meta <- raw_data_targeted %>%     
  filter(!(Type %in% c("Other"))) %>%
  select(Age, Gender, APOE, Type, GBAStatus, GBA_T369M, Id, Batch) %>%
  distinct() %>%
  mutate(across(c(Gender, APOE, Type, GBAStatus, GBA_T369M, Batch), ~as_factor(.x))) %>%
  mutate(Type = fct_relevel(Type, "CY", "CM", "CO", "AD", "PD"))
  

# sex
sex_batch_balance <- ggplot(subject_data_meta) + 
  geom_bar(aes(Gender)) + 
  facet_wrap(~Batch) + 
  labs(title = 'Sex-Batch Balance',
       x = 'Sex at birth')

# Age
age_batch_balance <- ggplot(subject_data_meta) + 
  geom_boxplot(aes(Batch, Age), fill = 'gray35') + 
  labs(title = 'Age-Batch Balance')

# phenotype
type_batch_balance <- ggplot(subject_data_meta) +
  geom_bar(aes(Type)) + 
  facet_wrap(~Batch) + 
  labs(title = 'Phenotype-Batch Balance',
       x = 'Phenotype')

ggsave("sex_batch_balance.png",
       plot = sex_batch_balance,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)
ggsave("age_batch_balance.png",
       plot = age_batch_balance,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)
ggsave("type_batch_balance.png",
       plot = type_batch_balance,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)


## PCA --------------------------

# subset the data to get rid of metadata and columns with > 10% NA
# note: untargeted_processed has only features with < 10% missingness, and scaled/centered features
# note: nipals handles missing values by using null weights, so no need to impute first.
pca_c_untargeted_x <- untargeted_c_processed %>% 
  dplyr::select(-any_of(metadata_cols))

pca_c_untargeted <-  pca_c_untargeted_x %>%
  opls(algoC = "nipals",
       parAsColFcVn = untargeted_c_processed$Gender)

pca_c_untargeted <- readRDS(here("aging_output_files", "pca_c_untargeted.Rds"))
# saveRDS(pca_c_untargeted, here("aging_output_files", "pca_c_untargeted.Rds"))


pca_scores_c_untargeted <- pca_c_untargeted@scoreMN %>%
  as_tibble() %>%
  bind_cols(untargeted_c_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id))
  
pc1_var_explained <- pca_c_untargeted@modelDF["p1","R2X"] %>%
  scales::percent()
pc2_var_explained <- pca_c_untargeted@modelDF["p2","R2X"] %>%
  scales::percent()
pc3_var_explained <- pca_c_untargeted@modelDF["p3","R2X"] %>%
  scales::percent()
pc4_var_explained <- pca_c_untargeted@modelDF["p4","R2X"] %>%
  scales::percent()

# look at median age for pc1 > 0 and pc1 <0
pca_scores_c_untargeted %>%
  transmute(Age, p1, 
            p1_greater_than_zero = p1 > 0) %>%
  group_by(p1_greater_than_zero) %>%
  summarize(p1_median_age = median(Age))

# by gender
pca_gender_untargeted <- ggplot(pca_scores_c_untargeted, aes(p1, p2, color = Sex)) + 
  geom_point(size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "Principal Component Analysis",
       subtitle = "Colored by Sex",
       x = str_glue("PC1 ({pc1_var_explained})"),
       y = str_glue("PC2 ({pc2_var_explained})")) + 
  scale_color_viridis_d()
ggsave("pca_gender_untar.png",
       plot = pca_gender_untargeted,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)

# by age
pca_age_untargeted <- ggplot(pca_scores_c_untargeted, aes(Age, p1, color = Age)) + 
  geom_point(size = 5) +
  scale_color_viridis_c() + 
  #stat_ellipse() +
  labs(title = "Principal Component Analysis",
       subtitle = str_glue("First PC by Age"),
       y = str_glue("PC1 ({pc1_var_explained})"),
       x = str_glue("Chronological Age (Years)")) 

ggsave("pca_age_untar.png",
       plot = pca_age_untargeted,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)





## PLS-DA  (using gender as response) -------------------------

# R2X / R2Y: % variance in X/Y explained by component
# Q2 : 1 - (predicted RSS / sum squares Y). tiny q2y means pretty bad predictive performance
plsda_c_untargeted <- untargeted_c_processed %>%
  dplyr::select(-any_of(metadata_cols)) %>%
  opls(y = untargeted_c_processed$Gender, algoC = "nipals", predI = 2)

plsda_c_untargeted <- readRDS(here("aging_output_files", "plsda_c_untargeted.Rds"))
# saveRDS(plsda_c_untargeted, here("aging_output_files", "plsda_c_untargeted.Rds"))



plsda_pc1_varx_explained <- plsda_c_untargeted@modelDF["p1","R2X"] %>%
  scales::percent()
plsda_pc2_varx_explained <- plsda_c_untargeted@modelDF["p2","R2X"] %>%
  scales::percent()


plsda_scores_c_untargeted <- plsda_c_untargeted@scoreMN %>%
  as_tibble() %>%
  bind_cols(untargeted_c_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id))

plsda_gender_untargeted <- ggplot(plsda_scores_c_untargeted, aes(p1, p2, color = Sex)) + 
  geom_point(size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "Partial Least Squares Discriminant Analysis",
       subtitle = "By Sex",
       x = str_glue("PC1 ({plsda_pc1_varx_explained})"),
       y = str_glue("PC2 ({plsda_pc2_varx_explained})")) + 
  scale_color_viridis_d()

ggsave("plsda_sex_untar.png",
       plot = plsda_gender_untargeted,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)



## compute first pc, do plsda in null space of transofmration matrix ------

x_c_nullpc1 <- pca_c_untargeted_x - pca_c_untargeted@scoreMN[,'p1'] %*% t(pca_c_untargeted@loadingMN[,'p1'])
plsda_c_untargeted_nullpc1 <- x_c_nullpc1 %>%
  opls(y = untargeted_c_processed$Gender, algoC = "nipals", predI = 2)

plsda_c_untargeted_nullpc1 <- readRDS(here("aging_output_files", "plsda_c_untargeted_nullpc1.Rds"))
# saveRDS(plsda_c_untargeted_nullpc1, here("aging_output_files", "plsda_c_untargeted_nullpc1.Rds"))



plsda_null_pc1_varx_explained <- plsda_c_untargeted_nullpc1@modelDF["p1","R2X"] %>%
  scales::percent()
plsda_null_pc2_varx_explained <- plsda_c_untargeted_nullpc1@modelDF["p2","R2X"] %>%
  scales::percent()


plsda_scores_c_untargeted_nullpc1 <- plsda_c_untargeted_nullpc1@scoreMN %>%
  as_tibble() %>%
  bind_cols(untargeted_c_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id))

pca_plsda_scores_c_untargeted <- pca_scores_c_untargeted %>%
  select("p1_pca" = p1, Id) %>%
  inner_join(plsda_scores_c_untargeted_nullpc1, by = "Id") %>%
  select(-p2) %>%
  rename("p1_plsda" = p1)


pca_plsda_untar <- ggplot(pca_plsda_scores_c_untargeted, aes(p1_pca, p1_plsda, color = Sex)) + 
  geom_point(aes(alpha = Age), size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "PCA/PLS-DA",
       x = str_glue("PCA: PC1 ({pc1_var_explained})"),
       y = str_glue("PLSDA: PC2 ({plsda_null_pc1_varx_explained})")) + 
  scale_color_viridis_d() +
  scale_fill_viridis_c()

ggsave("pca_plsda_untar.png",
       plot = pca_plsda_untar,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)


######################################################################



### Clock (linear glmnet age) ###
## Note: For the Amelia imputations, we keep empri = ~ 10% of the number of observations
##        This is a loose upper bound sugggested by the amelia paper




######################################################################








########

### GOT

########
message("GOT -------------------------------------------")

imputed_c_got5 <- filter_and_impute_multi(wide_data, c('CO', 'CY', 'CM'), empri = 10, scale = T)
got_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_got5, name = "GOT", color = NULL, imp_num = .x, map_got = FALSE))

imputed_c_got5 <- readRDS(here("aging_output_files", "imputed_c_got5.Rds"))
got_age_analysis <- readRDS(here("aging_output_files", "got_age_analysis.Rds"))
# saveRDS(imputed_c_got5, here("aging_output_files", "imputed_c_got5.Rds"))
# saveRDS(got_age_analysis, here("aging_output_files", "got_age_analysis.Rds"))

got_pred_df <- imputation_df(got_age_analysis)

# final plot
got_age_predtruth <- predtruth_plot(got_pred_df, name = "GOT")
ggsave(filename = "age_clock_GOT_pred.png", 
       plot = got_age_predtruth,
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)


got_avg_retained <- got_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

got_in_all <- got_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

got_in_at_least_one <- got_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")

got_retained_table <- tibble(
  data = "GOT",
  avg_retained = got_avg_retained,
  num_in_all = length(got_in_all),
  num_in_any = length(got_in_at_least_one)
)


## Create null model as baseline (it's dataset agnositic, since the observations are same)
null_df <- got_pred_df %>%
  ungroup() %>%
  select(truth, gender, apoe, apoe4, type) %>%
  mutate(pred = mean(truth))

(null_age_predtruth <- predtruth_plot(null_df, pred_name = "pred", name = "Mean model", errorbar = FALSE))








### Lipids

########
message("Lipids (amelia) -------------------------------------------")

imputed_c_lipids5 <- filter_and_impute_multi(wide_data_lipids, c('CO', 'CY', 'CM'), empri = 10)
lipids_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_lipids5, name = "Lipids", color = NULL, imp_num = .x))

imputed_c_lipids5 <- readRDS(here("aging_output_files", "imputed_c_lipids5.Rds"))
lipids_age_analysis <- readRDS(here("aging_output_files", "lipids_age_analysis.Rds"))
# saveRDS(imputed_c_lipids5, here("aging_output_files", "imputed_c_lipids5.Rds"))
# saveRDS(lipids_age_analysis, here("aging_output_files", "lipids_age_analysis.Rds"))


lipids_pred_df <- imputation_df(lipids_age_analysis)
  
lipids_age_predtruth <- predtruth_plot(lipids_pred_df, name = "Lipids")
ggsave(filename = "age_clock_lipids_pred.png", 
       plot = lipids_age_predtruth,
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)

lipids_avg_retained <- lipids_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_in_all <- lipids_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

lipids_in_at_least_one <- lipids_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")

lipids_retained_table <- tibble(
  data = "Lipids",
  avg_retained = lipids_avg_retained,
  num_in_all = length(lipids_in_all),
  num_in_any = length(lipids_in_at_least_one)
)

lipids_coef_table <- lipids_age_analysis %>%
  purrr::imap(~enframe(.x[[2]], value = str_glue("coef{.y}"))) %>%
  reduce(~full_join(.x, .y, by = "name")) %>%
  rowwise() %>%
  transmute(name,
            avgCoef = mean(c(coef1, coef2, coef3, coef4, coef5), na.rm = TRUE) %>% round(1)
  ) %>%
  ungroup() %>%
  arrange(desc(abs(avgCoef))) %>%
  filter(!(name %in% c("GenderM", "(Intercept)")))

# write separate tables for positive and negative coefs
lipids_coef_table %>%
  filter(avgCoef > 0) %>%
  write_csv(here("aging_tables", "pos_coef_table_lipids.csv"))

lipids_coef_table %>%
  filter(avgCoef < 0) %>%
  write_csv(here("aging_tables", "neg_coef_table_lipids.csv"))






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

### Targeted

########

message("Targeted -------------------------------------------")

imputed_c_targeted5 <- filter_and_impute_multi(wide_data_targeted, c('CO', 'CY', 'CM'), scale = T)
targeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_targeted5, name = "Targeted", color = NULL, imp_num = .x))

imputed_c_targeted5 <- readRDS(here("aging_output_files", "imputed_c_targeted5.Rds"))
targeted_age_analysis <- readRDS(here("aging_output_files", "targeted_age_analysis.Rds"))
# saveRDS(imputed_c_targeted5, here("aging_output_files", "imputed_c_targeted5.Rds"))
# saveRDS(targeted_age_analysis, here("aging_output_files", "targeted_age_analysis.Rds"))


targeted_pred_df <- imputation_df(targeted_age_analysis)


targeted_age_predtruth <- predtruth_plot(targeted_pred_df, name = "Targeted")
ggsave(filename = "age_clock_targeted_pred.png", 
       plot = targeted_age_predtruth,
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)

targeted_avg_retained <- targeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

targeted_in_all <- targeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

targeted_in_at_least_one <- targeted_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")

targeted_retained_table <- tibble(
  data = "Targeted",
  avg_retained = targeted_avg_retained,
  num_in_all = length(targeted_in_all),
  num_in_any = length(targeted_in_at_least_one)
)

targeted_coef_table <- targeted_age_analysis %>%
  purrr::imap(~enframe(.x[[2]], value = str_glue("coef{.y}"))) %>%
  reduce(~full_join(.x, .y, by = "name")) %>%
  rowwise() %>%
  transmute(name = str_replace_all(name, "_pos|_neg", ""),
            avgCoef = mean(c(coef1, coef2, coef3, coef4, coef5), na.rm = TRUE) %>% round(1)
            ) %>%
  ungroup() %>%
  arrange(desc(abs(avgCoef))) %>%
  filter(name != "(Intercept)")

# write separate tables for positive and negative coefs
targeted_coef_table %>%
  filter(avgCoef > 0) %>%
  write_csv(here("aging_tables", "pos_coef_table_targeted.csv"))

targeted_coef_table %>%
  filter(avgCoef < 0) %>%
  write_csv(here("aging_tables", "neg_coef_table_targeted.csv"))




# #### wild test
# m_targeted <- wide_data_targeted %>% 
#   select_if(is.numeric) %>%
#   select(-Age)
# 
# 
# m_age <- wide_data_untargeted$Age
# 
# m_targeted %>%
#   select_if(function(x) abs(cor(x, m_age)) > 0.3)
# 
# cor(m_age, m_targeted[,1])
# for(i in 1:ncol(m_targeted)){
#   if(abs(cor(m_targeted[,i], m_age)) > 0.3){
#     m_
#   }
# }
# 
# mt_untargeted <- wide_data_untargeted %>%
#   select_if(is.numeric) %>%
#   select(-Age) %>%
#   transmute_all(function(x) scale(x, center = T, scale = T) %>% Hmisc::impute()) 
# 
# # mt_imp <- amelia(t(mt_untargeted), empri = 10000000)
# mt_imp <- mice(t(mt_untargeted))
# # 1:5 %>%
# #   purrr::map(~ complete(mt_imp, .x) %>%
# #                t %>%
# #                as_tibble() %>%
# mt_untargeted %>%
#   transmute_all(function(x) abs(cor(x, m_age))) %>% 
#   slice(1) %>% 
#   gather() %>% 
#   ggplot() + 
#   geom_histogram(aes(value)) #%>%
#   #cowplot::plot_grid(plotlist = .)
# 
# 
# 
# 
# 
# untar_tiny <- wide_data_untargeted %>% select(setdiff(untargeted_in_at_least_one, "GenderM"), Age, Gender, Type, APOE, Id, GBA_T369M, GBAStatus)
# imp_untar_tiny <- filter_and_impute_multi(untar_tiny, c("CO", "CY", "CM"), transpose = F, empri = 200)
# untar_tiny_age_analysis <- purrr::map(1:5, ~age_control_analysis(imp_untar_tiny, name = "untargeted tiny", color = NULL, imp_num = .x, nlambda = 200))
# untar_tiny_pred_df <- imputation_df(untar_tiny_age_analysis)
# (untar_tiny_age_predtruth <- predtruth_plot(untar_tiny_pred_df, name = "untar tiny"))
# untar_tiny_avg_retained <- untar_tiny_age_analysis %>% 
#   purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
#   unlist %>% mean
# 
# untar_tiny_in_all <- untar_tiny_age_analysis %>% 
#   purrr::map(~.x[[2]] %>% names) %>% 
#   reduce(intersect) %>%
#   reduce(c) %>%
#   setdiff("(Intercept)")
# 
# 
# 
# untar_tiny_ad_amelia5 <- filter_and_impute_multi(untar_tiny, c('AD'), empri = 100)
# untar_tiny_pd_amelia5 <- filter_and_impute_multi(untar_tiny, c('PD'), empri = 100)
# untar_tiny_adpd_separate_amelia5 <- purrr::map2(untar_tiny_ad_amelia5, untar_tiny_pd_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE))
# untar_tiny_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(imp_untar_tiny[[.x]], new_data = untar_tiny_adpd_separate_amelia5[[.x]], nlambda = 200))
# untar_tiny_adpd_pred_df_separate <- imputation_df(untar_tiny_adpd_age_analysis_separate)
# (untar_tiny_adpd_age_predtruth <- predtruth_plot(untar_tiny_adpd_pred_df_separate, name = "untar tiny adpd"))
# 
# # match
# untar_tiny_adpd_matched_controls <- purrr::map(1:nrow(untar_tiny), ~find_control(.x, data = filter(untar_tiny, Type %in% c("AD", "PD")), 
#                                                                                            data_control = filter(untar_tiny, Type %in% c("CO", "CY", "CM")))) %>%
#   bind_rows(filter(untar_tiny, Type %in% c("AD", "PD"))) %>%
#   distinct(.keep_all = T)
# 
# # imputation_df() expects a list of lists, but only uses the first element of the inner list. so we add a random "1" to the output. doesn't do anything.
# untar_tiny_c_matched_age_analysis_table <- purrr::map(untar_tiny_age_analysis, ~list(.x[[1]] %>% filter(id %in% untar_tiny_adpd_matched_controls$Id), 1))
# untar_tiny_c_matched_age_pred_df <- imputation_df(untar_tiny_c_matched_age_analysis_table)
# 
# (untar_tiny_c_matched_age_predtruth <- predtruth_plot(untar_tiny_c_matched_age_pred_df, name = "untar_tiny (matched)", color = NULL))
# 
# ##### random n features
# 
# random_n_features_analysis <- function(n){
#   untar_random_100_features <- wide_data_untargeted %>%
#     names %>%
#     setdiff(c("Age", "Gender", "Type", "APOE", "Id", "GBA_T369M", "GBAStatus")) %>%
#     base::sample(size = n, replace = F)
#   
#   untar_tiny_random <- wide_data_untargeted %>% 
#     select(untar_random_100_features, Age, Gender, Type, APOE, Id, GBA_T369M, GBAStatus)
#   
#   imp_untar_tiny_random <- filter_and_impute_multi(untar_tiny_random, c("CO", "CY", "CM"), transpose = T, empri = 200)
#   untar_tiny_random_age_analysis <- purrr::map(1:5, ~age_control_analysis(imp_untar_tiny_random, name = "untargeted tiny random 100", color = NULL, imp_num = .x, nlambda = 200))
#   untar_tiny_random_pred_df <- imputation_df(untar_tiny_random_age_analysis)
#   untar_tiny_random_age_predtruth <- predtruth_plot(untar_tiny_random_pred_df, name = paste("untar tiny", n))
#   
#   untar_tiny_random_age_predtruth
#   
# }
# 
# 
# random_100_features <- purrr::map(1:3, ~random_n_features_analysis(100))
# random_200_features <- purrr::map(1:3, ~random_n_features_analysis(200))
# random_500_features <- purrr::map(1:3, ~random_n_features_analysis(500))
# random_1000_features <- purrr::map(1:3, ~random_n_features_analysis(1000))
# random_2000_features <- purrr::map(1:3, ~random_n_features_analysis(2000))
# random_3000_features <- purrr::map(1:3, ~random_n_features_analysis(3000))
# 
# 
# cowplot::plot_grid(plotlist = random_100_features)
# cowplot::plot_grid(plotlist = random_200_features)
# cowplot::plot_grid(plotlist = random_500_features)
# cowplot::plot_grid(plotlist = random_1000_features)
# cowplot::plot_grid(plotlist = random_2000_features)
# cowplot::plot_grid(plotlist = random_3000_features)
# 
# 
# # match
# untar_tiny_random_adpd_matched_controls <- purrr::map(1:nrow(untar_tiny_random), ~find_control(.x, data = filter(untar_tiny_random, Type %in% c("AD", "PD")), 
#                                                                                  data_control = filter(untar_tiny_random, Type %in% c("CO", "CY", "CM")))) %>%
#   bind_rows(filter(untar_tiny_random, Type %in% c("AD", "PD"))) %>%
#   distinct(.keep_all = T)
# 
# # imputation_df() expects a list of lists, but only uses the first element of the inner list. so we add a random "1" to the output. doesn't do anything.
# untar_tiny_random_c_matched_age_analysis_table <- purrr::map(untar_tiny_random_age_analysis, ~list(.x[[1]] %>% filter(id %in% untar_tiny_random_adpd_matched_controls$Id), 1))
# untar_tiny_random_c_matched_age_pred_df <- imputation_df(untar_tiny_random_c_matched_age_analysis_table)
# 
# (untar_tiny_random_c_matched_age_predtruth <- predtruth_plot(untar_tiny_random_c_matched_age_pred_df, name = "untar_tiny_random (matched)", color = NULL))
# 
# 


#######

### Unargeted

########

message("untargeted -------------------------------------------")

imputed_c_untargeted5 <- filter_and_impute_multi(wide_data_untargeted, c('CO', 'CY', 'CM'), empri = 10)
#save(imputed_c_untargeted5, file = "imputed_c_untargeted5_03022020")
# lambdas were reaching the end of their sequence due to hugeness of data, so we increase nlambda from 100 to 200
untargeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted5, name = "untargeted", color = NULL, imp_num = .x, nlambda = 200))
#save(untargeted_age_analysis, file = "untargeted_age_analysis_02032020.RData")


# saveRDS(imputed_c_untargeted5, here("aging_output_files", "imputed_c_untargeted5.Rds"))
# saveRDS(untargeted_age_analysis, here("aging_output_files", "untargeted_age_analysis.Rds"))
imputed_c_untargeted5 <- readRDS(here("aging_output_files", "imputed_c_untargeted5.Rds"))
untargeted_age_analysis <- readRDS(here("aging_output_files", "untargeted_age_analysis.Rds"))

untargeted_pred_df <- imputation_df(untargeted_age_analysis)

untargeted_age_predtruth <- predtruth_plot(untargeted_pred_df, name = "Untargeted")
ggsave(filename = "age_clock_untargeted_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = untargeted_age_predtruth,
       width = 14, height = 10)

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

untargeted_retained_table <- tibble(
  data = "Untargeted",
  avg_retained = untargeted_avg_retained,
  num_in_all = length(untargeted_in_all),
  num_in_any = length(untargeted_in_at_least_one)
)





#plot colored by different things
predtruth_plot(untargeted_pred_df, color = "apoe", errorbar = FALSE, name = "Untargeted")
predtruth_plot(untargeted_pred_df, color = "apoe4", errorbar = FALSE, name = "Untargeted")



########

### Combined Profiles

########

message("Combined GOT/lipids -------------------------------------------")

#' Helper functions to combine the imputed lipids/GOT/targeted/untargeted data
#' @param got_imp is a single element of filter_and_impute_multi(data = wide_data)
#' @param lipid_imp is a single element of filter_and_impute_multi(data = wide_data_lipids). this is the only one required
combine_got_lipids <- function(lipid_imp, got_imp = NULL, targeted_imp = NULL, untargeted_imp = NULL){
  
  
  # we only need to track apoe/gender/age once
  imputed_c_lipids_separate_withid <- lipid_imp[[1]] %>%
    as_tibble() %>%
    dplyr::select(-c(Age, GenderM, `(Intercept)`)) %>%
    dplyr::mutate(id = lipid_imp[[4]], 
                  apoe = lipid_imp[[3]],
                  type = lipid_imp[[2]])
  
  if(!is.null(got_imp)){
    # add the id columns to make the join easier
    imputed_c_got_separate_withid <- got_imp[[1]] %>% 
      as_tibble() %>%
      mutate(id = got_imp[[4]])  
    
    imputed_c_combined_separate_df <- imputed_c_lipids_separate_withid %>%
      inner_join(imputed_c_got_separate_withid, by = 'id') 
  }
  
  if(!is.null(targeted_imp)){
    imputed_c_targeted_separate_withid <- targeted_imp[[1]] %>% 
      as_tibble() %>%
      dplyr::select(-c(Age, GenderM, `(Intercept)`)) %>%
      mutate(id = targeted_imp[[4]])  
    
    imputed_c_combined_separate_df <- imputed_c_combined_separate_df %>%
      inner_join(imputed_c_targeted_separate_withid, by = "id")
  }
  
  if(!is.null(untargeted_imp)){
    imputed_c_untargeted_separate_withid <- untargeted_imp[[1]] %>%
      as_tibble() %>%
      dplyr::select(-c(Age, GenderM, `(Intercept)`)) %>%
      mutate(id = untargeted_imp[[4]])
    
    # join
    imputed_c_combined_separate_df <- imputed_c_combined_separate_df %>%
      inner_join(imputed_c_untargeted_separate_withid, by = "id")
    
  }
  
  imputed_c_combined_separate_df <- imputed_c_combined_separate_df %>%
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



imputed_c_combined_amelia5 <- purrr::pmap(list(imputed_c_got5, imputed_c_lipids5, imputed_c_targeted5, imputed_c_untargeted5), ~combine_got_lipids(..1, ..2, ..3, ..4))
combined_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_combined_amelia5, name = "Combined Profiles", color = NULL, imp_num = .x))

saveRDS(imputed_c_combined_amelia5, here("aging_output_files", "imputed_c_combined_amelia5.Rds"))
saveRDS(combined_age_analysis, here("aging_output_files", "combined_age_analysis.Rds"))

combined_pred_df <- imputation_df(combined_age_analysis)

(combined_age_predtruth <- predtruth_plot(combined_pred_df, name = "Combined Profiles"))
ggsave(filename = "age_clock_combined_pred.png", 
       plot = combined_age_predtruth,
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)

combined_avg_retained <- combined_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

combined_in_all <- combined_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")


## Correlation table suggests that Untargeted and Lipids are pretty far apart. try a model combining only these two
imputed_c_untar_lipids_amelia5 <- purrr::map2(imputed_c_lipids5, imputed_c_untargeted5, ~combine_got_lipids(.x, .y))
untar_lipids_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untar_lipids_amelia5, name = "Untargeted & Lipid Profiles", color = NULL, imp_num = .x))

untar_lipids_pred_df <- imputation_df(untar_lipids_age_analysis)

(untar_lipids_age_predtruth <- predtruth_plot(untar_lipids_pred_df, name = "Untargeted & Lipid Profiles"))
ggsave(filename = "age_clock_untar_lipids_pred.png", 
       plot = untar_lipids_age_predtruth,
       path = here("plots", "aging_figs"), 
       width = 14, height = 10)







######################

### Post selection inference

######################


# post selection
untargeted_post <- post_select(imputed_c_untargeted5, untargeted_age_analysis, "untargeted-post")
post_untargeted_pred_df <- imputation_df(untargeted_post)
post_untargeted_age_predtruth <- predtruth_plot(post_untargeted_pred_df, name = "Untargeted Post")
ggsave(filename = "post_age_clock_untargeted_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = post_untargeted_age_predtruth,
       width = 14, height = 10)


targeted_post <- post_select(imputed_c_targeted5, targeted_age_analysis, "targeted-post")
post_targeted_pred_df <- imputation_df(targeted_post)
post_targeted_age_predtruth <- predtruth_plot(post_targeted_pred_df, name = "Targeted Post")
ggsave(filename = "post_age_clock_targeted_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = post_targeted_age_predtruth,
       width = 14, height = 10)


lipids_post <- post_select(imputed_c_lipids5, lipids_age_analysis, "lipids-post")
post_lipids_pred_df <- imputation_df(lipids_post)
post_lipids_age_predtruth <- predtruth_plot(post_lipids_pred_df, name = "Lipids Post")
ggsave(filename = "post_age_clock_lipids_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = post_lipids_age_predtruth,
       width = 14, height = 10)

got_post <- post_select(imputed_c_got5, got_age_analysis, "got-post")
post_got_pred_df <- imputation_df(got_post)
post_got_age_predtruth <- predtruth_plot(post_got_pred_df, name = "GOT Post")
ggsave(filename = "post_age_clock_GOT_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = post_got_age_predtruth,
       width = 14, height = 10)




# ########
# 
# ### untargeted Raw
# 
# ########
# 
# imputed_c_untargeted_raw5 <- filter_and_impute_multi(wide_data_untargeted_raw, c('CO', 'CY', 'CM'), empri = 10)
# #save(imputed_c_untargeted_raw5, file = "imputed_c_untargeted_raw5_03022020")
# # lambdas were reaching the end of their sequence due to hugeness of data, so we increase nlambda from 100 to 200
# untargeted_raw_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted_raw5, name = "Untargeted, Raw", color = NULL, imp_num = .x, nlambda = 200))
# #save(untargeted_raw_age_analysis, file = "untargeted_raw_age_analysis_02032020.RData")
# 
# untargeted_raw_pred_df <- imputation_df(untargeted_raw_age_analysis)
# 
# (untargeted_raw_age_predtruth <- predtruth_plot(untargeted_raw_pred_df, name = "Untargeted, Raw")) 
# 

########

### untargeted Raw scaled

########

imputed_c_untargeted_rawscaled5 <- filter_and_impute_multi(wide_data_untargeted_rawscaled, c('CO', 'CY', 'CM'), empri = 10)
#save(imputed_c_untargeted_rawscaled5, file = "imputed_c_untargeted_rawscaled5_03022020")
# lambdas were reaching the end of their sequence due to hugeness of data, so we increase nlambda from 100 to 200
untargeted_rawscaled_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted_rawscaled5, name = "Untargeted, Raw Scaled", color = NULL, imp_num = .x, nlambda = 200))
#save(untargeted_rawscaled_age_analysis, file = "untargeted_rawscaled_age_analysis_02032020.RData")

imputed_c_untargeted_rawscaled5 <- readRDS(here("aging_output_files", "imputed_c_untargeted_rawscaled5.Rds"))
untargeted_rawscaled_age_analysis <- readRDS(here("aging_output_files", "untargeted_rawscaled_age_analysis.Rds"))


# saveRDS(imputed_c_untargeted_rawscaled5, here("aging_output_files", "imputed_c_untargeted_rawscaled5.Rds"))
# saveRDS(untargeted_rawscaled_age_analysis, here("aging_output_files", "untargeted_rawscaled_age_analysis.Rds"))


untargeted_rawscaled_pred_df <- imputation_df(untargeted_rawscaled_age_analysis)

untargeted_rawscaled_age_predtruth <- predtruth_plot(untargeted_rawscaled_pred_df, 
                                                     name = "Untargeted, Without Drift Correction")
ggsave(filename = "no_drift_correction_age_clock_untargeted_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = untargeted_rawscaled_age_predtruth,
       width = 14, height = 10)










# 
# ########
# 
# ### Untargeted, select columns with <90% missingness
# 
# ########
# 
# 
# message("Untargeted  subset-------------------------------------------")
# 
# 
# 
# # What about doing untargeted analysis on <90% missing?
# untargeted_less_90perc_missing <-wide_data_untargeted %>%
#   map_dbl(~ sum(is.na(.x))/length(.x)) %>% 
#   enframe(value = "perct_missing") %>% filter(perct_missing < .9)
# 
# wide_data_untargeted_less90perc_missing <- wide_data_untargeted %>% select_at(vars(untargeted_less_90perc_missing$name))    
# 
# 
# # analysis
# imputed_c_untargeted_90subset5 <- filter_and_impute_multi(wide_data_untargeted_less90perc_missing , c('CO', 'CY', 'CM'))
# untargeted_90subset_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted_90subset5, name = "untargeted_90subset", color = NULL, imp_num = .x))
# 
# untargeted_90subset_pred_df <- imputation_df(untargeted_90subset_age_analysis)
# 
# (untargeted_90subset_age_predtruth <- predtruth_plot(untargeted_90subset_pred_df, name = "untargeted_90subset")) 
# ggsave(filename = "age_clock_untargeted_90subset_pred.png") #14.9 x 8.21
# 
# untargeted_90subset_avg_retained <- untargeted_90subset_age_analysis %>% 
#   purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
#   unlist %>% mean
# 
# untargeted_90subset_in_all <- untargeted_90subset_age_analysis %>% 
#   purrr::map(~.x[[2]] %>% names) %>% 
#   reduce(intersect) %>%
#   setdiff("(Intercept)")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ########
# 
# ### Untargeted, select columns with <10% missingness
# 
# ########
# 
# 
# message("Untargeted  subset-------------------------------------------")
# 
# 
# 
# # What about doing untargeted analysis on <10% missing?
# untargeted_less_10perc_missing <-wide_data_untargeted %>%
#   map_dbl(~ sum(is.na(.x))/length(.x)) %>% 
#   enframe(value = "perct_missing") %>% filter(perct_missing < .1)
# 
# wide_data_untargeted_less10perc_missing <- wide_data_untargeted %>% select_at(vars(untargeted_less_10perc_missing$name))    
# 
# 
# # analysis
# imputed_c_untargeted_10subset5 <- filter_and_impute_multi(wide_data_untargeted_less10perc_missing , c('CO', 'CY', 'CM'))
# untargeted_10subset_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_untargeted_10subset5, name = "untargeted_10subset", color = NULL, imp_num = .x))
# 
# untargeted_10subset_pred_df <- imputation_df(untargeted_10subset_age_analysis)
# 
# (untargeted_10subset_age_predtruth <- predtruth_plot(untargeted_10subset_pred_df, name = "Untargeted (<10% missing)")) 
# ggsave(filename = "age_clock_untargeted_10subset_pred.png") #14.9 x 8.21
# 
# untargeted_10subset_avg_retained <- untargeted_10subset_age_analysis %>% 
#   purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
#   unlist %>% mean
# 
# untargeted_10subset_in_all <- untargeted_10subset_age_analysis %>% 
#   purrr::map(~.x[[2]] %>% names) %>% 
#   reduce(intersect) %>%
#   setdiff("(Intercept)")
# 

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


#############################

#### knockoff

#############################

my_glmnet_stat <- function(X, X_k, y){
  stat.glmnet_coefdiff(X, X_k, y, nlambda = 200, alpha = 0.5)
}

testx <- imputed_c_targeted5[[1]]$Y[,-which(colnames(imputed_c_targeted5[[1]]$Y) %in% c("(Intercept)", "Age"))]
testy <- as.matrix(imputed_c_targeted5[[1]]$Y[, 'Age'])
knockoff_tar <- knockoff.filter(testx, testy, knockoffs = create.second_order, statistic = my_glmnet_stat, fdr = 0.05)

print(knockoff_tar$selected)







########################

### comparing retained metabolite numbers

#######################

# join all the retained summary tables and make it long
retained_table <- bind_rows(got_retained_table, lipids_retained_table, targeted_retained_table, untargeted_retained_table) %>%
  pivot_longer(cols = -data, names_to = "stat", values_to = "num_metabolites") %>%
  mutate(num_features = 
           case_when(
             data == 'GOT' ~ length(setdiff(colnames(imputed_c_got5[[1]][[1]]), c('(Intercept)', 'Age', 'GenderM'))),
             data == 'Lipids' ~ length(setdiff(colnames(imputed_c_lipids5[[1]][[1]]), c('(Intercept)', 'Age', 'GenderM'))),
             data == 'Targeted' ~ length(setdiff(colnames(imputed_c_targeted5[[1]][[1]]), c('(Intercept)', 'Age', 'GenderM'))),
             data == 'Untargeted' ~ length(setdiff(colnames(imputed_c_untargeted5[[1]][[1]]), c('(Intercept)', 'Age', 'GenderM')))
           ),
         pct_of_total = scales::percent(num_metabolites / num_features, accuracy = 0.1)
  ) %>%
  mutate(stat = fct_reorder(stat, num_metabolites))

# write_csv(retained_table, here("aging_output_files", "retained_table.csv"))
retained_table <- read_csv(here("aging_output_files", "retained_table.csv"))
retained_metabolites_plot <- ggplot(retained_table, 
                                    aes(data, num_metabolites, fill = stat, label = pct_of_total)) + 
  geom_col(position = "dodge") + 
  ggrepel::geom_text_repel(position = position_dodge(width = 0.9), size = 8, 
                           vjust = -0.5, direction = "y",
                           segment.size = 0,
                           segment.alpha = 0) + 
  labs(title = "Nonzero Features over Imputations",
       x = "Data Source",
       y = "# of Metabolites") + 
  theme(legend.title = element_blank()) + 
  scale_fill_viridis_d(labels = c("avg_retained" = "Average retained", 
                                 "num_in_all" = "# in all", 
                                 "num_in_any" = "# in any"))
  

ggsave("retained_metabolites_bar.png",
       plot = retained_metabolites_plot, 
       path = here("plots", "aging_figs"),
       width = 18,
       height = 10)

 ########################

### Comparing redisualds bw models

#########################

# the residuals between models will be positively correlated because of regularization bias (not mean 0)
# here, we detrend the residuals against true age to better understand the correlation.

# note: these lm's are equiv to replacing the lm with lm(imp_avg ~ truth) directly
resid_corr_df <- select(got_pred_df, id, truth, "got" = imp_avg) %>%
  inner_join(select(lipids_pred_df, id, "lipids" = imp_avg), by = "id") %>%
  inner_join(select(targeted_pred_df, id, "targeted" =  imp_avg), by = "id") %>%
  inner_join(select(untargeted_pred_df, id, "untargeted" = imp_avg), by = "id") %>%
  transmute(Targeted = resid(lm((truth - targeted) ~ truth, data = .)),
            GOT = resid(lm((truth - got) ~ truth, data = .)),
            Untargeted = resid(lm((truth - untargeted) ~ truth, data = .)),
            Lipids = resid(lm((truth - lipids) ~ truth, data = .))
            ) 

# correlation test
resid_corr_df %>%
  summarize(lipid_got = (cor.test(GOT, Lipids, method = 'spearman'))$p.value,
         lipid_targeted = (cor.test(Lipids, Targeted, method = 'spearman'))$p.value,
         lipid_untargeted = (cor.test(Lipids, Untargeted, method = 'spearman'))$p.value,
         got_targeted = (cor.test(GOT, Targeted, method = 'spearman'))$p.value,
         got_untargeted = (cor.test(GOT, Untargeted, method = 'spearman'))$p.value,
         targeted_untargeted = (cor.test(Targeted, Untargeted, method = 'spearman'))$p.value
         ) %>%
  pivot_longer(everything(), 
               names_to = "pair", 
               values_to = "p_raw") %>%
  mutate(p_bh = p.adjust(p_raw, method = 'BH'))

resid_corr_matrix <- resid_corr_df %>%
  cor(method = "spearman") %>%
  round(2) 

write.csv(resid_corr_matrix, here("aging_tables", "resid_corr_matrix.csv"))

bind_rows("GOT" = got_pred_df, "Lipids" = lipids_pred_df, 
          "Targeted" = targeted_pred_df, "Untargeted" = untargeted_pred_df, 
          .id = "data") %>%
  select(data, id, imp_avg) %>%
  ggplot() + 
  geom_point(aes(y = id, x = imp_avg, color = data))



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



############################

# full models for coefs (targeted/lipids)

###########################

targeted_age_full <- purrr::map(1:5, ~age_control_analysis(imputed_c_targeted5, name = "Targeted, full", color = NULL, imp_num = .x, full_model = T))
lipids_age_full <- purrr::map(1:5, ~age_control_analysis(imputed_c_lipids5, name = "Lipids, full", color = NULL, imp_num = .x, full_model = T))

full_targeted_age_coefs <- get_importance_tables(targeted_age_full)

# write separate tables for positive and negative coefs
full_targeted_age_coefs %>%
  filter(`Avg Coef` > 0) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(str_glue(here("aging_tables", "{today()}-pos_coef_table_targeted.csv")))

full_targeted_age_coefs %>%
  filter(`Avg Coef` < 0) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(str_glue(here("aging_tables", "{today()}-neg_coef_table_targeted.csv")))

# full_targeted_ad_coefs %>%
#   filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
#   mutate(name = Name %>% str_to_lower() %>% str_trim(),
#          mapped_name = case_when(
#            #https://pubchem.ncbi.nlm.nih.gov/compound/12-Ketolithocholic-acid
#            name == "3\\?-hydroxy-12 ketolithocholic acid" ~ "12-Ketodeoxycholic acid",
#            #https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:27726 / https://www.genome.jp/dbget-bin/www_bget?cpd:C06560
#            name == "2-chloro-4,6-diamino-1,3,5-triazine" ~ "Deisopropyldeethylatrazine",
#            #http://www.hmdb.ca/metabolites/HMDB0006029
#            name == "acetyl-l-glutamine" ~ "N-Acetylglutamine",
#            #https://pubchem.ncbi.nlm.nih.gov/compound/Acetylornithine
#            name == "acetylornithine" ~  "N-Acetylornithine",
#            #https://pubchem.ncbi.nlm.nih.gov/compound/5-Aminovaleric-acid
#            name == "amino valerate" ~ "5-Aminopentanoic acid",
#            #https://pubchem.ncbi.nlm.nih.gov/compound/N_N-dimethylarginine
#            name == "dimethylarginine" ~ "Asymmetric dimethylarginine",
#            #https://pubchem.ncbi.nlm.nih.gov/compound/1826
#            name == "hiaa" ~ 	"5-Hydroxyindoleacetic acid",
#            TRUE ~ name)) %>%
#   select(mapped_name) %>%
#   write_tsv(here("ad_pd", "msea_ad_names_multivar_sig.txt"))


full_lipids_age_coefs <- get_importance_tables(lipids_age_full)
# write separate tables for positive and negative coefs
full_lipids_age_coefs %>%
  filter(`Avg Coef` > 0) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(str_glue(here("aging_tables", "{today()}-pos_coef_table_lipids.csv")))

full_lipids_age_coefs %>%
  filter(`Avg Coef` < 0) %>%
  filter(!(Name %in% c("(Intercept)", "Age", "GenderM"))) %>%
  write_csv(str_glue(here("aging_tables", "{today()}-neg_coef_table_lipids.csv")))



####################

## Mummichog

###################

## Untargeted --------------------------------------

imputed_c_untargeted_df <- imputed_c_untargeted5[[1]][[1]] %>%
  as_tibble() %>%
  dplyr::select(-any_of(c('(Intercept)', "GenderM")))

age_untargeted_p_values <- imputed_c_untargeted_df %>%
  names %>%
  setdiff('Age') %>%
  purrr::map(~age_metabolite_p(imputed_c_untargeted_df, .x)) %>%
  purrr::reduce(rbind) %>%
  as_tibble() %>%
  dplyr::rename('og_p_value' =  `Pr(>|t|)`) %>%
  dplyr::mutate(name = str_replace_all(name, "`", "")) %>%
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2)

age_untargeted_p_table <- cbind(age_untargeted_p_values, 
                                'bh_p_value' = p.adjust(age_untargeted_p_values$og_p_value, method = 'BH'))


mz_retention_untargeted <- raw_data_untargeted %>%
  group_by(Metabolite, Mode) %>%
  slice(1) %>%
  dplyr::select(Metabolite, Mode, `m/z`, `Retention time (min)`)

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_untargeted <- age_untargeted_p_table %>% 
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  #dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 't-score' = `t value`, Metabolite, Mode)

mummichog_untargeted %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("aging_output_files", str_glue('{lubridate::today()}_mummichog_untargeted_neg.txt')))

mummichog_untargeted %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("aging_output_files", str_glue('{lubridate::today()}_mummichog_untargeted_pos.txt')))


# mummichog 2 command: `mummichog -f {input_file} -o {output_name} -c 0.05 -m {neg/pos}`



mummichog_neg_plot <- mummichog_plot(here("aging_output_files", "1591847203.643466.06-10-2020_negResults", "tables", "mcg_pathwayanalysis_06-10-2020_negResults.tsv")) +
  labs(title = "Negative Mode")
ggsave("mummichog_neg.png",
       plot = mummichog_neg_plot, 
       path = here("plots", "aging_figs"),
       width = 24,
       height = 14)

mummichog_pos_plot <- mummichog_plot(here("aging_output_files", "1591847290.3476467.06-10-2020_posResults", "tables", "mcg_pathwayanalysis_06-10-2020_posResults.tsv")) +
  labs(title = "Positive Mode")
ggsave("mummichog_pos.png",
       plot = mummichog_pos_plot, 
       path = here("plots", "aging_figs"),
       width = 24,
       height = 14)


## GOT-MS ----------------------------------------

# using first imputation.
imputed_c_got_df <- imputed_c_got5[[1]][[1]] %>%
  as_tibble() %>%
  dplyr::select(-c('(Intercept)', GenderM))

age_got_p_values <- imputed_c_got_df %>%
  names %>%
  setdiff('Age') %>%
  purrr::map(~age_metabolite_p(imputed_c_got_df, .x)) %>%
  purrr::reduce(rbind) %>%
  as_tibble() %>%
  dplyr::rename('og_p_value' =  `Pr(>|t|)`) %>%
  mutate(name = str_replace_all(name, "Results_", "")) %>%
  tidyr::separate(col = name, into = c('Name','Mode'), sep = ' ')

age_got_p_table <- cbind(age_got_p_values, 
                                'bh_p_value' = p.adjust(age_got_p_values$og_p_value, method = 'BH'))


# need to find GOT mz / retention time beforew e can do this step

# mz_retention_got <- raw_data_got %>%
#   group_by(Metabolite, Mode) %>%
#   slice(1) %>%
#   dplyr::select(Metabolite, Mode, `m/z`, `Retention time (min)`)
# 
# # Left join because age_got has less columns. pull out only significant (.05 level)
# mummichog_got <- age_got_p_table %>% 
#   left_join(mz_retention_got, by = c('Metabolite', 'Mode')) %>%
#   dplyr::filter(bh_p_value < 0.05) %>%
#   dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 't-score' = `t value`, Metabolite, Mode)
# 


########################

# fit on full data for coefficients

########################
full_model_untargeted <- purrr::map(1:5, ~full_age_model_from_imputed(.x, imputed_c_untargeted5))




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

untargeted_adpd_separate_amelia5 <- readRDS(here("aging_output_files", "untargeted_adpd_separate_amelia5.Rds"))
untargeted_adpd_age_analysis_separate <- readRDS(here("aging_output_files", "untargeted_adpd_age_analysis_separate.Rds"))

# saveRDS(untargeted_adpd_separate_amelia5, here("aging_output_files", "untargeted_adpd_separate_amelia5.Rds"))
# saveRDS(untargeted_adpd_age_analysis_separate, here("aging_output_files", "untargeted_adpd_age_analysis_separate.Rds"))



untargeted_adpd_pred_df_separate <- imputation_df(untargeted_adpd_age_analysis_separate)

untargeted_adpd_avg_retained <- untargeted_adpd_age_analysis_separate %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

untargeted_adpd_in_all <- untargeted_adpd_age_analysis_separate %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

untargeted_adpd_in_at_least_one <- untargeted_adpd_age_analysis_separate %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")

untargeted_adpd_retained_table <- tibble(
  data = "untargeted_adpd",
  avg_retained = untargeted_adpd_avg_retained,
  num_in_all = length(untargeted_adpd_in_all),
  num_in_any = length(untargeted_adpd_in_at_least_one)
)


# will use for matched
adpd_null_pred <- mean(untargeted_adpd_pred_df_separate$truth) %>% rep(times = length(untargeted_adpd_pred_df_separate$truth))
adpd_rmse_null <- (untargeted_adpd_pred_df_separate$truth - adpd_null_pred)^2 %>% mean %>% sqrt %>% round(2)
adpd_mae_null <- (untargeted_adpd_pred_df_separate$truth - adpd_null_pred) %>% abs %>% mean %>% round(2)


untargeted_adpd_age_predtruth_separate <- predtruth_plot(untargeted_adpd_pred_df_separate, name = "Untargeted", data_name = "AD/PD, separate", color = "apoe4")
#ggsave(filename = "age_c_adpd_pred_untargeted_separate.png") #14.9 x 8.21

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

# saveRDS(untargeted_adpd_matched_controls, here("aging_output_files", "untargeted_adpd_matched_controls.Rds"))
untargeted_adpd_matched_controls <- readRDS(here("aging_output_files", "untargeted_adpd_matched_controls.Rds"))

untargeted_reduced <- readRDS(here('aging_output_files', 'untargeted_reduced.Rds'))
untargeted_c_matched_age_reduced <- untargeted_reduced$pred_df %>%
  filter(id %in% untargeted_adpd_matched_controls$Id)
# deprecated by reduced_anlaysis

# # imputation_df() expects a list of lists, but only uses the first element of the inner list. so we add a random "1" to the output. doesn't do anything.
# untargeted_c_matched_age_analysis_table <- purrr::map(untargeted_age_analysis, ~list(.x[[1]] %>% filter(id %in% untargeted_adpd_matched_controls$Id), 1))
# untargeted_c_matched_age_reduced <- imputation_df(untargeted_c_matched_age_analysis_table)
# 
# 
# untargeted_c_matched_age_predtruth <- predtruth_plot(untargeted_c_matched_age_reduced, name = "Untargeted (matched)", color = NULL)
# ggsave(filename = "age_clock_untargeted_pred_matched.png",
#        plot = untargeted_c_matched_age_predtruth,
#        width = 14,
#        height = 10,
#        path = here("plots", "aging_figs")) #14.9 x 8.21


adpd_matched_rsq <- cor(untargeted_c_matched_age_reduced$truth, untargeted_c_matched_age_reduced$imp_avg)^2 %>% round(2)
adpd_matched_rmse <- (untargeted_c_matched_age_reduced$truth - untargeted_c_matched_age_reduced$imp_avg)^2 %>% mean %>% sqrt %>% round(2)
adpd_matched_mae <- abs(untargeted_c_matched_age_reduced$truth - untargeted_c_matched_age_reduced$imp_avg) %>% mean %>% round(2)

# add indicator for points we want to highlight
untargeted_adpd_pred_df_separate <- untargeted_adpd_pred_df_separate %>%
  mutate(mse = (truth - imp_avg)^2) %>%
  group_by(type) %>%
  mutate(biggest_error = ifelse(mse == max(mse), 1,0)) %>%
  ungroup() %>%
  # for now, we aren't super interested in the AD outlier, so manually unhighlight that one
  mutate(biggest_error = ifelse(type == "AD", 0, biggest_error))

## Update regular ADPD plot with the stats for this one:
untargeted_adpd_age_predtruth_separate <- predtruth_plot(untargeted_adpd_pred_df_separate, name = "Untargeted AD/PD", data_name = "AD/PD, separate", color = "type") +
    geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, imp_avg, method = "pearson")^2 %>% round(2), 
                                                                              "<br>RMSE: ", (truth - imp_avg)^2 %>% mean %>% sqrt %>% round(2), " [",adpd_rmse_null,"]", 
                                                                              "<br>MAE: ", (truth - imp_avg) %>% abs %>% mean %>% round(2), " [",adpd_mae_null,"]")),
                  size =12) + 
    ggforce::geom_mark_circle(aes(truth, imp_avg, filter = biggest_error ==1), size = 1.5)
  

ggsave(filename = "age_adpd_pred_untargeted.png", 
       path = here("plots", "aging_figs"), 
       plot = untargeted_adpd_age_predtruth_separate,
       width = 16, height = 10) #14.9 x 8.21




### residuals plot, faceted by type

# first, we join the control preds with the ad/pd
untargeted_pred_df_with_adpd_oos <- untargeted_c_matched_age_reduced %>% 
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
targeted_ad_amelia5 <- filter_and_impute_multi(wide_data_targeted, c('AD'), empri = 5)
targeted_pd_amelia5 <- filter_and_impute_multi(wide_data_targeted, c('PD'), empri = 5)

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





###########

## Using targeted raw data, look at different levels of retained metabolites between matched controls and ad/pd

###########

# find matched controls
targeted_adpd_matched_controls <- purrr::map(1:nrow(wide_data_targeted), ~find_control(.x, data = filter(wide_data_targeted, Type %in% c("AD", "PD")), 
                                                                                           data_control = filter(wide_data_targeted, Type %in% c("CO", "CY", "CM")))) %>%
  bind_rows() %>%
  distinct(.keep_all = T)


targeted_adpd_retained <- wide_data_targeted %>%
  filter(Type %in% c("AD", "PD")) %>%
  select_at(vars(setdiff(targeted_coef_table$name, c("(Intercept)", "GenderM"))))

targeted_matched_c_retained <- targeted_adpd_matched_controls %>%
  select_at(vars(setdiff(targeted_coef_table$name, c("(Intercept)", "GenderM"))))
  
two_df_t_test <- function(df1, df2, name){
  first_col <- df1 %>% pull(name)
  second_col <- df2 %>% pull(name)

  # run a Welsh's t test (maybe better for unequal sample sizes since unequal var)
  result <- t.test(first_col, second_col, var.equal = FALSE)
  
  tibble(metabolite = name,
         t_stat = result$statistic,
         t_df = result$parameter,
         t_pval = result$p.value)
  
}

targeted_retained_t <- purrr::map_df(names(targeted_adpd_retained), 
                                        ~two_df_t_test(targeted_adpd_retained, targeted_matched_c_retained, .x)) %>%
  mutate(t_BH_pval = p.adjust(t_pval, method = "BH")) %>%
  arrange(t_BH_pval)

t_values_age_coefs <- left_join(targeted_retained_t, targeted_coef_table, by = c("metabolite" = "name")) %>%
  rowwise() %>%
  mutate(avg_coef = mean(c(coef1, coef2, coef3, coef4, coef5), na.rm = T)) %>%
  ungroup()

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



################

## Using targeted, separate imputation for AD/PD
## FIT ON ADPD data, pred on matched control

################

message("targeted Age model, fit on ADPD, predict on matched control -------------------------------------------")

targeted_reverse_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(targeted_adpd_separate_amelia5[[.x]], new_data = imputed_matched_c_targeted5[[.x]], nlambda = 200))

targeted_reverse_adpd_pred_df_separate <- imputation_df(targeted_reverse_adpd_age_analysis_separate)

(targeted_reverse_adpd_age_predtruth_separate <- predtruth_plot(targeted_reverse_adpd_pred_df_separate, name = "targeted", data_name = "Control, separate (fit on AD/PD)", color = NULL)) 
ggsave(filename = "age_c_reverse_adpd_pred_targeted_separate.png") #14.9 x 8.21



################################################################

### Fitting separate models for Males/Females

##############################################################

# Male first ------------------------

untar_male_rowindex <- which(imputed_c_untargeted5[[1]][[1]][,'GenderM'] == 1)
imputed_c_male_untargeted5 <- imputed_c_untargeted5 %>%
  purrr::map(function(x){
    x[[1]] <- x[[1]][untar_male_rowindex,]
    x[[2]] <- x[[2]][untar_male_rowindex]
    x[[3]] <- x[[3]][untar_male_rowindex]
    x[[4]] <- x[[4]][untar_male_rowindex]
    x
  })
male_untargeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_male_untargeted5, name = "male_untargeted", color = NULL, imp_num = .x, nlambda = 200))
male_untargeted_age_analysis <- readRDS( here("aging_output_files", "male_untargeted_age_analysis.Rds"))
# saveRDS(male_untargeted_age_analysis, here("aging_output_files", "male_untargeted_age_analysis.Rds"))


male_untargeted_pred_df <- imputation_df(male_untargeted_age_analysis)

male_untargeted_age_predtruth <- predtruth_plot(male_untargeted_pred_df, name = "Untargeted (Male)")
ggsave(filename = "age_clock_male_untargeted_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = male_untargeted_age_predtruth,
       width = 14, height = 10)

male_untargeted_avg_retained <- male_untargeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

male_untargeted_in_all <- male_untargeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

male_untargeted_in_at_least_one <- male_untargeted_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")

male_untargeted_retained_table <- tibble(
  data = "Untargeted (male)",
  avg_retained = male_untargeted_avg_retained,
  num_in_all = length(male_untargeted_in_all),
  num_in_any = length(male_untargeted_in_at_least_one)
)


# Female model -----------------------------------------

untar_female_rowindex <- which(imputed_c_untargeted5[[1]][[1]][,'GenderM'] == 0)
imputed_c_female_untargeted5 <- imputed_c_untargeted5 %>%
  purrr::map(function(x){
    x[[1]] <- x[[1]][untar_female_rowindex,]
    x[[2]] <- x[[2]][untar_female_rowindex]
    x[[3]] <- x[[3]][untar_female_rowindex]
    x[[4]] <- x[[4]][untar_female_rowindex]
    x
  })
female_untargeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_female_untargeted5, name = "female_untargeted", color = NULL, imp_num = .x, nlambda = 200))
female_untargeted_age_analysis <- readRDS(here("aging_output_files", "female_untargeted_age_analysis.Rds"))
# saveRDS(female_untargeted_age_analysis, here("aging_output_files", "female_untargeted_age_analysis.Rds"))

female_untargeted_pred_df <- imputation_df(female_untargeted_age_analysis)

female_untargeted_age_predtruth <- predtruth_plot(female_untargeted_pred_df, name = "Untargeted (female)")
ggsave(filename = "age_clock_female_untargeted_pred.png", 
       path = here("plots", "aging_figs"), 
       plot = female_untargeted_age_predtruth,
       width = 14, height = 10)

female_untargeted_avg_retained <- female_untargeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

female_untargeted_in_all <- female_untargeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff("(Intercept)")

female_untargeted_in_at_least_one <- female_untargeted_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff("(Intercept)")

female_untargeted_retained_table <- tibble(
  data = "Untargeted (Female)",
  avg_retained = female_untargeted_avg_retained,
  num_in_all = length(female_untargeted_in_all),
  num_in_any = length(female_untargeted_in_at_least_one)
)

# compare with random sample of same size (41 female, 44 males). let's do 41 first


#' quick function to do control age analysis on untargeted data using controls
#' with a random sample of size (a positive number smaller than 85 bc there are 85 controls)
#' 
age_model_c_random_sample <- function(size, seed){
  set.seed(seed)
  untar_random_rowindex <- sample(1:85, size, replace = F)
  imputed_c_random_untargeted5 <- imputed_c_untargeted5 %>%
    purrr::map(function(x){
      x[[1]] <- x[[1]][untar_random_rowindex,]
      x[[2]] <- x[[2]][untar_random_rowindex]
      x[[3]] <- x[[3]][untar_random_rowindex]
      x[[4]] <- x[[4]][untar_random_rowindex]
      x
    })
  random_untargeted_age_analysis <- purrr::map(1:5, ~age_control_analysis(imputed_c_random_untargeted5, name = str_glue("Untargeted (n = {size})"), color = NULL, imp_num = .x, nlambda = 200))
  
  random_untargeted_pred_df <- imputation_df(random_untargeted_age_analysis)
  
  random_untargeted_age_predtruth <- predtruth_plot(random_untargeted_pred_df, name = str_glue("Untargeted (n = {size})"))
  
  list(analysis = random_untargeted_age_analysis, 
       pred_df = random_untargeted_pred_df,
       pred_plot = random_untargeted_age_predtruth,
       rmse = with(random_untargeted_pred_df, sqrt(mean((truth - imp_avg)^2))))
}

random_sample_modeling <- tibble(
  sex = c(rep("Male", times = 10), rep("Female", times = 10))
)
random_sample_modeling <- random_sample_modeling %>%
  mutate(n = ifelse(sex == "Male", 44, 41),
         output = purrr::imap(n, ~age_model_c_random_sample(.x, .y))
         ) %>%
  mutate(output_plot = purrr::map(output, ~list(.x[['pred_plot']])),
         output_rmse = purrr::map(output, ~yardstick::rmse(.x[['pred_df']], truth = truth, estimate = imp_avg)))

# saveRDS(random_sample_modeling, here("aging_output_files", "random_sample_modeling.Rds"))

random_sample_modeling %>% 
  mutate(rmse = purrr::map(output, ~.x$rmse)) %>% 
  select(sex, rmse) %>% 
  unnest(cols = c(rmse)) %>% 
  group_by(sex) %>%
  summarize(mean(rmse))

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

c_gender_logistic <- purrr::map(1:5, ~csf_analysis(imputed_c_untargeted5, varname = "GenderM", imp_num = .x, nlambda = 200))

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

## The reference list must be IDs, so we need to map the names to KEGG/HMDB names.
# Post mortem, I'm coming back and renaming ones that I can identify (but that atempt 2 couldn't)
msea_table <-c_targeted_univariate_table %>%
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

## We need a reference list (ie include names on insignificant variables)
univar_targeted_names_all <-msea_table %>%
  select(mapped_name)
write_tsv(univar_targeted_names_all, here("aging_tables", "msea_targeted_names_ref.txt"))

# the list of significant ones
univar_targeted_names_sig <- msea_table %>%
  filter(bh_p_value < 0.05) %>%
  select(mapped_name)
write_tsv(univar_targeted_names_sig, here("aging_tables", "msea_targeted_names_sig.txt"))



# trying the matching tool on metaboAnalyst website -- based on univar_targeted_names_all
new_msea_ref <- read_csv(here("aging_tables", "msea_name_map.csv")) %>%
  mutate(Name = ifelse(is.na(Match), Query, Match)) %>%
  select(Name)
write_tsv(new_msea_ref, here("aging_tables", "msea_targeted_names_ref2.txt"))


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






### once we have the msea result table, recreate the enrichment bar plot in ggplot
msea_smpdb <- read_csv(here("aging_tables", "msea_ora_SMPDB_result.csv")) %>%
  mutate(bh_p = p.adjust(`Raw p`, method = "BH"),
         set = fct_reorder(X1, hits/total)) %>%
  filter(total > 2) %>%
  top_n(10, hits/total) #%>%
  #select("Set" = X1, "Total" = total, "Expected" = expected, "Hits" = hits)

write_csv(msea_smpdb, here("aging_tables", "msea_SMPDB_result_top10.csv"))

msea_smpdb_plot <- ggplot(msea_smpdb) + 
  geom_col(aes(set, hits/total)) + 
  geom_text(aes(set, hits/total, label = hits), size = 8, hjust = -1) + 
  # flip default scale
  #scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  coord_flip() + 
  labs(x = "",
       y = "# Hits / # Total",
       title = "SMPDB"
       )

ggsave("msea_smpdb.png",
       plot = msea_smpdb_plot,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)

### once we have the msea result table, recreate the enrichment bar plot in ggplot
msea_csf <- read_csv(here("aging_tables", "msea_ora_csf_result.csv")) %>%
  mutate(bh_p = p.adjust(`Raw p`, method = "BH"),
         set = fct_reorder(X1, hits/total)) %>%
  filter(total >= 2) %>%
  top_n(10, hits/total) #%>%
  # select("Set" = X1, "Total" = total, "Expected" = expected, "Hits" = hits)
write_csv(msea_csf, here("aging_tables", "msea_csf_result_top10.csv"))

msea_csf_plot <- ggplot(msea_csf) + 
  geom_col(aes(set, hits/total)) + 
  geom_text(aes(set, hits/total, label = hits), size = 8, hjust = -0.5) + 
  # flip default scale
  #scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  coord_flip() + 
  labs(x = "",
       y = "# Hits / # Total",
       title = "CSF Disease Library"
  )

ggsave("msea_csf.png",
       plot = msea_csf_plot,
       path = here("plots", "aging_figs"),
       width = 20,
       height = 10)




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
    labs(title = "Fitted Concentration as a Function of Age",
         subtitle = str_glue("Untargeted, Random Forest: Top {n}"),
         y = "Predicted Concentration") +
    ylim(c(-2.5,2.5))
  
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
c_conc_rf_top50 <- concentration_rf_plot(untargeted_univariate_c_conc_age_scale_rf, n = 50)
ggsave("fitted_conc_by_age_rf.png",
       plot = c_conc_rf_top50$flip_shift,
       width = 14,
       height = 10,
       path = here("plots", "aging_figs")
       )

# using same metabolites as full in c
full_conc_names <- full_conc_rf_top10$data %>% pull(name) %>% unique
c_conc_rf_fullnames <- concentration_rf_plot(untargeted_univariate_c_conc_age_scale_rf, names = full_conc_names)



######

## Missingness

#####


data_sources <- list(wide_data, wide_data_lipids, wide_data_targeted, wide_data_untargeted)
data_names <- list("GOT", "Lipids", "Targeted", "Untargeted")
pct_missing_data <- purrr::map2(data_sources, data_names, ~percent_missing_by_age(.x, .y)) %>%
  bind_rows %>%
  gather(key = "age_group", value = "pct_missing", -source) %>%
  mutate(age_group = fct_relevel(age_group, "[20,43]", "(43,65]", "(65,88]"))
  #mutate(Type = factor(Type, levels = c("CY", "CM", "CO", "AD", "PD"))) %>%
  #filter(Type %in% c("CY", "CM", "CO"))


missingness_overview_plot <- ggplot(pct_missing_data) + 
  geom_col(aes(x = source, y= pct_missing, fill = age_group), position = 'dodge') +
  labs(title = "Percent Missingness By Dataset",
       x = "Dataset",
       y = "Percent Missing",
       fill = "Age") + 
  scale_fill_viridis_d() + 
  theme(legend.text = element_text(size = 30))
ggsave("missingness_overview_bar.png",
       plot = missingness_overview_plot,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)



# Run chi squared test on missingness between 3 groups of controls (equal range)
#' @param data is one of the wide_data_* variants
#' @param dataset_name is what we call the dataset for the plot title, ie "Lipids"
#' @param cutoff is a value between 0 and 1, the BH cutoff p values to show
missingness_by_type_plot <- function(data, dataset_name, cutoff = 0.05){
  
  #create 3 equal length age groups for controls
  c_data <- data %>%
    filter(Type %in% c("CY", "CM", "CO")) %>%
    mutate(age_group = cut_interval(Age, n = 3, dig.lab = 2))
  
  n_by_type_untargeted <- c_data %>%
    group_by(age_group) %>%
    tally() %>%
    spread(key = 'age_group', value = 'n') %>%
    rename_all(function(x) paste('n', x, sep = '_'))
  
  missingness_by_type_untargeted_counts <- c_data %>%
    group_by(age_group) %>%
    group_map(~ map_int(.x, function(y) sum(is.na(y)))) %>%
    set_names(c_data$age_group %>% droplevels %>% levels) %>%
    lapply(function(x) enframe(x, name = 'name', value = 'num_missing')) 
  
  missingness_by_type_untargeted_pct <- c_data %>% 
    group_by(age_group) %>%
    group_map(~ map_dbl(.x, function(y) round(sum(is.na(y))/length(y), digits = 3))) %>%
    set_names(c_data$age_group %>% droplevels %>% levels) %>%
    lapply(function(x) enframe(x, name = 'name', value = 'pct_missing'))
  
  
  
  missingness_by_type_all <- purrr::reduce(missingness_by_type_untargeted_counts, inner_join, by = 'name') %>%
    #set the names to show type
    set_names(c('name', paste('num_missing',levels(droplevels(c_data$age_group)), sep = '_'))) %>%
    dplyr::filter(!(name %in% c('GBAStatus', 'cognitive_status', 'GBA_T369M'))) %>%
    cbind(n_by_type_untargeted) %>%
    dplyr::mutate(
      `pct_missing_[20,42]` = `num_missing_[20,42]`/`n_[20,42]`,
      `pct_missing_(42,64]` = `num_missing_(42,64]`/`n_(42,64]`,
      `pct_missing_(64,86]` = `num_missing_(64,86]`/`n_(64,86]`
    ) %>%
    #filter(reduce(list(pct_missing_AD, pct_missing_(42,64],pct_missing_(64,86],pct_missing_[20,42],pct_missing_PD), `==`)) %>%
    rowwise() %>%
    #mutate(p_value = (prop.test(x =  str_subset(names(.), 'num_missing'), n = str_subset(names(.), 'n_')))$p.value)
    dplyr::mutate(p_value = (prop.test(x = c(`num_missing_(42,64]`, `num_missing_(64,86]`, `num_missing_[20,42]`), 
                                       n = c(`n_(42,64]`,`n_(64,86]`, `n_[20,42]`)))$p.value) %>%
    cbind('bh_q_value' = p.adjust(.$p_value, method = 'BH')) %>%
    dplyr::filter(bh_q_value < cutoff) %>%
    gather('age_group', 'pct_missing', contains('pct_missing')) %>%
    mutate(age_group = factor(age_group, levels = c("pct_missing_[20,42]","pct_missing_(42,64]","pct_missing_(64,86]")))
  
  
  ggplot(missingness_by_type_all, aes(name, pct_missing)) +
    geom_col(aes(fill = age_group), position = "dodge") + 
    #geom_line() + 
    #geom_point(aes(color = Type), size = 3, position = position_jitter(width = 0.01, height = 0,seed = 1)) + 
    scale_fill_viridis_d(labels = c("[20,42]", "(42,64]", "(64,86]")) + 
    scale_y_continuous(labels = percent) + 
    theme(axis.text.x = element_text(angle= 90, hjust =1, size = 24)) + 
    labs(title = str_glue('Percent Missingness by Age: {dataset_name}'),
         subtitle = str_glue('Filtered BH q < {cutoff} (2-tailed Pearson Chi Squared Test)'),
         x = "Name",
         y = "Percent Missing",
         fill = "Age")
}


sig_missingness_lipids_plot <- missingness_by_type_plot(wide_data_lipids, "Lipids", cutoff = 0.05)
sig_missingness_untargeted_plot <- missingness_by_type_plot(wide_data_untargeted, "Untargeted", cutoff = 0.01)
ggsave("sig_missingness_lipids_bar.png",
       plot = sig_missingness_lipids_plot,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)

ggsave("sig_missingness_untar_bar.png",
       plot = sig_missingness_untargeted_plot,
       path = here("plots", "aging_figs"),
       width = 20,
       height = 10)



ggplot(missingness_by_type_all, aes(Type, pct_missing)) + 
  geom_line(color = name)



### we see that all 6 significant lipids are triglycerides. if 6 lipids were randomly selected, how likely is this outcome?
# number of triglycerides in the dataset
num_tri <- names(wide_data_lipids) %>%
  str_sub(1,3) %>%
  table() %>%
  magrittr::extract2('TAG')

# number of lipids in the dataset
num_lipids <- wide_data_lipids %>%
  select(-any_of(metadata_cols)) %>%
  ncol()
  
# if we treat this as a hypergeom experiment with 6 trials, the expected value is
expected_tri <- 6 * num_tri / num_lipids

# to get the x-fold enrichment, take the observed value / expected value. we observed all 6 tri
tri_fold_enrich <- 6 / expected_tri



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

# get just the age of controls. 
  # it's a little lazy of me to split up the features and response without an id col
  # I'm basically relying on the fact that filter() will act in a consistent way, so ordering is preserved
untargeted_c_age <- wide_data_untargeted %>% 
  filter(Type %in% c("CO", "CM", "CY")) %>%
  pull(Age)

# get the features
untargeted_c_na <- wide_data_untargeted %>%
  filter(Type %in% c("CO", "CM", "CY")) %>%
  select(-any_of(metadata_cols)) %>% 
  mutate_all(~ifelse(is.na(.x), 1, 0)) %>%
  as.matrix()
  

untargeted_c_na_full_model <- get_full_model(features= untargeted_c_na, untargeted_c_age, 
                                             alpha = 0.5, family = "gaussian", 
                                             penalize_AD_PD = FALSE, penalize_age_gender = FALSE)
# Note: If AD_ind, PD_ind are missing from the dataset, the flag penalize_AD_Pd doesn't do anything
untargeted_c_na_fitpred <- lapply(1:nrow(untargeted_c_na), function(x) loo_cvfit_glmnet(x, untargeted_c_na, untargeted_c_age, 
                                                                                        lambda = untargeted_c_na_full_model[[2]], 
                                                                                        full_fit = untargeted_c_na_full_model[[1]],
                                                                                        alpha = 0.5, family = 'gaussian', 
                                                                                        penalize_age_gender = FALSE, penalize_AD_PD = FALSE)
                                  )

untargeted_c_na_fitpred <- readRDS(here("aging_output_files", "untargeted_c_na_fitpred.Rds"))
# saveRDS(untargeted_c_na_fitpred, here("aging_output_files", "untargeted_c_na_fitpred.Rds"))

untargeted_c_na_pred <- lapply(untargeted_c_na_fitpred, function(x) x[[2]]) %>% unlist

# get mean metrics for the plot
untar_null_pred <- mean(untargeted_c_age) %>% rep(times = length(untargeted_c_age))
untar_rmse_null <- (untargeted_c_age - untar_null_pred)^2 %>% mean %>% sqrt %>% round(2)
untar_mae_null <- (untargeted_c_age - untar_null_pred) %>% abs %>% mean %>% round(2)

missingness_matrix_predtruth <- ggplot(tibble(truth = untargeted_c_age, 
              pred = untargeted_c_na_pred)
       ) +
  geom_point(aes(truth, pred), size = 4.5, position = position_dodge(width = .3)) +
  labs(title = paste0('Untargeted Missingness Matrix'),
       x = 'Chronological Age',
       y = 'Predicted Age') +
  geom_abline(intercept = 0, slope = 1) +
  expand_limits(x = 15, y = c(15, 100)) +
  geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, pred, method = "pearson")^2 %>% round(2),
                                                                            "<br>RMSE: ", (truth - pred)^2 %>% mean %>% sqrt %>% round(2), " (",untar_rmse_null,")",
                                                                            "<br>MAE: ", (truth - pred) %>% abs %>% mean %>% round(2), " (",untar_mae_null,")")),
                size = 12)

ggsave("missingness_matrix_predtruth.png",
       plot = missingness_matrix_predtruth,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)



#### Lipids missingness matrix ---------------------------------------

# get just the age of controls. 
# it's a little lazy of me to split up the features and response without an id col
# I'm basically relying on the fact that filter() will act in a consistent way, so ordering is preserved
lipids_c_age <- wide_data_lipids %>% 
  filter(Type %in% c("CO", "CM", "CY")) %>%
  pull(Age)

# get the features
lipids_c_na <- wide_data_lipids %>%
  filter(Type %in% c("CO", "CM", "CY")) %>%
  select(-any_of(metadata_cols)) %>% 
  mutate_all(~ifelse(is.na(.x), 1, 0)) %>%
  as.matrix()

# look at the distribution of columsn
count(enframe(colMeans(lipids_c_na)), value, sort = T)

#lipids_c_na_full_model <- get_full_model(features= lipids_c_na, lipids_c_age, alpha = 0.5, family = "gaussian", penalize_AD_PD = FALSE, penalize_age_gender = FALSE, nlambda = 200)
# Note: If AD_ind, PD_ind are missing from the dataset, the flag penalize_AD_Pd doesn't do anything
lipids_c_na_fitpred <- lapply(1:nrow(lipids_c_na), function(x) loo_cvfit_glmnet(x, lipids_c_na, lipids_c_age,
                                                                                alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE, penalize_AD_PD = FALSE, nlambda = 200))

lipids_c_na_fitpred <- readRDS(here("aging_output_files", "lipids_c_na_fitpred.Rds"))
# saveRDS(lipids_c_na_fitpred, here("aging_output_files", "lipids_c_na_fitpred.Rds"))

lipids_c_na_pred <- lapply(lipids_c_na_fitpred, function(x) x[[2]]) %>% unlist

# get mean metrics for the plot
lipids_null_pred <- mean(lipids_c_age) %>% rep(times = length(lipids_c_age))
lipids_rmse_null <- (lipids_c_age - lipids_null_pred)^2 %>% mean %>% sqrt %>% round(2)
lipids_mae_null <- (lipids_c_age - lipids_null_pred) %>% abs %>% mean %>% round(2)

missingness_matrix_lipids_predtruth <- ggplot(tibble(truth = lipids_c_age, 
                                              pred = lipids_c_na_pred)
) +
  geom_point(aes(truth, pred), size = 2.5, position = position_dodge(width = .2)) +
  labs(title = paste0('Control: True vs Predicted Age'),
       subtitle = paste0("Missingness Matrix"),
       x = 'True Age',
       y = 'Predicted Age') +
  geom_abline(intercept = 0, slope = 1) +
  expand_limits(x = 15, y = c(15, 100)) +
  geom_richtext(aes(x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste0("R^2: ", cor(truth, pred, method = "pearson")^2 %>% round(2),
                                                                            "<br>RMSE: ", (truth - pred)^2 %>% mean %>% sqrt %>% round(2), " (",lipids_rmse_null,")",
                                                                            "<br>MAE: ", (truth - pred) %>% abs %>% mean %>% round(2), " (",lipids_mae_null,")")),
                size = 12)

ggsave("missingness_matrix_lipids_predtruth.png",
       plot = missingness_matrix_lipids_predtruth,
       path = here("plots", "aging_figs"),
       width = 14,
       height = 10)









#sink(type="message")
#close(error_log)




