source(here::here("analysis", "starter.R"))



## Before we get into it, here's a quick distribution of age. could've picked any dataset -- same subjects
age_dist_ad_pd <- ggplot(wide_data_targeted %>% collapse_controls(), aes(Age)) +
  geom_histogram(bins = 10) +
  facet_wrap(~Type) + 
  labs(x = "Age",
       y = "Count",
       title = "Distribution of Age")

ggsave("age_dist.png",
       plot = age_dist_ad_pd,
       path = here("plots", "adpd_figs"),
       width= 14,
       height = 10)




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
pca_all_untargeted <- untargeted_all_processed %>% 
  dplyr::select(-any_of(metadata_cols)) %>%
  opls(algoC = "nipals",
       parAsColFcVn = untargeted_c_processed$Type)


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


############################

### AD/PD Logisitic Regression on untargeted ###

############################

message("Untargeted ADPD logistic -------------------------------------------")

#### AD First -----------------------------
# we want an AD indicator, but not a PD one
# untargeted_all_amelia5_ad_ind <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = TRUE, add_PD_ind = FALSE))
# untargeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(untargeted_all_amelia5_ad_ind, varname ="AD_ind", imp_num = .x, nlambda = 200))
untargeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))
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
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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



#### PD -------------------------------------
# we want a PD indicator, but not an AD one
# untargeted_all_amelia5_pd_ind <- purrr::map2(imputed_c_untargeted5, untargeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                       add_AD_ind = FALSE, add_PD_ind = TRUE))
# untargeted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(untargeted_all_amelia5_pd_ind, varname ="PD_ind", imp_num = .x, nlambda = 200))

untargeted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_untargeted5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T))
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
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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
targeted_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))
targeted_ad_logistic <- readRDS(here("ad_pd", "targeted_ad_logistic.Rds"))
# saveRDS(targeted_ad_logistic, here("ad_pd", "targeted_ad_logistic.Rds"))

# # compare with all
# targeted_ad_logistic_all <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))

targeted_ad_roc_table <- targeted_ad_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

targeted_ad_roc_plot <- ggplot(targeted_ad_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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
targeted_ad_coefs <- targeted_ad_logistic %>% 
  purrr::map(~.x[[1]] %>% enframe(value = "coef")) %>%
  reduce(full_join, by = "name") %>%
  rename("imp1" = 2, "imp2" = 3, "imp3" = 4, "imp4" = 5, "imp5" = 6) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.x), 0, .x)) %>%
  rowwise() %>%
  mutate(mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef)))

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

targeted_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_targeted5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T))
targeted_pd_logistic <- readRDS(here("ad_pd", "targeted_pd_logistic.Rds"))
# saveRDS(targeted_pd_logistic, here("ad_pd", "targeted_pd_logistic.Rds"))

# targeted_pd_logistic %>%
#   purrr::map(~.x[[2]]) %>%
#   cowplot::plot_grid(plotlist = .)

targeted_pd_roc_table <- targeted_pd_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

targeted_pd_roc_plot <- ggplot(targeted_pd_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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
  mutate(mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef)))


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


lipids_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))
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
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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
  mutate(mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef)))

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

lipids_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_lipids5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T))
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
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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
  mutate(mean_coef = mean(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(mean_coef)))


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


got_ad_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_got5, varname ="AD_ind", imp_num = .x, nlambda = 200, AD_ind = T))
got_ad_logistic <- readRDS(here("ad_pd", "got_ad_logistic.Rds"))
# saveRDS(got_ad_logistic, here("ad_pd", "got_ad_logistic.Rds"))


got_ad_roc_table <- got_ad_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

got_ad_roc_plot <- ggplot(got_ad_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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

got_pd_logistic <- purrr::map(1:5, ~logistic_control_analysis(imputed_all_got5, varname ="PD_ind", imp_num = .x, nlambda = 200, PD_ind = T))
got_pd_logistic <- readRDS(here("ad_pd", "got_pd_logistic.Rds"))
# saveRDS(got_pd_logistic, here("ad_pd", "got_pd_logistic.Rds"))

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
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2) + 
  labs(title = "ROC: PD vs C",
       subtitle = TeX('got'),
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

### Univariate AD/PD Logistic Regresion on untargeted

### For use with mummichog

##########################

mz_retention_untargeted <- raw_data_untargeted %>%
  group_by(Metabolite, Mode) %>%
  slice(1) %>%
  dplyr::select(Metabolite, Mode, `m/z`, `Retention time (min)`)



######## AD ----------------------------------------------------------

# add ad indicator
untargeted_all_amelia5_ad_ind <- imputed_all_untargeted5 %>%
  purrr::map(~cbind(.x[[1]], "AD_ind" = ifelse(imputed_all_untargeted5[[1]][[2]] == "AD", 1, 0)) %>% list())


untargeted_univar_ad_logistic <- bh_univariate_age(untargeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE)
untargeted_univar_ad_logistic <- readRDS(here("ad_pd", "untargeted_univar_ad_logistic.Rds"))
# saveRDS(untargeted_univar_ad_logistic, here("ad_pd", "untargeted_univar_ad_logistic.Rds"))

# Left join because age_untargeted has less columns. pull out only significant (.05 level)
mummichog_ad_untargeted <- untargeted_univar_ad_logistic %>% 
  tidyr::separate(col = name, into = c('Metabolite1','Metabolite2', 'Mode'), sep = '_') %>%
  tidyr::unite(col = "Metabolite", Metabolite1, Metabolite2) %>%
  left_join(mz_retention_untargeted, by = c('Metabolite', 'Mode')) %>%
  #dplyr::filter(bh_p_value < 0.05) %>%
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = og_p_value, 'z-score' = `z value`, Metabolite, Mode)

mummichog_ad_untargeted %>%
  dplyr::filter(Mode == 'neg') %>%
  write_tsv(here("ad_pd", 'mummichog_ad_neg.txt'))

mummichog_ad_untargeted %>%
  dplyr::filter(Mode == 'pos') %>%
  write_tsv(here("ad_pd", 'mummichog_ad_pos.txt'))



######## PD ----------------------------------------------------------

# add pd indicator
untargeted_all_amelia5_pd_ind <- imputed_all_untargeted5 %>%
  purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(imputed_all_untargeted5[[1]][[2]] == "PD", 1, 0)) %>% list())


untargeted_univar_pd_logistic <- bh_univariate_age(untargeted_all_amelia5_pd_ind, var = "PD_ind", family = "binomial", conc = FALSE)
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

############################

### Univariate AD/PD Logisitic Regression on targeted ###
### for use with MSEA

############################

message("Targeted Univariate ADPD Logistic -------------------------------------------")

#### AD First -----------------------------
# we want an AD indicator, but not a PD one

#empri = .1(num PD/AD)
# targeted_ad_amelia5 <- filter_and_impute_multi(wide_data_targeted, c('AD'), empri = 6)
# targeted_pd_amelia5 <- filter_and_impute_multi(wide_data_targeted, c('PD'), empri = 5)
# 
# targeted_adpd_separate_amelia5 <- purrr::map2(targeted_ad_amelia5, targeted_pd_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE))
# 
# targeted_all_amelia5_ad_ind <- purrr::map2(imputed_c_targeted5, targeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
#                                                                                                                 add_AD_ind = TRUE, add_PD_ind = FALSE))
# removing all columns with >10% missingness (to help imputation)
# colnames_less10perct_targeted <- wide_data_targeted %>% 
#   map_dbl(~ sum(is.na(.x))/nrow(wide_data_targeted) < .1) %>%
#   names() %>%
#   append(c("GBAStatus", "GBA_T369M"))
# 
# #### TODO: NEED TO PLAY WITH EMPRI!!! IT SHOULDN'T BE THIS HIGH
# imputed_less10perct_targeted <- wide_data_targeted %>%
#   select_at(vars(colnames_less10perct_targeted)) %>%
#   filter_and_impute_multi(types = c("CO", "CY", "CM", "AD", "PD"), empri = 500) 
# 
# # add ad_indicator. note that we have to add the extra "list()" at the end because this removes metadata, and bh_univariate_age expects a list of lists
# targeted_all_amelia5_ad_ind <- imputed_less10perct_targeted %>%
#   purrr::map(~cbind(.x[[1]], "AD_ind" = ifelse(imputed_less10perct_targeted[[1]][[2]] == "AD", 1, 0)) %>% list())

targeted_all_amelia5_ad_ind <- imputed_all_targeted5 %>%
  purrr::map(~cbind(.x[[1]], "AD_ind" = ifelse(imputed_all_targeted5[[1]][[2]] == "AD", 1, 0)) %>% list())


# Create table with bh-corrected p values
targeted_univar_ad_logistic <- bh_univariate_age(targeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE)
targeted_univar_ad_logistic <- readRDS(here("ad_pd", "targeted_univar_ad_logistic.Rds"))
# saveRDS(targeted_univar_ad_logistic, here("ad_pd", "targeted_univar_ad_logistic.Rds"))


targeted_univar_ad_logistic %>%
  filter(bh_p_value < 0.05)




#### PD -------------------------------------
# # we want a PD indicator, but not an AD one
# targeted_all_amelia5_pd_ind <- purrr::map2(imputed_c_targeted5, targeted_adpd_separate_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE, 
# add_AD_ind = FALSE, add_PD_ind = TRUE))

# add ad_indicator. note that we have to add the extra "list()" at the end because this removes metadata, and bh_univariate_age expects a list of lists
# targeted_all_amelia5_pd_ind <- imputed_less10perct_targeted %>%
#   purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(imputed_less10perct_targeted[[1]][[2]] == "PD", 1, 0)) %>% list())

targeted_all_amelia5_pd_ind <- imputed_all_targeted5 %>%
  purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(imputed_all_targeted5[[1]][[2]] == "PD", 1, 0)) %>% list())


# Create table with bh-corrected p values
targeted_univar_pd_logistic <- bh_univariate_age(targeted_all_amelia5_pd_ind, var = "PD_ind", family = "binomial", conc = FALSE)
targeted_univar_pd_logistic <- readRDS(here("ad_pd", "targeted_univar_pd_logistic.Rds"))
# saveRDS(targeted_univar_pd_logistic, here("ad_pd", "targeted_univar_pd_logistic.Rds"))

targeted_univar_pd_logistic %>%
  filter(bh_p_value < 0.05)


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
# isolating just the PD subjects
targeted_imputed_all_pd_index <- which(imputed_all_targeted5[[1]][[2]] == "PD")
targeted_imputed_all_pd_features <- purrr::map(imputed_all_targeted5, ~.x[[1]][targeted_imputed_all_pd_index,])
targeted_imputed_all_pd_metadata <- purrr::map(1:5, function(x) purrr::map(2:length(imputed_all_targeted5), function(y) imputed_all_targeted5[[x]][[y]][targeted_imputed_all_pd_index]))

# get data in form of filter_and_impute_multi output
targeted_imputed_all_pd <- purrr::map(1:5, ~append(targeted_imputed_all_pd_metadata[[.x]], list(targeted_imputed_all_pd_features[[.x]]), after = 0))


targeted_gba_logistic <- purrr::map(1:5, ~logistic_control_analysis(targeted_imputed_all_pd, varname ="GBA Status", imp_num = .x, nlambda = 200, GBA_ind = T))
targeted_gba_logistic <- readRDS(here("ad_pd", "targeted_gba_logistic.Rds"))
# saveRDS(targeted_gba_logistic, here("ad_pd", "targeted_gba_logistic.Rds"))


targeted_gba_roc_table <- targeted_gba_logistic %>%
  purrr::map(~.x$roc_df) %>%
  bind_rows(.id = "imp")

targeted_gba_roc_plot <- ggplot(targeted_gba_roc_table) + 
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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
  mutate(median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(median_coef)))



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
  geom_line(mapping = aes(fpr, tpr, group = imp), lwd = 1.5) + 
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
  mutate(median_coef = median(c(imp1, imp2, imp3, imp4, imp5))) %>%
  arrange(desc(abs(median_coef)))


lipids_gba_avg_retained <- lipids_gba_logistic %>% 
   purrr::map(~.x[[1]] %>% setdiff("(Intercept)") %>% length) %>% 
  unlist %>% mean

lipids_gba_in_all <- lipids_gba_logistic %>% 
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






