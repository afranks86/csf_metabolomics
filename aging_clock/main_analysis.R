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
       path = here("plots", "aging_figs", 'main_figs'),
       width= 83,
       height = 50,
       dpi = 300, units = 'mm')

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

pca_c_combined_x <- combined_c_processed %>% 
  dplyr::select(-any_of(metadata_cols))

pca_c_combined <-  pca_c_combined_x %>%
  opls(algoC = "nipals",
       parAsColFcVn = combined_c_processed$Gender)

pca_c_combined <- readRDS(here("aging_output_files", "pca_c_combined.Rds"))
# saveRDS(pca_c_combined, here("aging_output_files", "pca_c_combined.Rds"))


pca_scores_c_combined <- pca_c_combined@scoreMN %>%
  as_tibble() %>%
  bind_cols(combined_c_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id))

pc1_var_explained_combined <- pca_c_combined@modelDF["p1","R2X"] %>%
  scales::percent()
pc2_var_explained_combined <- pca_c_combined@modelDF["p2","R2X"] %>%
  scales::percent()

# look at median age for pc1 > 0 and pc1 <0
pca_scores_c_combined %>%
  transmute(Age, p1, 
            p1_greater_than_zero = p1 > 0) %>%
  group_by(p1_greater_than_zero) %>%
  summarize(p1_median_age = median(Age))

# by gender
pca_gender_combined <- ggplot(pca_scores_c_combined, aes(p1, p2, color = Sex)) + 
  geom_point(size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "Principal Component Analysis",
       subtitle = "Colored by Sex",
       x = str_glue("PC1 ({pc1_var_explained_combined})"),
       y = str_glue("PC2 ({pc2_var_explained_combined})")) + 
  scale_color_viridis_d()

# by age
pca_age_combined <- ggplot(pca_scores_c_combined, aes(Age, p1, color = Age)) + 
  geom_point(size = 5) +
  scale_color_viridis_c() + 
  #stat_ellipse() +
  labs(title = "Principal Component Analysis",
       subtitle = str_glue("First PC by Age"),
       y = str_glue("PC1 ({pc1_var_explained_combined})"),
       x = str_glue("Chronological Age (Years)")) 

## PLS-DA  (using gender as response) -------------------------

# R2X / R2Y: % variance in X/Y explained by component
# Q2 : 1 - (predicted RSS / sum squares Y). tiny q2y means pretty bad predictive performance
plsda_c_combined <- combined_c_processed %>%
  dplyr::select(-any_of(metadata_cols)) %>%
  opls(y = combined_c_processed$Gender, algoC = "nipals", predI = 2)

plsda_c_combined <- readRDS(here("aging_output_files", "plsda_c_combined.Rds"))
# saveRDS(plsda_c_combined, here("aging_output_files", "plsda_c_combined.Rds"))



plsda_pc1_varx_explained_combined <- plsda_c_combined@modelDF["p1","R2X"] %>%
  scales::percent()
plsda_pc2_varx_explained_combined <- plsda_c_combined@modelDF["p2","R2X"] %>%
  scales::percent()


plsda_scores_c_combined <- plsda_c_combined@scoreMN %>%
  as_tibble() %>%
  bind_cols(combined_c_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id))

plsda_gender_combined <- ggplot(plsda_scores_c_combined, aes(p1, p2, color = Sex)) + 
  geom_point(size = 5) +
  stat_ellipse(size = 2) +
  labs(title = "Partial Least Squares Discriminant Analysis",
       subtitle = "By Sex",
       x = str_glue("PC1 ({plsda_pc1_varx_explained_combined})"),
       y = str_glue("PC2 ({plsda_pc2_varx_explained_combined})")) + 
  scale_color_viridis_d()

## compute first pc, do plsda in null space of transofmration matrix ------

x_c_nullpc1_combined <- pca_c_combined_x - pca_c_combined@scoreMN[,'p1'] %*% t(pca_c_combined@loadingMN[,'p1'])
plsda_c_combined_nullpc1 <- x_c_nullpc1_combined %>%
  opls(y = combined_c_processed$Gender, algoC = "nipals", predI = 2)

plsda_c_combined_nullpc1 <- readRDS(here("aging_output_files", "plsda_c_combined_nullpc1.Rds"))
# saveRDS(plsda_c_combined_nullpc1, here("aging_output_files", "plsda_c_combined_nullpc1.Rds"))

plsda_null_pc1_varx_explained_combined <- plsda_c_combined_nullpc1@modelDF["p1","R2X"] %>%
  scales::percent()
plsda_null_pc2_varx_explained_combined <- plsda_c_combined_nullpc1@modelDF["p2","R2X"] %>%
  scales::percent()


plsda_scores_c_combined_nullpc1 <- plsda_c_combined_nullpc1@scoreMN %>%
  as_tibble() %>%
  bind_cols(combined_c_processed %>% select(Age, "Sex" = Gender, APOE, Type, Id))

pca_plsda_scores_c_combined <- pca_scores_c_combined %>%
  select("p1_pca" = p1, Id) %>%
  inner_join(plsda_scores_c_combined_nullpc1, by = "Id") %>%
  select(-p2) %>%
  rename("p1_plsda" = p1)


pca_plsda_combined <- ggplot(pca_plsda_scores_c_combined, aes(p1_pca, p1_plsda, shape = Sex)) + 
  geom_point(aes(color = Age), size = 1) +
  stat_ellipse(size = 0.5) +
  labs(title = "PCA/PLS-DA",
       x = str_glue("PCA: PC1 ({pc1_var_explained_combined})"),
       y = str_glue("PLSDA: PC2 ({plsda_null_pc1_varx_explained_combined})")) + 
  scale_color_viridis_c() #+
#scale_fill_viridis_c()

ggsave("pca_plsda_combined.png",
       plot = pca_plsda_combined,
       path = here("plots", "aging_figs"),
       width = 83,
       height = 60,
       units = 'mm', dpi = 300)


######################################################################

### Clock (linear glmnet age prediction models) ###

######################################################################


targeted_reduced <- reduced_age_analysis(wide_data_targeted)
untargeted_reduced <- reduced_age_analysis(wide_data_untargeted)
got_reduced <- reduced_age_analysis(wide_data)
lipids_reduced <- reduced_age_analysis(wide_data_lipids)

### combined profile
combined_reduced <- reduced_age_analysis(wide_data_combined)

figure_3 <- predtruth_plot(targeted_reduced$pred_df, name = "Targeted") /
  predtruth_plot(untargeted_reduced$pred_df, name = "Untargeted") /
  predtruth_plot(got_reduced$pred_df, name = "GOT") /
  predtruth_plot(lipids_reduced$pred_df, name = "Lipids")

ggsave('figure_3.png',
       plot = figure_3,
       path = here('plots', 'aging_figs', 'main_figs'),
       width = 83, height = 200, dpi = 300, units = 'mm')


####### Untargeted Variations #################

####### rawscaled (to compare performance to drift-corrected)
untargeted_reduced_rawscaled <- reduced_age_analysis(wide_data_untargeted_rawscaled)
saveRDS(untargeted_reduced_rawscaled, file = here('aging_output_files', 'untargeted_reduced_rawscaled.Rds'))

####### models for each sex (to verify decision not to penalize sex)
wide_data_untargeted_male <- wide_data_untargeted %>%
  filter(Gender == 'M')
wide_data_untargeted_female <- wide_data_untargeted %>%
  filter(Gender == 'F')

untargeted_reduced_male <- reduced_age_analysis(wide_data_untargeted_male)
untargeted_reduced_female <- reduced_age_analysis(wide_data_untargeted_female)

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

cowplot::plot_grid(plotlist = purrr::map(wide_data_untargeted_n_male, ~.x$figure))
cowplot::plot_grid(plotlist = purrr::map(wide_data_untargeted_n_female, ~.x$figure))


####### post selection

targeted_reduced_post <- reduced_age_analysis(wide_data_targeted, post_select = T)
untargeted_reduced_post <- reduced_age_analysis(wide_data_untargeted, post_select = T)
got_reduced_post <- reduced_age_analysis(wide_data, post_select = T)
lipids_reduced_post <- reduced_age_analysis(wide_data_lipids, post_select = T)


##### cor > .3 (following Thompson et al's clock)

targeted_reduced_cor <- reduced_age_analysis(wide_data_targeted, remove_cor_age = T)
untargeted_reduced_cor <- reduced_age_analysis(wide_data_untargeted, remove_cor_age = T)
got_reduced_cor <- reduced_age_analysis(wide_data, remove_cor_age = T)
lipids_reduced_cor <- reduced_age_analysis(wide_data_lipids, remove_cor_age = T)



############################

# full models
# (targeted/lipids because they are directly interpretable)

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


################

## Full model fit on untargeted
## -> prediction on AD/PD

################

untargeted_ad_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('AD'), empri = 5)
untargeted_pd_amelia5 <- filter_and_impute_multi(wide_data_untargeted, c('PD'), empri = 5)

untargeted_adpd_separate_amelia5 <- purrr::map2(untargeted_ad_amelia5, untargeted_pd_amelia5, ~merge_datasets(.x, .y, include_metadata = TRUE, include_age = TRUE))

untargeted_adpd_age_analysis_separate <- 1:5 %>% purrr::map(~full_model_new_data(imputed_c_untargeted5[[.x]], new_data = untargeted_adpd_separate_amelia5[[.x]], nlambda = 200))


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



## Comparing metrics for matched controls/ ADPD


# match controls for AD/PD, join with the ADPD data, and remove duplicates
untargeted_adpd_matched_controls <- purrr::map(1:nrow(wide_data_untargeted), ~find_control(.x, data = filter(wide_data_untargeted, Type %in% c("AD", "PD")), 
                                                                                           data_control = filter(wide_data_untargeted, Type %in% c("CO", "CY", "CM")))) %>%
  bind_rows(filter(wide_data_untargeted, Type %in% c("AD", "PD"))) %>%
  distinct(.keep_all = T)

untargeted_c_matched_age_reduced <- untargeted_reduced$pred_df %>%
  filter(id %in% untargeted_adpd_matched_controls$Id)

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
                size =rel(2)) + 
  ggforce::geom_mark_circle(aes(truth, imp_avg, filter = biggest_error ==1), expand = unit(3, 'mm'), size = 0.25)


figure_4a <- predtruth_plot(untargeted_c_matched_age_reduced, name = "Untargeted (matched)")

figure_4 <- figure_4a / 
  untargeted_adpd_age_predtruth_separate

####################

## Mummichog (pathway analysis on untargeted using univariate lm p-values)

###################

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
mummichog_pos_plot <- mummichog_plot(here("aging_output_files", "1591847290.3476467.06-10-2020_posResults", "tables", "mcg_pathwayanalysis_06-10-2020_posResults.tsv")) +
  labs(title = "Positive Mode")


##########

## MSEA (set enrichment analysis for targeted)

#######

# univariate regression to get "significant" metabolite list
c_targeted_univariate_table <- bh_univariate_age(imputed_c_targeted5)

## The reference list must be IDs, so we need to map the names to KEGG/HMDB names.
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

# the list of significant ones
univar_targeted_names_sig <- msea_table %>%
  filter(bh_p_value < 0.05) %>%
  select(mapped_name)


# get sig metabolite names from the univariate analysis
# plug this into MetaboAnalystR
targeted_sig_names <- univar_targeted_names %>%
  str_replace_all("_neg", "") %>%
  str_replace_all("_pos", "") %>%
  # an alternate name for 12 ketolithocholic acid according to https://pubchem.ncbi.nlm.nih.gov/compound/12-Ketolithocholic-acid
  str_replace("3\\?-Hydroxy-12 Ketolithocholic Acid", "12-Ketodeoxycholic acid")

### once we have the msea result table from MetaboAnalystR, recreate the enrichment bar plot in ggplot
msea_smpdb <- read_csv(here("aging_tables", "msea_ora_SMPDB_result.csv")) %>%
  mutate(bh_p = p.adjust(`Raw p`, method = "BH"),
         set = fct_reorder(X1, hits/total)) %>%
  filter(total > 2) %>%
  top_n(10, hits/total) #%>%
#select("Set" = X1, "Total" = total, "Expected" = expected, "Hits" = hits)

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


### once we have the msea result table, recreate the enrichment bar plot in ggplot
msea_csf <- read_csv(here("aging_tables", "msea_ora_csf_result.csv")) %>%
  mutate(bh_p = p.adjust(`Raw p`, method = "BH"),
         set = fct_reorder(X1, hits/total)) %>%
  filter(total >= 2) %>%
  top_n(10, hits/total) #%>%
# select("Set" = X1, "Total" = total, "Expected" = expected, "Hits" = hits)

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


############################

### univariate RF on metabolite concentration ~ Age ###
### To see if there's a natural plateu in metabolite concentration as a function of age

############################


# get full untargeted dataset
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


full_conc_rf_top10 <- concentration_rf_plot(untargeted_univariate_full_conc_age_scale_rf, n = 50) 
c_conc_rf_top50 <- concentration_rf_plot(untargeted_univariate_c_conc_age_scale_rf, n = 50)



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
  scale_y_continuous(labels = scales::percent_format())

figure_1 <- age_dist_c / missingness_overview_plot

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


#### Fitting a missingness model with a 1/0 matrix ---------------------------------------

# get just the age of controls. 
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




#### Lipids missingness matrix ---------------------------------------

# get just the age of controls. 
lipids_c_age <- wide_data_lipids %>% 
  filter(Type %in% c("CO", "CM", "CY")) %>%
  pull(Age)

# get the features
lipids_c_na <- wide_data_lipids %>%
  filter(Type %in% c("CO", "CM", "CY")) %>%
  select(-any_of(metadata_cols)) %>% 
  mutate_all(~ifelse(is.na(.x), 1, 0)) %>%
  as.matrix()

#lipids_c_na_full_model <- get_full_model(features= lipids_c_na, lipids_c_age, alpha = 0.5, family = "gaussian", penalize_AD_PD = FALSE, penalize_age_gender = FALSE, nlambda = 200)
# Note: If AD_ind, PD_ind are missing from the dataset, the flag penalize_AD_Pd doesn't do anything
lipids_c_na_fitpred <- lapply(1:nrow(lipids_c_na), function(x) loo_cvfit_glmnet(x, lipids_c_na, lipids_c_age,
                                                                                alpha = 0.5, family = 'gaussian', penalize_age_gender = FALSE, penalize_AD_PD = FALSE, nlambda = 200))


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


