source(here::here("analysis", "starter.R"))



## descriptive summary table
panuc_data <- read_csv('E:/Projects/metabolomics/ND_Metabolomics/data/PANUC/panuc.csv')
ad_c_tracking <- read_csv('E:/Projects/metabolomics/ND_Metabolomics/data/ADTracking.csv') %>%
  select(LPAge, Sex = Gender, Type = PrimaryDx, OnsetAge) %>%
  mutate(Type = fct_collapse(Type,'Control' = c('CO', 'CY', 'CM')))
ad_c_tracking_wrace <- readxl::read_xlsx('E:/Projects/metabolomics/ND_Metabolomics/data/ADTracking 2_20210621.xlsx') %>%
  transmute(type2 = fct_collapse(PrimaryDx, "C" = c("CO", "CM", "CY")),
            Gender, Race,
            is_white = Race == "White")
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





# figure 1

## Missing data overview plot

data_sources <- list(wide_data, wide_data_lipids, wide_data_targeted, wide_data_untargeted)
data_names <- list("GOT", "Lipids", "Targeted", "Untargeted")
pct_missing_data_adpd <- purrr::map2(data_sources, data_names, ~percent_missing_by_type(.x, .y)) %>%
  bind_rows %>%
  gather(key = "Type", value = "pct_missing", -source) %>%
  collapse_controls()

age_dist_ad_pd <- ggplot(wide_data_targeted %>% collapse_controls(), aes(Age)) +
  geom_histogram(bins = 10) +
  facet_wrap(~Type) + 
  labs(x = "Age",
       y = "Count",
       title = "Distribution of Age")


missingness_overview_adpd_plot <- ggplot(pct_missing_data_adpd) + 
  geom_col(aes(x = source, y= pct_missing, fill = Type), position = 'dodge') +
  labs(title = "Proportion of Missingness By Dataset",
       x = "Dataset",
       y = "Proportion Missing")
ggsave("missingness_overview_adpd_bar.png",
       plot = missingness_overview_adpd_plot,
       path = here("plots", "adpd_figs"),
       width = 14,
       height = 10)

figure_1 <- age_dist_ad_pd + missingness_overview_adpd_plot +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
ggsave("figure_1.jpg",
       type = 'cairo',
       device = 'jpg',
       plot = figure_1,
       path = here("plots", "adpd_figs", "main_figs"),
       width = 83 * 3,
       height = 100,
       units = 'mm', dpi = 300)


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

pca_type_untargeted <- ggplot(pca_scores_all_untargeted, aes(p1, p2, color = Type)) + 
  geom_point(size = 1) +
  stat_ellipse(size = 0.5) +
  labs(title = "Principal Component Analysis",
       subtitle = "Colored by Phenotype",
       x = str_glue("PC1 ({pc1_all_var_explained})"),
       y = str_glue("PC2 ({pc2_all_var_explained})"))

ggsave("figure_2.jpg",
       type = 'cairo',
       device = 'jpg',
       plot = pca_type_untargeted,
       path = here("plots", "adpd_figs", "main_figs"),
       width = 125,
       height = 100,
       units = 'mm', dpi = 300)


###########

# fitting iterative classification models. 
# warning: these take a very long time to run on a local machine.
  # I recommend validating the methods and checking code using a tiny fraction of the data.

###########

### Models for AD vs c------------------
targeted_ad_reduced_detrend <- reduced_logistic_analysis(data = wide_data_targeted_detrend, 
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


### Models for PD vs c ------------------

targeted_pd_reduced_detrend <- reduced_logistic_analysis(data = wide_data_targeted_detrend, 
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



### Models for AD vs PD ---------------------

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

figure_3 <- untargeted_ad_c_roc_detrend + untargeted_pd_c_roc_detrend + untargeted_adpd_c_roc_detrend +
  plot_annotation(title = "ROC curves using the untargeted profile", theme = theme(title = element_text(size = rel(2))),
                  tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

ggsave("figure_3.jpg",
       type = 'cairo',
       device = 'jpg',
       plot = figure_3,
       path = here("plots", "adpd_figs", "main_figs"),
       width = 83*5,
       height = 150,
       units = 'mm', dpi = 300)



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

untargeted_univar_ad_logistic_list <- furrr::future_map(1:5, ~bh_univariate_age(untargeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = ad_c_index, get_coefs = T))
# saveRDS(untargeted_univar_ad_logistic_list, here("ad_pd", "untargeted_univar_ad_logistic_list.Rds"))


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
  dplyr::select('mz' = `m/z`, 'rtime' = `Retention time (min)`, 'p-value' = bh_p_value, 'z-score' = `z value`, Metabolite, Mode, Estimate)

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

pd_c_index <- which(imputed_all_untargeted5[[1]][[2]] %in% c('PD', 'CO', 'CY', 'CM'))

untargeted_univar_pd_logistic_list <- 
  furrr::future_map(1:5, ~bh_univariate_age(untargeted_all_amelia5_pd_ind, 
                                            types_index = pd_c_index, 
                                            var = "PD_ind", family = "binomial", conc = FALSE, imp_num = .x, get_coefs = T))

# saveRDS(untargeted_univar_pd_logistic_list, here("ad_pd", "untargeted_univar_pd_logistic_list.Rds"))

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


######## PD vs AD ----------------------------------------------------------


pd_ad_index <- which(imputed_all_untargeted5[[1]][[2]] %in% c('PD', 'AD'))

untargeted_univar_adpd_logistic_list <- furrr::future_map(1:5, ~bh_univariate_age(untargeted_all_amelia5_pd_ind, types_index = pd_ad_index, var = "PD_ind", family = "binomial", conc = FALSE, imp_num = .x, get_coefs = T))
# saveRDS(untargeted_univar_adpd_logistic_list, here("ad_pd", "untargeted_univar_adpd_logistic_list.Rds"))


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


# reporting large univar coefs
untargeted_univar_adpd_median_coef <- untargeted_univar_adpd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(Estimate == median(Estimate)) %>%
  # break ties
  slice(1) %>%
  ungroup() %>%
  arrange(desc(abs(Estimate))) %>%
  mutate(OR_est = exp(Estimate))




mummichog_pd_pos_plot <- mummichog_plot(here("ad_pd", "mummichog_pd_pos", "tables", "mcg_pathwayanalysis_mummichog_pd_pos.tsv")) +
  labs(title = "PD vs Controls") +#, subtitle = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title.x = element_text(size = rel(1.5)),
        plot.title = element_text(size = rel(1.75)))

mummichog_ad_pos_plot <- mummichog_plot(here("ad_pd", "mummichog_ad_pos", "tables", "mcg_pathwayanalysis_mummichog_ad_pos.tsv")) +
  labs(title = "AD vs Controls")+#,
  #subtitle = "Positive Mode") + 
  theme(axis.text.y = element_text(size = rel(2)),
        axis.title.x = element_text(size = rel(1.5)),
        plot.title = element_text(size = rel(1.75)))

figure_4 <- mummichog_pd_pos_plot + mummichog_ad_pos_plot +
  plot_annotation(title = "Mummichog Identified Pathways (Positive Mode)",
                  theme = theme(plot.title = element_text(size = rel(2))))


ggsave("figure_4.jpg",
       type = 'cairo',
       device = 'jpg',
       plot = figure_4,
       path = here("plots", "adpd_figs", "main_figs"),
       width = 83*6,
       height = 250,
       units = 'mm', dpi = 300)



############################

### Univariate AD/PD Logisitic Regression on targeted ###
### for use with MSEA

############################

#### AD vs controls -----------------------------
targeted_all_amelia5_ad_ind <- imputed_all_targeted5 %>%
  purrr::map(~cbind(.x[[1]], "AD_ind" = ifelse(imputed_all_targeted5[[1]][[2]] == "AD", 1, 0)) %>% list())

ad_c_index_targeted <- which(imputed_all_targeted5[[1]][[2]] %in% c('AD', "CO", "CM", "CY"))

# Create table with bh-corrected p values
targeted_univar_ad_logistic_list <- purrr::map(1:5, ~bh_univariate_age(targeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = ad_c_index_targeted, get_coefs = T))
targeted_univar_ad_logistic <- targeted_univar_ad_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()


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


# reporting large univar coefs
targeted_univar_ad_median_coef <- targeted_univar_ad_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(Estimate == median(Estimate)) %>%
  # break ties
  slice(1) %>%
  ungroup() %>%
  arrange(desc(abs(Estimate))) %>%
  mutate(OR_est = exp(Estimate))



#### PD -------------------------------------
# # we want a PD indicator, but not an AD one
targeted_all_amelia5_pd_ind <- imputed_all_targeted5 %>%
  purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(imputed_all_targeted5[[1]][[2]] == "PD", 1, 0)) %>% list())

pd_c_index_targeted <- which(imputed_all_targeted5[[1]][[2]] %in% c('PD', "CO", "CM", "CY"))
# Create table with bh-corrected p values
targeted_univar_pd_logistic_list <- purrr::map(1:5, ~bh_univariate_age(targeted_all_amelia5_pd_ind, var = "PD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = pd_c_index_targeted, get_coefs = T))

# for each metabolite, use the smallest p value out of the 5 imputations. this is very generous
targeted_univar_pd_logistic <- targeted_univar_pd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  slice(1) %>%
  ungroup()

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



# reporting large univar coefs
targeted_univar_pd_median_coef <- targeted_univar_pd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(Estimate == median(Estimate)) %>%
  # break ties
  slice(1) %>%
  ungroup() %>%
  arrange(desc(abs(Estimate))) %>%
  mutate(OR_est = exp(Estimate))


#### AD vs PD -----------------------------
# classify AD against PD (ad is positive class)

ad_pd_index_targeted <- which(imputed_all_targeted5[[1]][[2]] %in% c('AD', "PD"))

# Create table with bh-corrected p values
targeted_univar_adpd_logistic_list <- purrr::map(1:5, ~bh_univariate_age(targeted_all_amelia5_ad_ind, var = "AD_ind", family = "binomial", conc = FALSE, imp_num = .x, types_index = ad_pd_index_targeted, get_coefs = T))
targeted_univar_adpd_logistic <- targeted_univar_adpd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()

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




# reporting large univar coefs
targeted_univar_adpd_median_coef <- targeted_univar_adpd_logistic_list %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(Estimate == median(Estimate)) %>%
  # break ties
  slice(1) %>%
  ungroup() %>%
  arrange(desc(abs(Estimate))) %>%
  mutate(OR_est = exp(Estimate))



###################

# Missingness

###################

untar_missingness_adpd_plot <- ggplot(missingness_by_type_untargeted_sig_pct, aes(fct_reorder(name, pct_missing), pct_missing)) +
  geom_col(aes(fill = fct_relevel(Type, "AD", after = 1)), position = "dodge", width = 0.75) + 
  #scale_color_manual(labels = c("CY", "CM", "CO", "AD", "PD"), values = c("lightskyblue", "dodgerblue", "blue", "darkgreen", "purple")) +
  theme(axis.text.x = element_text(angle= 90, hjust =1)) + 
  #theme(axis.text.x = element_blank()) + 
  labs(title = 'Percent Missingness by Type',
       subtitle = 'Untargeted, FDR < 0.05',
       y = "Percent Missing",
       x = "Metabolite",
       fill = "Type")  


ggsave("figure_a1.jpg",
       type = 'cairo',
       device = 'jpg',
       plot = untar_missingness_adpd_plot,
       path = here("plots", "adpd_figs", "appendix_figs"),
       width = 83*2,
       height = 100,
       units = 'mm', dpi = 300)



#### fitting models using only missingness indicators (and age and sex).
#    (again, warning that these models take a very long time to run)


targeted_ad_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
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

# PD vs controls


targeted_pd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
                                                                  target ="PD_ind", PD_ind = T, 
                                                                  types = c("PD", "CO", "CM", "CY"),
                                                                  penalize_age_gender = T
)


untargeted_pd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_untargeted_naind, 
                                                                    target ="PD_ind", PD_ind = T, 
                                                                    types = c("PD", "CO", "CM", "CY"),
                                                                    penalize_age_gender = T
)


lipids_pd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_lipids_naind, 
                                                                target ="PD_ind", PD_ind = T, 
                                                                types = c("PD", "CO", "CM", "CY"),
                                                                penalize_age_gender = T
)

#AD vs PD

targeted_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'targeted_adpd_reduced_naind_nopenalize.Rds'))
untargeted_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'untargeted_adpd_reduced_naind_nopenalize.Rds'))
lipids_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'lipids_adpd_reduced_naind_nopenalize.Rds'))
got_adpd_reduced_naind_nopenalize <- readRDS(file = here('aging_output_files', 'got_adpd_reduced_naind_nopenalize.Rds'))


targeted_adpd_reduced_naind_nopenalize <- reduced_logistic_analysis(data = wide_data_targeted_naind, 
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







###########################

## Appendix C:  (figure A4) 


############################

### fitting models (again, warning that this takes a very long time) ----


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


######## looking at univariate p values after orthogonalization

t_adpd_mod <- logistic_control_analysis(imputed_all_targeted5, varname ="AD_ind", imp_num = 1, nlambda = 200, AD_ind = T, types = c("PD", "AD"), full_model = T, include_age_gender = FALSE,
                                        alpha = 0.5)

t_adpd_coef <- coef(t_adpd_mod, s = "lambda.min") %>% 
  as.matrix() %>%
  as_tibble(rownames= "id") %>%
  set_names(c("id", "beta")) %>%
  filter(beta != 0)


t_impute5_orthog <- furrr::future_map(imputed_all_targeted5, function(lst){
  lst[[1]]  <- adpd_orthog(lst[[1]], coef = t_adpd_coef)
  lst
})


t_all_amelia5_pd_ind_orthog <- t_impute5_orthog %>%
  purrr::map(~cbind(.x[[1]], "PD_ind" = ifelse(t_impute5_orthog[[1]][[2]] == "PD", 1, 0)) %>% list())

pd_c_index <- which(t_impute5_orthog[[1]][[2]] %in% c('PD', 'CO', 'CY', 'CM'))

t_univar_pd_list_orthog <- 
  furrr::future_map(seq_along(t_all_amelia5_pd_ind_orthog), 
                    ~bh_univariate_age(t_all_amelia5_pd_ind_orthog, 
                                       types_index = pd_c_index, 
                                       var = "PD_ind", family = "binomial", conc = FALSE, imp_num = .x))

t_univar_pd_orthog <- t_univar_pd_list_orthog %>%
  bind_rows(.id = "imp") %>%
  group_by(name) %>%
  filter(bh_p_value == median(bh_p_value)) %>%
  # break ties
  slice(1) %>%
  ungroup()

t_univar_pd_orthog  <- readRDS(here("ad_pd", "t_univar_pd_orthog.Rds"))
# saveRDS(t_univar_pd_orthog, here("ad_pd", "t_univar_pd_orthog.Rds"))


tar_pval_compare <- targeted_univar_pd_logistic %>%
  select(name, imp, "full" = bh_p_value) %>%
  inner_join(t_univar_pd_orthog, by = c("imp", "name")) %>%
  rename("orthog" = bh_p_value) %>%
  mutate(is_sig = case_when(
    (full < 0.05 & orthog > 0.05) | (full > 0.05 & orthog < 0.05) ~ "one",
    full < 0.05 &  orthog < 0.05 ~ "both",
    TRUE ~ "neither") %>%
      fct_relevel("both", "one")
  )
#pivot_longer(c("full", "bh_p_value"), names_to = "var", values_to = "bh_p_val")  %>%
#mutate(var = ifelse(var == "bh_p_value", "orthog", var))



figure_a4 <- tar_pval_compare %>%
  mutate(across(all_of(c("full", "orthog")), ~ -log10(.x))) %>%
  ggplot() +
  geom_point(aes(full, orthog, color = is_sig), size = 5) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    title = "Comparison of univariate p-values before and after orthogonalization",
    subtitle = "Targeted profile, logistic regression classifying PD vs control",
    x = "BH-corrected -log10 p-value before orthogonalization",
    y = "BH-corrected -log10 p-value after orthogonalization",
    color = "p < 0.05?"
  ) #+
#scale_color_viridis_d()

ggsave("figure_a4.jpg",
       type = 'cairo',
       device = 'jpg',
       plot = figure_a4,
       path = here("plots", "adpd_figs", "appendix_figs"),
       width = 83*2,
       height = 100,
       units = 'mm', dpi = 300)


######################


# appendix:  not penalizing age and sex

######################

untargeted_ad_c_roc <- roc_fig(untargeted_ad_reduced$roc_df, untargeted_ad_reduced$predtruth_df) + 
  labs(title = 'Untargeted')
targeted_ad_c_roc <- roc_fig(targeted_ad_reduced$roc_df, targeted_ad_reduced$predtruth_df) + 
  labs(title = 'Targeted')
lipids_ad_c_roc <- roc_fig(lipids_ad_reduced$roc_df, lipids_ad_reduced$predtruth_df) + 
  labs(title = 'Lipids')
got_ad_c_roc <- roc_fig(got_ad_reduced$roc_df, got_ad_reduced$predtruth_df) + 
  labs(title = 'GOT')


roc_ad_reduced_grid <- untargeted_ad_c_roc + targeted_ad_c_roc + lipids_ad_c_roc + 
  plot_annotation(title = "AD vs C",theme = theme(title = element_text(size = rel(1.5)))#,
                  #subtitle = str_glue("(Age, sex only AUC: {meta_ad_c_auc})")) +
  ) +
  plot_layout(ncol = 3, nrow = 1)



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

# get null model from adpd_glmnet.R (base_ad_c...)
meta_pd_c_auc <- round(base_pd_c_roc$auc[1], 3)

roc_pd_reduced_grid <- untargeted_pd_c_roc + targeted_pd_c_roc + lipids_pd_c_roc + 
  plot_annotation(title = "PD vs C", theme = theme(title = element_text(size = rel(1.5)))#,
                  #subtitle = str_glue("(Age, sex only: {meta_pd_c_auc})")) +
  )+
  plot_layout(ncol = 3, nrow = 1)



untargeted_adpd_roc <- roc_fig(untargeted_adpd_reduced$roc_df, untargeted_adpd_reduced$predtruth_df) + 
  labs(title = 'Untargeted')
targeted_adpd_roc <- roc_fig(targeted_adpd_reduced$roc_df, targeted_adpd_reduced$predtruth_df) + 
  labs(title = 'Targeted')
lipids_adpd_roc <- roc_fig(lipids_adpd_reduced$roc_df, lipids_adpd_reduced$predtruth_df) + 
  labs(title = 'Lipids')
roc_adpd_reduced_grid <- untargeted_adpd_roc + targeted_adpd_roc + lipids_adpd_roc + 
  plot_annotation(title = "AD vs PD", theme = theme(title = element_text(size = rel(1.5)))#,
                  #subtitle = str_glue("(Age, sex only: {meta_adpd_auc})")) +
  )+
  plot_layout(ncol = 3, nrow = 1)


roc_ad_reduced_grid

ggsave("figure_a5a.tiff",
       type = 'cairo',
       device = 'tiff',
       plot = roc_ad_reduced_grid,
       path = here("plots", "adpd_figs", "appendix_figs"),
       width = 83*3,
       height = 100,
       units = 'mm', dpi = 300)

ggsave("figure_a5b.tiff",
       type = 'cairo',
       device = 'tiff',
       plot = roc_pd_reduced_grid,
       path = here("plots", "adpd_figs", "appendix_figs"),
       width = 83*3,
       height = 100,
       units = 'mm', dpi = 300)

ggsave("figure_a5c.jpg",
       type = 'cairo',
       device = 'jpeg',
       plot = roc_adpd_reduced_grid,
       path = here("plots", "adpd_figs", "appendix_figs"),
       width = 83*3,
       height = 100,
       units = 'mm', dpi = 300)





#############################

# Appendix. models for LED, and GBA

##############################

### GBA

targeted_GBA_reduced <- reduced_logistic_analysis(data = wide_data_targeted, 
                                                  target ="GBA_ind", GBA_ind = T, 
                                                  types = c('PD')
)
untargeted_GBA_reduced <- reduced_logistic_analysis(data = wide_data_untargeted, 
                                                    target ="GBA_ind", GBA_ind = T, 
                                                    types = c("PD")
)

lipids_GBA_reduced <- reduced_logistic_analysis(data = wide_data_lipids, 
                                                target ="GBA_ind", GBA_ind = T, 
                                                types = c("PD")
)


### LED

targeted_led_reduced <- reduced_age_analysis(wide_data_targeted_panuc, 
                                             target = 'led',
                                             types = c('PD'),
                                             include_led = T)
untargeted_led_reduced <- reduced_age_analysis(wide_data_untargeted_panuc, 
                                               target = 'led',
                                               types = c('PD'),
                                               include_led = T)
lipids_led_reduced <- reduced_age_analysis(wide_data_lipids_panuc, 
                                           target = 'led',
                                           types = c('PD'),
                                           include_led = T)
