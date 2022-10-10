library(tidyverse)
library(here)

### correlations over age


#### Targeted ----------------------------------------
#pick out some significant metablolites
load(here("targeted_age_analysis_02032020.RData"))

# get all sig metabolites in all imputations
targeted_in_any <- targeted_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff(c("GenderM","(Intercept)"))


#' Function to compute spearman correlation bewteen two metabolites for a given type (CO, CM, ..)
#' 
cor_by_type <- function(data, m1, m2, type){
  data %>%
    filter(Type == type) %>%
    select(m1, m2) %>%
    mutate_all(~scale(.x, center = T)) %>%
    # NA means that there were not complete cases
    summarize(cor = cor(!!sym(m1), !!sym(m2), use = "na.or.complete",method = "spearman")) %>%
    as.numeric()
  
}


# is this ok to use on spearman?
fisher_z_transform <- function(rho){
  0.5 * log((1 + rho)/(1 - rho))
}


pairwise_comparisons_type <- crossing(m1 = targeted_in_any,
                                 m2 = targeted_in_any,
                                 type = c("CY", "CM", "CO", "AD", "PD")) %>%
  rowwise() %>%
  mutate(cor = cor_by_type(wide_data_targeted, m1, m2, type)) %>%
  ungroup() %>%
  mutate(type = fct_relevel(type, "CY", "CM", "CO", "AD", "PD")) %>%
  filter(m1 != m2)
  

ggplot(pairwise_comparisons_type) +
  geom_boxplot(aes(x = type, y = cor))

ggplot(pairwise_comparisons_type) + 
  geom_density_ridges(aes(cor, type, fill = type),
                      quantile_lines = T, quantiles = 2,
                      jittered_points = T,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7) + 
  labs(title = "spearman correlations of significant targeted metabolites")



### sliding age window, fit lm on cor ~ age within each age window for each metabolite
# pull out significant slopes, and plot 

cor_by_age <- function(data, m1, m2, age_min, age_max){
  data %>%
    filter(Age >= age_min & Age < age_max) %>%
    select(m1, m2, Age) %>%
    mutate_at(vars(-Age), ~scale(.x, center = T)) %>%
    summarize(cor = cor(!!sym(m1), !!sym(m2), use = "na.or.complete",method = "spearman")) %>%
    as.numeric()
}

# Age ranges from 20 to 88 in our data. We increment by 10
age_ranges <- seq(20, 90, by = 10)

pairwise_comparisons_age <- crossing(m1 = targeted_in_any,
                                     m2 = targeted_in_any,
                                     age_min = head(age_ranges, -1)) %>%
  mutate(age_max = age_min + 10) %>%
  # above creates all combos of age min and age max. we only want a few
  # we could for sure speed this up...
  rowwise() %>%
  mutate(cor = cor_by_age(wide_data_targeted, m1, m2, age_min, age_max),
         fcor = fisher_z_transform(cor)) %>%
  ungroup() %>%
  # cor is trivially 1 when m1 = m2
  filter(m1 != m2) %>%
  # we also need to get rid of duplicates (ie m1 = m2, and then m2 = m1)
  mutate(key = paste0(pmin(m1, m2), pmax(m1, m2), sep = "")) %>%
  distinct(key, age_min, .keep_all = T)



get_lm_p <- function(meta1, meta2){
  # can use age_min or age_max. doesn't matter because they're 1-1
  data <- pairwise_comparisons_age %>%
    filter(m1 == meta1, m2 == meta2)
  fit <- lm(fcor ~ age_min, data = data)
  coef(summary(fit))["age_min", "Pr(>|t|)"]
}

cor_by_age_lm <- pairwise_comparisons_age %>%
  group_by(key) %>%
  summarize(p_val = get_lm_p(m1, m2)) %>%
  ungroup() %>%
  mutate(p_val = p.adjust(p_val, method = "BH")) %>%
  arrange(p_val)


most_int_cor <- cor_by_age_lm %>%
  head(10) %>%
  pull(key)

pairwise_comparisons_age %>%
  filter(key %in% most_int_cor) %>%
  ggplot() +
  geom_line(aes(age_min, fcor, color = key)) + 
  labs(title = "Fisher Correlation in 10 Age Bins",
       subtitle = "10 most significant metabolites")







## fisher transform -> z test on two groups ---------------
median_age <- median(wide_data_targeted$Age)

cor_2groups <- function(data, m1, m2, group){
  if(group == "Young"){
    data <- data %>%
      filter(Age <= median_age)
  } else{
    data <- data %>%
      filter(Age > median_age)
  }
  data %>%
    select(m1, m2, Age) %>%
    mutate_at(vars(-Age), ~scale(.x, center = T)) %>%
    summarize(cor = cor(!!sym(m1), !!sym(m2), use = "na.or.complete",method = "spearman")) %>%
    as.numeric()
}


pairwise_comparisons_2groups <- crossing(m1 = targeted_in_any,
                                     m2 = targeted_in_any,
                                     age_group = c("Young","Old")) %>%
  # above creates all combos of age min and age max. we only want a few
  # we could for sure speed this up...
  rowwise() %>%
  mutate(cor = cor_2groups(wide_data_targeted, m1, m2, age_group),
         fcor = fisher_z_transform(cor)) %>%
  ungroup() %>%
  # cor is trivially 1 when m1 = m2
  filter(m1 != m2) %>%
  # we also need to get rid of duplicates (ie m1 = m2, and then m2 = m1)
  mutate(key = paste0(pmin(m1, m2), pmax(m1, m2), sep = "")) %>%
  distinct(key, age_group, .keep_all = T) %>%
  arrange(m1, m2)


youngcor <- pairwise_comparisons_2groups %>%
  filter(age_group == "Young") %>% 
  pull(fcor)
oldcor <- pairwise_comparisons_2groups %>%
  filter(age_group == "Old") %>% 
  pull(fcor)

t.test(youngcor, oldcor)



#### Untargeted --------------------------------
#pick out some significant metablolites
load(here("untargeted_age_analysis_02032020.RData"))

# get all sig metabolites in all imputations
untargeted_in_any <- untargeted_age_analysis %>%
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(c) %>%
  setdiff(c("GenderM","(Intercept)"))



untar_pairwise_comparisons <- crossing(m1 = untargeted_in_all,
                                 m2 = untargeted_in_all,
                                 type = c("CY", "CM", "CO")) %>%
  filter(m1 != m2) %>%
  rowwise() %>%
  mutate(cor = cor_by_type(wide_data_untargeted, m1, m2, type)) %>%
  ungroup() %>%
  mutate(type = fct_relevel(type, "CY", "CM", "CO"))


ggplot(untar_pairwise_comparisons) +
  geom_boxplot(aes(x = type, y = cor))

ggplot(untar_pairwise_comparisons) + 
  geom_density_ridges(aes(cor, type, fill = type),
                      quantile_lines = T, quantiles = 2,
                      jittered_points = T,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7) + 
  labs(title = "spearman correlations of significant untargeted metabolites")


crossing(m1 = targeted_in_all,
         m2 = targeted_in_all,
         type = c("CY", "CM", "CO")) %>%
  filter(m1 != m2) %>%
  rowwise() %>%
  mutate(cor = cor_by_type(wide_data_targeted, m1, m2, type)) %>%
  ungroup() %>%
  mutate(type = fct_relevel(type, "CY", "CM", "CO"))



# sub bullets for intro
# motivation

