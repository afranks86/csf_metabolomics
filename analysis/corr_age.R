library(tidyverse)
library(here)

### correlations over age


#### Targeted ----------------------------------------
#pick out some significant metablolites
load(here("targeted_age_analysis_02032020.RData"))

# get all sig metabolites in all imputations
targeted_in_all <- targeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff(c("GenderM", "(Intercept)"))


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

pairwise_comparisons <- crossing(m1 = targeted_in_all,
                                 m2 = targeted_in_all,
                                 type = c("CY", "CM", "CO")) %>%
  rowwise() %>%
  mutate(cor = cor_by_type(wide_data_targeted, m1, m2, type)) %>%
  ungroup() %>%
  mutate(type = fct_relevel(type, "CY", "CM", "CO")) %>%
  filter(m1 != m2)
  

ggplot(pairwise_comparisons) +
  geom_boxplot(aes(x = type, y = cor))

ggplot(pairwise_comparisons) + 
  geom_density_ridges(aes(cor, type, fill = type)) + 
  labs(title = "spearman correlations of significant targeted metabolites")




#### Untargeted --------------------------------
#pick out some significant metablolites
load(here("untargeted_age_analysis_02032020.RData"))

# get all sig metabolites in all imputations
untargeted_in_all <- untargeted_age_analysis %>% 
  purrr::map(~.x[[2]] %>% names) %>% 
  reduce(intersect) %>%
  setdiff(c("GenderM", "(Intercept)"))



untar_pairwise_comparisons <- crossing(m1 = untargeted_in_all,
                                 m2 = untargeted_in_all,
                                 type = c("CY", "CM", "CO")) %>%
  rowwise() %>%
  mutate(cor = cor_by_type(wide_data_untargeted, m1, m2, type)) %>%
  ungroup() %>%
  mutate(type = fct_relevel(type, "CY", "CM", "CO")) %>%
  filter(m1 != m2)


ggplot(untar_pairwise_comparisons) +
  geom_boxplot(aes(x = type, y = cor))

ggplot(untar_pairwise_comparisons) + 
  geom_density_ridges(aes(cor, type, fill = type)) + 
  labs(title = "spearman correlations of significant untargeted metabolites")


