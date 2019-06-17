## envelope model for CSF data
## Envelope models do much worse when sighat is a poor estimator of sigma
##
library(tidyverse)
library(ggExtra)
library(magrittr)
library(made4)
library(gbm3)
library(huge)
library(microbenchmark)
library(mvtnorm)
library(mvnfast)
library(rstiefel)
library(Amelia)
library(modelr)
library(robust)
library(ggridges)



source("envelope_functions.R")
source("utility.R")

## load("preprocessed_gotms_data.RData")
load("preprocessed_gotms_data-2018-03-21.RData")

all_ids <- subject_data$Metabolite %>% sub(" Results", "", .) %>% unique

pos_trans <- read_csv("../data/got-ms/standards/transition_data_pos.csv")
neg_trans <- read_csv("../data/got-ms/standards/transition_data_neg.csv")
neg_trans <- neg_trans %>% mutate(Name = as.character(Name))
all_trans <- read_csv("../data/got-ms/standards/all_transitions.csv")

all_features <- left_join(tibble(Name=all_ids), all_trans, by="Name")
all_features <- all_features %>% mutate(`Precursor Ion` = round(`Precursor Ion`),
                                        `Product Ion` = round(`Product Ion`))
all_features <- all_features %>% filter(!is.na(`Precursor Ion`))

standards_pos <- read_csv("../data/got-ms/standards/standards_pos.csv")
standards_neg <- read_csv("../data/got-ms/standards/standards_neg.csv")


find_match <- function(df, mode="pos") {

    if(mode == "pos") {
        standards <- standards_pos
    } else {
        standards <- standards_neg
    }

    if(nrow(df) > 1) {
       df <- df[1, ]
    }
    
    pre <- df$`Precursor Ion`
    prod <- df$`Product Ion`

    precursor_diffs <- abs(standards$`Precursor Ion` - pre)
    indices <- which(precursor_diffs == min(precursor_diffs, na.rm=TRUE))

    standards_sub <- standards[indices, ]
    product_diffs <- abs(standards_sub$`Product Ion` - prod)
    indices <- which(product_diffs == min(product_diffs, na.rm=TRUE))

    best_match <- standards_sub[indices, ]
    best_match$PreOrig <- as.numeric(pre)
    best_match$ProdOrig <- as.numeric(prod)
    best_match
}

pos_matches <- pos_trans %>% group_by(Name) %>% do(find_match(.))
pos_matches %>% dplyr::select(Name, Metabolite, PreOrig, `Precursor Ion`,
                              ProdOrig, `Product Ion`)

neg_matches <- neg_trans %>% group_by(Name) %>% do(find_match(., mode="neg"))
neg_matches %>% dplyr::select(Name, Metabolite, PreOrig, `Precursor Ion`,
                              ProdOrig, `Product Ion`)

pos_matches$Mode <- "pos"
neg_matches$Mode <- "neg"

all_matches <- rbind(pos_matches, neg_matches)

pos_matches %>% dplyr::select(Name) %>% as.character

pos_matches[403,"Name"]

all_matches %>% filter(Name == "254-5.55") %>% 
    dplyr::select(Name, Metabolite, PreOrig,
                  `Precursor Ion`, ProdOrig, `Product Ion`)


all_matches %>% filter(Name == "10")
all_matches %>% filter(Name == "171-2.779")
all_matches %>% filter(Name == "196")

save(all_matches, file="../data/got-ms/identification_map.RData")
