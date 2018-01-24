library(tidyverse)
library(broom)
library(magrittr)
library(made4)
library(ggjoy)
library(gbm3)

positive_mode_data <- read_csv("~/course/ND_Metabolomics/data/csf_positive.csv")
negative_mode_data <- read_csv("~/course/ND_Metabolomics/data/csf_negative.csv")
positive_mode_data$RunIndex <- negative_mode_data$RunIndex <- 1:nrow(positive_mode_data)

################################
## Detrend data
################################

## Get positive QC data
positive_mode_data$Mode <- ifelse(grepl("QC", positive_mode_data$Code), "QC_pos", "pos")
negative_mode_data$Mode <- ifelse(grepl("QC", negative_mode_data$Code), "QC_neg", "neg")

pos_long <- positive_mode_data %>% gather(key = Metabolite, value = Abundance, -one_of("Code", "RunIndex", "Mode"))
neg_long <- negative_mode_data %>% gather(key = Metabolite, value = Abundance, -one_of("Code", "RunIndex", "Mode"))
pos_long <- pos_long %>% mutate_at("Abundance", funs(log(.)))
neg_long <- neg_long %>% mutate_at("Abundance", funs(log(.)))


## Combine positive and negative mode
csf_long <- rbind(pos_long, neg_long)

## Split QC and sample data
QC_indices <- grepl("QC", csf_long$Code)
QC_long <- csf_long %>% filter(QC_indices)
nd_data <- csf_long %>% filter(!QC_indices)
nd_data$Index <- nd_data$Code %>% gsub(pattern = "(NEG-CSF-06202017-)|(CSF-06122017-)", replacement = "") %>% as.numeric()

## 0 abundace peps are NA for now
nd_data$Abundance[nd_data$Abundance == -Inf] <- NA

fit_boosted_model <- function(df, ntrees=10000, cv.folds=10, min_node_size=10) {
  print(df$Metabolite[1])

  not_na_indices <- which(!is.na(df$Abundance))
  df.omitted <- df %>% na.omit(Abundance)

  boost_fit_cv <- gbm(
    Abundance ~ RunIndex, data = df.omitted, distribution = "laplace",
    n.trees = ntrees, cv.folds = 10, n.minobsinnode = min_node_size, verbose = FALSE
  )
  opt_iters <- gbm.perf(boost_fit_cv, method = "cv")
  print(opt_iters)

  best_fit <- predict(boost_fit_cv, df.omitted, opt_iters)

  abund <- df$Abundance
  bf <- rep(NA, length(abund))

  abund[not_na_indices] <- df.omitted$Abundance - best_fit
  bf[not_na_indices] <- best_fit

  data.frame(Raw = df$Abundance, Abundance = abund, Trend = bf, RunIndex = df$RunIndex, Code = df$Code, Index = df$Index, Mode = df$Mode)
}

detrended <- nd_data %>% group_by(Metabolite, Mode) %>% do(fit_boosted_model(.))

detrended_QC <- QC_long %>% mutate(Index = RunIndex) %>% group_by(Metabolite, Mode) %>% do(fit_boosted_model(., cv.folds = 22, min_node_size = 2))


## Join with Tracking data
tracking_data <- read_csv("~/course/ND_Metabolomics/data/NDTracking.csv")
subject_data <- inner_join(tracking_data, detrended, by = "Index")

## Type2 collapse (CO, CM, CY) to C
subject_data %<>% mutate(Type2 = ifelse(Type %in% c("CY", "CM", "CO"), "C", as.character(.$Type))) %>% mutate(Type2 = factor(Type2))

## APOE NA's to "Unknown"
subject_data <- subject_data %>% mutate(APOE=ifelse(is.na(APOE), "Unknown", as.character(APOE))) %>% mutate(APOE=factor(APOE))

## Fit linear models and save residuals

subject_data$Residual <- NA
for(met in unique(subject_data$Metabolite)) {
    met_indices <- (subject_data$Metabolite == met & !is.na(subject_data$Raw))
    lmod <- subject_data %>% filter(met_indices) %>% lm(Abundance ~ Age + Gender + Type2 + APOE, data=.)
    
    subject_data$Residual[met_indices] <- lmod$residuals

}

lmod_info <- tibble(pred=rep(c("Age", "GenderM", "AD", "Other", "PD", "APOE24", "APOE33", "APOE34", "APOE44"), 2),
                    stat=rep(c("Est", "Pval"), each=9))
for(met in unique(subject_data$Metabolite)) {
    met_indices <- (subject_data$Metabolite == met & !is.na(subject_data$Raw))
    lmod <- subject_data %>% filter(met_indices) %>% lm(Abundance ~ Age + Gender + relevel(Type2, 2) + APOE, data=.)

    lmod_info[[met]] <- as.vector(summary(lmod)$coef[2:10, c(1, 4)])
}

linear_model_info <- lmod_info %>% gather(key="Metabolite", value="Value", -one_of("pred", "stat"))



## biocLite("KEGGREST")
library(KEGGREST)
listDatabases()

metabolites <- unique(subject_data$Metabolite)

matches <- sapply(metabolites, function(met) {
    kegg_matches <- keggFind("compound", met)
    pat1 <- paste0("^", met)
    pat2 <- paste0("; ", met, ";")
    match_idx1 <- grep(pat1, kegg_matches)
    match_idx2 <- grep(pat2, kegg_matches)    
    if(length(match_idx1) > 1) {
        return(kegg_matches[match_idx1][1])
    } else if(length(match_idx2) > 1) {
        return(kegg_matches[match_idx2][1])
    } else {
        return(kegg_matches[1])
    }

    
})
matches_vec <- unlist(matches)
metabolite_kegg_map <- tibble(metabolite = metabolites, kegg_id = sapply(strsplit(names(matches_vec), ":", fixed=TRUE), function(x) x[2]), description=matches_vec)

metabolite_kegg_map[100, ]

pathways <- sapply(metabolite_kegg_map$kegg_id, function(cid) {
    print(cid)
    if(is.na(cid)) {
        NULL
    } else {
        keggGet(cid)[[1]]$PATHWAY
    }
})

metabolite_kegg_map$pathways <- pathways

save(subject_data, QC_long, linear_model_info, metabolite_kegg_map, file = "preprocessed_csf_data.RData")
