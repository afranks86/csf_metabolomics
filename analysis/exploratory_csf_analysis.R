library(tidyverse)
library(ggExtra)
library(magrittr)
library(made4)
library(ggridges)
library(gbm3)
library(huge)

source("utility.R")

load("../data/preprocessed_csf_data.RData")


processed_files <- dir(pattern="^preprocessed_csf_data*")

## Most recent file
load(max(processed_files[grep("-20+", processed_files)]))

## Notes:
## Raw: raw intensities
## Abundance: detrended intensiteis (raw intensity - boosted fit on RunIndex)
## Trend: boosted fit on RunIndex

subject_data %<>% mutate_at(c("Type", "Gender", "APOE"), factor)

## Look at average raw data
subject_data %>% group_by(RunIndex, Mode) %>%
    summarize(tot=sum(Raw)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() + ggtitle("Positive") +
    geom_smooth(method="lm", se=FALSE)

## Look at quality control
QC_long %>% group_by(RunIndex, Mode) %>%
    summarize(tot=sum(Abundance)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() + ggtitle("Positive") +
    geom_smooth(method="lm", se=FALSE)

## Focus on raw neg
subject_data %>% filter(Mode=="neg") %>%
    group_by(RunIndex, Mode) %>%
    summarize(tot=sum(Raw)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)

## Focus on raw pos
subject_data %>% filter(Mode=="pos") %>%
    group_by(RunIndex, Mode) %>%
    summarize(tot=sum(Raw)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() +
    ggtitle("Negative") + geom_smooth(method="lm", se=FALSE)


## check that detrended data looks sensible (for a random sample of metabolites)
subject_data %>%
    filter(Metabolite %in% sample(unique(subject_data$Metabolite), 4)) %>%
    ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) +
    geom_line() +
    geom_point() +
    theme(legend.position="none")

## Which metabolites had the worst drift?
subject_data %>% group_by(Metabolite) %>%
    summarise(var = (max(Trend, na.rm=TRUE) -
                     min(Trend, na.rm=TRUE)) / sd(Raw, na.rm=TRUE)) %>%
    arrange(desc(var))

## Look at some trends fits for a few metabolites with large trends
intercepts <- (subject_data %>% group_by(Batch) %>%
               summarise(mx = max(RunIndex)))$mx

subject_data %>% filter(Metabolite %in%
                        c("Creatinine", "Fructose", "Xanthine", "Indole")) %>%
    ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) +
    geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=RawScaled)) +
    geom_vline(xintercept = intercepts)

## Plot raw intensity versus log variance of the abundance (detrended)
subject_data  %>% group_by(Metabolite) %>%
    summarise(med=median(Raw),
              var=log(var(Abundance))) %>%
    ggplot(aes(x=med, y=var)) +
    geom_point() +
    geom_smooth(method="lm")

subject_data  %>% group_by(Metabolite) %>%
    summarise(med=median(Raw),
              var=log(var(Abundance))) %>%
    cor.test(.$med, .$var, data=.)

subject_data  %>% group_by(Metabolite) %>%
    summarise(med=median(Raw), var=log(var(Abundance))) %>%
    lm(.$var ~ .$med, data=.) %>% summary

## Pairwise scatter plots
wide_data <- subject_data %>%
    dplyr::select(-one_of("Code", "Mode", "RunIndex",
                          "Raw", "RawScaled", "Trend", "Type2")) %>%
    spread(key=Metabolite, value=Abundance)

p <- wide_data %>% dplyr::select(sort(sample(8:ncol(wide_data), 2))) %>%
    rename_all(make.names) %>%
    ggplot(aes_string(x=colnames(.)[1], y=colnames(.)[2])) +
    geom_point() + geom_density_2d()
ggMarginal(p, size=5) 

## Look what happens after normalization
norm_subj <- subject_data %>%
    group_by(Metabolite) %>%
    mutate(normalized_abundance = normalize_density(Abundance))

norm_wide <- norm_subj %>%
    dplyr::select(-one_of("Code", "Mode", "RunIndex", "Raw", "RawScaled",
                          "Trend", "Abundance", "Residual")) %>%
    spread(key=Metabolite, value=normalized_abundance)
p <- norm_wide %>% dplyr::select(sort(sample(8:ncol(norm_wide), 2))) %>%
    rename_all(make.names) %>%
    ggplot(aes_string(x=colnames(.)[1], colnames(.)[2])) +
    geom_point() + geom_density2d()
ggMarginal(p, sizef=5)

#######################
## Compute all pairwise correlations and look at some scatter plots
#######################

cormat <- cor(wide_data %>% dplyr::select(-(1:7)), use="pairwise.complete.obs")
covmat <- cov(wide_data %>% dplyr::select(-(1:7)), use="pairwise.complete.obs")
wide_data %>% dplyr::select(-(1:7)) %>% colnames
## Find most correlated
cor_tri <- cormat * upper.tri(diag(108))
indices <- order(abs(cor_tri), decreasing=TRUE)
cor_vec <- cor_tri[indices]
row_idx <- indices %% nrow(cormat)
col_idx <- floor(indices / nrow(cormat)) + 1


length(row_idx)
length(col_idx)
length(cor_vec)
head(names(cor_vec))
length(rownames(cormat)[row_idx])
head(cbind(rownames(cormat)[row_idx], rownames(cormat)[col_idx], cor_vec), n=20)


cormat <- cor(wide_data %>% filter(Type=="AD") %>%
              dplyr::select(-(1:7)),
              use="pairwise.complete.obs", method="spearman")

## Find most correlated
cor_tri <- cormat * upper.tri(diag(108))
indices <- order(abs(cor_tri), decreasing=TRUE)
cor_vec <- cor_tri[indices]
row_idx <- indices %% nrow(cormat)
col_idx <- floor(indices / nrow(cormat)) + 1

head(cor_vec, n=300)

met_id <- 919
p <- wide_data_res %>% dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>%
    rename_all(make.names) %>%
    dplyr::select(1, row_idx[met_id] + 1 , col_idx[met_id] + 1) %>%
    ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) +
    geom_point()
ggMarginal(p, size=5)

head(cormat["Creatine",][order(abs(cormat["Creatine", ]),
                               decreasing=TRUE)], n=50)

met_id <- 919
p <- wide_data_res %>%
    dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>%
    rename_all(make.names) %>%
    dplyr::select(1,  "Kynurenine", "Cystine") %>%
    ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) +
    geom_point()
ggMarginal(p, size=5)

p <- wide_data_res %>% filter(Type2 %in% c("C", "AD")) %>%
    dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>%
    rename_all(make.names) %>%
    dplyr::select(1:2,  "Creatine", "Homocysteine") %>%
    ggplot(aes_string(x=colnames(.)[3], colnames(.)[4], col="Type2")) +
    geom_point()
ggMarginal(p, size=5)

cor.test(as.numeric(unlist(wide_data_res[, "Homocysteine"])),
         as.numeric(unlist(wide_data_res[, "Creatine"])), method="spearman")

#######################
## QC and variance analysis
#######################
colnames(QC_long)

covmat <- QC_long %>% spread(key=Metabolite, value=Abundance) %>%
    filter(Mode == "QC_neg") %>%
    dplyr::select(-one_of(c("Code", "RunIndex", "Mode"))) %>%
    cov(use="pairwise.complete.obs") %>% as.matrix
dim(covmat)

covmat <- QC_long %>% filter(Mode == "QC_neg") %>%
    spread(key=Metabolite, value=Abundance) %>%
    dplyr::select(-one_of(c("Code", "RunIndex", "Mode"))) %>%
    cov(use="pairwise.complete.obs") %>% as.matrix


QC_long %>% spread(key=Metabolite, value=Abundance) %>%
    dplyr::select(Fructose) %>% as.matrix


## Do QC variances correlate with raw sample variances?

## QC variance decreases with increasing intensity!
met_var_qc <- QC_long %>% group_by(Metabolite) %>%
    summarize(qc_var = var(Abundance),
              qc_mad = median(abs(Abundance -
                                  median(Abundance, na.rm=TRUE)), na.rm=TRUE),
              qc_mean=mean(Abundance))

met_var_qc %>% ggplot(aes(x=qc_mean, y=log(qc_var))) +
    geom_point() + geom_smooth(method="lm")

met_var_qc %>% cor.test(.$qc_mean, log(.$qc_var), data=.)

## Subject data variance decreases with increasing intensity!
met_var_subj <-  subject_data %>% group_by(Metabolite) %>%
    summarize(subject_var = var(Abundance, na.rm=TRUE),
              subject_mad = median(abs(Raw - median(Raw, na.rm=TRUE)), na.rm=TRUE),
              subject_mean = mean(Raw, na.rm=TRUE))

met_var_subj %>% ggplot(aes(x=subject_mean, y=log(subject_mad))) +
    geom_point() + geom_smooth(method="lm")

met_var_subj %>% cor.test(.$subject_mean, log(.$subject_var), data=.)

met_var <- left_join(met_var_qc, met_var_subj, by="Metabolite")


## log variances 
met_var %>% ggplot(aes(x=log(qc_var), y=log(subject_var))) +
    geom_point() + geom_smooth(method="lm") + geom_abline()

met_var %>% ggplot(aes(x=log(qc_mad), y=log(subject_mad))) +
    geom_point() + geom_smooth(method="lm") + geom_abline()

met_var %>% ggplot(aes(x=qc_mean, y=subject_mean)) +
    geom_point() + geom_smooth(method="lm")

met_var %>% cor.test(log(.$qc_mad),
                     log(.$subject_mad),
                     data=., method="spearman")

## QC mean intensity and subject means matchup
met_var %>% ggplot(aes(x=qc_mean, y=subject_mean)) +
    geom_point() + geom_smooth(method="lm") + geom_abline()

met_var %>% ggplot(aes(x=qc_mean, y=subject_mean)) +
    geom_point() + geom_smooth(method="lm") + geom_abline()

met_var %>% dplyr::select(2:7) %>% as.matrix %>% cor(method="spearman")

###############################
## Metabolites that show effects with age
###############################


fit_lm <- function(df, type="lm", predictors) {
    form <- as.formula(paste("Abundance", paste(predictors, collapse = " + "), 
                             sep = " ~ "))

    if(sum(is.na(df$Abundance)) > 50)
        return(NULL)
    else {
        if(type=="lm")
            return(lm(form, data=df))
        else if(type == "rq")
            return(rq(form, data=df, tau=0.5))
        
    }
}

pdf("../figs/targeted/Age_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")
age_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "Gender", "APOE"))) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["Age", 1],
              pval=summary(lmres)$coefficients["Age", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(abs(pval))

mets <- age_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=types),
               fill=Type)) +  geom_density_ridges() +
    facet_wrap(~ Metabolite, nrow=4, scales="free") + ylab("")

dev.off()

###############################
## Metabolites that show effects with Gender
###############################

pdf("../figs/targeted/Gender_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")
gender_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "Gender", "APOE"))) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["GenderM", 1],
              pval=summary(lmres)$coefficients["GenderM", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- gender_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Gender), fill=Gender)) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=2, scales="free")
dev.off()

###############################
## AD vs Controls (old)
###############################

pdf("../figs/targeted/AD_comparison.pdf", width=14)
types <- c("CO", "AD")
## Check if control is differen than AD
ad_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "Type", "APOE", "Gender"))) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["TypeCO", 1],
              pval=summary(lmres)$coefficients["TypeCO", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(abs(pval))

mets <- ad_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type), fill=Type)) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=2, scales="free")
dev.off()

###############################
## PD vs Controls (old)
###############################

pdf("../figs/targeted/PD_comparison.pdf", width=14)
types <- c("CO", "PD")
## Check if control is differen than PD
pd_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "Type", "APOE", "Gender"))) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["TypePD", 1],
              pval=summary(lmres)$coefficients["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- pd_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type), fill=Type)) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=2, scales="free")
dev.off()


############################################

pdf("../figs/targeted/APOE_comparison.pdf", width=14)
## APOE
types <- c("AD", "CO")
## Check if control is differen than AD
apoe_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "Type", "APOE", "Gender"))) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["APOE44", 1],
              pval=summary(lmres)$coefficients["APOE44", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- apoe_diff$Met
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(APOE), fill=factor(APOE))) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=2, scales="free")
dev.off()


############ Fit linear model and compute residuals########################



linear_model_info %>% filter(stat == "Pval" & pred == "AD" & Value < 0.05)

linear_model_info %>% filter(stat == "Pval" & pred == "APOE33" & Value < 0.05)

linear_model_info %>% filter(stat == "Pval" & pred == "GenderM" & Value < 0.05)

linear_model_info %>% filter(stat == "Pval" & pred == "PD" & Value < 0.05)
linear_model_info %>% filter(stat == "Pval" & pred == "AD" & Value < 0.05)





tmp_data <- subject_data
tmp_data$Abundance[is.na(tmp_data$Abundance)] <- 0
lm_data <- tmp_data %>% group_by(Metabolite) %>% do(res=lm(Abundance ~ Gender + Type, data=.)$residuals)
tmp_data$residuals <- unlist(lm_data$res)
tst <- tmp_data %>% dplyr::select(-c(Abundance, Raw, Trend, Mode, Code, RunIndex)) %>% spread(key=Metabolite, value=residuals)

glimpse(tst[, 10:30])
glimpse(tst[10, ])

met_var_qc <- QC_long %>% group_by(Metabolite) %>% summarize(qc_var = var(Abundance))
met_var_subj <-  tst %>% gather(key=Metabolite, value=Residual, -(1:7)) %>% group_by(Metabolite) %>% summarize(subject_var = var(Residual, na.rm=TRUE))
met_var <- left_join(met_var_qc, met_var_subj, by="Metabolite")

met_var %>% cor.test(log(.$qc_var), log(.$subject_var),  data=., method="spearman")
met_var %>% ggplot(aes(x=qc_var, y=subject_var)) + geom_point() + geom_smooth(method="lm") + geom_abline()


## Principle component analysis for control subjects at different ages
comparison_groups <- c("CY", "CM", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0

sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Gender", "Age", "APOE", "Batch", "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)
type_vec <- sub_mat$Type
pcres <- sub_mat %>%  select(-one_of("Id", "Type")) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=2, classvec=type_vec)

resid <- sub_mat %>%  select(-one_of("Id", "Type"))
colors = rainbow(length(unique(type_vec)))
names(colors) = unique(type_vec)
ecb = function(x,y){ plot(x, t='n'); text(x, labels=type_vec, col=colors[type_vec]) }
tsne(resid, epoch_callback = ecb)

## M v F
comparison_groups <- c("CY", "CM", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0

sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Type", "Age", "APOE", "Batch", "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)
type_vec <- sub_mat$Gender
pcres <- sub_mat %>%  select(-one_of("Id", "Gender")) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=5, classvec=type_vec)

## Principle component analysis for AD and old controls

comparison_groups <- c("PD", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0
sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Gender", "Age", "APOE", "Batch", "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)
type_vec <- sub_mat$Type
pcres <- sub_mat %>%  select(-one_of("Id", "Type")) %>% as.matrix() %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=2, classvec=type_vec)

comparison_groups <- c("AD", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0
sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Gender", "Age", "APOE", "Batch", "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)
type_vec <- sub_mat$Type
pcres <- sub_mat %>%  select(-one_of("Id", "Type")) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=5, classvec=type_vec)

resid <- sub_mat %>%  select(-one_of("Id", "Type"))
colors = rainbow(length(unique(type_vec)))
names(colors) = unique(type_vec)
ecb = function(x,y){ plot(x, t='n'); text(x, labels=type_vec, col=colors[type_vec]) }
tsne(resid, epoch_callback = ecb)
## Test



CoefficientMat <- longQC %>% group_by(Metabolite) %>% summarize(coef=lm(log(Abundance) ~ Index)$coefficients[2], pval=(summary(lm(log(Abundance) ~ Index)))$coefficients[2, 4])

subject_data[subject_data==0] <- NA
tmp <- sapply(colnames(subject_data)[9:90], function(metabolite) {
    t.test(log(filter(subject_data, Type=="PD")[[metabolite]]),
           log(filter(subject_data, Type=="CO")[[metabolite]]))$p.value
})
names(tmp) <- colnames(subject_data)[9:90]
sort(tmp, decreasing=TRUE)
