library(tidyverse)
library(ggExtra)
library(magrittr)
library(made4)
library(ggridges)
library(gbm3)
library(huge)
library(patchwork)
library(scales)
library(quantreg)
library(CCA)
library(glmnet)

source("utility.R")

processed_files <- dir(pattern="^preprocessed_lipid_data*")

## load("preprocessed_gotms_data.RData")

## Most recent file
load(max(processed_files[grep("-20+", processed_files)]))

## Notes:
## RawScaled: raw intensities
## Abundance: detrended intensiteis (raw intensity - boosted fit on RunIndex)
## Trend: boosted fit on RunIndex



######################################
## Run Index uncorrelated with other features
#####################################

wide_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend", "Batch", "Name",
                          "Id", "Index", "GBAStatus", "GBA_T369M",
                          "cognitive_status")) %>%
    spread(key=Lipid, value=Abundance)

library(rsample)
boots <- bootstraps(wide_data, times = 100)

compute_cancor <- function(split) {

    df <- split$data[split$in_id, ]
    Xmat <- df %>% dplyr::select(Gender, Age, APOE, Type) %>%
        mutate_all(as.numeric) %>%
        as.matrix
    cancor(df$RunIndex, Xmat)$cor
}

names(boots$splits[[1]])
res <- boots %>% mutate(cancor = as.double(purrr::map(splits, compute_cancor)))

summary(res$cancor)

purrr::map(splits, function(x) compute_cancor(x))


compute_cancor2 <- function(split) {

    df <- split$data[split$in_id, ]
    Xmat <- df %>% dplyr::select(Gender, Age, APOE, Type) %>%
        mutate_all(as.numeric) %>% as.matrix
    cancor(df$RunIndex, Xmat)$cor
}



######################################
## Look at median abundance by mode
#####################################

## Raw
subject_data %>%
    group_by(Mode, RunIndex) %>%
    summarize(tot=median(Raw, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap(~Mode)

## Corrected
subject_data %>%
    group_by(Mode, RunIndex) %>%
    summarize(tot=median(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap(~Mode)

## Check ACF for samples
subject_data %>%
    group_by(Mode, RunIndex) %>%
    summarize(tot=median(Raw, na.rm=TRUE)) %>%
    filter(Mode == "Sample") %>% ungroup %>%
    dplyr::select(tot) %>%
    acf(lag.max=100)

subject_data %>% group_by(Mode, RunIndex) %>%
    summarize(tot=median(Abundance, na.rm=TRUE)) %>%
    filter(Mode == "Sample") %>%
    ungroup %>%
        dplyr::select(tot) %>%
        acf(lag.max=100)

## compare total abundance to total raw)
subject_data %>%
    group_by(Mode, RunIndex) %>%
    summarize(SummedAbundance=sum(Abundance, na.rm=TRUE),
              SummedRaw=sum(Raw, na.rm=TRUE)) %>%
    ungroup %>%
    ggplot(aes(x=SummedAbundance, y=SummedRaw, col=Mode)) +
    geom_point()

subject_data %>% group_by(RunIndex) %>%
    summarize(tot=median(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) +
    geom_line() +
    ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)


subject_data %>% filter(RunIndex == 92) %>% dplyr::select(Raw) %>% mutate(exp(Raw)) %>% summary
subject_data %>% filterRunIndex == 93) %>% dplyr::select(Raw) %>% mutate(exp(Raw)) %>% summary
subject_data %>% filter(RunIndex %in% c(94, 93)) %>% dplyr::select(one_of(c("Lipid", "Raw", "RunIndex"))) %>%
    spread(key=RunIndex, value=Raw) %>% ggplot(aes(x=`94`, y=`93`)) + geom_point() + geom_abline()

## check that detrended data looks sensible (for a random sample of metabolites)
lpds <- subject_data %>%
    group_by(Lipid) %>%
    summarize(na_counts = sum(is.na(Abundance))) %>%
    filter(na_counts < 20) %>% dplyr::select(Lipid)

subject_data %>% filter(Lipid %in% sample(lpds$Lipid, 2)) %>%
    ggplot(aes(x=RunIndex, y=Abundance, colour=Lipid)) +
    geom_line() +
    geom_point(aes(shape=Type))# + theme(legend.position="none")

## Which metabolites had the worst drift?
drift_lipids <- subject_data %>%
    group_by(Lipid) %>%
    summarise(var = (quantile(Trend, 0.8, na.rm=TRUE) -
                     quantile(Trend, 0.2, na.rm=TRUE)) /
                  sd(RawScaled, na.rm=TRUE)) %>%
    arrange(desc(var))

## Look at some trends fits for a few metabolites with large trends
lipid_id <- 2
lipid_name <- drift_lipids$Lipid[lipid_id]

subject_data %>% filter(Lipid == lipid_name) %>%
    ggplot(aes(x=RunIndex, y=Batch)) + geom_line()


cols <- hue_pal()(2)
trend_plot <- subject_data %>%
    filter(Lipid %in% lipid_name) %>%
    ggplot(aes(x=RunIndex, y=RawScaled, col="blue")) +
    geom_point(aes(shape=Type), col=cols[1]) +
    geom_line(col=cols[1]) +
    geom_line(aes(x=RunIndex, y=Trend), col=cols[2], size=1) 

ymax <- ggplot_build(trend_plot)$layout$panel_params[[1]]$y.range[2]
trend_plot + geom_segment(aes(x=RunIndex-1/2,
                              xend=RunIndex+1/2,
                              y=ymax,
                              yend=ymax, col=factor(Batch), size=2)) +
    ggtitle(lipid_name)


subject_data %>% filter(Lipid %in% lipid_name) %>%
    ggplot(aes(x=RunIndex, y=Abundance)) +
    geom_point(aes(shape=Type), col=cols[1]) +
    geom_line(col=cols[1]) +
    geom_smooth(method="lm", col=cols[2]) + ggtitle(lipid_name)


## Plot raw intensity versus log variance of the abundance (detrended)
subject_data  %>% group_by(Lipid) %>%
    summarise(MedianIntensity=median(RawScaled), LogVar=log(var(Abundance))) %>%
    ggplot(aes(x=MedianIntensity, y=LogVar)) +
    geom_point() +
    geom_smooth(method="lm") +
    ggtitle("Raw abundance vs log variance")

## Compute correlation
subject_data  %>% group_by(Lipid) %>%
    summarise(med=median(Raw), var=log(var(Abundance))) %>%
    cor.test(.$med, .$var, data=.)

## linear regression
subject_data  %>% group_by(Lipid) %>%
    summarise(med=median(Raw), var=log(var(Abundance))) %>%
    lm(.$var ~ .$med, data=.) %>%
    summary

p <- wide_data %>% dplyr::select(sort(sample(8:ncol(wide_data), 2))) %>%
    rename_all(make.names) %>%
    ggplot(aes_string(x=colnames(.)[1], y=colnames(.)[2])) +
    geom_point() +
    geom_density_2d()
ggMarginal(p, size=5) 

## Look what happens after normalization
norm_subj <- subject_data %>%
    group_by(Lipid) %>%
    mutate(normalized_abundance = normalize_density(Abundance))

norm_wide <- norm_subj %>%
    dplyr::select(-one_of("Code", "Mode", "RunIndex", "Raw",
                          "Trend", "Abundance", "Residual")) %>%
    spread(key=Lipid, value=normalized_abundance)

p <- norm_wide %>% dplyr::select(sort(sample(8:ncol(norm_wide), 2))) %>%
    rename_all(make.names) %>%
    ggplot(aes_string(x=colnames(.)[1], colnames(.)[2])) +
    geom_point() +
    geom_density2d()
ggMarginal(p, size=5)

#######################
## Compute all pairwise correlations and look at some scatter plots
#######################

wide_data2 <- subject_data %>% unite(tmp, Lipid, Mode) %>%
    dplyr::select(-one_of("RunIndex", "Abundance", "Trend", "Name",
                          "Index", "Data File")) %>% spread(key=tmp, value=Raw)

cormat <- cor(wide_data2 %>% dplyr::select(-(1:7)), use="pairwise.complete.obs")

covmat <- cov(wide_data %>% dplyr::select(-(1:7)), use="pairwise.complete.obs")
wide_data %>% dplyr::select(-(1:7)) %>% colnames
## Find most correlated
cor_tri <- cormat * upper.tri(diag(nrow(cormat)))
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


met_id <- 12
p <- wide_data %>% dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>%
    rename_all(make.names) %>%
    dplyr::select(1, row_idx[met_id] + 1 , col_idx[met_id] + 1) %>%
    ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) +
    geom_point()  + abline()
ggMarginal(p, size=5)


p <- QC_long %>%
    spread(key=Lipid, value=Abundance) %>%
    dplyr::select(-c(1, 2, 3)) %>%
    rename_all(make.names) %>%
    dplyr::select(Leucine, Isoleucine) %>%
    ggplot(aes(x=Leucine, y=Isoleucine)) + geom_point()  + abline()
ggMarginal(p, size=5)

## Compute correlations on residuals (todo: correction for attenuation?)
cormat <- cor(wide_data_res %>%
              dplyr::select(-(1:7)),
              use="pairwise.complete.obs",
              method="spearman")

## Find most correlated
cor_tri <- cormat * upper.tri(diag(108))
indices <- order(abs(cor_tri), decreasing=TRUE)
cor_vec <- cor_tri[indices]
row_idx <- indices %% nrow(cormat)
col_idx <- floor(indices / nrow(cormat)) + 1

head(cor_vec, n=300)

met_id <- 919
p <- wide_data_res %>%
    dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>%
    rename_all(make.names) %>%
    dplyr::select(1, row_idx[met_id] + 1 , col_idx[met_id] + 1) %>%
    ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) +
    geom_point()
ggMarginal(p, size=5)


#######################
## QC and variance analysis
#######################

## Do QC variances correlate with raw sample variances?

## QC variance decreases with increasing intensity!
met_var_qc <- QC_long %>% group_by(Lipid) %>%
    summarize(qc_var = var(Abundance),
              qc_mad = median(abs(Abundance - median(Abundance, na.rm=TRUE)), na.rm=TRUE),
              qc_median=median(Abundance))

met_var_qc %>% ggplot(aes(x=qc_median, y=log(qc_var))) +
    geom_point() + geom_smooth(method="lm")

met_var_qc %>% cor.test(.$qc_mean, log(.$qc_var), data=.)

## Subject data variance decreases with increasing intensity!
met_var_subj <-  subject_data %>% group_by(Lipid) %>%
    summarize(subject_var = var(Abundance, na.rm=TRUE),
              subject_mad = median(abs(Raw - median(Raw, na.rm=TRUE)), na.rm=TRUE),
              subject_mean = mean(Raw, na.rm=TRUE))

met_var_subj %>% ggplot(aes(x=subject_mean, y=log(subject_var))) +
    geom_point() + geom_smooth(method="lm")

met_var_subj %>% cor.test(.$subject_mean, log(.$subject_var), data=.)

met_var <- left_join(met_var_qc, met_var_subj, by="Lipid")


## log variances 
met_var %>% ggplot(aes(x=log(qc_var), y=log(subject_var))) +
    geom_point() + geom_smooth(method="lm")
met_var %>% ggplot(aes(x=log(subject_var), y=log(qc_var))) +
    geom_point() + geom_smooth(method="lm") + geom_abline()

geom_abline()
    
met_var %>% cor.test(log(.$qc_var), log(.$subject_var), data=., method="spearman")

## QC mean intensity and subject means matchup
met_var %>% ggplot(aes(x=qc_mean, y=subject_mean)) +
    geom_point() + geom_smooth(method="lm") + geom_abline()

met_var %>% ggplot(aes(x=qc_median, y=subject_mean)) +
    geom_point() + geom_smooth(method="lm") + geom_abline()

###############################
## Lipids that show effects with age
###############################


##  do(lmres=rq(Abundance ~ Age + Type + Gender + APOE, data=., tau=0.5)) %>%

## Null model test
pdf("../figs/lipids/ridx_comparison.pdf", width=14)
types <- c("CO", "CY", "CM")
ridx_diff <- subject_data %>% group_by(Lipid) %>%
    filter(Type %in% types) %>% 
    do(lmres = fit_lm(.)) %>%
    summarize(Lip=Lipid,
              est=ifelse(is.null(lmres),
                         NA,
                         summary(lmres)$coefficients["RunIndex", 1]),
              pval=ifelse(is.null(lmres),
                          NA,
                          summary(lmres)$coefficients["RunIndex", 4])) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

lips <- ridx_diff$Lip[1:16]
subject_data %>% filter(Lipid %in% lips, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(cut(RunIndex, 3)),
               fill=factor(cut(RunIndex, 3)))) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Lipid, nrow=4)
dev.off()

## Group by age
pdf("../figs/lipids/Age_comparison.pdf", width=14)

age_diff <- subject_data %>% group_by(Lipid, Mode) %>%
    filter(Type %in% types) %>%
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender"), type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Lip=Lipid,
              Mode=Mode,
              est=lmres["Age", 1],
              pval=lmres["Age", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

write_csv(age_diff, path="../results/lipid_features/age_features.csv")

lips <- age_diff$Lip[1:16]
subject_data %>% filter(Lipid %in% lips, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=types),
               fill=Type)) +  geom_density_ridges(scale = 4) + theme_bw() +
    facet_wrap(~ Lipid, nrow=4, scales="free")

dev.off()

Yage <- wide_data %>% dplyr::select(Age) %>% as.matrix
X <- wide_data %>% dplyr::select(-one_of(c("Age", "Gender", "Type", "APOE",
                                           "Mode", "RunIndex", "Type2"))) %>%
    as.matrix

X <- X[, colMeans(is.na(X)) < 0.2]
colnames(X) <- make.names(colnames(X))
mice_X <- mice(X)

Xcomp <- complete(mice_X) %>% as.matrix

res <- cv.glmnet(x=Xcomp, y=Yage, family="gaussian", alpha=0.5)
coef(res, res$lambda.min)
plot(res)


###############################
## Lipids that show effects with Gender
###############################


pdf("../figs/lipids/Gender_comparison.pdf", width=14)


gender_diff <- subject_data %>% group_by(Lipid, Mode) %>%
    filter(Type %in% types) %>%
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender"), type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Lip=Lipid,
              Mode=Mode,
              est=lmres["GenderM", 1],
              pval=lmres["GenderM", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

write_csv(gender_diff, path="../results/lipid_features/gender_features.csv")

lips <- gender_diff$Lip[1:16]
subject_data %>% filter(Lipid %in% lips, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Gender, levels=types),
               fill=Gender)) +  geom_density_ridges(scale = 4) + theme_bw() +
    facet_wrap(~ Lipid, nrow=4, scales="free")

dev.off()

Ygender <- wide_data %>% .$Gender

res <- cv.glmnet(x=Xcomp, y=Ygender, family="binomial", alpha=1)
coef(res, res$lambda.min)
plot(res)
predict(res, Xcomp)


###############################
## AD vs Controls (old)
###############################

pdf("../figs/lipids/AD_comparison.pdf", width=14)
types <- c("CO", "AD")

ad_diff <- subject_data %>% group_by(Lipid, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres = fit_lm(., type="lm",
                      predictors=c("Age", "Type", "APOE", "Gender"))) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Lip=Lipid,
              Mode=Mode,
              est=lmres["TypeCO", 1],
              pval=lmres["TypeCO", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

write_csv(ad_diff, path="../results/lipid_features/ad_features.csv")

lips <- ad_diff$Lip[1:16]
subject_data %>% filter(Lipid %in% lips, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = Type,
               fill=Type)) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Lipid, nrow=4)

dev.off()

Ytype <- wide_data %>% filter(Type %in% c("AD", "CO")) %>% .$Type2 %>% droplevels
Xcomp_ad <- Xcomp[wide_data$Type %in% c("AD", "CO"),]
res <- cv.glmnet(x=Xcomp_ad, y=Ytype, family="binomial", alpha=1)
res2 <- glmnet(x=Xcomp_ad, y=Ytype, family="binomial", alpha=1, lambda=res$lambda.min)
coef(res, res$lambda.min)
par(mfrow=c(1,2))
preds <- predict(res2, Xcomp_ad, type="response")
table(Ytype)

library(broom)
library(AUC)
r <- roc(preds, Ytype)
auc_val <- auc(r)

td <- tidy(r)
td

library(ggplot2)
  
ggplot(td, aes(fpr, tpr)) +  geom_line() + theme_bw() + geom_abline(lty=2)


hist(Xcomp_ad %*%  coef(res, res$lambda.min)[-1,]  + coef(res, res$lambda.min)[1])
Xcomp_ad[, which(coef(res, res$lambda.min) != 0)]


###############################
## PD vs Controls (old)
###############################

pdf("../figs/lipids/PD_comparison.pdf", width=14)

types <- c("CO", "PD")

pd_diff <- subject_data %>% group_by(Lipid, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres = fit_lm(., type="lm",
                      predictors=c("Age", "Type", "APOE", "Gender"))) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Lip=Lipid,
              Mode=Mode,
              est=lmres["TypePD", 1],
              pval=lmres["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

write_csv(pd_diff, path="../results/lipid_features/pd_features.csv")

lips <- pd_diff$Lip[1:16]
subject_data %>% filter(Lipid %in% lips, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = Type,
               fill=Type)) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Lipid, nrow=4)

dev.off()

Ytype <- wide_data %>%
    filter(Type %in% c("PD", "CO")) %>% .$Type2 %>% droplevels

Xcomp_pd <- Xcomp[wide_data$Type %in% c("PD", "CO"),]
res <- cv.glmnet(x=Xcomp_pd, y=Ytype, family="binomial", alpha=0.5)
coef(res, res$lambda.min)
hist(predict(res, Xcomp_pd, type="response"), breaks=50)

res2 <-glmnet(x=Xcomp_pd, y=Ytype, family="binomial", lambda=res$lambda.min)




############################################

pdf("../figs/lipids/APOE_comparison.pdf", width=14)
## APOE
types <- c("AD", "CO")
## Check if control is differen than AD
apoe_diff <- subject_data %>%
    group_by(Lipid, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "Type", "APOE", "Gender"))) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Lip=Lipid,
              Mode=Mode,
              est=lmres["APOE44", 1],
              pval=lmres["APOE44", 4]) %>%
    filter(pval < 0.05) 

write_csv(apoe_diff, path="../results/lipid_features/apoe_features.csv")

lips <- apoe_diff$Lip[1:16]
meds <- subject_data %>%
    filter(Lipid %in% lips, Type %in% types) %>%
    group_by(APOE, Lipid) %>%
    summarise(med = quantile(Abundance, 0.25, na.rm=TRUE)) %>%
    ungroup
p <- subject_data %>%
    filter(Lipid %in% lips, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(APOE),
               fill=factor(APOE))) +
    geom_density_ridges(scale = 4) + facet_wrap(~ Lipid, nrow=4)
p + geom_segment(data=meds, aes(x=med, xend=med, y=0, yend=0.75, color=APOE),
                 size=1)
dev.off()


############ Look at correlation in quality control data ########################

QCcor <- QC_long %>% unite(MetMode, Lipid, Mode) %>% spread(key=MetMode, value=Abundance) %>%
    dplyr::select(-one_of(c("Name", "Data File", "RunIndex", "SampleId"))) %>%
    cor(use="pairwise.complete.obs") %>% as.matrix

cor_tri <- QCcor * upper.tri(diag(nrow(QCcor)))
plot(sort(cor_tri[upper.tri(diag(nrow(QCcor)))], decreasing=TRUE))
hist(cor_tri[upper.tri(diag(nrow(QCcor)))], breaks=100)



p <- QC_long %>% spread(key=Lipid, value=Abundance) %>% dplyr::select(-c(1, 2, 3)) %>% rename_all(make.names) %>% dplyr::select(Leucine, Isoleucine) %>% ggplot(aes(x=Leucine, y=Isoleucine)) + geom_point()  + abline()
ggMarginal(p, size=5)

## Compute correlations on residuals (todo: correction for attenuation?)
cormat <- cor(wide_data_res %>% dplyr::select(-(1:7)), use="pairwise.complete.obs", method="spearman")


