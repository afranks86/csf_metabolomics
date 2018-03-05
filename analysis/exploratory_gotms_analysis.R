library(tidyverse)
library(ggExtra)
library(magrittr)
library(made4)
library(ggjoy)
library(gbm3)
library(huge)
library(patchwork)

source("utility.R")

load("preprocessed_gotms_data.RData")


## Notes:
## Raw: raw intensities
## Abundance: detrended intensiteis (raw intensity - boosted fit on RunIndex)
## Trend: boosted fit on RunIndex

## Look at average raw data
subject_data %>% group_by(Mode, RunIndex) %>% summarize(tot=sum(exp(Raw), na.rm=TRUE)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Positive") + geom_smooth(method="lm", se=FALSE)

subject_data %>% group_by(Mode, RunIndex) %>% summarize(tot=sum(exp(Raw), na.rm=TRUE)) %>%  filter(Mode == "pos02") %>% ungroup %>% dplyr::select(tot) %>% acf(lag.max=100)

subject_data %>% group_by(Mode, RunIndex) %>% summarize(tot=sum(Abundance, na.rm=TRUE)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Positive") + geom_smooth(method="lm", se=FALSE)

subject_data %>% group_by(Mode, RunIndex) %>% summarize(tot=sum(Abundance, na.rm=TRUE), tot2=log(sum(exp(Raw), na.rm=TRUE))) %>% ungroup %>% ggplot(aes(x=tot, y=tot2, col=Mode)) + geom_point()

## Focus on raw neg
subject_data %>% filter(Mode=="neg") %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Raw, na.rm=TRUE)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Negative") + geom_smooth(method="lm", se=FALSE)

## Focus on raw pos
subject_data %>% filter(Mode=="pos01") %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Raw, na.rm=TRUE)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Negative") + geom_smooth(method="lm", se=FALSE)

subject_data %>% filter(Mode=="pos01") %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Abundance, na.rm=TRUE)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Negative") + geom_smooth(method="lm", se=FALSE)

## check that detrended data looks sensible (for a random sample of metabolites)
subject_data %>% filter(Metabolite %in% sample(unique(subject_data$Metabolite), 2)) %>% ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) + geom_line() + geom_point(aes(shape=Type))# + theme(legend.position="none")

## Which metabolites had the worst drift?
drift_mets <- subject_data %>% group_by(Metabolite, Mode) %>% summarise(var = (max(Trend, na.rm=TRUE) - min(Trend, na.rm=TRUE))/sd(Raw, na.rm=TRUE)) %>% arrange(desc(var))

## Look at some trends fits for a few metabolites with large trends

batch_positions <- subject_data %>% group_by(Batch, Mode) %>% summarise(mx = max(RunIndex))

metid <- 50
intercepts <- batch_positions %>% ungroup() %>% filter(Mode == batch_positions$Mode[5]) %>% dplyr::select(mx) %>% unlist %>% as.numeric
subject_data %>% filter(Metabolite %in% drift_mets$Metabolite[metid] & Mode == drift_mets$Mode[metid]) %>% ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) + geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=Raw)) + geom_vline(xintercept=intercepts)

## Plot raw intensity versus log variance of the abundance (detrended)
subject_data  %>% group_by(Metabolite) %>% summarise(med=median(Raw), var=log(var(Abundance))) %>% ggplot(aes(x=med, y=var)) + geom_point() + geom_smooth(method="lm")
subject_data  %>% group_by(Metabolite) %>% summarise(med=median(Raw), var=log(var(Abundance))) %>% cor.test(.$med, .$var, data=.)
subject_data  %>% group_by(Metabolite) %>% summarise(med=median(Raw), var=log(var(Abundance))) %>% lm(.$var ~ .$med, data=.) %>% summary

## Pairwise scatter plots
wide_data <- subject_data %>% unite(tmp, Metabolite, Mode) %>% dplyr::select(-one_of("RunIndex", "Raw", "Trend", "Name", "Index", "Data File")) %>% spread(key=tmp, value=Abundance)

p <- wide_data %>% dplyr::select(sort(sample(8:ncol(wide_data), 2))) %>% rename_all(make.names) %>% ggplot(aes_string(x=colnames(.)[1], y=colnames(.)[2])) + geom_point() + geom_density_2d()
ggMarginal(p, size=5) 

## Look what happens after normalization
norm_subj <- subject_data %>% group_by(Metabolite) %>% mutate(normalized_abundance = normalize_density(Abundance))
norm_wide <- norm_subj %>% dplyr::select(-one_of("Code", "Mode", "RunIndex", "Raw", "Trend", "Abundance", "Residual")) %>% spread(key=Metabolite, value=normalized_abundance)
p <- norm_wide %>% dplyr::select(sort(sample(8:ncol(norm_wide), 2))) %>% rename_all(make.names) %>% ggplot(aes_string(x=colnames(.)[1], colnames(.)[2])) + geom_point() + geom_density2d()
ggMarginal(p, size=5)

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


met_id <- 12
p <- wide_data %>% dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>% rename_all(make.names) %>% dplyr::select(1, row_idx[met_id] + 1 , col_idx[met_id] + 1) %>% ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) + geom_point()  + abline()
ggMarginal(p, size=5)



p <- QC_long %>% spread(key=Metabolite, value=Abundance) %>% dplyr::select(-c(1, 2, 3)) %>% rename_all(make.names) %>% dplyr::select(Leucine, Isoleucine) %>% ggplot(aes(x=Leucine, y=Isoleucine)) + geom_point()  + abline()
ggMarginal(p, size=5)

## Compute correlations on residuals (todo: correction for attenuation?)
cormat <- cor(wide_data_res %>% dplyr::select(-(1:7)), use="pairwise.complete.obs", method="spearman")

## cormat <- cor(wide_data_res %>% filter(Type=="AD") %>% dplyr::select(-(1:7)), use="pairwise.complete.obs", method="spearman")

## Find most correlated
cor_tri <- cormat * upper.tri(diag(108))
indices <- order(abs(cor_tri), decreasing=TRUE)
cor_vec <- cor_tri[indices]
row_idx <- indices %% nrow(cormat)
col_idx <- floor(indices / nrow(cormat)) + 1

head(cor_vec, n=300)

met_id <- 919
p <- wide_data_res %>% dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>% rename_all(make.names) %>% dplyr::select(1, row_idx[met_id] + 1 , col_idx[met_id] + 1) %>% ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) + geom_point()
ggMarginal(p, size=5)

head(cormat["Creatine",][order(abs(cormat["Creatine", ]), decreasing=TRUE)], n=50)
met_id <- 919
p <- wide_data_res %>% dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>% rename_all(make.names) %>% dplyr::select(1,  "Kynurenine", "Cystine") %>% ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) + geom_point()
ggMarginal(p, size=5)

p <- wide_data_res %>% filter(Type2 %in% c("C", "AD")) %>% dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>% rename_all(make.names) %>% dplyr::select(1:2,  "Creatine", "Homocysteine") %>% ggplot(aes_string(x=colnames(.)[3], colnames(.)[4], col="Type2")) + geom_point()
ggMarginal(p, size=5)



#######################
## QC and variance analysis
#######################

## Do QC variances correlate with raw sample variances?

## QC variance decreases with increasing intensity!
met_var_qc <- QC_long %>% group_by(Metabolite) %>%
    summarize(qc_var = var(Abundance),
              qc_mad = median(abs(Abundance - median(Abundance, na.rm=TRUE)), na.rm=TRUE),
              qc_mean=mean(Abundance))

met_var_qc %>% ggplot(aes(x=qc_mean, y=log(qc_var))) + geom_point() + geom_smooth(method="lm")

met_var_qc %>% cor.test(.$qc_mean, log(.$qc_var), data=.)

## Subject data variance decreases with increasing intensity!
met_var_subj <-  subject_data %>% group_by(Metabolite) %>%
    summarize(subject_var = var(Abundance, na.rm=TRUE),
              subject_mad = median(abs(Raw - median(Raw, na.rm=TRUE)), na.rm=TRUE),
              subject_mean = mean(Raw, na.rm=TRUE))

met_var_subj %>% ggplot(aes(x=subject_mean, y=log(subject_var))) + geom_point() + geom_smooth(method="lm")

met_var_subj %>% cor.test(.$subject_mean, log(.$subject_var), data=.)

met_var <- left_join(met_var_qc, met_var_subj, by="Metabolite")


## log variances 
met_var %>% ggplot(aes(x=log(qc_var), y=log(subject_var))) + geom_point() + geom_smooth(method="lm") + geom_abline() 
met_var %>% cor.test(log(.$qc_var), log(.$subject_var), data=.)

## QC mean intensity and subject means matchup
met_var %>% ggplot(aes(x=qc_mean, y=subject_mean)) + geom_point() + geom_smooth(method="lm") + geom_abline()

met_var %>% ggplot(aes(x=qc_mean, y=subject_mean)) + geom_point() + geom_smooth(method="lm") + geom_abline()

###############################
## Metabolites that show effects with age
###############################

pdf("../figs/gotms/Age_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")
age_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              Mode=Mode,
              est=summary(lmres)$coefficients["Age", 1],
              pval=summary(lmres)$coefficients["Age", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- age_diff$Met[1:32]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=types),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=8)

subject_data %>% group_by(Metabolite) %>%
    do(lmres=lm(Abundance ~ Age, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients[2, 1],
              pval=summary(lmres)$coefficients[2, 4]) %>%
    filter(pval < 0.05) %>%
    arrange(est)
dev.off()

metid <- 200
intercepts <- batch_positions %>% ungroup() %>% filter(Mode == batch_positions$Mode[5]) %>% dplyr::select(mx) %>% unlist %>% as.numeric

subject_data %>% filter(Metabolite == "1061 Results") %>%
    filter(Type %in% types) %>% ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) + geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=Raw)) + geom_point(aes(x=RunIndex, y=Raw, shape=Type)) + geom_line(aes(x=RunIndex, y=Abundance)) + geom_point(aes(x=RunIndex, y=Abundance, shape=Type)) 

subject_data %>% filter(Metabolite == "734 Results") %>% ggplot(aes(x=RunIndex, y=Abundance, colour=Type)) + geom_point()
subject_data %>% filter(Metabolite == "734 Results") %>% arrange(RunIndex) %>% dplyr::select(Trend) %>% acf(lag.max=50)

###############################
## Metabolites that show effects with Gender
###############################

pdf("../figs/gotms/Gender_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")
gender_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              Mode=Mode,
              est=summary(lmres)$coefficients["GenderM", 1],
              pval=summary(lmres)$coefficients["GenderM", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- gender_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Gender),
               fill=Gender)) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=2)
dev.off()

###############################
## AD vs Controls (old)
###############################

pdf("../figs/gotms/AD_comparison.pdf", width=14)
types <- c("CO", "AD")
## Check if control is differen than AD
ad_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=summary(lmres)$coefficients["TypeCO", 1],
              pval=summary(lmres)$coefficients["TypeCO", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- ad_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=2)
dev.off()

###############################
## PD vs Controls (old)
###############################

pdf("../figs/gotms/PD_comparison.pdf", width=14)
types <- c("CO", "PD")
## Check if control is differen than PD
pd_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=summary(lmres)$coefficients["TypePD", 1],
              pval=summary(lmres)$coefficients["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- pd_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=2)
dev.off()

###############################
## AD vs PD
###############################

pdf("../figs/gotms/AD_PD_comparison.pdf", width=14)
types <- c("AD", "PD")
## Check if control is differen than PD
ad_pd_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=summary(lmres)$coefficients["TypePD", 1],
              pval=summary(lmres)$coefficients["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- ad_pd_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=2)
dev.off()



############################################

pdf("../figs/gotms/APOE_comparison.pdf", width=14)
## APOE
types <- c("AD", "CO")
## Check if control is differen than AD
apoe_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Gender + factor(APOE), data=.)) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=summary(lmres)$coefficients["factor(APOE)44", 1],
              pval=summary(lmres)$coefficients["factor(APOE)44", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- apoe_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(APOE),
               fill=factor(APOE))) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=4)
dev.off()

## Some potentially interesting candidates

lst1 <- Reduce(intersect, list(ad_diff$Met, apoe_diff$Met))

lst2 <- Reduce(intersect , list(ad_pd_diff$Met, pd_diff$Met, gender_diff$Met))

lst3 <- Reduce(intersect, list(age_diff$Met, gender_diff$Met))

ij1 <- inner_join(age_diff, gender_diff, by=c("Met", "Mode"))
ij2 <- inner_join(ad_diff, apoe_diff, by=c("Met", "Mode"))
ij3 <- inner_join(ad_pd_diff, pd_diff, by=c("Met", "Mode"))

interesting_features <- bind_rows(ij1, ij2, ij3) %>% dplyr::select(Met, Mode) %>% distinct()


write.csv(interesting_features, file="~/Desktop/csf_ids.csv")


############ Look at correlation in quality control data ########################

QCcor <- QC_long %>% unite(MetMode, Metabolite, Mode) %>% spread(key=MetMode, value=Abundance) %>%
    dplyr::select(-one_of(c("Name", "Data File", "RunIndex", "SampleId"))) %>%
    cor(use="pairwise.complete.obs") %>% as.matrix

cor_tri <- QCcor * upper.tri(diag(nrow(QCcor)))
plot(sort(cor_tri[upper.tri(diag(nrow(QCcor)))], decreasing=TRUE))
hist(cor_tri[upper.tri(diag(nrow(QCcor)))], breaks=100)



p <- QC_long %>% spread(key=Metabolite, value=Abundance) %>% dplyr::select(-c(1, 2, 3)) %>% rename_all(make.names) %>% dplyr::select(Leucine, Isoleucine) %>% ggplot(aes(x=Leucine, y=Isoleucine)) + geom_point()  + abline()
ggMarginal(p, size=5)

## Compute correlations on residuals (todo: correction for attenuation?)
cormat <- cor(wide_data_res %>% dplyr::select(-(1:7)), use="pairwise.complete.obs", method="spearman")


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
pcres <- sub_mat %>%  select(-one_of("Id", "Type")) %>% prcomp
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
