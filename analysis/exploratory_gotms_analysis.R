library(tidyverse)
library(ggExtra)
library(magrittr)
library(made4)
library(ggridges)
library(quantreg)
library(gbm3)
library(huge)
library(patchwork)
library(scales)

source("utility.R")

processed_files <- dir(pattern="^preprocessed_gotms_data*")

## load("preprocessed_gotms_data.RData")

## Most recent file
nms <- load(max(processed_files[grep("-20+", processed_files)]))

load("../data/got-ms/identification_map.RData")


## Notes:
## RawScaled: raw intensities
## Abundance: detrended intensiteis (raw intensity - boosted fit on RunIndex)
## Trend: boosted fit on RunIndex

##############################
## Look at median abundance by mode
##############################

## Plot Raw
subject_data %>% group_by(Mode, RunIndex) %>%
    summarize(tot=median(Raw, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap(~Mode)

## Plot Corrected
subject_data %>% group_by(Mode, RunIndex) %>%
    summarize(tot=median(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot, colour=Mode)) +
    geom_line() +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap(~Mode)

## Check ACF for pos03
subject_data %>% group_by(Mode, RunIndex) %>%
    summarize(tot=median(Raw, na.rm=TRUE)) %>%
    filter(Mode == "pos03") %>% ungroup %>%
        dplyr::select(tot) %>% acf(lag.max=100)

## compare total abundance to total raw)
subject_data %>% group_by(Mode, RunIndex) %>%
    summarize(SummedAbundance=sum(Abundance, na.rm=TRUE),
              SummedRaw=sum(Raw, na.rm=TRUE)) %>%
    ungroup %>%
        ggplot(aes(x=SummedAbundance, y=SummedRaw, col=Mode)) + geom_point()

## Focus on raw neg
subject_data %>%
    filter(Mode=="neg") %>%
    group_by(RunIndex) %>%
    summarize(tot=sum(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) + geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)

## Focus on raw pos 01
subject_data %>% filter(Mode=="pos01") %>%
    group_by(RunIndex) %>%
    summarize(tot=sum(Raw, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) +
    geom_line() + ggtitle("Negative") + geom_smooth(method="lm", se=FALSE)

subject_data %>% filter(Mode=="pos01") %>%
    group_by(RunIndex) %>%
    summarize(tot=sum(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) + geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)

## Focus on raw pos 02
subject_data %>%
    filter(Mode=="pos02") %>%
    group_by(RunIndex) %>%
    summarize(tot=sum(Raw, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) + geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)

subject_data %>% filter(Mode=="pos02") %>%
    group_by(RunIndex) %>%
    summarize(tot =sum(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) + geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)

## Focus on raw pos 03
subject_data %>% filter(Mode=="pos03") %>%
    group_by(RunIndex) %>%
    summarize(tot=sum(Raw, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) + geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)

subject_data %>% filter(Mode=="pos03") %>%
    group_by(RunIndex) %>%
    summarize(tot=median(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) + geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)

subject_data %>% filter(Mode=="pos03") %>%
    group_by(RunIndex) %>%
    summarize(min=sum(Abundance, na.rm=TRUE), ri = unique(RunIndex)) %>%
    dplyr::select(min) %>% unlist %>% which.min

subject_data %>% filter(Mode=="pos03", RunIndex == 92) %>%
    dplyr::select(Raw) %>% mutate(exp(Raw)) %>% summary
subject_data %>% filter(Mode=="pos03", RunIndex == 93) %>%
    dplyr::select(Raw) %>% mutate(exp(Raw)) %>% summary

## Plot metabolites against one another
subject_data %>%
    filter(Mode=="pos03", RunIndex %in% c(94, 93)) %>%
    dplyr::select(one_of(c("Metabolite", "Raw", "RunIndex"))) %>%
    spread(key=RunIndex, value=Raw) %>%
    ggplot(aes(x=`94`, y=`93`)) + geom_point() + geom_abline()

## check that detrended data looks sensible (for a random sample of metabolites)
subject_data %>%
    filter(Metabolite %in% sample(unique(subject_data$Metabolite), 2)) %>%
    ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) +
    geom_line() +
    geom_point(aes(shape=Type))# + theme(legend.position="none")

## Which metabolites had the worst drift?
drift_mets <- subject_data %>%
    group_by(Metabolite, Mode) %>%
    summarise(var = (max(Trend, na.rm=TRUE) - min(Trend, na.rm=TRUE)) /
                  sd(RawScaled, na.rm=TRUE)) %>%
    arrange(desc(var))

## Look at some trends fits for a few metabolites with large trends

metid <- 2
met_name <- drift_mets$Metabolite[metid]
mode <- drift_mets$Mode[metid]


subject_data %>%
    filter(Metabolite == met_name & Mode == mode) %>%
    ggplot(aes(x=RunIndex, y=Batch)) + geom_line()

qc_indices <- QC_long %>%
    filter(grepl(mode, .$Mode)) %>%
    dplyr::select(RunIndex) %>%
    unique %>% unlist

QC_long %>% filter(Metabolite == drift_mets$Metabolite[metid]) %>%
    ggplot(aes(x=RunIndex, y=Abundance)) + geom_line()


## Large drift metabolite
cols <- hue_pal()(2)
trend_plot <- subject_data %>%
    filter(Metabolite %in% met_name & Mode == mode) %>%
    ggplot(aes(x=RunIndex, y=RawScaled, col="blue")) +
    geom_point(aes(shape=Type), col=cols[1]) +
    geom_line(col=cols[1]) +
    geom_line(aes(x=RunIndex, y=Trend), col=cols[2], size=1)


subject_data %>% filter(Metabolite %in% met_name & Mode == mode) %>%
    ggplot(aes(x=RunIndex, y=Abundance)) +
    geom_point(aes(shape=Type), col=cols[1]) +
    geom_line(col=cols[1]) +
    geom_smooth(method="lm", col=cols[2]) + ggtitle(met_name)


## Plot raw intensity versus log variance of the abundance (detrended)
subject_data  %>% group_by(Metabolite) %>%
    summarise(MedianIntensity=median(RawScaled), LogVar=log(var(Abundance))) %>%
    ggplot(aes(x=MedianIntensity, y=LogVar)) +
    geom_point() +
    geom_smooth(method="lm") +
    ggtitle("Raw abundance vs log variance")

## Compute correlation
subject_data  %>% group_by(Metabolite) %>%
    summarise(med=median(Raw), var=log(var(Abundance))) %>%
    cor.test(.$med, .$var, data=.)

## linear regression
subject_data  %>% group_by(Metabolite) %>%
    summarise(med=median(Raw), var=log(var(Abundance))) %>%
    lm(.$var ~ .$med, data=.) %>%
    summary

## Pairwise scatter plots
wide_data <- subject_data %>%
    unite(tmp, Metabolite, Mode) %>%
    dplyr::select(-one_of("RunIndex", "Raw", "Trend",
                          "Name", "Index", "Data File")) %>%
    spread(key=tmp, value=Abundance)

p <- wide_data %>%
    dplyr::select(sort(sample(8:ncol(wide_data), 2))) %>%
    rename_all(make.names) %>%
    ggplot(aes_string(x=colnames(.)[1], y=colnames(.)[2])) +
    geom_point() +
    geom_density_2d()
ggMarginal(p, size=5) 

## Look what happens after normalization
norm_subj <- subject_data %>%
    group_by(Metabolite) %>%
    mutate(normalized_abundance = normalize_density(Abundance))
norm_wide <- norm_subj %>%
    dplyr::select(-one_of("Code", "Mode", "RunIndex",
                          "Raw", "Trend", "Abundance", "Residual")) %>%
    spread(key=Metabolite, value=normalized_abundance)
p <- norm_wide %>% dplyr::select(sort(sample(8:ncol(norm_wide), 2))) %>%
    rename_all(make.names) %>%
    ggplot(aes_string(x=colnames(.)[1], colnames(.)[2])) +
    geom_point() +
    geom_density2d()
ggMarginal(p, size=5)

#######################
## Compute all pairwise correlations and look at some scatter plots
#######################

wide_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend",
                          "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=Abundance)
dim(wide_data)

cormat <- cor(wide_data %>% dplyr::select(-(1:11)), use="pairwise.complete.obs")

covmat <- cov(wide_data %>% dplyr::select(-(1:11)), use="pairwise.complete.obs")

## Find most correlated
cor_tri <- cormat * upper.tri(diag(nrow(cormat)))
indices <- order(abs(cor_tri), decreasing=TRUE)
cor_vec <- cor_tri[indices]
row_idx <- indices %% nrow(cormat)
col_idx <- floor(indices / nrow(cormat)) + 1



wide_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend",
                          "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=Abundance)
cormat <- cor(wide_data %>% dplyr::select(-(1:11)), use="pairwise.complete.obs")

## Find most negatively correlated
cor_tri <- cormat * upper.tri(diag(nrow(cormat)))
indices <- order(cor_tri)
cor_vec <- cor_tri[indices]
row_idx <- indices %% nrow(cormat)
col_idx <- floor(indices / nrow(cormat)) + 1

row_idx[1:50]
col_idx[1:50]
lapply(1:100, function(i) cormat[row_idx[i], col_idx[i], drop=FALSE])


sub_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "Abundance", "Trend",
                          "Name","Data File")) %>%
    spread(key=Metabolite, value=RawScaled) %>%
    dplyr::select(`RunIndex`, APOE,
                  `1388-7.043 Results_pos03`,
                  `325 Results_pos03`) %>%
    na.omit


p1 <- sub_data %>%
    ggplot() +
    geom_point(aes(x=RunIndex,
                   y=`1388-7.043 Results_pos03`,
                   col=APOE,
                   shape=factor(abs(`1388-7.043 Results_pos03`) > 2)), size=2) +
    scale_color_discrete_sequential(palette="Viridis")
p2 <- sub_data %>%
    ggplot() +
    geom_point(aes(x=RunIndex, y=`325 Results_pos03`,
                   col=APOE,
                   shape=factor(abs(`1388-7.043 Results_pos03`) > 2)), size=2) +
    theme(legend.pos='None') + scale_color_discrete_sequential(palette="Viridis")
p3 <- sub_data %>%
    ggplot() +
    geom_point(aes(x=RunIndex, y= scale(`1388-7.043 Results_pos03`) - scale(`325 Results_pos03`),
                   col=APOE,
                   shape=factor(abs(`1388-7.043 Results_pos03`) > 2))) +
    theme(legend.pos='None') + scale_color_discrete_sequential(palette="Viridis")
p1 + p2 + p3

wide_data %>% dplyr::select(`1388-7.043 Results_pos03`, `325 Results_pos03`) %>%
    ggplot() + geom_point(aes(x=1:198,
                              y=`1388-7.043 Results_pos03`+`325 Results_pos03`))

wide_data %>% dplyr::select(`18 Results_neg`, `35 Results_neg`) %>%
    ggplot() + geom_point(aes(x=`18 Results_neg`, y=`35 Results_neg`))

lapply(101:200, function(z) {
    print(z)
    cormat[row_idx[z], col_idx[1], drop=FALSE]
    })

plot(eigen(cormat)$values)
evec <- eigen(cormat)$vector[, 2]
names(evec) <- rownames(cormat)
sort(abs(evec), decreasing=TRUE)

length(row_idx)
length(col_idx)
length(cor_vec)
head(names(cor_vec))
length(rownames(cormat)[row_idx])
head(cbind(rownames(cormat)[row_idx], rownames(cormat)[col_idx], cor_vec), n=20)


met_id <- 12
p <- wide_data %>%
    dplyr::select(-c(1, 3, 4, 5, 6, 7)) %>%
    rename_all(make.names) %>%
    dplyr::select(1, row_idx[met_id] + 1 , col_idx[met_id] + 1) %>%
    ggplot(aes_string(x=colnames(.)[2], colnames(.)[3], col="Type")) +
    geom_point()  +
    abline()

ggMarginal(p, size=5)









#######################
## QC and variance analysis
#######################

## Do QC variances correlate with raw sample variances?

## QC variance decreases with increasing intensity!
met_var_qc <- QC_long %>%
    group_by(Metabolite) %>%
    summarize(qc_var = var(Abundance),
              qc_mad = median(abs(Abundance - median(Abundance, na.rm=TRUE)),
                              na.rm=TRUE),
              qc_median=median(Abundance))

met_var_qc %>% ggplot(aes(x=qc_median, y=log(qc_var))) +
    geom_point() + geom_smooth(method="lm")

met_var_qc %>% cor.test(.$qc_mean, log(.$qc_var), data=.)

## Subject data variance decreases with increasing intensity!
met_var_subj <-  subject_data %>% group_by(Metabolite) %>%
    summarize(subject_var = var(Abundance, na.rm=TRUE),
              subject_mad = median(abs(Raw - median(Raw, na.rm=TRUE)), na.rm=TRUE),
              subject_mean = mean(Raw, na.rm=TRUE))

met_var_subj %>% ggplot(aes(x=subject_mean, y=log(subject_var))) +
    geom_point() + geom_smooth(method="lm")

met_var_subj %>% cor.test(.$subject_mean, log(.$subject_var), data=.)

met_var <- left_join(met_var_qc, met_var_subj, by="Metabolite")


## log variances 
met_var %>% ggplot(aes(x=log(qc_var), y=log(subject_var))) +
    geom_point() +
    geom_smooth(method="lm")
met_var %>% ggplot(aes(x=log(subject_var), y=log(qc_var))) +
    geom_point() +
    geom_smooth(method="lm") + geom_abline()

geom_abline()
    
met_var %>% cor.test(log(.$qc_var), log(.$subject_var), data=., method="spearman")

## QC mean intensity and subject means matchup
met_var %>% ggplot(aes(x=qc_mean, y=subject_mean)) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_abline()

met_var %>% ggplot(aes(x=qc_median, y=subject_mean)) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_abline()



###### Look at variance by type
subject_data %>% group_by(Metabolite, Type) %>%
    summarise(var = var(Abundance), abund=mean(Raw)) %>%
    filter(Type %in% c("CO", "AD", "PD")) %>%
    ggplot() + geom_point(aes(x=abund, y=var, col=Type)) 


###############################
## Null Test
###############################

types <- c("CY", "CM", "CO")
## Null model test
pdf("../figs/gotms/ridx_comparison.pdf", width=14)
ridx_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ RunIndex + Age + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              Mode=Mode,
              est=summary(lmres)$coefficients["RunIndex", 1],
              pval=summary(lmres)$coefficients["RunIndex", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- ridx_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(cut(RunIndex, 3)),
               fill=factor(cut(RunIndex, 3)))) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=4)
dev.off()

###############################
## Metabolites that show effects with Gender
###############################

pdf("../figs/gotms/Gender_comparison.pdf", width=14)

types <- c("CO", "CM", "CY")
gender_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., type="rq", predictors = c("Age", "APOE", "Gender"))) %>%
    summarize(Met=Metabolite,
              Mode=Mode,
              est=lmres["GenderM", 1],
              pval=lmres["GenderM", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)


gender_diff_matches <- gender_diff %>%
    mutate(Mode2 = str_sub(Mode, 1, 3),
           Metabolite = str_replace(Met, " Results", "")) %>% 
    left_join(all_matches %>% rename(MetaboliteName = Metabolite),
              by=c("Mode2"="Mode", "Metabolite" = "Name"))

write_csv(gender_diff_matches,
          path="../results/gotms_features/gender_features.csv")

gender_diff_matches <- gender_diff_matches %>% slice(1:16) 


subject_data %>%
    left_join(gender_diff_matches, .,
              by=c("Met"="Metabolite", "Mode"="Mode")) %>%
    filter(Type %in% types) %>%
    ggplot(aes(x = Abundance, y = Gender, fill=Gender)) +
    geom_density_ridges(scale = 4, alpha=0.75, na.rm=FALSE) +
    facet_wrap(~ MetaboliteName + Metabolite, nrow=2, scales="free") +
    theme_bw()


dev.off()

###############################
## Metabolites that show effects with Age
###############################

pdf("../figs/gotms/Age_comparison.pdf", width=14)

types <- c("CO", "CM", "CY")
age_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>%
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender"), type="rq")) %>%
    summarize(Met=Metabolite,
              Mode=Mode,
              est=lmres["Age", 1],
              pval=lmres["Age", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

age_diff_matches <- age_diff %>%
    mutate(Mode2 = str_sub(Mode, 1, 3),
           Metabolite = str_replace(Met, " Results", "")) %>% 
    left_join(all_matches %>% rename(MetaboliteName = Metabolite),
              by=c("Mode2"="Mode", "Metabolite" = "Name"))

write_csv(age_diff_matches, path="../results/gotms_features/age_features.csv")

age_diff_matches <- age_diff_matches %>% slice(1:16) 


subject_data %>%
    left_join(age_diff_matches, .,
              by=c("Met"="Metabolite", "Mode"="Mode")) %>%
    filter(Type %in% types) %>%
    mutate(Type = factor(Type, levels=c("CY", "CM", "CO"))) %>% 
    ggplot(aes(x = Abundance, y = Type, fill=Type)) +
    geom_density_ridges(scale = 4, alpha=0.75, na.rm=FALSE) +
    facet_wrap(~ MetaboliteName + Metabolite, nrow=2, scales="free") +
    theme_bw()

dev.off()

###############################
## AD vs Controls (old)
###############################

pdf("../figs/gotms/AD_comparison.pdf", width=14)
types <- c("CO", "AD")
## Check if control is differen than AD
ad_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender", "Type"),
                    type="rq")) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["TypeCO", 1],
              pval=lmres["TypeCO", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)
    ## arrange(desc(abs(est)))

ad_diff_matches <- ad_diff %>%
    mutate(Mode2 = str_sub(Mode, 1, 3),
           Metabolite = str_replace(Met, " Results", "")) %>% 
    left_join(all_matches %>% rename(MetaboliteName = Metabolite),
              by=c("Mode2"="Mode", "Metabolite" = "Name"))

write_csv(ad_diff_matches,
          path="../results/gotms_features/ad_features.csv")

ad_diff_matches <- ad_diff_matches %>% slice(1:16) 

subject_data %>%
    left_join(ad_diff_matches, .,
              by=c("Met"="Metabolite", "Mode"="Mode")) %>%
    filter(Type %in% types) %>%
    ggplot(aes(x = Abundance, y = Type, fill=Type)) +
    geom_density_ridges(scale = 4, alpha=0.75, na.rm=FALSE) +
    facet_wrap(~ MetaboliteName + Metabolite, nrow=2, scales="free") +
    theme_bw()



dev.off()

###############################
## PD vs Controls (old)
###############################

pdf("../figs/gotms/PD_comparison.pdf", width=14)
types <- c("CO", "PD")
## Check if control is differen than PD
pd_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender", "Type"),
                    type="rq")) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["TypePD", 1],
              pval=lmres["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

pd_diff_matches <- pd_diff %>%
    mutate(Mode2 = str_sub(Mode, 1, 3),
           Metabolite = str_replace(Met, " Results", "")) %>% 
    left_join(all_matches %>% rename(MetaboliteName = Metabolite),
              by=c("Mode2"="Mode", "Metabolite" = "Name"))

write_csv(pd_diff_matches, path="../results/gotms_features/pd_features.csv")

pd_diff_matches <- pd_diff_matches %>% slice(1:16) 

subject_data %>%
    left_join(pd_diff_matches, .,
              by=c("Met"="Metabolite", "Mode"="Mode")) %>%
    filter(Type %in% types) %>%
    ggplot(aes(x = Abundance, y = Type, fill=Type)) +
    geom_density_ridges(scale = 4, alpha=0.75, na.rm=FALSE) +
    facet_wrap(~ MetaboliteName + Metabolite, nrow=2, scales="free") +
    theme_bw()
dev.off()


############################################

pdf("../figs/gotms/APOE_comparison.pdf", width=14)
## APOE
types <- c("AD", "CO")
## Check if control is differen than AD
apoe_diff <- subject_data %>%
    group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender", "Type"),
                    type="rq")) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["APOE44", 1],
              pval=lmres["APOE44", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

apoe_diff_matches <- apoe_diff %>%
    mutate(Mode2 = str_sub(Mode, 1, 3),
           Metabolite = str_replace(Met, " Results", "")) %>% 
    left_join(all_matches %>% rename(MetaboliteName = Metabolite),
              by=c("Mode2"="Mode", "Metabolite" = "Name"))

write_csv(apoe_diff_matches, path="../results/gotms_features/apoe_features.csv")

apoe_diff_matches <- apoe_diff_matches %>% slice(1:16) 

p <- subject_data %>%
    left_join(apoe_diff_matches, .,
              by=c("Met"="Metabolite", "Mode"="Mode")) %>%
    filter(Type %in% types) %>%
    ggplot(aes(x = Abundance, y = APOE, fill=APOE)) +
    geom_density_ridges(scale = 4, alpha=0.75, na.rm=FALSE) +
    facet_wrap(~ MetaboliteName + Metabolite, nrow=2, scales="free") +
    theme_bw()

meds <- subject_data %>%
    left_join(apoe_diff_matches, .,
              by=c("Met"="Metabolite", "Mode"="Mode")) %>%
    filter(Type %in% types) %>%
    group_by(APOE, Metabolite, MetaboliteName) %>%
    summarise(med = quantile(Abundance, 0.5, na.rm=TRUE)) %>%
    ungroup

p + geom_segment(data=meds,
                 aes(x=med, xend=med, y=0, yend=0.75, color=factor(APOE)),
                 size=1) 
dev.off()


###############################
## AD vs PD
###############################

pdf("../figs/gotms/AD_PD_comparison.pdf", width=14)
types <- c("AD", "PD")
## Check if control is differen than PD
ad_pd_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres= lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              Mode = Mode,
              est=summary(lmres)$coefficients["TypePD", 1],
              pval=summary(lmres)$coefficients["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- ad_pd_diff$Met[1:8]
subject_data %>%
    left_join(., age_diff,
              by=c("Metabolite"="Met", "Mode"="Mode")) %>%
    filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type),
               fill=Type)) +  geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite + MetaboliteName, nrow=2)
dev.off()




apoe_diff <- apoe_diff %>% mutate(sgn = sign(est))
tmp <- apoe_diff %>% dplyr::select(Met, sgn, est)
subject_data %>% 
    filter(Type %in% types, Metabolite %in% apoe_diff$Met) %>%
    left_join(.,
              tmp,
              by=c("Metabolite"="Met")) %>% 
    group_by(Metabolite) %>%    
    ## mutate(MetAbund = Abundance / sd(Abundance)) %>%
    mutate(MetAbund = Abundance) %>%
    ungroup %>%
    group_by(Id) %>%
    summarise(sum = sum(MetAbund * sgn, na.rm=TRUE), APOE=APOE[1]) %>%
    ggplot(aes(x = sum, y = factor(APOE),
               fill=factor(APOE))) +  geom_density_ridges(scale=4) +
    theme_bw()
    



mets <- apoe_diff$Met[1:16]
meds <- subject_data %>% filter(Metabolite %in% mets) %>% group_by(APOE, Metabolite) %>% summarise(med = quantile(Abundance, 0.25, na.rm=TRUE)) %>% ungroup
p <- subject_data %>% filter(Metabolite %in% mets, Age > 60) %>%
    ggplot(aes(x = Abundance, y = factor(APOE),
               fill=factor(APOE))) +  geom_density_ridges(scale = 4) + facet_wrap(~ Metabolite, nrow=4)
p + geom_segment(data=meds, aes(x=med, xend=med, y=0, yend=0.75, color=APOE), size=1)


#########################################################

Y <- wide_data %>%
    dplyr::select(-one_of("Type", "Type2", "Gender", "Age", "APOE", "Batch",
                          "Data File", "Index", "GBAStatus", "Id",
                          "GBA_T369M", "cognitive_status")) %>%
    as.matrix()
Y[Y==0] <- NA
dim(Y)

Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
Y <- t(Yt) %>% as.matrix

## APOE 33 is the baseline
X <- wide_data %>% 
    dplyr::select(one_of("Gender", "Age", "Type2", "APOE")) %>%
    mutate(Type2 = factor(Type2, c("C", "AD", "PD"))) %>%
    mutate(APOE = factor(APOE, c("33", "22", "23", "24", "34", "44"))) %>%
    model_matrix(~ Age + Type2 + Gender + APOE - 1) %>% as.matrix
colnames(X)
p <- nrow(Y)

## This overfits
betaHat_MLE <- solve(t(X) %*% X) %*% t(X) %*% (Y %*% D)

cvmfit <- cv.glmnet(X, as.matrix(Y %*% D), family = "mgaussian",
                    type.measure="mae", nfolds=10)

glmres <- glmres <- glmnet(X, Y %*% D, family="mgaussian", alpha=0, lambda=cvmfit$lambda.min)

betaHat <- do.call(cbind, glmres$beta)

V <- svd(betaHat)$v
dim(V)

U <- svd(betaHat)$u
dim(U)

tib <- tibble(YV = Y %*% betaHat[1, ],
              Type = factor(wide_data$Type, levels=c("CY", "CM", "CO")),
              Gender = factor(wide_data$Gender),
              APOE = wide_data$APOE)

tib %>% filter(Type %in% c("CY", "CM", "CO")) %>% ggplot() +
    geom_density_ridges(aes(y = Type, x=YV, fill=Type))

tib %>% filter( %in% c("CY", "CM", "CO")) %>% ggplot() +
    geom_density_ridges(aes(y = Type, x=YV, fill=Type))

tib <- tibble(YV = Y %*% svd(betaHat[6:10, ])$v[, 1],
              Type = factor(wide_data$Type), 
              Gender = factor(wide_data$Gender),
              APOE = wide_data$APOE)

tib %>%  filter(Type %in% c("CO", "AD", "PD")) %>% ggplot() +
    geom_density_ridges(aes(y = APOE, x=YV, fill=APOE))












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
lm_data <- tmp_data %>%
    group_by(Metabolite) %>%
    do(res=lm(Abundance ~ Gender + Type, data=.)$residuals)

tmp_data$residuals <- unlist(lm_data$res)
tst <- tmp_data %>%
    dplyr::select(-c(Abundance, Raw, RawScaled,
                     Trend, Mode, Code, RunIndex)) %>%
    spread(key=Metabolite, value=residuals)

glimpse(tst[, 10:30])
glimpse(tst[10, ])

## Principle component analysis for control subjects at different ages
comparison_groups <- c("CY", "CM", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0

sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Gender", "Age", "APOE", "Batch",
                   "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)

type_vec <- sub_mat$Type
pcres <- sub_mat %>%  select(-one_of("Id", "Type")) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=2, classvec=type_vec)

resid <- sub_mat %>%  select(-one_of("Id", "Type"))
colors = rainbow(length(unique(type_vec)))
names(colors) = unique(type_vec)
ecb = function(x,y){
    plot(x, t='n')
    text(x, labels=type_vec, col=colors[type_vec])
}
tsne(resid, epoch_callback = ecb)

## M v F
comparison_groups <- c("CY", "CM", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0

sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Type", "Age", "APOE",
                   "Batch", "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)
type_vec <- sub_mat$Gender
pcres <- sub_mat %>%  select(-one_of("Id", "Gender")) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=5, classvec=type_vec)

## Principle component analysis for AD and old controls

comparison_groups <- c("PD", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0
sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Gender", "Age", "APOE",
                   "Batch", "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)
type_vec <- sub_mat$Type
pcres <- sub_mat %>%  select(-one_of("Id", "Type")) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=2, classvec=type_vec)

comparison_groups <- c("AD", "CO")
tmp <- subject_data
tmp$Abundance[is.na(tmp$Abundance)] <- 0
sub_mat  <- tmp %>% filter(Type %in% comparison_groups) %>%
    select(-one_of("Code", "Gender", "Age", "APOE",
                   "Batch", "Index", "Raw", "Trend", "RunIndex", "Mode")) %>%
    spread(key=Metabolite, value=Abundance)

type_vec <- sub_mat$Type
pcres <- sub_mat %>%  select(-one_of("Id", "Type")) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=5, classvec=type_vec)

resid <- sub_mat %>%  select(-one_of("Id", "Type"))
colors = rainbow(length(unique(type_vec)))
names(colors) = unique(type_vec)
ecb = function(x,y){
    plot(x, t='n')
    text(x, labels=type_vec, col=colors[type_vec])
}
tsne(resid, epoch_callback = ecb)
## Test



CoefficientMat <- longQC %>%
    group_by(Metabolite) %>%
    summarize(coef=lm(log(Abundance) ~ Index)$coefficients[2],
              pval=(summary(lm(log(Abundance) ~ Index)))$coefficients[2, 4])

subject_data[subject_data==0] <- NA
tmp <- sapply(colnames(subject_data)[9:90], function(metabolite) {
    t.test(log(filter(subject_data, Type=="PD")[[metabolite]]),
           log(filter(subject_data, Type=="CO")[[metabolite]]))$p.value
})
names(tmp) <- colnames(subject_data)[9:90]
sort(tmp, decreasing=TRUE)


