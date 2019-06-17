library(tidyverse)
library(ggExtra)
library(magrittr)
library(made4)
library(ggridges)
library(gbm3)
library(huge)
library(patchwork)
library(scales)

source("utility.R")

processed_files <- dir(path="../preprocessed_data",
                       pattern="^preprocessed_untargeted_data*")

## Most recent file
load(max(processed_files[grep("-20+", processed_files)]))

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
    filter(Mode == "pos") %>% ungroup %>%
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
subject_data %>% filter(Mode=="pos") %>%
    group_by(RunIndex) %>%
    summarize(tot=sum(Raw, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) +
    geom_line() + ggtitle("Positive") + geom_smooth(method="lm", se=FALSE)

subject_data %>% filter(Mode=="neg") %>%
    group_by(RunIndex) %>%
    summarize(tot=sum(Abundance, na.rm=TRUE)) %>%
    ggplot(aes(x=RunIndex, y=tot)) + geom_line() + ggtitle("Negative") +
    geom_smooth(method="lm", se=FALSE)


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

metid <- 1
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
    unite("Metabolite", c("Metabolite", "Mode")) %>%
    dplyr::select(1:11, 31:37) %>% 
    dplyr::select(-one_of("RunIndex", "Raw", "RawScaled", "Trend",
                          "Name", "Index")) %>%
    spread(key=Metabolite, value=Abundance)

wide_data$Type
dim(wide_data)

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




fit_lm <- function(df, predictors, type="lm") {
    form <- as.formula(paste("Abundance", paste(predictors, collapse = " + "), 
                             sep = " ~ "))

    if(sum(is.na(df$Abundance)) > 50)
        return(NULL)
    else {
        if(type=="lm") {
            return(summary(lm(form, data=df))$coefficients)
        }
        else if(type == "rq") {

            rq_fit <- rq(form, data=df, tau=0.5)
            coefs <- summary(rq_fit, se="boot")$coefficients
            colnames(coefs)[1] <- "Estimate" ## to match lm
            coefs

        }
    
        
    }
}


###############################
## Metabolites that show effects with age
###############################

types <- c("CY", "CM", "CO")
## Null model test
pdf("../figs/untargeted/ridx_comparison.pdf", width=14)
ridx_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(.,
                    predictors = c("RunIndex", "Age", "APOE", "Gender"),
                    type="rq")) %>%
    summarize(Met=Metabolite,
              Mode=Mode,
              est=lmres["RunIndex", 1],
              pval=lmres["RunIndex", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- ridx_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(cut(RunIndex, 3)),
               fill=factor(cut(RunIndex, 3)))) +
    geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=4)
dev.off()

pdf("../figs/untargeted/Age_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")
age_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender"), type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Met=Metabolite,
              Mode=Mode,
              est=lmres["Age", "Estimate"],
              pval=lmres["Age", "Pr(>|t|)"]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

mets <- age_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=types),
               fill=Type)) +  geom_density_ridges(scale = 4) + theme_bw() +
    facet_wrap(~ Metabolite, nrow=4, scales="free")

dev.off()

metid <- 200
intercepts <- batch_positions %>% ungroup() %>%
    filter(Mode == batch_positions$Mode[5]) %>%
    dplyr::select(mx) %>% unlist %>% as.numeric

subject_data %>% filter(Metabolite == "1061 Results") %>%
    filter(Type %in% types) %>%
    ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) +
    geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=Raw)) +
    geom_point(aes(x=RunIndex, y=Raw, shape=Type)) +
    geom_line(aes(x=RunIndex, y=Abundance)) +
    geom_point(aes(x=RunIndex, y=Abundance, shape=Type)) 

###############################
## Metabolites that show effects with Gender
###############################

pdf("../figs/untargeted/Gender_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")

gender_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "APOE", "Gender"), type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Met=Metabolite,
              Mode=Mode,
              est=lmres["GenderM", "Estimate"],
              pval=lmres["GenderM", "Pr(>|t|)"]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

mets <- gender_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = Gender,
               fill=Gender)) +  geom_density_ridges(scale = 4) + theme_bw() +
    facet_wrap(~ Metabolite, nrow=4, scales="free")

dev.off()

###############################
## AD vs Controls (old)
###############################

pdf("../figs/untargeted/AD_comparison.pdf", width=14)
types <- c("CO", "AD")
## Check if control is differen than AD
ad_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(., predictors = c("Age", "Type", "Gender", "APOE"),
                    type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["TypeCO", "Estimate"],
              pval=lmres["TypeCO", "Pr(>|t|)"]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)


mets <- ad_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type),
               fill=Type)) +  geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=4, scales="free")
dev.off()

###############################
## PD vs Controls (old)
###############################

pdf("../figs/untargeted/PD_comparison.pdf", width=14)
types <- c("CO", "PD")
## Check if control is differen than PD
pd_diff <- subject_data %>% group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(.,
                    predictors = c("Age", "Type", "Gender", "APOE"),
                    type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["TypePD", "Estimate"],
              pval=lmres["TypePD", "Pr(>|t|)"]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)


mets <- pd_diff$Met[1:16]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type),
               fill=Type)) +  geom_density_ridges(scale = 4) +
    facet_wrap(~ Metabolite, nrow=4, scales="free")
dev.off()

############################################

pdf("../figs/untargeted/APOE_comparison.pdf", width=14)
## APOE
types <- c("AD", "CO")
## Check if control is differen than AD
apoe_diff <- subject_data %>%
    group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(.,
                    predictors = c("Age", "Gender", "Type", "APOE"),
                    type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["APOE44", "Estimate"],
              pval=lmres["APOE44", "Pr(>|t|)"]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)


mets <- apoe_diff$Met[1:16]
meds <- subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    group_by(APOE, Metabolite) %>%
    summarise(med = quantile(Abundance, 0.5, na.rm=TRUE)) %>%
    ungroup 
p <- subject_data %>%
    filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(APOE),
               fill=factor(APOE))) +  geom_density_ridges(scale=4) +
    facet_wrap(~ Metabolite, nrow=4, scales="free") + theme_bw()
p + geom_segment(data=meds,
                 aes(x=med, xend=med, y=0, yend=0.75, color=factor(APOE)),
                 size=1) 
dev.off()


pdf("../figs/untargeted/PDCO_APOE_comparison.pdf", width=14)
## APOE
types <- c("PD", "CO")
## Check if control is differen than AD
pd_apoe_diff <- subject_data %>%
    group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>% 
    do(lmres=fit_lm(.,
                    predictors = c("Age", "Gender", "Type", "APOE"),
                    type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["APOE34", "Estimate"],
              pval=lmres["APOE34", "Pr(>|t|)"]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

mets <- pd_apoe_diff$Met[1:16]
meds <- subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    group_by(APOE, Metabolite) %>%
    summarise(med = quantile(Abundance, 0.5, na.rm=TRUE)) %>%
    ungroup 
p <- subject_data %>%
    filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(APOE),
               fill=factor(APOE))) +  geom_density_ridges(scale=4) +
    facet_wrap(~ Metabolite, nrow=4, scales="free") + theme_bw()
p + geom_segment(data=meds,
                 aes(x=med, xend=med, y=0, yend=0.75, color=factor(APOE)),
                 size=1) 
dev.off()

intersect(apoe_diff$Met, pd_apoe_diff$Met)

pdf("../figs/untargeted/PD_APOE_comparison.pdf", width=14)
## PD only
types <- c("PD")
## Check if control is differen than AD
cog_diff <- subject_data %>%
    group_by(Metabolite, Mode) %>%
    filter(Type %in% types) %>%
    mutate(cogstat = ifelse(cognitive_status == "No cognitive impairment",
                            "N",
                            "Y")) %>% 
    do(lmres=fit_lm(.,
                    predictors = c("Age", "cogstat"),
                    type="rq")) %>%
    filter(!is.null(lmres)) %>% 
    summarize(Met=Metabolite,
              Mode = Mode,
              est=lmres["cogstatN", "Estimate"],
              pval=lmres["cogstatN", "Pr(>|t|)"]) %>%
    filter(pval < 0.05) %>%
    arrange(pval)

mets <- apoe_diff$Met[1:16]
meds <- subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    group_by(APOE, Metabolite) %>%
    summarise(med = quantile(Abundance, 0.5, na.rm=TRUE)) %>%
    ungroup 
p <- subject_data %>%
    filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(APOE),
               fill=factor(APOE))) +  geom_density_ridges(scale=4) +
    facet_wrap(~ Metabolite, nrow=4, scales="free") + theme_bw()
p + geom_segment(data=meds,
                 aes(x=med, xend=med, y=0, yend=0.75, color=factor(APOE)),
                 size=1) 
dev.off()
