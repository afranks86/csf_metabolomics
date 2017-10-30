library(tidyverse)
library(magrittr)
library(made4)
library(ggjoy)
library(gbm3)
library(huge)

load("normalized_csf_data.RData")

## Look at average

subject_data %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Raw)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Positive") + geom_smooth(method="lm", se=FALSE)

QC_long %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Abundance)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Positive") + geom_smooth(method="lm", se=FALSE)


## Look at Raw
subject_data %>% filter(Mode=="neg") %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Raw)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Negative") + geom_smooth(method="lm", se=FALSE)

## check that detrended data looks sensible (for a random sample of metabolites)
subject_data %>% filter(Metabolite %in% sample(unique(subject_data$Metabolite), 4)) %>% ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) + geom_line() + geom_point() + theme(legend.position="none")

## Which metabolites had the worst drift?
subject_data %>% group_by(Metabolite) %>% summarise(var = (max(Trend, na.rm=TRUE) - min(Trend, na.rm=TRUE))/sd(Raw, na.rm=TRUE)) %>% arrange(desc(var))

## Look at some trends fits for a few metabolites
intercepts <- (subject_data %>% group_by(Batch) %>% summarise(mx = max(RunIndex)))$mx
subject_data %>% filter(Metabolite %in% c("Creatinine", "Fructose", "Xanthine")) %>% ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) + geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=Raw)) + geom_vline(xintercept = intercepts)

## Plot intensity versus log variance
subject_data  %>% group_by(Metabolite) %>% summarise(med=median(Raw), var=log(var(Abundance))) %>% ggplot(aes(x=med, y=var)) + geom_point() + geom_smooth(method="lm")
subject_data  %>% group_by(Metabolite) %>% summarise(med=median(Raw), var=log(var(Abundance))) %>% cor.test(.$med, .$var, data=.)


## and residuals..
subject_data %>% filter(Metabolite %in% c("Creatinine", "Tryptophan", "Methionine")) %>% ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) + geom_line() + geom_point()

## Do QC variances correlate with raw sample variances?
met_var_qc <- QC_long %>% group_by(Metabolite) %>% summarize(var = var(Abundance))
met_var_subj <-  subject_data %>% group_by(Metabolite) %>% summarize(var = var(Abundance, na.rm=TRUE))
met_var <- left_join(met_var_qc, met_var_subj, by="Metabolite")

met_var %>% ggplot(aes(x=log(var.x), y=log(var.y))) + geom_point() + geom_smooth(method="lm")



###############################
## Metabolites that show effects with age
###############################

pdf("../figs/Age_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")
age_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["Age", 1],
              pval=summary(lmres)$coefficients["Age", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- age_diff$Met[1:8]
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=types),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=2)

subject_data %>% group_by(Metabolite) %>%
    do(lmres=lm(Abundance ~ Age, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients[2, 1],
              pval=summary(lmres)$coefficients[2, 4]) %>%
    filter(pval < 0.05) %>%
    arrange(est)
dev.off()

###############################
## Metabolites that show effects with Gender
###############################

pdf("../figs/Gender_comparison.pdf", width=14)
types <- c("CO", "CM", "CY")
gender_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Gender, data=.)) %>%
    summarize(Met=Metabolite,
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

pdf("../figs/AD_comparison.pdf", width=14)
types <- c("CO", "AD")
## Check if control is differen than AD
ad_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
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

pdf("../figs/PD_comparison.pdf", width=14)
types <- c("CO", "PD")
## Check if control is differen than PD
pd_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
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

pdf("../figs/AD_PD_comparison.pdf", width=14)
types <- c("AD", "PD")
## Check if control is differen than PD
ad_pd_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
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

pdf("../figs/APOE_comparison.pdf", width=14)
## APOE
types <- c("AD", "CO")
## Check if control is differen than AD
apoe_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% types) %>% 
    do(lmres=lm(Abundance ~ Age + Gender + factor(APOE), data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["factor(APOE)44", 1],
              pval=summary(lmres)$coefficients["factor(APOE)44", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

mets <- apoe_diff$Met
subject_data %>% filter(Metabolite %in% mets, Type %in% types) %>%
    ggplot(aes(x = Abundance, y = factor(APOE),
               fill=factor(APOE))) +  geom_joy(scale = 4) + theme_joy() + facet_wrap(~ Metabolite, nrow=2)
dev.off()


############################################
tmp_data <- subject_data
tmp_data$Abundance[is.na(tmp_data$Abundance)] <- 0
lm_data <- tmp_data %>% group_by(Metabolite) %>% do(res=lm(Abundance ~ Age + Gender + Type, data=.)$residuals)
tmp_data$residuals <- unlist(lm_data$res)
tst <- tmp_data %>% spread(key=Metabolite, value=residuals, vars(-Abundaance, Raw)




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
