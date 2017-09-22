library(tidyverse)
library(magrittr)
library(made4)
library(ggjoy)
library(gbm3)

positive_mode_data <- read_csv("~/course/ND_Metabolomics/data/csf_positive.csv")
negative_mode_data <- read_csv("~/course/ND_Metabolomics/data/csf_negative.csv")
positive_mode_data$RunIndex <- negative_mode_data$RunIndex <- 1:nrow(positive_mode_data)

################################
## Normalize mode data
################################

## Get positive QC data
positive_mode_data$Mode <- ifelse(grepl("QC", positive_mode_data$Code), "QC_pos", "pos")
negative_mode_data$Mode <- ifelse(grepl("QC", negative_mode_data$Code), "QC_neg", "neg")

QC_pos <- positive_mode_data %>% filter(Mode == "QC_pos")

pos_long <- positive_mode_data %>% gather(key=Metabolite, value=Abundance, -one_of("Code", "RunIndex", "Mode")) 
neg_long <- negative_mode_data %>% gather(key=Metabolite, value=Abundance, -one_of("Code", "RunIndex", "Mode"))
pos_long <- pos_long %>% mutate_at("Abundance", funs(log(.)))
neg_long <- neg_long %>% mutate_at("Abundance", funs(log(.)))

## Look at average
pdf("~/Desktop/positive.pdf")
pos_long %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Abundance)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Positive") + geom_smooth(method="lm", se=FALSE)
dev.off()

pdf("~/Desktop/negative.pdf")
neg_long %>% group_by(RunIndex, Mode) %>% summarize(tot=sum(Abundance)) %>%  ggplot(aes(x=RunIndex, y=tot, colour=Mode)) + geom_line() + ggtitle("Negative") + geom_smooth(method="lm", se=FALSE)
dev.off()

#####################################33

csf_long <- rbind(pos_long, neg_long)

QC_indices <- grepl("QC", csf_long$Code)
QC_long <- csf_long %>% filter(QC_indices)
nd_data <- csf_long %>% filter(!QC_indices)
nd_data$Index <- nd_data$Code %>% gsub(pattern="(NEG-CSF-06202017-)|(CSF-06122017-)", replacement="") %>% as.numeric

nd_data$Abundance[nd_data$Abundance==-Inf] <- NA

fit_boosted_model <- function(df, ntrees=10000, cv.folds=10){
    print(df$Metabolite[1])
    not_na_indices <- which(!is.na(df$Abundance))
    df.omitted <- df %>% na.omit(Abundance)

    boost_fit_cv <- gbm(Abundance ~ RunIndex, data=df.omitted,
                        distribution="laplace", n.trees=ntrees, cv.folds=10, verbose=FALSE)
    opt_iters <- gbm.perf(boost_fit_cv, method="cv")
    print(opt_iters)

    best_fit <- predict(boost_fit_cv, df.omitted, opt_iters)


    
    abund <- df$Abundance
    bf <- rep(NA, length(abund))
    
    abund[not_na_indices] <- df.omitted$Abundance - best_fit
    bf[not_na_indices] <- best_fit
    
    data.frame(Raw = df$Abundance, Abundance = abund, Trend = bf, RunIndex = df$RunIndex, Code=df$Code, Index=df$Index, Mode=df$Mode)

}

detrended <- nd_data %>% group_by(Metabolite, Mode) %>% do(fit_boosted_model(.))
save(detrended, file="normalized_csf_data.RData")

## some plots
detrended %>% filter(Metabolite %in% sample(unique(detrended$Metabolite), 4)) %>% ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) + geom_line() + geom_point() + theme(legend.position="none")

detrended %>% filter(Metabolite %in% sample(unique(detrended$Metabolite), 4)) %>% ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) + geom_line() + geom_point() + theme(legend.position="none")

detrended %>% filter(Metabolite %in% c("Creatinine", "Creatine", "Methionine")) %>% ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) + geom_line() + geom_point()
detrended %>% filter(Metabolite %in% c("Creatinine", "Creatine", "Methionine")) %>% ggplot(aes(x=RunIndex, y=Trend, colour=Metabolite)) + geom_line() + geom_point() + geom_line(aes(x=RunIndex, y=Raw))

###############################
## Join with Tracking data
###############################

tracking_data <- read_csv("~/course/ND_Metabolomics/data/NDTracking.csv")

subject_data <- inner_join(tracking_data, detrended, by="Index")

sj_wide <- subject_data %>% spread(Metabolite, Abundance)

subject_data %>% group_by(Metabolite) %>%
    do(lmres=lm(Abundance ~ Age, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients[2, 1],
              pval=summary(lmres)$coefficients[2, 4]) %>%
    filter(pval < 0.05) %>%
    arrange(est)

## Not right, need to look at residuals
cor.test(sj_wide[["Kynurenine"]], sj_wide[["Anthranilic acid"]])
cor.test(sj_wide[["Kynurenine"]], sj_wide[["Decanoylcarnitine"]])
cor.test(sj_wide[["Kynurenine"]], sj_wide[["Serotonin"]])


met <- "Kynurenine"
met <- "Serotonin"
met <- "Caffeine"
met <- "Methionine"
met <- "Creatinine"
met <- "DOPA"
met <- "Sorbitol"
met <- "Mannose"
x <- c("CY", "CM", "CO")

subject_data %>% group_by(Metabolite) %>% summarise(var = (max(Trend, na.rm=TRUE) - min(Trend, na.rm=TRUE))/sd(Raw, na.rm=TRUE)) %>% arrange(desc(var))

subject_data %>% filter(Metabolite == met, Type %in% x) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=x),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy()

subject_data %>% filter(Metabolite == met, Type %in% x) %>%
    ggplot(aes(x = Raw, y = factor(Type, levels=x),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy()

x <- c("PD", "CO", "AD")
met <- "Kynurenine"
met <- "Serotonin"
met <- "Caffeine"
met <- "Methionine"
met <- "Decanoylcarnitine"
met <- "Anthranilic acid"
met <-  "HIAA"
met <-  "Leucine"
met <-  "Xanthosine"
met <-  "Hydroxyproline"

met <-  "Sarcosine"
met <-  "Homocysteine"
met <-  "Indole-3-acetic acid"

met <-  "Nonadecanoic acid"
met <-  "Methylhistamine"



## Type
subject_data %>% filter(Metabolite == met, Type %in% x) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=x),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy() + ggtitle(met)

## Gender
subject_data %>% filter(Metabolite == met, Type %in% c("CO", "CM", "CY")) %>%
    ggplot(aes(x = Abundance, y = Gender,
               fill=Gender)) +  geom_joy(scale = 4) + theme_joy() + ggtitle(met)c


gender_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("CO", "CM", "CY")) %>% 
    do(lmres=lm(Abundance ~ Age + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["GenderM", 1],
              pval=summary(lmres)$coefficients["GenderM", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

age_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("CO", "CM", "CY")) %>% 
    do(lmres=lm(Abundance ~ Age, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["Age", 1],
              pval=summary(lmres)$coefficients["Age", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

## Check if control is differen than AD
ad_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("CO", "AD")) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["TypeCO", 1],
              pval=summary(lmres)$coefficients["TypeCO", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

## Check if control is differen than PD
pd_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("CO", "PD")) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["TypePD", 1],
              pval=summary(lmres)$coefficients["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))


pd_ad_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("PD", "AD")) %>% 
    do(lmres=lm(Abundance ~ Age + Type + Gender, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients["TypePD", 1],
              pval=summary(lmres)$coefficients["TypePD", 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

## Check differences in Gender
gender_diff <- subject_data %>% group_by(Metabolite) %>%
    do(lmres=lm(Abundance ~ Gender, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients[2, 1],
              pval=summary(lmres)$coefficients[2, 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

## Check that global abundances are different
subject_data %>% group_by(Type) %>% summarize(avg = mean(Abundance))


## Check if conrol different than PD
pd_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("CO", "PD")) %>% 
    do(lmres=lm(Abundance ~ Age + Type, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients[2, 1],
              pval=summary(lmres)$coefficients[2, 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))

    


## Principle component analysis for control subjects at different ages
subject_data$Abundance[subject_data$Abundance == -Inf] <- -3


comparison_groups <- c("CY", "CM", "CO")

sub_mat  <- subject_data %>% filter(Type %in% comparison_groups) %>%
    filter(Mode=="Positive") %>% 
    spread(key=Metabolite, value=Abundance)
type_vec <- sub_mat$Type
pcres <- sub_mat %>% select(11:ncol(.)) %>% prcomp
plotarrays(pcres$x, axis1=1, axis2=2, classvec=type_vec)

## M v F
type_vec <- sub_mat$Gender
plotarrays(pcres$x, axis1=3, axis2=2, classvec=type_vec)

## Principle component analysis for AD and old controls
comparison_groups <- c("AD", "CO")
type_vec <- subject_data %>% filter(Type %in% comparison_groups) %>%
    na_if(0) %>% na.omit %>% select("Type")
pcres <- subject_data %>% filter(Type %in% comparison_groups) %>% select(metabolite_indices) %>% na_if(0) %>% na.omit %>% log %>% prcomp()
plotarrays(pcres$x, axis1=2, axis2=4, classvec=type_vec$Type)

## Principle component analysis for AD and old controls
comparison_groups <- c("PD", "CO")
type_vec <- subject_data %>% filter(Type %in% comparison_groups) %>%
    na_if(0) %>% select("Type")
pcres <- subject_data %>% filter(Type %in% comparison_groups) %>% select(metabolite_indices) %>%mutate_all(funs(ifelse(.==0, 1, .))) %>% log %>% prcomp()
plotarrays(pcres$x, axis1=2, axis2=3, classvec=type_vec$Type)



CoefficientMat <- longQC %>% group_by(Metabolite) %>% summarize(coef=lm(log(Abundance) ~ Index)$coefficients[2], pval=(summary(lm(log(Abundance) ~ Index)))$coefficients[2, 4])

subject_data[subject_data==0] <- NA
tmp <- sapply(colnames(subject_data)[9:90], function(metabolite) {
    t.test(log(filter(subject_data, Type=="PD")[[metabolite]]),
           log(filter(subject_data, Type=="CO")[[metabolite]]))$p.value
})
names(tmp) <- colnames(subject_data)[9:90]
sort(tmp, decreasing=TRUE)
