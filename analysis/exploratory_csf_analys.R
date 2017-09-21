library(tidyverse)
library(magrittr)
library(made4)
library(ggjoy)
library(gbm)

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


fit_boosted_model <- function(df, ntrees=100000, cv.folds=10){

    boost_fit_cv <- gbm(Abundance ~ RunIndex + as.factor(Mode), data=df,
                         distribution="laplace", n.trees=ntrees, cv.folds=10, verbose=FALSE)
    opt_iters <- gbm.perf(boost_fit_cv, method="cv")
    print(opt_iters)
    boost_fit <- gbm(Abundance ~ RunIndex, data=df, distribution="laplace", n.trees=opt_iters)
    df$Abundance - boost_fit$fit
}

detrended <- nd_data %>% group_by(Metabolite, Mode) %>% do(data.frame(Abundance_Detrend=fit_boosted_model(.), RunIndex=.$RunIndex, Code=.$Code))

## some plots
detrended %>% ggplot(aes(x=RunIndex, y=Abundance_Detrend, colour=Metabolite)) + geom_line() + geom_point()
csf_long %>% filter(Metabolite %in% c("Creatinine", "Creatine", "Methionine")) %>% ggplot(aes(x=RunIndex, y=Abundance, colour=Metabolite)) + geom_line() + geom_point()

nd_data <- inner_join(nd_data, detrended, by="Code")




















met <- "Serotonin"
met <- "Caffeine"
met <- "Methionine"
met <- "Decanoylcarnitine"
met <- "Anthranilic acid"
met <-  "HIAA"
met <-  "Leucine"
met <-  "Homocysteine"
met <-  "4-Aminobutyric acid"
met <- "Agmatine"
met <- "2-Chloro-4,6-diamino-1,3,5-triazine"
met <- "Hydrocortisone"
met <- "Nicotinamide"
met <- "Creatinine"

cdf <- pos_long %>% filter(Metabolite %in% c(met)) %>% filter(Mode != "QC_pos")

widths_vec <- seq(2, 200, by=2)
err_vec <- c()
for(width in widths_vec) {
    print(width)
    err <- 0
    for(i  in 1:nrow(cdf)) {
        ridx <- cdf$RunIndex[-i]
        abund <- cdf$Abundance[-i]
        res <- sapply(1:(length(ridx) - width), function(i) median(abund[i:(i+width)]))
        res <- c(rep(res[1], width/2), res, rep(res[length(res)], width/2))
        if( i == nrow(cdf)) {
            err <- err + abs(cdf$Abundance[i] - res[length(ridx)])
        } else {
            err <- err + abs(cdf$Abundance[i] - res[match(TRUE, i < ridx)])
        }

    }
    err <- err / nrow(cdf)
    err_vec <- c(err_vec, err)
}

width <- widths_vec[which.min(err_vec)]
print(width)
plot(cdf$RunIndex, cdf$Abundance, type="l", col="black")
res <- sapply(1:(length(cdf$RunIndex) - width), function(i) median(cdf$Abundance[i:(i+width)]))
lines((width/2):(length(cdf$RunIndex)-(width/2)-1), res, type="l", lwd=3, col="red")

rf <- randomForest(cdf$Abundance ~ cdf$RunIndex, maxnodes=12, keep.inbag=TRUE)

boosted_trees <- gbm(Abundance ~ RunIndex + as.factor(Mode), data=cdf, distribution="laplace", n.trees=100000, cv.folds=10, verbose=FALSE)
opt_iters <- gbm.perf(boosted_trees, method="cv")

boosted_trees <- gbm(Abundance ~ RunIndex, data=cdf, distribution="laplace", n.trees=opt_iters)
plot(cdf$RunIndex, cdf$Abundance, type="l", col="black")
lines(cdf$RunIndex, boosted_trees$fit, col="red", lwd=3)

res1 <- randomForest(cdf$Abundance ~ cdf$RunIndex, maxnodes = 5)
res2 <- randomForest(cdf$Abundance ~ cdf$RunIndex, maxnodes = 3)
res3 <- randomForest(cdf$Abundance ~ cdf$RunIndex, maxnodes = 1)

for(met in unique(pos_long$Metabolite)) {
    print(met)
    cdf <- pos_long %>% filter(Metabolite %in% c(met))
    print(which.min(sapply(1:30, function(i) randomForest(cdf$Abundance ~ cdf$RunIndex, maxnodes = i)$mse[500])))
    print("---------")
}

plot(cdf$RunIndex, cdf$Abundance, type="l", col="black")
lines(cdf$RunIndex, rf$predicted, type="l", lwd=2, col="green")
lines(cdf$RunIndex, res1$predicted, type="l", lwd=2, col="red")
lines(cdf$RunIndex, res2$predicted, type="l", lwd=2, col="blue")
lines(cdf$RunIndex, res3$predicted, type="l", lwd=2, col="orange")

plot(res$mse, col="green", type="l", ylim=c(0, max(res$mse)))
lines(res2$mse, col="blue", type="l")
lines(res1$mse, col="red", type="l")
lines(res3$mse, col="orange", type="l")

pos_dat[pos_dat==0] <- 1
for(met in colnames(pos_dat)[-c(1:20, 84)]) {
    pdf(sprintf("../figs/%s.pdf", met))
    plot(pos_dat$RunIndex, log(pos_dat[[met]]), type="l", main=met, xlab="Index", ylab="Abundance")
    coefs <- lm(log(pos_dat[[met]]) ~ pos_dat$RunIndex)$coefficients
    abline(a=coefs[1], b=coefs[2], lty=2, col="black")
    lines(QC_pos$RunIndex, log(QC_pos[[met]]), col="blue")
    coefs <- lm(log(QC_pos[[met]]) ~ QC_pos$RunIndex)$coefficients
    abline(a=coefs[1], b=coefs[2], lty=2, col="blue")
    dev.off()

    pdf(sprintf("../figs/%s_bootstrap.pdf", met))
    rm <- runmed(log(pos_dat[[met]]), k=51)
    ## bootstrap
    
    boot_diff <- sapply(1:1000, function(i) {
        rm_boot <- runmed(sample(log(pos_dat[[met]]), nrow(pos_dat), replace=TRUE), k=51)
        rm_boot[170] - rm_boot[20]
    })
   
    hist(boot_diff, breaks=50, main= mean((rm[170]-rm[20]) > boot_diff))
    abline(v=(rm[170]-rm[20]), lwd=2, col="red")

    dev.off()
}



## Look at QC data by metabolite
ggplot(data=QC_long, aes(x=RunIndex, y=Abundance, colour=Metabolite)) +
    geom_line() + coord_trans(y="log10") + theme(legend.position="none")

## Plot variance in QC vs median by metabolite
QC_long %>% group_by(Metabolite, Mode) %>%
    summarize(med=median(log(Abundance), na.rm=TRUE), var=log(var(log(Abundance)))) %>%
    ggplot(data=., aes(x=med, y=var, color=Mode)) + geom_point()

## Regression coefficients
QC_long <- QC_pos %>% gather(key=Metabolite, value=Abundance, -one_of("Code", "RunIndex"))

out <- QC_long %>% group_by(Metabolite) %>% summarize(coef=lm(log(Abundance) ~ RunIndex)$coefficients[2], pval=(summary(lm(log(Abundance) ~ RunIndex)))$coefficients[2, 4])

posout <- pos_long %>% group_by(Metabolite) %>% summarize(coef=lm(Abundance ~ RunIndex)$coefficients[2], pval=(summary(lm(Abundance ~ RunIndex)))$coefficients[2, 4])

tracking_data <- read_csv("~/course/ND_Metabolomics/data/NDTracking.csv")

subject_data <- inner_join(tracking_data, nd_data, by="Index")


metabolite_columns <- -one_of("Code", "RunIndex", vars=colnames(QC))
normalizing_constants <-
    sapply(subject_data$RunIndex,
           function(x) {
               min_index <- max(which(QC$RunIndex < x))
               max_index <- min(which(QC$RunIndex > x))
               min_val <- QC$RunIndex[min_index]
               max_val <- QC$RunIndex[max_index]
               frac <- (max_val - x) / (max_val - min_val)
               unlist(frac * QC[min_index, metabolite_columns] +
                   (1-frac) * QC[max_index, metabolite_columns])
           }) %>% t %>% as.matrix

metabolite_indices <- -one_of("Id", "Type", "Gender", "Age", "APOE",
                              "Batch", "Code", "Index",
                              vars=colnames(subject_data))

## Normalization is weird don't do for now
## subject_data[, metabolite_indices] <- subject_data %>% select(metabolite_indices) / normalizing_constants


hist(log(filter(subject_data, Type=="AD", Metabolite=="Homocysteine")$Abundance), breaks=30, col="red")
hist(log(filter(subject_data, Type=="PD", Metabolite=="Homocysteine")$Abundance), breaks=30, col="blue", add=TRUE)



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
x <- c("CY", "CM", "CO")
subject_data %>% filter(Metabolite == met, Type %in% x) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=x),
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
## Type
subject_data %>% filter(Metabolite == met, Type %in% x) %>%
    ggplot(aes(x = Abundance, y = factor(Type, levels=x),
               fill=Type)) +  geom_joy(scale = 4) + theme_joy() + ggtitle(met)

## Gender
subject_data %>% filter(Metabolite == met, Type %in% c("CO", "CM", "CY")) %>%
    ggplot(aes(x = Abundance, y = Gender,
               fill=Gender)) +  geom_joy(scale = 4) + theme_joy() + ggtitle(met)

met <- "Kynurenine"
subject_data %>% filter(Metabolite == met) %>%
    mutate(Type=ifelse(Type %in% c("CO", "CM", "CY"), "C", Type)) %>% 
    ggplot(aes(x=Age, y=Abundance, color=Type)) + geom_smooth(method="lm") + ggtitle(met)

met <- "HIAA"
subject_data %>% filter(Metabolite == met) %>%
    mutate(Type=ifelse(Type %in% c("CO", "CM", "CY"), "C", Type)) %>% 
    ggplot(aes(x=Age, y=Abundance, color=Type)) + geom_smooth(method="lm") + ggtitle(met)

met <- "Methionine"
subject_data %>% filter(Metabolite == met) %>%
    mutate(Type=ifelse(Type %in% c("CO", "CM", "CY"), "C", Type)) %>% 
    ggplot(aes(x=Age, y=Abundance, color=Type)) + geom_smooth(method="lm") + ggtitle(met)

met <- "Dopamine"

subject_data %>% filter(Metabolite == met) %>%
    mutate(Type=ifelse(Type %in% c("CO", "CM", "CY"), "C", Type)) %>% 
    ggplot(aes(x=Age, y=Abundance, color=Type)) + geom_smooth(method="lm") + ggtitle(met)

subject_data %>% filter(Metabolite == met) %>%
    mutate(Type=ifelse(Type %in% c("CO", "CM", "CY"), "C", Type)) %>%
    summarize(pval=lm(Abundance ~ Age)$coefficients)


subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("CO", "CM", "CY")) %>% 
    do(lmres=lm(Abundance ~ Age, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients[2, 1],
              pval=summary(lmres)$coefficients[2, 4]) %>%
    filter(pval < 0.05) %>%
    arrange(abs(est))

## Check if control is differen than AD
ad_diff <- subject_data %>% group_by(Metabolite) %>%
    filter(Type %in% c("CO", "AD")) %>% 
    do(lmres=lm(Abundance ~ Age + Type, data=.)) %>%
    summarize(Met=Metabolite,
              est=summary(lmres)$coefficients[2, 1],
              pval=summary(lmres)$coefficients[2, 4]) %>%
    filter(pval < 0.05) %>%
    arrange(desc(abs(est)))]

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
