## envelope model for CSF data
## Envelope models do much worse when sighat is a poor estimator of sigma
##
library(tidyverse)
library(ggExtra)
library(magrittr)
library(microbenchmark)
library(mvtnorm)
library(mvnfast)
library(rstiefel)
library(mgCov)
library(Amelia)
library(modelr)
library(robust)
library(ggridges)
library(Matrix)
library(car)
library(patchwork)
library(glmnet)

source("utility.R")

processed_files <- dir(pattern="^preprocessed_gotms_data*")
## Most recent file
load(max(processed_files[grep("-20+", processed_files)]))
load("../data/got-ms/identification_map.RData")

wide_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "Trend",
                          "RunIndex", "Name","Data File")) %>%
    spread(key=Metabolite, value=Abundance)
dim(wide_data)

## Impute missing values
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

met_var_qc <- QC_long %>%
    mutate(Abundance=replace(Abundance, Abundance == -Inf, NA)) %>% 
    group_by(Metabolite, Mode) %>%
    summarize(qc_var = var(Abundance, na.rm=TRUE),
              qc_mad = median(abs(Abundance - median(Abundance, na.rm=TRUE)), na.rm=TRUE),
              qc_mean=mean(Abundance))

met_var_qc <- met_var_qc %>% ungroup %>%
    mutate(Mode = str_replace(met_var_qc$Mode, "QC_", "")) %>%
    unite(Metabolite, Metabolite, Mode)

## Match QC variances with Y columns
met_var_qc %<>% arrange(colnames(Y))


D <- Diagonal(ncol(Y), 1/sqrt(apply(Y, 2, function(x) var(x, na.rm=TRUE))))

## This overfits
betaHat_MLE <- solve(t(X) %*% X) %*% t(X) %*% (Y %*% D)

cvmfit <- cv.glmnet(X, as.matrix(Y %*% D), family = "mgaussian",
                    type.measure="mae", nfolds=10)

glmres <- glmres <- glmnet(X, Y %*% D, family="mgaussian", alpha=0, lambda=cvmfit$lambda.min)

betaHat <- do.call(cbind, glmres$beta)


residual <- Y -  X %*% betaHat
s <- getRank(residual)

## create residual list
## types <- c("CO", "PD", "AD")
## types <- c("CY", "CM", "CO")
types <- unique(wide_data$Type)

ngroups <- length(types)
residual_tibble <- as_tibble(as.matrix(residual))
residual_tibble$Type <- wide_data$Type
nvec_table <- residual_tibble %>%
    filter(Type %in% types) %>%
    group_by(Type) %>%
    summarise(n = n())
nvec <- nvec_table$n
names(nvec) <- nvec_table$Type
nvec <- nvec[types]

residuals_list <- Slist <- list()
for (type in types) {
    group_resid <- residual_tibble %>%
        filter(wide_data$Type == type) %>%
        select(-Type) %>% as.matrix()
  residuals_list[[type]] <- group_resid
}

p <- ncol(residuals_list[[1]])

## Weighted initialization
isoVar <- sapply(1:ngroups, function(i)
    median(apply(residuals_list[[i]], 2, var)))

weightsList <- lapply(1:ngroups, function(i) {
    evals <- svd(residuals_list[[i]]/sqrt(nvec[i]))$d^2
    dp <- ifelse(evals/isoVar[i] <= (1+sqrt(p/nvec[i])),
                 0,
                 sqrt((1-p/nvec[i]/(evals/isoVar[i]-1)^2) /
                      (1+p/nvec[i]/(evals/isoVar[i]-1))))

    weights <- 1/(1-dp) - 1
    weights[1:min(nvec[i], s)]
})


combined <- do.call(cbind, lapply(1:ngroups, function(k) {
    svdk <- svd(t(residuals_list[[k]]))$u[, 1:min(nvec[k], s)] %*%
                                      diag(weightsList[[k]][1:min(nvec[k], s)])
}))

Vinit <- svd(combined)$u[, 1:s]
Vinit <- Vinit %*% rustiefel(s, s)

EMFit <- subspaceEM(residuals_list, S=s, Vstart=Vinit, verbose=TRUE)

Vfit <- EMFit$V

compute_variance_explained(Vfit, residuals_list, 1/EMFit$PrecVec)

samples <- fitBayesianSpike(V=Vfit, Ylist=residuals_list, R=s, Q=0,
                                 niters=10000, nskip=100)

save(EMFit, Vfit, s, p, nvec, residuals_list,
     file=sprintf("../results/EM_ND-%s.RData", format(Sys.Date(), "%m-%d")))
save(samples,
     file=sprintf("../results/samples_ND-%s.RData", format(Sys.Date(), "%m-%d")))

level_key <- all_matches$Metabolite
names(level_key) <- paste(all_matches$Name, all_matches$Mode, sep="_")

metabolite_matches <- dplyr::recode(
    str_replace(rownames(Vfit), " Results_(pos|neg)[0-9]*", "_\\1"),
    !!!level_key)

rownames(Vstar) <- metabolite_matches
rownames(Vfit) <- metabolite_matches

metabolite_matches[order(abs(betaHat["Age", ]), decreasing=TRUE)] %>% head(50)
sort(abs(betaHat["Age", ]), decreasing=TRUE)
metabolite_matches[order(abs(betaHat["APOE44", ]), decreasing=TRUE)] %>% head(50)
sort(abs(betaHat["Age", ]), decreasing=TRUE)
metabolite_matches[order(abs(betaHat["Type2C", ]), decreasing=TRUE)] %>% head(50)
sort(abs(betaHat["Age", ]), decreasing=TRUE)
metabolite_matches[order(abs(betaHat["APOE", ]), decreasing=TRUE)] %>% head(20)


mag <- apply(Vstar[, 1, drop=FALSE], 1, function(x) sqrt(sum(x^2)))
pos_mag <- apply(Vstar[, 1, drop=FALSE], 1, function(x) sqrt(sum(x^2)))


## CO vs AD
groups_to_plot <- c(1, 2)
names(groups_to_plot) <- types[groups_to_plot]

create_plots(V=Vfit, samples, group1=1, group2=2,
             to_plot = groups_to_plot, view=c(3, 4),
             plot_type="biplot", label_size=2.5)

## CO vs PD
groups_to_plot <- c(2, 5)
names(groups_to_plot) <- types[groups_to_plot]

pdf("../figs/mgCov/pd_biplot1.pdf")
create_plots(V=Vfit, samples, group1=2, group2=5,
             to_plot = groups_to_plot, view=c(1, 2),
             plot_type = "biplot", label_size=3)
dev.off()

pdf("../figs/mgCov/pd_biplot2.pdf")
create_plots(V=Vfit, samples, group1=2, group2=5,
             to_plot = groups_to_plot, view=c(3, 4),
             plot_type="biplot", label_size=2.5)
dev.off()


## Age analysis
types
groups_to_plot <- c(4, 3, 2)
names(groups_to_plot) <- types[groups_to_plot]

pdf("../figs/mgCov/aging_biplot.pdf")
create_plots(V=Vfit, samples, group1=2, group2=4,
             to_plot = groups_to_plot, view=c(1, 2),
             plot_type="biplot")
dev.off()

pdf("../figs/mgCov/aging_posterior.pdf")
create_plots(V=Vfit, samples, group1=2, group2=4,
             to_plot = groups_to_plot, view=c(1, 2),
             plot_type="posterior")
dev.off()






pooledResiduals <- do.call(rbind, residuals_list)
preds <- makePrediction(pooledResiduals, Vfit,
                        samples$Osamps,
                        samples$omegaSamps,
                        samples$s2samps,
                        ngroups=length(nvec),
                        nsamps=100,
                        numToAvg=50)




## For how many
all_types <- rep(types, nvec)
sapply(unique(all_types), function(type) {
    max_idx <- apply(preds[all_types==type, ], 1, which.max)
    type_id <- which(unique(all_types)==type)
    mean(max_idx == type_id)
})



proj_dat <- as_tibble(pooledResiduals %*% Vfit %*% O)
colnames(proj_dat) <- c("V1", "V2")
proj_dat$Type <- all_types
proj_dat %>% filter(Type %in% types[to_plot]) %>% ggplot() + geom_point(aes(x=V1, y=V2, col=Type))

proj_dat %>% ggplot() + geom_density_ridges(aes(x=V1, y=Type))
proj_dat %>% ggplot() + geom_density_ridges(aes(x=V2, y=Type))

