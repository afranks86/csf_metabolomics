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
library(Amelia)
library(modelr)
library(robust)
library(ggridges)
library(envelopeR)
library(Matrix)

source("utility.R")

## load("preprocessed_gotms_data.RData")
##load("preprocessed_gotms_data-2018-03-21.RData")

processed_files <- dir(pattern="^preprocessed_gotms_data*")
## Most recent file
load(max(processed_files[grep("-20+", processed_files)]))


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

X <- wide_data %>% 
  dplyr::select(one_of("Gender", "Age", "Type2", "APOE")) %>%
  mutate(APOE = factor(APOE, c("22", "23", "24", "33", "34", "44"))) %>%
  model_matrix(~ Age + Type2 + Gender + APOE - 1) %>% as.matrix

W <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select("RunIndex", "Metabolite", "Abundance") %>% 
    spread(key=Metabolite, value=Abundance)

W <- W %>% dplyr::select("RunIndex")


## Try TSNE
tsne_fit <- Rtsne(Y) 
tsne_plot <- tibble(x=tsne_fit$Y[, 1], y=tsne_fit$Y[, 2],
                    Type=wide_data$Type2, Age=wide_data$Age) %>%
    ggplot(aes(x=x, y=y)) + geom_point(aes(col=Age, shape=Type), size=3) +
    scale_colour_distiller(palette="Spectral") + xlab("tsne 1") + ylab("tsne 2")

save_plot("~/Desktop/gotms_tsne.png", tsne_plot)

Age_comparison <- wide_data %>% filter(Type2 == "C") %>% 
    mutate(Type2 = droplevels(Type2), Type = droplevels(Type))

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

s <- ncol(X)
r <- getRank(Y) - ncol(X)
q <- ncol(X)

betaHat <- solve(t(X) %*% X) %*% t(X) %*% Y
residual <- Y -  X %*% betaHat
getRank(residual)

Y[1:5, 1:5]
ncol(X)
dim(svd(residual)$v)
tr((t(svd(betaHat)$v) %*% svd(residual)$v) %*%
   (t(svd(residual)$v) %*%  svd(betaHat)$v)) / 10


betav <- svd(betaHat)$v
Vinit <- cbind(betav, NullC(betav)[, 1:(s+r-q)])

D <- Diagonal(ncol(Y), 1/sqrt(met_var_qc$qc_var))
D <- Diagonal(ncol(Y), sqrt(met_var_qc$qc_var))
D <- Diagonal(ncol(Y), 1/sqrt(apply(Y, 2, function(x) var(x, na.rm=TRUE))))
D <- Diagonal(ncol(Y), 1)

## works well with large s and small r??

prior_counts <- 10

res <- envelopeR::fit_envelope(Y=Y, X=X, D=D, s=s, r=r,
                               prior_counts=prior_counts,
                               maxIters=200, U1=diag(s), U0=diag(r),
                               Vinit="OLS", use_py=FALSE)
Vfit <- res$V

plot(res$beta_env , betaHat)

## colnames(Y)[order(rowSums(Vfit[, 1:2]^2), decreasing=TRUE)] %>% head(n=20)

eta <- res$eta_hat
Yproj <- as.matrix(Y) %*% Vfit[, 1:s]
YprojPerp <- as.matrix(Y) %*% Vfit[, (s+1):(s+r)]


beta_env <- res$beta_env
colnames(beta_env) <- colnames(Y)

sort(abs(beta_env["APOE23", ]- beta_env["APOE44", ]), decreasing=TRUE) %>%
    head(n=100)

sort(abs(beta_env["APOE44", ]), decreasing=TRUE) %>% head(n=50)

eta
YprojGender <- as.matrix(Y) %*% Vfit[, 1:s] %*% t(eta)[, "GenderM"]
YprojPD <- as.matrix(Y) %*% Vfit[, 1:s] %*% t(eta)[, "Type2PD"]
YprojAD <- as.matrix(Y) %*% Vfit[, 1:s] %*% t(eta)[, c("Type2AD", "GenderM")]
APOEproj <- as.matrix(Y) %*%
    rowMeans((Vfit[, 1:s] %*%
              svd(eta[grepl("APOE", rownames(eta)), ])$v[, 1:2]))


Ysub <- as.matrix(Y) %*% Vfit[, 1:s]
dim(Ysub)


##############################
## age plot
##############################
head(X)
res2 <- fit_envelope(Yproj, X[, "Age", drop=FALSE], s=2, r=6,
                     maxIters=1000, use_py=FALSE,
                     searchParams=list(rho=0.1, eta=0.9))
YP_age <- Yproj %*% res2$V[, 1:2]
tib <- tibble(x=YP_age[, 1], y=YP_age[, 2], Gender=factor(X[, "GenderM"]), Age=as.numeric(X[, "Age"]), Type=factor(wide_data$Type, c("CY", "CM", "CO")), Type2=wide_data$Type2)
age_plot <- tib %>% ggplot(aes(x=x, y=y)) + geom_point(aes(colour=Age), size=3) + scale_colour_distiller(palette="Spectral") + xlab("Dimension 1") + ylab("Dimension 2")
save_plot("../Figs/age_plot_gotms.png", age_plot)

pdf("~/Desktop/age_ridges.pdf")
res <- fit_envelope(YP_age, X[, "Age", drop=FALSE], distn="normal", s=1, r=1,
                    searchParams=list(rho=0.1, eta=0.9))
YP_age <- Yproj %*% res$V[, 1]
tib <- tibble(x=YP_age[, 1], y=wide_data$Type)
tib %>%  filter(y %in% c("CY", "CM", "CO")) %>%  ggplot() + geom_density_ridges(aes(x=x, y=y, fill=y)) + xlab("")
dev.off()

tib %>% filter(y > -2)  %>% ggplot() + geom_point(aes(x=x, y=y, col=Age, shape=Gender), size=2) + scale_color_gradientn(colors=rainbow(3))
tib %>% ggplot() + geom_point(aes(x=x, y=y, col=Age, shape=Gender), size=2) + scale_color_gradientn(colors=rainbow(3))

age_loadings <- Vfit[, 1:s] %*% eta[1, ]
age_loadings <- age_loadings / sqrt(sum(age_loadings^2))
hist(age_loadings, breaks=100)
colnames(Y)[abs(age_loadings) > 0.05]
write.csv(colnames(Y)[abs(age_loadings) > 0.05],
          file="../results/age_features.csv",
          row.names=FALSE, quote=FALSE)

##############################
## PD plot
##############################
Yold <- Y[!(wide_data$Type %in% c("CY", "CM")),]
Xold <- X[!(wide_data$Type %in% c("CY", "CM")),]
res <- fit_envelope(Yold, Xold[, , drop=FALSE], s=2, r=40, use_py=TRUE)

YP_pd <- Yold %*% res$V[, 1:2]
tib <- tibble(x1=YP_pd[,1], x2=YP_pd[, 2], y=Xold[, "Type2PD"])
tib %>%  ggplot(aes(x=x1, y=factor(y))) + geom_density_ridges()
tib %>%  ggplot(aes(x=x2, y=factor(y))) + geom_density_ridges()
tib %>%  ggplot() + geom_point(aes(x=x1, y=x2, col=y)) + xlab("")

tib <- tibble(x=Yproj[, "Type2PD"], y=Yproj[, "Type2C"], Gender=factor(X[, "GenderM"]), Age=as.numeric(X[, "Age"]), Type=wide_data$Type2)
tib %>% filter(Type %in% c("C", "PD")) %>% 
    ggplot() + geom_point(aes(x=x, y=y, col=Type, shape=Type))

pd_loadings <- (Vfit[, 1:s] %*% t(eta))[, c("TypePD", "TypeCO")] %*%  prcomp(Yproj[, c("TypePD", "TypeCO")])$rotation[, 2]

tibble(x=as.numeric(Y %*% pd_loadings), Type=wide_data$Type) %>%
    filter(Type %in% c("CO", "PD")) %>% 
    ggplot() + geom_density_ridges(aes(x=x, y=Type, fill=Type))

names(pd_loadings) <- colnames(Y)
hist(pd_loadings , breaks=100)
colnames(Y)[abs(pd_loadings) > 0.5]

pd_loadings <- pd_loadings[abs(pd_loadings) > quantile(abs(pd_loadings), 0.95)]

write.csv(colnames(Y)[abs(pd_loadings) > 0.1],
          file="../results/pd_features.csv",
          row.names=FALSE, quote=FALSE)

##############################
## AD plot
##############################
res <- fit_envelope(Yproj[, c("Type2AD", "Type2C")], X[, c("Type2AD"), drop=FALSE], s=1, r=1, prior_counts=10)

dim(Yproj)

YP_ad <- Yproj[, c("Type2AD", "Type2C")] %*% res$V[, 1]
t.test(YP_ad[wide_data$Type2=="AD",],YP_ad[wide_data$Type2=="C",])
tib <- tibble(x=YP_ad, y=wide_data$Type2)
tib %>%  ggplot() + geom_density_ridges(aes(x=x, y=y, fill=y)) + xlab("")

tib <- tibble(x=Yproj[, "Type2AD"], y=Yproj[, "Type2C"], Gender=factor(X[, "GenderM"]), Age=as.numeric(X[, "Age"]), Type=wide_data$Type2)
tib %>% filter(Type %in% c("C", "AD")) %>% 
    ggplot() + geom_point(aes(x=x, y=y, col=Type, shape=Type))




tib <- tibble(x=Yproj[, "TypeAD"], y=Yproj[, "TypeCO"], APOE=factor(wide_data$APOE), Age=as.numeric(X[, "Age"]), Gender=as.numeric(X[, "GenderM"]), Type=wide_data$Type)
ad_plot <- tib %>% filter(Type %in% c("CO", "AD")) %>% 
    ggplot() + geom_point(aes(x=x, y=y, col=Type, shape=as.factor(Type)))

ad_plot


tib <- tibble(x=Yproj[, 2], y=Yproj[, 3], APOE=factor(wide_data$APOE), Age=as.numeric(X[, "Age"]), Gender=as.numeric(X[, "GenderM"]), Type=wide_data$Type)
ad_plot <- tib %>% filter(Type %in% c("CO", "AD")) %>% 
    ggplot() + geom_point(aes(x=x, y=y, col=Type, shape=as.factor(APOE)))

ad_plot


ad_loadings <- (Vfit[, 1:s] %*% t(eta)[, 1:9])[, 2:3] %*%  prcomp(Yproj[, 2:3])$rotation[, 2]

tibble(x=as.numeric(Y %*% ad_loadings), Type=wide_data$Type) %>%
    filter(Type %in% c("CO", "AD")) %>% 
    ggplot() + geom_density_ridges(aes(x=x, y=Type, fill=Type))

hist(ad_loadings , breaks=100)
colnames(Y)[abs(ad_loadings) > 0.05]
write.csv(colnames(Y)[abs(ad_loadings) > 0.05],
          file="../results/ad_features.csv",
          row.names=FALSE, quote=FALSE)

library(slipper)
xsamps <- tib %>% filter(Type %in% c("CO", "AD")) %>% group_by(Type) %>%
    do(slipper(mean(x), B=100, df=.))
ysamps <- tib %>% filter(Type %in% c("CO", "AD")) %>% group_by(Type) %>%
    do(slipper(mean(y), B=100, df=.))

bind_cols(xsamps, ysamps[, -c(1:2)]) %>% ggplot() + geom_point(aes(x=value, y=value1, col=Type))


ad_means <- tib %>% filter(Type %in% c("CO", "AD")) %>% group_by(Type) %>%
    summarise(meanx=median(x), meany=median(y))    




##############################
## Gender plot
##############################
tibble(GenderAxis=as.numeric(YprojGender), Gender=factor(X[, "GenderM"]), Age=X[, "Age"]) %>% 
    ggplot() + geom_density_ridges(aes(x=GenderAxis, y=Gender, fill=Gender))


gender_loadings <- Vfit[, 1:s] %*% eta["GenderM", ]
gender_loadings <- gender_loadings / sqrt(sum(gender_loadings^2))
hist(gender_loadings, breaks=100)
colnames(Y)[abs(gender_loadings) > 0.05]
write.csv(colnames(Y)[abs(gender_loadings) > 0.05],
          file="../results/gender_features.csv",
          row.names=FALSE, quote=FALSE)

##############################
## APOE plot
##############################
colnames(Y)
res <- fit_envelope(Y[wide_data$Type != "PD",],
                    X[wide_data$Type != "PD", grep("APOE", colnames(X)), drop=FALSE],
                    s=4, r=30, maxIters=250, use_py=FALSE)

YP_apoe <- Y[wide_data$Type != "PD",] %*% res$V[, 1:2]
tib <- tibble(x=YP_apoe[, 1], y=YP_apoe[, 2], col=wide_data$APOE[wide_data$Type != "PD"])
tib %>% ggplot(aes(x=x, y=y)) + geom_point(aes(col=col))
tib %>% ggplot() + geom_density_ridges(aes(x=x, y=col, fill=col)) + xlab("")



apoe_plot <- tibble(x=as.numeric(APOEproj),
                    APOE=wide_data$APOE,
                    Gender=wide_data$Gender) %>%
    ggplot() + geom_density_ridges(aes(x=x, y=APOE, fill=APOE)) +
    scale_fill_brewer(type="qual", palette=6)
apoe_plot

apoe_plot2 <- tibble(x=as.numeric(APOEproj), APOE=wide_data$APOE, Gender=wide_data$Gender) %>% unite("combo", c("APOE", "Gender"), remove=FALSE) %>%
    filter(APOE != "Unknown") %>% 
    ggplot() + geom_density_ridges(aes(x=x, y=combo, fill=APOE, col=Gender), alpha=0.5)
apoe_plot2

apoe_plot3 <- tibble(x=as.numeric(APOEproj), APOE=wide_data$APOE, Gender=wide_data$Gender, Type2=wide_data$Type2) %>% unite("combo", c("APOE", "Gender"), remove=FALSE) %>%
    filter(APOE != "Unknown") %>%
    filter(Type2 == "C") %>% 
    ggplot() + geom_density_ridges(aes(x=x, y=combo, fill=APOE, col=Gender), alpha=0.5)
apoe_plot3


tibble(x=as.numeric(APOEproj), APOE=wide_data$APOE, Gender=wide_data$Gender) %>% unite("combo", c("APOE", "Gender"), remove=FALSE) %>% select(combo) %>% table

meds <- tibble(x=as.numeric(APOEproj), APOE=factor(wide_data$APOE)) %>%
    filter(APOE != "Unknown") %>% 
    group_by(APOE) %>%
    summarise(med = median(x, na.rm=TRUE)) %>% ungroup()

apoe_plot + geom_segment(data=meds, aes(x=med, xend=med, y=0, yend=0.75, color=APOE), size=1.5)

apoe_loadings <- svd(Vfit[, 1:s] %*% t(eta[c("APOE24", "APOE33", "APOE34", "APOE44"), ]))$u[, 1]
apoe_loadings <- apoe_loadings / sqrt(sum(apoe_loadings^2))
hist(apoe_loadings, breaks=100)
colnames(Y)[abs(apoe_loadings) > 0.05]
write.csv(colnames(Y)[abs(apoe_loadings) > 0.05],
          file="../results/apoe_features.csv",
          row.names=FALSE, quote=FALSE)

bind_cols(as_tibble(Yproj), as_tibble(X)) %>% filter(Type %in% c("CO", "AD"))
AD_analysis <- fit_envelope(Yproj, X, s=s, r=r, prior_counts=1000, maxIters=2000, U1=diag(s), nu1=100, U0=diag(r), nu0 = 100, Vinit="OLS")




ggplot(tibble(x=Yproj[, 5], y=Yproj[, 4], Gender=factor(X[, "GenderM"]), Age=X[, "Age"], Type=wide_data$Type2)) + geom_point(aes(x=x, y=y, shape=Type, col=Type))

ggplot(tibble(Axis=Yproj[, 1]-Yproj[, 2], Gender=Xplot$Gender, Age=Xplot$Age)) + geom_density_ridges(aes(x=Axis, y=Age < 40, fill=Age < 40))

ggplot(tibble(Axis=Yproj[, 1], Gender=Xplot$Gender, Age=Xplot$Age)) + geom_density_ridges(aes(x=Axis, y=Gender, fill=Gender))

median(Yproj[Xplot$Gender=="F", 1]-Yproj[Xplot$Gender=="F", 2])
median(Yproj[Xplot$Gender=="M", 1]-Yproj[Xplot$Gender=="M", 2])


dim(X1[, 2] %*% tmp[2, ,drop=FALSE] %*% t(Vfit))

Xplot <- Age_comparison %>% 
  dplyr::select(one_of("Gender", "Age"))



proj_data <- bind_cols(as.tibble(as.matrix(Yproj)), Xplot)
summary_data <- proj_data %>% group_by(cut(Age, 3)) %>% summarise(med1=median(V1), med2=median(V2))
ggplot(proj_data)  + geom_point(aes(x = V1, y = V2, col = cut(Age, 3), shape=Gender, size=4)) +  geom_point(data=summary_data, aes(x=med1, y=med2, size=3,))



colnames(Y)[order(rowSums(Vfit^2), decreasing=TRUE)] %>% head(n=20)

plot(log(diag(D)), log(apply(Y, 2, function(x) var(x, na.rm=TRUE))))

plot(log(apply(Y, 2, function(x) var(x, na.rm=TRUE)))), log(apply(Y %*% D, 2, function(x) var(x, na.rm=TRUE))))

hist(log(apply(Y %*%  D, 2, function(x) var(x, na.rm=TRUE))))










folds <- crossv_kfold(wide_data, k = nrow(wide_data))
folds <- folds[1, ]

folds %<>% mutate(res = map(train, function(r) fit_resampled(r, S = S)))
folds %<>% mutate(prediction = map2(test, res, test_out_of_sample))

folds %<>% mutate(
  env = unlist(map(prediction, function(x) x[[1]])),
  mle = unlist(map(prediction, function(x) x[[2]])),
  mse_env = unlist(map(prediction, function(x) x[[3]])),
  mse_mle = unlist(map(prediction, function(x) x[[4]]))
)

bind_rows(
  summarise(folds, type = "ENV", mean = mean(env), med = median(env)),
  summarise(folds, type = "MLE", mean = mean(mle), med = median(mle)),
  summarise(folds, type = "MSE_ENV", mean = mean(mse_env), med = median(mse_env)),
  summarise(folds, type = "MSE_MLE", mean = mean(mse_mle), med = median(mse_mle)),
)

VfitLst <- lapply(folds$res, function(x) x$Vfit)
tst <- Reduce(function(y, x) {
  y %*% t(y) %*% x %*% t(x)
}, VfitLst[sample(length(VfitLst), 10)])
tr(tst) / S
tr(v1 %*% t(v1) %*% v2 %*% t(v2)) / S

met_var_qc <- QC_long %>%
  group_by(Metabolite) %>%
  summarize(
    qc_var = var(Abundance),
    qc_mad = median(abs(Abundance - median(Abundance, na.rm = TRUE)), na.rm = TRUE),
    qc_mean = mean(Abundance)
  )

fit_resampled <- function(resamp, S=2) {

  df <- resamp$data[resamp$idx, ]

  Y <- df %>%
    dplyr::select(-one_of("Index", "Type", "Type2", "Gender", "Age", "APOE", "Index"))

  X <- df %>%
    mutate(Type2 = relevel(Type2, 2)) %>%
    dplyr::select(one_of("Type2", "Gender", "Age", "APOE", "Index")) %>%
    model_matrix(~ Age + Type2 + Gender + APOE - 1) %>%
    dplyr::select(-one_of("APOEUnknown"))

  fit_envelope(Y, X, S)

}

test_out_of_sample <- function(test, params) {
  df <- test$data[test$idx, ]

  Y <- df %>%
    dplyr::select(-one_of("Id", "Type", "Type2", "Gender", "Age", "APOE", "Batch", "Index")) %>%
    as.matrix()

  X <- df %>%
    mutate(Type2 = relevel(Type2, 2)) %>%
    dplyr::select(one_of("Type2", "Gender", "Age", "APOE", "Batch", "Index")) %>%
    model_matrix(~ Age + Type2 + Gender + APOE - 1) %>%
    dplyr::select(-one_of("APOEUnknown")) %>%
    as.matrix()

  betaHatNew <- params$betaHatNew
  betaHat <- params$betaHat
  SigmaHatNew <- params$SigmaHatNew
  SigmaHat <- params$SigmaHat

  N <- nrow(Y)
  P <- ncol(Y)

  Z <- (Y - X %*% t(betaHat))
  ZNew <- (Y - X %*% t(betaHatNew))

  ## Out of sample log likelihood

  ll_env <- -N / 2 * (determinant(SigmaHatNew, log = TRUE)$modulus + ZNew %*% solve(SigmaHatNew) %*% t(ZNew))
  ll_mle <- -N / 2 * (determinant(SigmaHat, log = TRUE)$modulus + Z %*% solve(SigmaHat) %*% t(Z))
  mse_mle <- sum(Z ^ 2)
  mse_env <- sum(ZNew ^ 2)

  list(
    pred_env = ll_env,
    pred_mle = ll_mle,
    mse_mle = mse_mle,
    mse_env = mse_env
  )
}

##### Exploratory analysis on full data ################3

par(mfrow = c(1, 2))

## Plot data on subspace and random projection
Age_comparison <- wide_data %>%
  filter(Type %in% c("CO", "PD")) %>%
  mutate(Type2 = droplevels(Type2), Type = droplevels(Type))

Age_comparison <- wide_data %>% filter(Type2 == "C") %>% 
mutate(Type2 = droplevels(Type2), Type = droplevels(Type))

Y <- Age_comparison %>% dplyr::select(-one_of("Index", "Type", "Type2", "Gender", "Age", "APOE", "Index"))
Y <- as.matrix(Y)
n <- nrow(Y)
Y <- rmvnorm(n, mean=rep(0, ncol(Y)), sigma= ((t(Y) %*% Y) / n + diag(ncol(Y))))
Y <- rmvnorm(n, mean=rep(0, ncol(Y)), sigma= diag(ncol(Y)))

X <- Age_comparison %>%
  dplyr::select(one_of("Gender", "Age", "Type", "APOE")) %>%
  mutate(APOE = droplevels(APOE))

X <- Age_comparison %>%
  dplyr::select(one_of("Gender", "Age", "Type"))



res <- fit_envelope(Y, X, 3)


Vfit <- res$Vfit
Vfit <- rustiefel(ncol(Y), S)

Yproj <- as.matrix(Y) %*% Vfit[, 1:2] %>% as.tibble()
proj_data <- bind_cols(Yproj, X)

proj_data %>% ggplot(aes(x = V1, y = V2, col = Type, shape=Gender)) + geom_point()
proj_data %>% ggplot(aes(x = (V1-V2), y = V3, col = Age, shape=Gender)) + geom_point()

proj_data %>% ggplot(aes(x = V1, y = V2, col = Type, shape=Gender)) + geom_point()

proj2 <- proj_data %>% filter(Type2 == "C")
fit_envelope(proj2, 2)

Yboot %>% bootstrap(10, by_group = TRUE) %>% group_by(Type) %>% do(data.frame(mv1 = mean(.$V1), mv2 = mean(.$V2)))

Yboot %>%
  nest(-Type) %>%
  mutate(cors_boot = map(data, ~ bootstrap(., 10) %>% do(tidy(colMeans(.$V1, .$V2))))) %>%
  unnest(cors_boot)

ADsamps1 <- Yboot %>% filter(Type == "AD") %>% slipper(mean(V5), B = 100)
PDsamps1 <- Yboot %>% filter(Type == "PD") %>% slipper(mean(V5), B = 100)
COsamps1 <- Yboot %>% filter(Type == "CO") %>% slipper(mean(V5), B = 100)

ADsamps2 <- Yboot %>% filter(Type == "AD") %>% slipper(mean(V3), B = 100)
PDsamps2 <- Yboot %>% filter(Type == "PD") %>% slipper(mean(V3), B = 100)
COsamps2 <- Yboot %>% filter(Type == "CO") %>% slipper(mean(V3), B = 100)


par(mfrow = c(1, 1))
plot(ADsamps1$value, ADsamps2$value, col = "red", pch = 19, ylim = c(-0.5, 0.5), xlim = c(-0.5, 0.5), cex = 0.5)
points(PDsamps1$value, PDsamps2$value, col = "blue", pch = 19, ylim = c(-0.5, 0.5), xlim = c(-0.5, 0.5), cex = 0.5)
points(COsamps1$value, COsamps2$value, col = "green", pch = 19, ylim = c(-0.5, 0.5), xlim = c(-0.5, 0.5), cex = 0.5)

head(Vfit[order(abs(Vfit[, 5]), decreasing = TRUE), 1], n = 10)
head(Vfit[order(abs(Vfit[, 3]), decreasing = TRUE), 1], n = 10)

tmp <- Vfit[, c(5, 3)] %>% apply(1, function(x) sum(x ^ 2))
nms <- names(head(tmp[order(tmp, decreasing = TRUE)], n = 40))
text(Vfit[nms, c(5, 3)], labels = nms)

boostrap_proj <- sapply(1:100, function(i) colMeans(Yproj[sample(nrow(Yproj), replace = TRUE), ]))
plot(
  Yproj %*% rustiefel(S, 2), pch = 19,
  col = cols, main = "Subspace Projection"
)
plot(
  Y %*% rustiefel(P, 2),
  pch = 19, col = cols, main = "Random Projection"
)

## Envelope estimates
etaHat <- t(Yproj) %*% X %*% solve(t(X) %*% X)
betaHatNew <- Vfit %*% etaHat



resNew <- Y - betaHatNew %*% X
