## envelope model for CSF data
## Envelope models do much worse when sighat is a poor estimator of sigma
##
library(tidyverse)
library(ggExtra)
library(magrittr)
library(made4)
library(ggjoy)
library(gbm3)
library(huge)
library(microbenchmark)
library(mvtnorm)
library(mvnfast)
library(rstiefel)
library(Amelia)
library(modelr)
library(robust)

source("~/course/rstiefel/R/opt.stiefel.R")
source("envelope_functions.R")
source("utility.R")

## load("preprocessed_gotms_data.RData")
load("preprocessed_gotms_data-2018-03-21.RData")

## Dimension of envelope reduction
s <- 2

wide_data <- subject_data %>%     
    filter(!(Type %in% c("Other"))) %>%
    unite("Metabolite", c("Metabolite", "Mode")) %>% 
    mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Raw", "RawScaled", "RunIndex", "Trend", "Batch", "Name", "Id", "Data File")) %>%
    spread(key=Metabolite, value=Abundance)

## Impute missing values
Y <- wide_data %>%
  dplyr::select(-one_of("Type", "Type2", "Gender", "Age", "APOE", "Index")) %>%
  as.matrix()
Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
Y <- t(Yt) %>% as.matrix

X <- wide_data %>% 
  dplyr::select(one_of("Gender", "Age", "Type", "APOE")) %>%
  mutate(APOE = droplevels(APOE)) %>% 
  model_matrix(~ Age + Type + Gender + APOE - 1) %>%
  dplyr::select(-one_of("APOEUnknown")) %>% 
  as.matrix


Age_comparison <- wide_data %>% filter(Type2 == "C") %>% 
mutate(Type2 = droplevels(Type2), Type = droplevels(Type))

Y <- wide_data %>% mutate(Type2 = droplevels(Type2), Type = droplevels(Type)) %>% dplyr::select(-one_of("Index", "Type", "Type2", "Gender", "Age", "APOE", "Index")) 
Yt <- amelia(t(Y), m = 1, empri = 100)$imputations$imp1
Y <- t(Yt) %>% as.matrix

met_var_qc <- QC_long %>%
    mutate(Abundance=replace(Abundance, Abundance == -Inf, NA)) %>% 
    group_by(Metabolite, Mode) %>%
    summarize(qc_var = var(Abundance, na.rm=TRUE),
              qc_mad = median(abs(Abundance - median(Abundance, na.rm=TRUE)), na.rm=TRUE),
              qc_mean=mean(Abundance))

met_var_qc <- met_var_qc %>% ungroup %>% mutate(Mode = str_replace(met_var_qc$Mode, "QC_", "")) %>% unite(Metabolite, Metabolite, Mode)

## Match QC variances with Y columns
met_var_qc %<>% arrange(colnames(Y))

r <- 20

betaHat <- solve(t(X) %*% X) %*% t(X) %*% Y
betav <- svd(betaHat)$v
Vinit <- cbind(betav, NullC(betav)[, 1:r])

D <- Diagonal(ncol(Y), 1/sqrt(met_var_qc$qc_var))
D <- Diagonal(ncol(Y), 1/sqrt(apply(Y, 2, function(x) var(x, na.rm=TRUE))))


## works well with large s and small r??
s <- 10
r=10
Vfit <- fit_envelope(Y, X, s=s, r=r, prior_counts=1000)


## colnames(Y)[order(rowSums(Vfit[, 1:2]^2), decreasing=TRUE)] %>% head(n=20)

Yproj <- as.matrix(Y) %*% Vfit[, 1:s]
YprojPerp <- as.matrix(Y) %*% Vfit[, (s+1):(s+r)]

eta <- solve((1 + 10000/nrow(Y))*t(X) %*% X) %*% t(X) %*% Yproj
eta

## Vfit2 <- Vfit[, 1:s] %*% svd(eta)$v
Vfit2 <- Vfit[, 1:s] %*% t(eta)
YprojGender <- as.matrix(Y) %*% Vfit[, 1:s] %*% t(eta)[, 7]
YprojPD <- as.matrix(Y) %*% Vfit[, 1:s] %*% t(eta)[, 6]
APOEproj <- as.matrix(Y) %*% rowMeans((Vfit[, 1:s] %*% svd(eta[8:11, ])$v[, 1:2]))


Yproj <- as.matrix(Y) %*% Vfit2

## Age
tib <- tibble(x=Yproj[, 1], y=Yproj[, 2], Gender=factor(X[, "GenderM"]), Age=as.numeric(X[, "Age"]), Type=wide_data$Type2)
tib <- tib %>% filter(Type == "C")
ggplot(tib) + geom_point(aes(x=x, y=y, col=Age, shape=Type)) + scale_color_gradientn(colors=heat.colors(10))


tib <- tibble(x=Yproj[, 6], y=Yproj[, 4], Gender=factor(X[, "GenderM"]), Age=as.numeric(X[, "Age"]), Type=wide_data$Type2)
tib <- tib %>% filter(Type %in% c("C", "PD"))
ggplot(tib) + geom_point(aes(x=x, y=y, col=Type, shape=Type))



ggplot(tibble(GenderAxis=as.numeric(YprojGender), Gender=factor(X[, "GenderM"]), Age=X[, "Age"])) + geom_joy(aes(x=GenderAxis, y=Gender, fill=Gender))
ggplot(tibble(x=as.numeric(YprojPD), Type=wide_data$Type2)) + geom_joy(aes(x=x, y=Type, fill=Type))
ggplot(tibble(x=as.numeric(APOEproj), APOE=wide_data$APOE)) + geom_joy(aes(x=x, y=APOE, fill=APOE))





ggplot(tibble(x=Yproj[, 5], y=Yproj[, 4], Gender=factor(X[, "GenderM"]), Age=X[, "Age"], Type=wide_data$Type2)) + geom_point(aes(x=x, y=y, shape=Type, col=Type))

ggplot(tibble(Axis=Yproj[, 1]-Yproj[, 2], Gender=Xplot$Gender, Age=Xplot$Age)) + geom_joy(aes(x=Axis, y=Age < 40, fill=Age < 40))

ggplot(tibble(Axis=Yproj[, 1], Gender=Xplot$Gender, Age=Xplot$Age)) + geom_joy(aes(x=Axis, y=Gender, fill=Gender))

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
