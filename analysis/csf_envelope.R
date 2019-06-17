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
library(made4)
library(gbm3)
library(huge)
library(modelr)
library(robust)
library(cowplot)
library(colorspace)



source("utility.R")
load("../data/preprocessed_csf_data.RData")


## X <- matrix(rnorm(N, sd = 5), ncol=1)
wide_data <- subject_data %>% 
  filter(!(Type %in% c("Other"))) %>%
  mutate(Type = droplevels(Type), Type2 = droplevels(Type2)) %>%
    dplyr::select(-one_of("Code", "Mode", "RunIndex",
                          "Raw", "Trend", "RawScaled")) %>%
  spread(key = Metabolite, value = Abundance)

Y <- wide_data %>%
    dplyr::select(-one_of("Id", "Type", "Type2", "Gender",
                          "Age", "APOE", "Batch", "Index",
                          "GBAStatus", "GBA_T369M", "cognitive_status")) %>%
    as.matrix()

Y <- amelia(Y, m = 1, empri = 100)$imputations$imp1


X <- wide_data %>% 
  dplyr::select(one_of("Gender", "Age", "Type2", "APOE")) %>%
  mutate(APOE = factor(APOE, c("22", "23", "24", "33", "34", "44"))) %>%
  model_matrix(~ Age + Type2 + Gender + APOE) %>% as.matrix

s <- getRank(Y)
r <- 0

res <- fit_envelope(Y, X, s=s, r=r, Vinit="OLS",
                    prior_counts=10, use_py=FALSE, maxIters=1000,
                    U1=0.1*diag(s), U0=0.1*diag(r),
                    searchParams=list(rho=0.1, eta=0.9))


names(res)
plot(res$beta_ols, res$beta_env, xlim=c(-2, 2), ylim=c(-2, 2))

plot(res$beta_env["APOE44", ])

Vfit <- res$V
rownames(Vfit) <- colnames(Y)
eta <- res$eta_hat
dim(eta)

### Age Fit
Z <- eta["Age", ] %*% t(Vfit)
names(Z) <- rownames(Vfit)
Z[order(abs(Z), decreasing=TRUE)]

tibble(Yproj=as.numeric(eta["Age", ] %*% t(Vfit) %*% t(Y)),
       Age=wide_data$Age, Type=wide_data$Type2) %>%
    ggplot() + geom_point(aes(x=Age, y=Yproj, col=Type)) +
    geom_smooth(aes(x=Age, y=Yproj), method="lm")


### Age Fit
Z <- as.numeric(eta["GenderM", ] %*% t(Vfit))
names(Z) <- rownames(Vfit)
Z[order(abs(Z), decreasing=TRUE)]

tibble(Yproj= as.numeric(Z %*% t(Y)), Gender=wide_data$Gender,
       Type=wide_data$Type2) %>%
    ggplot() + geom_density_ridges(aes(y=Gender, x=Yproj))

t.test(as.numeric(Z %*% t(Y))[wide_data$Gender == "M"],
       as.numeric(Z %*% t(Y))[wide_data$Gender == "F"])

## Check APOE for controls vs AD
indices <- wide_data$Type != "PD"
res <- fit_envelope(Y[indices, ], X[indices, ], s=s, r=r, Vinit="OLS",
                    prior_counts=10, use_py=FALSE, maxIters=1000,
                    U1=0.1*diag(s), U0=0.1*diag(r),
                    searchParams=list(rho=0.1, eta=0.9))




res <- fit_envelope(Y[wide_data$Type == "AD",],
                    X[wide_data$Type == "AD",
                      grep("APOE", colnames(X)), drop=FALSE],
                    s=s, r=r, use_py=FALSE)

YP_apoe <- Y[wide_data$Type == "AD",] %*% res$V %*% svd(res$eta)$v[, 2:3]
tib <- tibble(x=YP_apoe[, 1], col=wide_data$APOE[wide_data$Type == "AD"])
tib %>% ggplot(aes(x=x, y=y, col=col)) + geom_point()
tib %>% ggplot() + geom_density_ridges(aes(x=x, y=col, fill=col)) + xlab("")

plot(apply(Y, 1, mean))






X %*% beta_env

Y$Methylhistamine

beta_env <- res$beta_env
colnames(beta_env) <- colnames(Y)

YprojGender <- as.matrix(Y) %*% Vfit[, 1:s] %*% t(eta)[, "GenderM"]
YprojPD <- as.matrix(Y) %*% Vfit[, 1:s] %*% t(eta)[, "Type2PD"]
APOEproj <- as.matrix(Y) %*% rowMeans((Vfit[, 1:s] %*% svd(eta[grepl("APOE", rownames(eta)), ])$v[, 1:2]))

apoe_plot <- tibble(x=as.numeric(APOEproj), APOE=wide_data$APOE, Gender=wide_data$Gender) %>%
    ggplot() + geom_density_ridges(aes(x=x, y=APOE, fill=APOE)) + scale_fill_brewer(type="qual", palette=6)
apoe_plot




Ysub <- as.matrix(Y) %*% Vfit[, 1:s]
dim(Ysub)


APOE_indicator <- wide_data %>% filter(Type2 != "PD") %>% dplyr::select(APOE)
APOE_indicator
tibble(x = as.numeric(unlist(Y2 %*% V[, 2])),
       y=as.factor(unlist(APOE_indicator))) %>%
    ggplot(aes(x=x, y=y)) +
    geom_density_ridges(aes(fill=y)) +
    scale_fill_brewer(type="qual", palette=6)


wide_data[, -(1:8)] <- Y

folds <- crossv_kfold(wide_data, k = nrow(wide_data))
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
    dplyr::select(-one_of("Id", "Type", "Type2", "Gender", "Age", "APOE", "Batch", "Index"))

  X <- df %>%
    mutate(Type2 = relevel(Type2, 2)) %>%
    dplyr::select(one_of("Type2", "Gender", "Age", "APOE", "Batch", "Index")) %>%
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
AD_comparison <- wide_data %>%
  filter(Type %in% c("CO", "PD")) %>%
  mutate(Type2 = droplevels(Type2), Type = droplevels(Type))

Y <- AD_comparison %>% dplyr::select(-one_of("Id", "Type", "Type2", "Gender", "Age", "APOE", "Batch", "Index"))
if(any(is.na(Y))) {
    Y <- amelia(Y, m = 1, empri = 1000)$imputations$imp1
}

X <- AD_comparison %>%
  mutate(Type2 = relevel(Type2, 2)) %>%
  dplyr::select(one_of("Type2", "Gender", "Age", "APOE")) %>%
  mutate(APOE = droplevels(APOE)) %>% model.matrix

X <- wide_data %>% 
  dplyr::select(one_of("Gender", "Age", "Type2", "APOE")) %>%
  mutate(APOE = droplevels(factor(APOE, c("Unknown", "23", "24", "33", "34", "44")))) %>% 
  model_matrix(~ Age + Type2 + Gender + APOE - 1)
X <- X %>% select(-"Type2C") %>% as.matrix

res <- fit_envelope(as.matrix(Y), as.matrix(X), s=20, use_py=FALSE)


Vfit <- res$Vfit

rbpal <- colorRampPalette(c("red", "white", "blue"))
cols <- rbpal(10)[as.numeric(cut(X[, "Age"], breaks = 10))]

cols <- case_when(
  X[, "Type2C"] == 1 ~ "red",
  X[, "Type2AD"] == 1 ~ "blue",
  X[, "Type2PD"] == 1 ~ "orange",
  TRUE ~ "green"
)

Yproj <- as.matrix(Y) %*% Vfit %>% as.tibble()
proj_data <- bind_cols(Yproj, X)

proj_data %>% ggplot(aes(x = V1, y = V2, col = Type2)) + geom_point()
proj_data %>% ggplot(aes(x = V1, y = V3, col = Age)) + geom_point()

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
