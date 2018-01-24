## Testing envelope models
## Envelope models do much worse when sighat is a poor estimator of sigma
##

library(microbenchmark)
source("~/course/rstiefel/R/opt.stiefel.R")
library(mvtnorm)
library(mvnfast)

## Optimal threshold from Gavish, Donoho 2014
getRank <- function(Y) {
  svals <- svd(Y)$d

  m <- max(nrow(Y), ncol(Y))
  n <- min(nrow(Y), ncol(Y))

  if (m == n) {
    rank <- sum(svals > 2.858 * median(svals))
  } else {
    beta <- n / m
    omeg <- 0.56 * beta ^ 3 - 0.95 * beta ^ 2 + 1.82 * beta + 1.43
    rank <- sum(svals > omeg * median(svals))
  }

  rank
}

## Generate some fake data
P <- 200
N <- 100
S <- 2

## Generate envelope
V <- rstiefel::rustiefel(P, S)
Vperp <- rstiefel::NullC(V)

## regression parameters

## eta small relative to evals gives based improvements in MSE

eta <- c(-1, 1)
Beta <- V %*% eta

## covariates
X <- rnorm(N, sd = 5)

## covariance matrices
psi <- diag(sort(rexp(S, 1 / 1), decreasing = TRUE))
## psi <- diag(c(10, 10))
## psi.0 <- diag(rep(100, P-S))
psi.0 <- diag(sort(rexp(P - S, 1 / 100), decreasing = TRUE))

Sigma <- V %*% psi %*% t(V) + Vperp %*% psi.0 %*% t(Vperp)

Y <- sapply(1:N, function(i) rmvn(1, Beta %*% X[i], Sigma))

Xtest <- rnorm(N, sd = 5)
Ytest <- sapply(1:N, function(i) rmvn(1, Beta %*% Xtest[i], Sigma))

rnk <- getRank(Y - Beta %*% X)

## Plot Y colored by X on VV^T and a random subspace UU^T
rbpal <- colorRampPalette(c("red", "white", "blue"))
cols <- rbpal(10)[as.numeric(cut(X, breaks = 10))]

par(mfrow = c(1, 2))
plot(
  t(Y) %*% V, pch = 19,
  col = cols, main = "Subspace Projection"
)
plot(
  t(Y) %*% rstiefel::rustiefel(P, S),
  pch = 19, col = cols, main = "Random Projection"
)


betaHat <- Y %*% X %*% solve(t(X) %*% X)
plot(Beta, betaHat)
sum((Beta - betaHat) ^ 2)

##################### Test BENV #############################


res <- BENV::BayesEnvelopeMC(t(Y), matrix(X, ncol=1), S, 5)

## res$beta are the coefficients
## res$eta are the projected coefficients

###############################################################

tr(Vfit %*% t(Vfit) %*% V %*% t(V))


Yproj <- t(Vfit) %*% Y
par(mfrow = c(1, 2))


## Plot data on subspace and random projeciton
plot(
  t(Y) %*% Vfit, pch = 19,
  col = cols, main = "Subspace Projection"
)
plot(
  t(Y) %*% rustiefel(P, S),
  pch = 19, col = cols, main = "Random Projection"
)

## Envelope estimates
etaHat <- Yproj %*% X %*% solve(t(X) %*% X)
betaHatNew <- Vfit %*% etaHat

## MAD coefficients
print(sprintf("MAD MLE: %f", median(abs(Beta - betaHat))))
print(sprintf("MAD ENV: %f", median(abs(Beta - betaHatNew))))

## MSE coefficients
print(sprintf("MSE MLE: %f", sum((Beta - betaHat) ^ 2)))
print(sprintf("MSE ENV: %f", sum((Beta - betaHatNew) ^ 2)))

## RMSE predictions
print(sprintf("Pred MLE: %f", sqrt(sum((Ytest - betaHat %*% Xtest) ^ 2)/ (N*P))))
print(sprintf("Pred ENV: %f", sqrt(sum((Ytest - betaHatNew %*% Xtest) ^ 2)/ (N*P))))

sum(Beta^2)

resNew <- Y - betaHatNew %*% X
Mnew <- resNew %*% t(resNew) / N


plot(Beta, betaHatNew)
plot(Beta, betaHat)
cor.test(Beta, betaHatNew)
cor.test(Beta, betaHat)

## Compare covariance estimates
SigmaInv <- V %*% solve(psi) %*% t(V) + Vperp %*% solve(psi.0) %*% t(Vperp)
SigmaHat <- Mnew 

steins_loss_env <- tr(SigmaHat %*% SigmaInv) - log(det(SigmaHat %*% SigmaInv)) - P

SigmaHat <- M

steins_loss_mle <- tr(SigmaHat %*% SigmaInv) - log(det(SigmaHat %*% SigmaInv)) - P

print(sprintf("Stein MLE: %f", steins_loss_mle))
print(sprintf("Stein ENV: %f", steins_loss_env))

