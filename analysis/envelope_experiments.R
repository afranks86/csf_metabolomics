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
P <- 50
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

res <- (Y - betaHat %*% X)
rankM <- getRank(res)
rankMU <- getRank(Y)

MU <- Y %*% t(Y) / N + 0.01*diag(P)
## MUeigen <- eigen(MU)
## MU <- MUeigen$vectors[, 1:rankMU] %*% diag(MUeigen$values[1:rankMU]) %*% t(MUeigen$vectors[, 1:rankMU])
## MUinv <- MUeigen$vectors[, 1:rankMU] %*% diag(1/MUeigen$values[1:rankMU]) %*% t(MUeigen$vectors[, 1:rankMU])

M <- res %*% t(res) / N
## Meigen <- eigen(M)
## M <- Meigen$vectors[, 1:rankM] %*% diag(Meigen$values[1:rankM]) %*% t(Meigen$vectors[, 1:rankM])

## M <- M + 0.01*diag(nrow(M))
## MUinv <- M + 0.01*diag(nrow(MUinv))

MUinv <- solve(MU)

F <- function(V, M, MUinv) {
  determinant(t(V) %*% M %*% V, logarithm = TRUE)$modulus +
    determinant(t(V) %*% MUinv %*% V, logarithm = TRUE)$modulus
}

dF <- function(V, M, MUinv) {
  2 * M %*% V %*% solve(t(V) %*% M %*% V) +
    2 * MUinv %*% V %*% solve(t(V) %*% MUinv %*% V)
}

## start <- Sys.time()
## Vfit <- optStiefel(function(V) F(V, M, MUinv),
##                    function(V) dF(V, M, MUinv),
##                    method="curvilinear",
##                    Vinit=rstiefel::rustiefel(P, S),
##                    verbose=TRUE,
##                    searchParams=list(tau=1/2, rho1=0.5, rho2=0.99),
##                    maxIters=10000,
##                    maxLineSearchIters=500)
## print(Sys.time()-start)

start <- Sys.time()
Vfit <- optStiefel(
  function(V) F(V, M, MUinv),
  function(V) dF(V, M, MUinv),
  method = "bb",
  Vinit = rustiefel(P, S),
  verbose = TRUE,
  maxIters = 1000,
  maxLineSearchIters = 20
)
print(Sys.time() - start)

## Vfit <- optStiefel(function(V) F(V, M, MUinv),
##                    function(V) dF(V, M, MUinv),
##                    method="curvilinear",
##                    Vinit=rustiefel(P, S),
##                    verbose=TRUE,
##                    maxIters=10000,
##                    searchParams=list(tau=1, rho1=0.1, rho2=0.99),
##                    maxLineSearchIters=200)

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

