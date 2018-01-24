## Testing envelope models
## Envelope models do much worse when sighat is a poor estimator of sigma
##

library(microbenchmark)
source("~/course/rstiefel/R/opt.stiefel.R")
library(mvtnorm)
library(mvnfast)
library(rstiefel)

## Generate some fake data
P <- 300
N <- 100
S <- 20

## Generate envelope
V <- rstiefel::rustiefel(P, S)
Vperp <- rstiefel::NullC(V)

## regression parameters

## eta small relative to evals gives based improvements in MSE

## eta <- matrix(c(-1, 1), ncol=1)

eta <- matrix(c(-10, 10, 10, 2), nrow=2, ncol=2)
eta <- matrix(rnorm(S^2, 0, 1), nrow=S)
Beta <- V %*% eta

## covariates

X <- matrix(rmvnorm(N, rep(0, S), sigma=diag(rep(5^2, S))), ncol=S)
#X <- rmvnorm(N, mean=c(0,0), sigma=diag(c(5,5)))

## covariance matrices
psi <- diag(sort(rexp(S, 1 / 100), decreasing = TRUE))
## psi <- diag(c(10, 10))
## psi.0 <- diag(rep(100, P-S))
psi.0 <- diag(sort(rexp(P - S, 1 / 100), decreasing = TRUE))

ecdf(diag(psi.0))(diag(psi))

Sigma <- V %*% psi %*% t(V) + Vperp %*% psi.0 %*% t(Vperp)

Y <- sapply(1:N, function(i) rmvn(1, Beta %*% X[i, ], Sigma))

betaHat <- Y %*% X %*% solve(t(X) %*% X)
plot(Beta, betaHat)
sum((Beta - betaHat) ^ 2)

res <- Y - betaHat %*% t(X)

gam <- svd(betaHat)$u
gam0 <- rstiefel::NullC(gam)

print(sprintf("Gamma similarity: %f", tr(V %*% t(V) %*% gam %*% t(gam)) / S))

Yproj <- t(gam) %*% res
psihat <- Yproj %*% t(Yproj) / N

Yproj0 <- t(gam0) %*% Y
psihat.0 <- Yproj0 %*% t(Yproj0) / N

SigmaHat <- Y %*% t(Y) / N + diag(rep(1, P))
SigmaHat2 <- gam %*% psihat %*% t(gam) + gam0 %*% psihat.0 %*% t(gam0) + diag(rep(1, P))

SigmaInv <- V %*% solve(psi) %*% t(V) + Vperp %*% solve(psi.0) %*% t(Vperp)

loss1 <- tr(SigmaHat %*% SigmaInv) - log(det(SigmaHat %*% SigmaInv)) - P

loss2 <- tr(SigmaHat2 %*% SigmaInv) - log(det(SigmaHat2 %*% SigmaInv)) - P

print(sprintf("Loss ratio: %f", loss1 / loss2))
