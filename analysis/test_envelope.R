library(microbenchmark)
source("~/course/rstiefel/R/opt.stiefel.R")
library(mvtnorm)
library(mvnfast)

## Generate some fake 
P <- 100
N <- 100
S <- 2

## Generate envelope
V <- rstiefel::rustiefel(P, S)
Vperp <- rstiefel::NullC(V)

## regression parameters
## eta <- rnorm(S, 0, 10)
eta <- c(-10, 10)
Beta <- V %*% eta

## covariates
X <- rnorm(N)

## covariance matrices
## psi <- diag(sort(rexp(S, 1/20), decreasing=TRUE))
psi <- diag(c(20, 10))
# psi.0 <- diag(rep(1000, P-S))
psi.0 <- diag(sort(rexp(P-S, 1/10), decreasing=TRUE))

Sigma <- V %*% psi %*% t(V) + Vperp %*% psi.0 %*% t(Vperp)

Y <- sapply(1:N, function(i) rmvn(1, Beta %*% X[i], Sigma))

## Plot Y colored by X on VV^T and a random subspace UU^T
rbpal <- colorRampPalette(c('red','white', 'blue'))
cols <- rbpal(10)[as.numeric(cut(X, breaks = 10))]

par(mfrow=c(1,2))
plot(t(Y) %*% V, pch=19,
     col=cols, main="Subspace Projection")
plot(t(Y) %*% rstiefel::rustiefel(P, S),
     pch=19, col=cols, main="Random Projection")


betaHat <- Y %*% X %*% solve(t(X) %*% X)
plot(Beta, betaHat)
sum((Beta - betaHat)^2)

res <- (Y - betaHat %*% X)

MU <- Y %*% t(Y)/N 
M <- res %*% t(res)/N
MUinv <- solve(MU)

F <- function(V, M, MUinv) {
    determinant(t(V) %*% M %*% V, logarithm=TRUE)$modulus +
    determinant(t(V) %*% MUinv %*% V, logarithm=TRUE)$modulus
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
Vfit <- optStiefel(function(V) F(V, M, MUinv),
                   function(V) dF(V, M, MUinv),
                   method="bb",
                   Vinit=rustiefel(P, S),
                   verbose=TRUE,
                   maxIters=1000,
                   maxLineSearchIters=20)
print(Sys.time()-start)

## Vfit <- optStiefel(function(V) F(V, M, MUinv),
##                    function(V) dF(V, M, MUinv),
##                    method="curvilinear",
##                    Vinit=rustiefel(P, S),
##                    verbose=TRUE,
##                    maxIters=10000,
##                    searchParams=list(tau=1, rho1=0.1, rho2=0.99),
##                    maxLineSearchIters=200)

Yproj <- t(Vfit) %*% Y
par(mfrow=c(1,2))

plot(t(Y) %*% Vfit, pch=19,
     col=cols, main="Subspace Projection")
plot(t(Y) %*% rustiefel(P, S),
     pch=19, col=cols, main="Random Projection")

etaHat <- Yproj %*% X %*% solve(t(X) %*% X)
betaHatNew <- Vfit %*% etaHat

resNew <- Y - betaHatNew %*% X
Mnew <- resNew %*% t(resNew)/N

median(abs(Beta-betaHat))
median(abs(Beta - betaHatNew))

sum((Beta-betaHat)^2)
sum((Beta - betaHatNew)^2)

plot(Beta, betaHatNew)
plot(Beta, betaHat)

SigmaInv <- V %*% solve(psi) %*% t(V) + Vperp %*% solve(psi.0) %*% t(Vperp)
SigmaHat <- Mnew

tr(SigmaHat %*% SigmaInv) - log(det(SigmaHat %*% SigmaInv)) - P

SigmaInv <- V %*% solve(psi) %*% t(V) + Vperp %*% solve(psi.0) %*% t(Vperp)
SigmaHat <- M

tr(SigmaHat %*% SigmaInv) - log(det(SigmaHat %*% SigmaInv)) - P

## Check covariance similarity
