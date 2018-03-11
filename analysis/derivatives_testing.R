library(madness)


## Shrinkage estimation
F <- function(V, etaHat) {

    Z <- Y %*% V
    LnInv <- solve(t(X) %*% X + L0)
    Bn <- LnInv %*% (t(X) %*% Z)
    Mn <- M0 + t(Z - X %*% Bn) %*% (Z - X %*% Bn) + t(etaHat - Bn) %*% Lam0 %*% (etaHat - Bn)

    return(determinant(Mn)$modulus)
    
}

dF <- function(V, etaHat) {

    madV <- madness(V)
    
    Z <- Y %*% madV
    LnInv <- solve(t(X) %*% X + L0)
    Bn <- LnInv %*% (t(X) %*% Z)
    Mn <- M0 + t(Z - X %*% Bn) %*% ( Z - X %*% Bn) + t(etaHat - Bn) %*% Lam0 %*% (etaHat - Bn)

    return(matrix(dvdx(log(det(Mn))), ncol=20))
    
}

X <- matrix(rnorm(32), nrow=2)
madX <- madness(X)
madX

tmp <- tr(solve(diag(16) - (t(madX) %*% madX)))

val(tmp)

matrix(dvdx(tmp), ncol=2)
tr <- function(X) { sum(diag(X)) }

F <- function(V, S) {
    log(tr(V %*% t(V) %*% S - S))
}

dF <- function(V, S) {
    U <- tr(t(V) %*% S %*% V) - tr(S)
    (2 / U) * S %*% V
}


V <- matrix(rnorm(128), ncol=4)
Y <- rmvnorm(100, sigma=diag(nrow(V)))
Vmad <- madness(V)

S <- t(Y) %*% Y
F(V, S)

start <- Sys.time()
res1 <- dF(V, S)
print(Sys.time()-start)

start <- Sys.time()
res2 <- matrix(dvdx(F(Vmad, S)), ncol=ncol(V))
print(Sys.time()-start)

U <- tr(S) - tr(t(V) %*% S %*% V)
