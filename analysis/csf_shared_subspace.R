library(tidyverse)
library(magrittr)
library(made4)
library(ggthemes)
library(ggridges)
library(gbm3)

residuals <- read_csv("residual_data.csv")

rlong <- residuals %>% gather(key = Metabolite, value = Abundance, -one_of("Id", "Type", "Gender", "Age", "APOE", "Batch", "Index"))

normalize_density <- function(x) {

  ## get winsorized estiamte of F
  n <- length(x)
  cumsum((1:n) / n)

  ## transform to uniform and then normalize
}

rlong %>% ggplot(aes(x = Abundance, y = factor(Type), fill = Type)) + geom_density_ridges()

types <- c("CO", "PD", "AD")

ngroups <- 3
nvec_table <- residuals %>% filter(Type %in% types) %>% group_by(Type) %>% summarise(n = n())
nvec <- nvec_table$n
names(nvec) <- nvec_table$Type
nvec <- nvec[types]



residuals_list <- Slist <- list()
for (type in types) {
  group_resid <- residuals %>% filter(Type == type) %>% select(-(1:7)) %>% as.matrix()
  residuals_list[[type]] <- group_resid
  Slist[[type]] <- t(group_resid) %*% group_resid
}

P <- ncol(residuals_list[[1]])


pooled_residuals <- do.call(rbind, residuals_list)
median(svd(pooled_residuals)$d ^ 2 / sum(nvec))
samp <- rmvnorm(sum(nvec), mean = rep(0, P), sigma = diag(median(pooled_evals), P))

plot(
  (1:nvec[1]) / nvec[1], svd(residuals_list[[1]] / sqrt(nvec[1]))$d ^ 2, pch = 19, cex = 0.75, type = "l", col = 1,
  ylim = c(0, 1)
)
lines((1:nvec[2]) / nvec[2], svd(residuals_list[[2]] / sqrt(nvec[2]))$d ^ 2, pch = 19, cex = 0.75, col = 2)
lines((1:nvec[3]) / nvec[3], svd(residuals_list[[3]] / sqrt(nvec[3]))$d ^ 2, pch = 19, cex = 0.75, col = 3)
lines((1:P) / P, (svd(samp / sqrt(sum(nvec)))$d ^ 2)[1:P], pch = 19, cex = 0.75, col = "blue")

pooled_evals <- svd(pooled_residuals)$d ^ 2 / sum(nvec)
samp <- rmvnorm(sum(nvec), mean = rep(0, P), sigma = diag(median(pooled_evals), P))
samp_evals <- svd(samp)$d ^ 2 / sum(nvec)

plot(pooled_evals)
lines(samp_evals)

ad_evals <- svd(residuals_list[[3]])$d ^ 2 / nvec[[3]]
samp <- rmvnorm(nvec[3], mean = rep(0, P), sigma = diag(median(ad_evals), P))
samp_evals <- svd(samp)$d ^ 2 / sum(nvec)

plot(ad_evals)
lines(samp_evals)

pd_evals <- svd(residuals_list[[2]])$d ^ 2 / nvec[[2]]
samp <- rmvnorm(nvec[2], mean = rep(0, P), sigma = diag(median(pd_evals), P))
samp_evals <- svd(samp)$d ^ 2 / sum(nvec)

plot(pd_evals)
lines(samp_evals)



## compute the measurement variance
isoVar <- sapply(1:ngroups, function(i) median(apply(residuals_list[[i]], 2, var)))
weightsList <- sapply(1:ngroups, function(i) {
  evals <- svd(residuals_list[[i]] / sqrt(nvec[i]))$d ^ 2
  dp <- ifelse(evals / isoVar[i] <= (1 + sqrt(P / nvec[i])), 0, sqrt((1 - P / nvec[i] / (evals / isoVar[i] - 1) ^ 2) / (1 + P / nvec[i] / (evals / isoVar[i] - 1))))
  weights <- 1 / (1 - dp) - 1
  weights
})

source("~/course/shared-subspace/helper.R")
source("~/course/shared-subspace/fit-subspace.R")
library(rstiefel)

S <- getRank(pooled_residuals)
R <- 2


Vinit <- svd(do.call(cbind, lapply(1:ngroups, function(k) svd(t(residuals_list[[k]]))$u[, 1:S])))$u[, 1:S]
Vinit <- svd(t(residuals_list[["CO"]]))$u[, 1:S]
PrecVec <- rep(5, ngroups)
Vinit <- Vinit[, sample(1:ncol(Vinit))]
EMFit <- subspaceEM(Slist, P = P, S = S, R = R, nvec = nvec, Vstart = Vinit, verbose = TRUE, rho1 = 1e-5, rho2 = 1 - 1e-5)

## EMFit <- subspaceEM(Slist, P=P, S=S, R=R, nvec=nvec, Vstart=rstiefel::rustiefel(P, S), verbose=TRUE, rho1=1e-5, rho2=1-1e-5)
Vinit <- EMFit$V

evalRatios <- sapply(1:ngroups, function(k) {
  numer <- sum(eigen(t(Vinit[, 1:S]) %*% Slist[[k]] %*% Vinit[, 1:S])$values) / nvec[k]
  denom <- sum(svd(residuals_list[[k]])$d[1:S] ^ 2) / nvec[k]
  correction <- (denom - 1 / EMFit$PrecVec[k] * S * P / nvec[k]) / denom
  (numer / denom) / correction
})

par(mar = c(6.1, 4.1, 4.1, 2.1))
barplot(evalRatios, ylim = c(0, 1), xlab = "", main = "", space = 0.2, cex.axis = 1.5, col = "#646464", col.axis = "black", col.main = "black", border = NA)

OmegaList <- Ulist <- list()
for (k in 1:ngroups) {
  eigK <- eigen(diag(S) - EMFit$PhiList[[k]] / EMFit$PrecVec[k])

  OmegaList[[k]] <- eigK$values
  Ulist[[k]] <- Vinit %*% eigK$vectors
}
s2vec <- 1 / EMFit$PrecVec

initSS <- list(V = Vinit, Ulist = Ulist, OmegaList = OmegaList, s2vec = s2vec)
samples <- fitSubspace(
  P, S, R, Q = S - R, Slist,
  nvec, ngroups, init = initSS,
  niters = 1000, draw = c(V = FALSE)
)

with(
  samples,
  posteriorPlot(
    Osamps[1:2, 1:2, , ], omegaSamps[1:2, , ],
    s2samps, nsamps = 100, groups = 1:ngroups,
    probRegion = 0.90, hline = NULL,
    col = 1:3, lty = rep(1, 3), logRatio = TRUE,
    plotPoints = FALSE
  )
)

legend("topright", legend = types, col = 1:3, lty = rep(1, 3), bty = "n", lwd = 3)


############### Biplot ##########################

contourColors <- c("blue", "red", "green")
mag <- apply(EMFit$V[, 1:2], 1, function(x) sqrt(sum(x ^ 2)))
indices <- which(mag > quantile(mag, 0.90))

Vsub <- EMFit$V[indices, ]
topright <- which(Vsub[, 1] > 0 & Vsub[, 2] > 0)
topleft <- which(Vsub[, 1] < 0 & Vsub[, 2] > 0)
bottom <- which(Vsub[, 2] < 0)

plotlim <- 1.1 * max(abs(Vsub[, 1:2]))
plot(
  EMFit$V[, 1:2], col = "light grey", pch = 19, cex = 0.5,
  xlim = c(-plotlim, plotlim), ylim = c(-plotlim, plotlim), xlab = "V1", ylab = "V2"
)
points(
  Vsub[topright, 1], Vsub[topright, 2], xlab = "V1", ylab = "V2", cex = 1.5,
  xlim = c(-plotlim, plotlim), ylim = c(-plotlim, plotlim), pch = 15, col = "black"
)
points(
  Vsub[topleft, 1], Vsub[topleft, 2], xlab = "V1", ylab = "V2", cex = 1.5,
  xlim = c(-plotlim, plotlim), ylim = c(-plotlim, plotlim), pch = 17, col = "black"
)
points(
  Vsub[bottom, 1], Vsub[bottom, 2], xlab = "V1", ylab = "V2", cex = 1.5,
  xlim = c(-plotlim, plotlim), ylim = c(-plotlim, plotlim), pch = 16, col = "black"
)
# text(EMFit$V[indices, 1], EMFit$V[indices, 2], labels=rownames(EMFit$V)[indices],
#     xlab="V1", ylab="V2", cex=0.5, pos=4)

abline(h = 0, v = 0)
legend(
  "topright", legend = unique(types)[1:3], lty = 1, lwd = 2, col = contourColors,
  cex = 1.3, bty = "n", ncol = 1, title = "Types"
)

library(car)
for (k in 1:3) {
  eigK <- eigen(solve(EMFit$PhiList[[k]][1:2, 1:2]))
  lambda <- eigK$values
  print(lambda)
  evecs <- eigK$vectors

  maxIndex <- which.max(lambda)
  lamRatio <- lambda[maxIndex] / lambda[-maxIndex]
  angle <- atan(evecs[2, maxIndex] / evecs[1, maxIndex])
  print(angle)

  ## arrows(0, 0, evecs[1, 1]*.05, evecs[2, 1]*.05, col=k, lwd=3)
  ## arrows(0, 0, -evecs[1, 1]*.05, -evecs[2, 1]*.05, col=k, lwd=3)
  ellipse(center = c(0, 0), shape = evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2]) / lambda[maxIndex], radius = 0.05, col = contourColors[k], center.cex = 000)
  ellipse(center = c(0, 0), shape = evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2]) / lambda[maxIndex], radius = 0.125, col = contourColors[k], center.cex = 000)
  ellipse(center = c(0, 0), shape = evecs[, 1:2] %*% diag(lambda[1:2]) %*% t(evecs[, 1:2]) / lambda[maxIndex], radius = 0.2, col = contourColors[k], center.cex = 000)
  points(angle, lamRatio, col = k, pch = 19, cex = 2)
}
