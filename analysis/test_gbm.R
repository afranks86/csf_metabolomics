library(gbm3)

n <-  1000
for(i in 1:100) {
    print(i)
    x <- 1:n
    y <- 3 * x + rexp(n) * ifelse(runif(n) < 0.5, -1, 1)
    df <- data.frame(x=x, y=y)
    fit_cv <- gbm(y ~ x, data=df, distribution="laplace", n.trees=10000, cv.folds=10, verbose=FALSE)
    opt_iters <- gbm.perf(fit_cv, method="cv")
    best_fit <- predict(fit_cv, df, opt_iters)

}
