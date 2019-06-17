kegg_map <- read_csv("~/course/ND_Metabolomics/data/kegg_map.csv", trim_ws=TRUE) %>% dplyr::select(1:6)

n2k <- function(nms) {
    kegg_map %>% filter(METABOLITE %in% nms)
}

k2n <- function(k) {

}

normalize_density <-  function(x) {
    n <- length(x)
    delta <- 1 / (4 * n^(1/4) * sqrt(pi*log(n)))
    
    uniform_transf <- n / (n+1) * ecdf(x)(x)

    windsorized <- ifelse(uniform_transf < delta, delta, uniform_transf)
    windsorized <- ifelse(windsorized > 1- delta, 1 - delta, windsorized)

    qnorm(windsorized)

}

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


fit_lm <- function(df, predictors, type="lm") {
    form <- as.formula(paste("Abundance", paste(predictors, collapse = " + "), 
                             sep = " ~ "))

    if(sum(is.na(df$Abundance)) > 50)
        return(NULL)
    else {
        if(type=="lm") {
            return(summary(lm(form, data=df))$coefficients)
        }
        else if(type == "rq") {

            rq_fit <- rq(form, data=df, tau=0.5)
            coefs <- summary(rq_fit, se="boot")$coefficients
            colnames(coefs)[1] <- "Estimate" ## to match lm
            coefs

        }
    
        
    }
}

