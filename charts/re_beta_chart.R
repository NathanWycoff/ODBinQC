#!/usr/bin/Rscript
#  re_beta_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.23.2018

## A chart to model binomial data with inherent batch-to-batch variation using the 
## beta-binomial model.

require('hypergeo')

#A numerically stable choose function

#PMF functions
bb.pmf <- function(k,n,a,b) {
    l <- lchoose(n,k) + lbeta(k + a, n - k + b) - lbeta(a, b)
    exp(l)
}


#Calculate the cdf using the generalized hypergeometric function
bb.cdf.hypergeom <- function(k, n, a, b) {
    ## first try to do it using the hypergeometric function
    num <- lbeta(b + n - k - 1, a + k + 1)

    #prepare the generalized hypergeomtric function
    u <- c(1, a + k + 1, k - n + 1)
    l <- c(k + 2, k + 2 - b - n)
    num <- num + log(abs(genhypergeo(u, l, 1)))

    denom <- lbeta(a, b) + lbeta(n - k, k + 2) * (n + 1)

    ret <- 1 - exp(num - denom)


    return(ret)
}

#Calculate the cdf by summing pmf's
bb.cdf.pmf <- function(k, n, a, b) {
    ret <- sum(sapply(0:k, function(i) bb.pmf(i, n, a, b)))
    return(ret)
}

bb.cdf <- function(k,n,a,b) {
    # First try to do it using the hypergeometric function
    num <- lbeta(b + n - k - 1, a + k + 1)

    #prepare the generalized hypergeomtric function
    u <- c(1, a + k + 1, k - n + 1)
    l <- c(k + 2, k + 2 - b - n)
    num <- num + log(abs(genhypergeo(u, l, 1)))

    denom <- lbeta(a, b) + lbeta(n - k, k + 2) + log(n + 1)

    ret <- 1 - exp(num - denom)

    # If that didn't work, let's sum up a bunch of pmf's
    if (is.na(ret) || is.nan(ret)) {
        ret <- sum(sapply(0:k, function(i) bb.pmf(i, n, a, b)))
    }

    return(ret)
}

bb.quantile <- function(p, n, a, b) {
    #cum_dists <- sapply(0:n, function(i) bb.cdf(i, n, alpha, beta))
    cum_dists <- cumsum(sapply(0:n, function(i) bb.pmf(i, n, alpha, beta)))
    ret <- sapply(p, function(i) sum(cum_dists < i, na.rm = TRUE))
    return(ret)
}


#Calculate quantiles for the beta-binomial distribution.


re.beta.chart <- function(alpha, beta) {
    ret <- list()
    ret$get.lims <- function(N, ps = c(0.025, 0.975)) {
        sapply(ps, function(p) bb.quantile(p, N, alpha, beta))
    }
    return(ret)
}
