#!/usr/bin/Rscript
#  lib.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.22.2018

## Functions useful for multiple files in this repository.

#' Get the mean of the beta-binom dist for fixed params
bb.mean <- function(alpha, beta, n) {
    alpha * n / (alpha + beta)
}

#' Get the variance of the beta dist for fixed params
bb.var <- function(alpha, beta, n) {
    num <- n * alpha * beta * (alpha + beta + n)
    denom <- (alpha + beta)^2 * (alpha + beta + 1)
    return(num / denom)
}

#' Get Control Chart Limits
#' 
#' Get the limits at which a control chart will signal.
#' @param chart The control chart the signalling limits of which are desired.
#' @param n The sample size. Not all charts need this.
#' @return A tuple of the form c(lower_lim, upper_lim), both between 0 and 1.
get_lims <- function(chart, n = NULL) UseMethod("get_lims")

#' Calibrate In Control ARL
#'
#' Tune a control chart such that it has approximately the desired in sample Average Run Length (ARL) (if the true parameters are known. If they are esimated, this will be even more approximate).
#' @param chart The control chart to tune.
#' @param target_arl A positive scalar, the desired in control ARL.
#' @return The chart object passed to this function suitably modified to achieve target ARL.
cal_arl <- function(chart, target_arl, n, alpha, beta) UseMethod('cal_arl')

#' Estimate Control-Chart Parameters.
#'
#' Estimate the parameters for a control chart given data.
est_params <- function(chart, X, N, ...) UseMethod('est_params')

#' Moment Matching for Various Distributions
#'
#' Moment match various distributions to the beta-binomial, and back.
bb_mm <- function(target, alpha = NULL, beta = NULL, N = NULL, mu = NULL, 
                  sig = NULL) {
    if (target == 'laney') {
        rho <- alpha / (alpha + beta)
        m <- bb.mean(alpha, beta, N)
        v <- bb.var(alpha, beta, N)
        sig_p <- sqrt(rho * (1 - rho))
        sig_pi <- sig_p / sqrt(N)
        sig_z <- sqrt(v) / (N * sig_pi)
        return(list(rho = rho, sig_p = sig_p, sig_z = sig_z))

    } else if (target == 'betabinom') {

        alpha <- N * (mu* (N - mu) - sig^2) / ((1 + (N-mu) / mu) * 
                                               (N*sig^2 - mu * (N - mu)))
        beta <- alpha * (N - mu) / mu

        if (alpha < 0 || beta < 0) {
            warning("alpha or beta is less than zero (desired mean/variance combination not possible at this sample size)")
        }

        return(list(alpha = alpha, beta = beta))
    } else {
        stop("Target not implemented.")
    }
}

#' Plot a Control Chart Given Data
#'
#' Plot a control chart and its limits given data
plot.control_chart <- function(chart, X, N, force_01 = TRUE) {
    m <- length(X)
    if ((length(X) != length(N)) && length(N) != 1) {
        stop("'N' should either be a scalar, indicating constant sample size, or a vector of the same length of 'X'")
    }
    
    #Get the chart's limits for all the data
    lims <- lapply(N, function(n) get_lims(chart, n))
    l_lims <- sapply(lims, function(lim) lim[1])
    u_lims <- sapply(lims, function(lim) lim[2])

    #Get sample proportions
    rho_hats <- X / N

    #Get the boundaries of our chart's y-axis
    if (force_01) {
        ylim <- c(0,1)
    } else {
        l <- min(c(l_lims, rho_hats))
        u <- max(c(u_lims, rho_hats))
        ylim <- c(l, u)
    }

    #Figure out which ones are out of bounds:
    is_oob <- (rho_hats < l_lims | rho_hats > u_lims)
    cols <- c('black', 'red')

    #Make our chart!
    plot(NA, NA, xlim = c(1, m), ylim = ylim, main = class(chart), 
         xlab = "Run", ylab = "Observed Proportion")
    points(1:m, X / N, col = cols[is_oob + 1])
    points(1:m, l_lims, type = 'l', col = 'red')
    points(1:m, u_lims, type = 'l', col = 'red')
}

#PMF functions
bb.pmf <- function(k,n,a,b, lspace = FALSE) {
    l <- lchoose(n,k) + lbeta(k + a, n - k + b) - lbeta(a, b)
    if (lspace) {
        return(l)
    } else {
        return(exp(l))
    }
}

# cdf func
#TODO: Closed form cdf is possible with use of special functions
bb.cdf <- function(k, n, a, b) {
    sum(sapply(0:k, function(i) bb.pmf(i, n, a, b)))
}

#' Maximum Likelihood Estimation for the Beta-Binomial Distribution.
#'
#' Numerically optimize the likelihood of the beta-binomial distribution given some observations. The inital guess is determined via method of moments.
#' @param X An integer vector, the observed counts
#' @param N Either an integer vector of length length(X), or a scalar, indicating the sample size for all trials or for each trial, repectively.
bb.mle <- function(X, N) {
    #Check inputs
    m <- length(X)
    if (length(N) == 1) {
        N <- rep(N, m)
    }
    if (length(N) != m) {
        stop("'N' should either be a scalar, indicating constant sample size, or a vector of the same length of 'X'")
    }

    # Compose our cost function, the negative log likelihood
    nllik <- function(params) -sum(sapply(1:m, function(i) 
                          bb.pmf(X[i], N[i], params[1], params[2], lspace = TRUE)))

    # Initialize using method of moments estimators.
    mu <- mean(X)
    sig <- sd(X)
    init <- unlist(bb_mm('betabinom', mu = mu, sig = sig, N = mean(N)))
    if (min(init) < 1e-3) {
        init <- pmax(1e-3, init)
    }

    # Do the optim
    ret <- optim(init, nllik)$par
    return(list(alpha = ret[1], beta = ret[2]))
}

#' Calculate Beta-Binomial Quantiles.
#'
#' Calculate quantiles of the Beta-Binomial distribution in a numerically stable manner.
#' @param p The desired probability level, a scalar in [0,1]
#' @param n The sample size
#' @param a The first positive shape parameter, also referred to as 'alpha'
#' @param b The second positive shape parameter, also referred to as 'beta'
#' @return The corresponding quantile, an integer in {0, ..., n}
bb.quantile <- function(p, n, a, b) {
    #cum_dists <- sapply(0:n, function(i) bb.cdf(i, n, alpha, beta))
    cum_dists <- cumsum(sapply(0:n, function(i) bb.pmf(i, n, a, b)))
    ret <- sapply(p, function(i) sum(cum_dists < i, na.rm = TRUE))
    return(ret)
}

