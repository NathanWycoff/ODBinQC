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

#' Get Beta RE Chart Limits
#' 
#' Get the limits at which a Beta random effects chart will signal.
#' @param chart The control chart the signalling limits of which are desired.
#' @param n The sample size. Required for this chart.
#' @return A tuple of the form c(lower_lim, upper_lim), both between 0 and 1.
get_lims.re_beta_chart <- function(chart, n) {
    l_lim <- bb.quantile(chart$l_p, n, chart$alpha, chart$beta) / n
    u_lim <- bb.quantile(chart$u_p, n, chart$alpha, chart$beta) / n
    return(c(l_lim, u_lim))
}

#' Calibrate In Control ARL
#'
#' Tune a control chart such that it has approximately the desired in sample Average Run Length (ARL) (if the true parameters are known. If they are esimated, this will be even more approximate).
#' @param chart The control chart to tune.
#' @param target_arl A positive scalar, the desired in control ARL.
#' @return The chart object passed to this function suitably modified to achieve target ARL.
cal_arl <- function(chart, target_arl, n, alpha, beta) UseMethod('cal_arl')

#' Calibrate In Control ARL
#'
#' Tune a control chart such that it has approximately the desired in sample Average Run Length (ARL) (if the true parameters are known. If they are esimated, this will be even more approximate).
#' @param chart The control chart to tune.
#' @param target_arl A positive scalar, the desired in control ARL.
#' @param ... Additional arguments ignored, but available for compatibility with siblings of this function.
#' @return The chart object passed to this function suitably modified to achieve target ARL.
cal_arl.re_beta_chart <- function(chart, target_arl, ...) {
    #Just ignores extra args
    target_p <- 1 / target_arl
    chart$l_p <- 0 + target_p / 2
    chart$u_p <- 1 - target_p / 2
    return(chart)
}

#' Estimate Control-Chart Parameters.
#'
#' Estimate the parameters for a control chart given data.
est_params <- function(chart, X, N) UseMethod('est_params')

#' Moment Matching for Various Distributions
#'
#' Moment match various distributions to the beta-binomial, and back.
bb_mm <- function(target, alpha = NULL, beta = NULL, N = NULL, mu = NULL, sig = NULL,
                  inits = 1e2) {
    if (target == 'laney') {
        rho <- alpha / (alpha + beta)
        m <- bb.mean(alpha, beta, N)
        v <- bb.var(alpha, beta, N)
        sig_p <- sqrt(rho * (1 - rho))
        sig_pi <- sig_p / sqrt(N)
        sig_z <- sqrt(v) / (N * sig_pi)
        return(list(rho = rho, sig_p = sig_p, sig_z = sig_z))

    } else if (target == 'betabinom') {

        c1 <- (N - mu) / mu
        alpha <- N * (mu* (N - mu) - sig^2) / ((1 + (N-mu) / mu) * 
                                               (sig^2*N - mu * (N - mu)))
        beta <- alpha * (N - mu) / mu

        return(list(alpha = alpha, beta = beta))
    } else {
        stop("Target not implemented.")
    }
}
