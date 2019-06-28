#!/usr/bin/Rscript
#  re_beta_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.23.2018

## A chart to model binomial data with inherent batch-to-batch variation using the 
## beta-binomial model.

library(hypergeo)

#' The Beta Random Effects Control Chart
#'
#' A control chart for monitoring proportions with legitimate batch-to-batch variation.
#'
#' Creates a control chart based on the binomial distribution with success probability assumed to vary within a reasonable range from batch to batch. If X and N are specified, this method conducts maximum likelihood estimation to determine suitable values of alpha and beta (typical useage will be in phase I monitoring). If alpha and beta are specified, these parameters will not be estimated.
#' @param X An integer vector representing the observed data.
#' @param N An integer vector representing the sample sizes. N and X should be of the same length.
#' @param alpha A positive scalar shape parameter.
#' @param beta A positive scalar shape parameter.
#' @param l_p A scalar indicating the lower quantile at which to signal. May be set to match a specified ARL using cal_arl.
#' @param u_p A scalar indicating the upper quantile at which to signal. May be set to match a specified ARL using cal_arl.
#' @param def_inf A string, one of {'mle', 'rmm'} specifying the default inference method. 'mle' is Maximum Likelihood, 'rmm' is Robust Method of Moments.
#' @return An object of class control.chart representing the chart.
re_beta_chart <- function(alpha = NA, beta = NA, l_p = 0.025, u_p = 0.975, def_inf = 'rmm') {
    ## Check inputs

    #Create the object.
    ret <- list()
    ret$alpha <- alpha
    ret$beta <- beta
    ret$l_p = l_p
    ret$u_p = u_p
    ret$def_inf <- def_inf
    class(ret) <- c('re_beta_chart', 'control_chart')
    return(ret)
}

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


#' Estimate Beta-Binomial Chart parameters
#'
#' Estimate the parameters for a control chart given data.
#' @param chart The chart object which estimation is to be performed for.
#' @param X A integer vector of observed counts.
#' @param N Either an integer vector of the same length as X, indicating sample sizes, or a scalar integer for constant sample size.
#' @param type A string specifying the kind of inference, 'mle' for maximum likelihood and 'rmm' for robust moment matching (see <paper name> for details). Defaults to the chart's default.
#' @return The chart with parameters modified. 
est_params.re_beta_chart <- function(chart, X, N, type) {
    if ((length(X) != length(N)) && length(N) != 1) {
        stop("'N' should either be a scalar, indicating constant sample size, or a vector of the same length of 'X'")
    }

    if (missing(type)) {
        type <- chart$def_inf
    }

    if (type == 'rmm') {
        mu <- mean(X)
        sig <- mean(abs(diff(X))) / 1.128

        ret <- bb_mm('betabinom', N = mean(N), mu = mu, sig = sig)
        chart$alpha <- ret$alpha
        chart$beta <- ret$beta
    } else if (type == 'mle') {
        ests <- bb.mle(X, N)
        chart$alpha <- ests$alpha
        chart$beta <- ests$beta
    } else {
        stop("Inference type not implemented, choose either 'rmm' or 'mle'")
    }

    return(chart)
}

#' Calibrate Beta RE Chart In Control ARL
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
