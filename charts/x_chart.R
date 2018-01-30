#!/usr/bin/Rscript
#  charts/x_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.23.2018

#' The X Chart.
#' 
#' Monitor proportion data approximated by a normal distribution, treated as individual observations.
#'
#' @param X An integer vector of observed values.
#' @param N An integer vector of sample sizes.
#' @param mu The true overall mean. Estimated using mean(X) if not specified.
#' @param sig The true overall standard deviation. Estimated using moving ranges if not specified.
#' @param k A positive scalar indicating how many standard deviation away from the mean the symmetric signalling limits should be. 
#' @return An object of class control.chart representing the chart.
x_chart <- function(mu, sig, k = 3) {
    ## Create the return object.
    ret <- list()
    ret$mu <- mu
    ret$sig <- sig
    ret$k <- k
    class(ret) <- c('x_chart', 'control_chart')
    return(ret)
}

#' Get X Chart Limits
#' 
#' Get the limits at which an X chart will signal.
#' @param chart The control chart the signalling limits of which are desired.
#' @param n The sample size. Not used for this chart.
#' @return A tuple of the form c(lower_lim, upper_lim), both between 0 and 1.
get_lims.x_chart <- function(chart, n = NULL) {
    l_lim <- pmax(0, chart$mu - chart$k * chart$sig)
    u_lim <- pmin(1, chart$mu + chart$k * chart$sig)
    return(c(l_lim, u_lim))
}

#' Estimate X Chart Parameters.
#'
#' Estimate the probability of success for a P chart
#' @param chart The chart object which estimation is to be performed for.
#' @param X A integer vector of observed counts.
#' @param N Either an integer vector of the same length as X, indicating sample sizes, or a scalar integer for constant sample size.
#' @return The chart with parameter modified. 
est_params.x_chart <- function(chart, X, N) {
    if ((length(X) != length(N)) && length(N) != 1) {
        stop("'N' should either be a scalar, indicating constant sample size, or a vector of the same length of 'X'")
    }

    mu <- mean(X / N)

    if (length(X) == 1) {
        stop("Cannot estimate standard deviation with only 1 observation")
    }
    #Calculate sd using range
    ranges <- abs(diff(X / N))
    r_bar <- mean(ranges)
    sig <- r_bar / 1.128

    #Store the vals
    chart$mu <- mu
    chart$sig <- sig

    return(chart)
}

#' Calibrate X Chart In Control ARL
#'
#' Tune an X chart such that it has approximately the desired in sample Average Run Length (ARL) (if the true parameters are known. If they are esimated, this will be even more approximate), under the assumption that data are being drawn form a Beta-Binomial model. 
#' @param chart The control chart to tune.
#' @param target_arl A positive scalar, the desired in control ARL.
#' @param n The sample size, or, in the case of varying sample sizes, the expected sample size. If an expected sampe size is passed, it will introduce more error and there will be a larger discrepency between desired and achieved ARL.
#' @param alpha A positive scalar, the first shape parameter of the beta-binom distribution.
#' @param beta A positive scalar, the second shape parameter of the beta-binom distibution.
#' @return The chart object passed to this function with k, the number of standard deviation from the mean after which to signal, suitably modified to achieve the desired in control ARL.
#TODO: It may be possible to avoid specifying alpha and beta, and instead infer them based on the mu and sig values from the chart. I'm not sure if this is desireable yet.
cal_arl.x_chart <- function(chart, target_arl, n, alpha, beta) {
    #Define a cost function returning descrepancy between target and specified.
    target_p <- 1 / target_arl
    cost <- function(k) {
        chart$k <- k
        lims <- get_lims(chart) * n
        p <- bb.cdf(floor(lims[1]), n, alpha, beta) + 
            (1-bb.cdf(ceiling(lims[2]), n, alpha, beta))
        return((p - target_p)^2)
    }
    chart$k <- optimize(cost, c(0, 10))$minimum
    return(chart)
}
