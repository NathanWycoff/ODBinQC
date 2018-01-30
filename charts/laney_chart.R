#!/usr/bin/Rscript
#  charts/laney_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2018

## Make a p' chart as described by Laney in Improved Control Charts for Attributes

#' The Laney Control Chart
#'
#' The Laney p' Chart to monitor binary valued quality characteristics with legitimate batch-to-batch variation.

#' Creates a control chart based on the binomial distribution with success probability assumed to vary within a reasonable range from batch to batch. If X and N are specified, this method estimates the mean and variance using the sample mean and the sample ranges.
#' @param rho The true overall probability of success
#' @param sig_z The true variation on the z-scale; the batch-to-batch variation.
#' @param sig_p The true variation due to binomial draws, this is sqrt(n_i) * sig_pi from Laney's paper.
#' @param k A positive scalar indicating how many standard deviation away from the mean the symmetric signalling limits should be. 
#' @return An object of class control.chart representing the chart.
laney_chart <- function(rho = NA, sig_p = NA, sig_z = NA, k = 3) {
    ret <- list()
    ret$rho <- rho
    ret$sig_p <- sig_p
    ret$sig_z <- sig_z
    ret$k <- k
    class(ret) <- c('laney_chart', 'control_chart')
    return(ret)
}

#' Get Laney Chart Limits
#' 
#' Get the limits at which laney chart will signal.
#' @param chart The control chart the signalling limits of which are desired.
#' @param n The sample size. 
#' @return A tuple of the form c(lower_lim, upper_lim), both between 0 and 1.
get_lims.laney_chart <- function(chart, n) {
    l_lim <- pmax(0, chart$rho - chart$k * chart$sig_z * chart$sig_p / sqrt(n))
    u_lim <- pmin(1, chart$rho + chart$k * chart$sig_z * chart$sig_p / sqrt(n))
    return(c(l_lim, u_lim))
}


#' Estimate Laney Chart Parameters.
#'
#' Estimate the parameters for a control chart given data.
#' @param chart The chart object which estimation is to be performed for.
#' @param X A integer vector of observed counts.
#' @param N Either an integer vector of the same length as X, indicating sample sizes, or a scalar integer for constant sample size.
#' @return The chart with parameters modified. 
est_params.laney_chart <- function(chart, X, N) {
    if ((length(X) != length(N)) && length(N) != 1) {
        stop("'N' should either be a scalar, indicating constant sample size, or a vector of the same length of 'X'")
    }
    rho <- mean(X / N)
    sig_p <- sqrt(rho * (1 - rho))

    Z <- (X/N - rho) / (sig_p / sqrt(N))

    ranges <- abs(diff(Z))
    r_bar <- mean(ranges)
    sig_z <- r_bar / 1.128

    #Store the vals
    chart$rho <- rho
    chart$sig_p <- sig_p
    chart$sig_z <- sig_z

    return(chart)
}

#' Calibrate Laney Chart In Control ARL
#'
#' Tune an Laney chart such that it has approximately the desired in sample Average Run Length (ARL) (if the true parameters are known. If they are esimated, this will be even more approximate), under the assumption that data are being drawn form a Beta-Binomial model. 
#' @param chart The control chart to tune.
#' @param target_arl A positive scalar, the desired in control ARL.
#' @param n The sample size, or, in the case of varying sample sizes, the expected sample size. If an expected sampe size is passed, it will introduce more error and there will be a larger discrepency between desired and achieved ARL.
#' @param alpha A positive scalar, the first shape parameter of the beta-binom distribution.
#' @param beta A positive scalar, the second shape parameter of the beta-binom distibution.
#' @return The chart object passed to this function with k, the number of standard deviation from the mean after which to signal, suitably modified to achieve the desired in control ARL.
#TODO: It may be possible to avoid specifying alpha and beta, and instead infer them based on the mu and sig values from the chart. I'm not sure if this is desireable yet.
cal_arl.laney_chart <- function(chart, target_arl, n, alpha, beta) {
    #Define a cost function returning descrepancy between target and specified.
    target_p <- 1 / target_arl
    cost <- function(k) {
        chart$k <- k
        lims <- get_lims(chart, n) * n
        p <- bb.cdf(floor(lims[1]), n, alpha, beta) + 
            (1-bb.cdf(ceiling(lims[2]), n, alpha, beta))
        return((p - target_p)^2)
    }
    chart$k <- optimize(cost, c(0, 10))$minimum
    return(chart)
}
