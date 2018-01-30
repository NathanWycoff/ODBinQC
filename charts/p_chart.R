#!/usr/bin/Rscript
#  charts/p_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.29.2018

#' A P Chart
#'
#' Monitor binomial proportions using a normality approximation to the binomial distribution.
#' @param rho A scalar numeric between 0 and 1, the true probability of success.
#' @param k A positive scalar numeric, how many standard deviations away to place the control limits.
p_chart <- function(rho = NULL, k = 3) {
    ## Create the return object.
    ret <- list()
    ret$rho <- rho
    ret$k <- k
    class(ret) <- c('p_chart', 'control_chart')
    return(ret)
}

#' Get X Chart Limits
#' 
#' Get the limits at which an X chart will signal.
#' @param chart The control chart the signalling limits of which are desired.
#' @param n The sample size. Not used for this chart.
#' @return A tuple of the form c(lower_lim, upper_lim), both between 0 and 1.
get_lims.p_chart <- function(chart, n) {
    l_lim <- pmax(0, chart$rho - chart$k * sqrt(chart$rho * (1 - chart$rho) / n))
    u_lim <- pmin(1, chart$rho + chart$k * sqrt(chart$rho * (1 - chart$rho) / n))
    return(c(l_lim, u_lim))
}

#' Estimate P Chart Parameter.
#'
#' Estimate the probability of success for a P chart
#' @param chart The chart object which estimation is to be performed for.
#' @param X A integer vector of observed counts.
#' @param N Either an integer vector of the same length as X, indicating sample sizes, or a scalar integer for constant sample size.
#' @return The chart with parameter modified. 
est_params.p_chart <- function(chart, X, N) {
    if ((length(X) != length(N)) && length(N) != 1) {
        stop("'N' should either be a scalar, indicating constant sample size, or a vector of the same length of 'X'")
    }
    rho <- mean(X / N)

    #Store the vals
    chart$rho <- rho

    return(chart)
}

#' Calibrate P Chart In Control ARL.
#'
#' Optimize k, the number of standard deviations from the mean at which to place the symmetric control limits to, as closely as possible, match a desired in control ARL.
#'
#' @param chart The control chart to tune.
#' @param target_arl A positive scalar, the desired in control ARL.
#' @param n The sample size, or, in the case of varying sample sizes, the expected sample size. If an expected sampe size is passed, it will introduce more error and there will be a larger discrepency between desired and achieved ARL.
#' @param beta A positive scalar, the second shape parameter of the beta-binom distibution.
#' @param How many times to evaluate the cost function along a 1D grid? Unfortunately, the cost surface for this function is so gross that the built in numerical optimization methods fail, so a manual evaluation over a grid is necessary, followed by a numerical refinement.
#' @return The chart object passed to this function with k, the number of standard deviation from the mean after which to signal, suitably modified to achieve the desired in control ARL.
cal_arl.p_chart <- function(chart, target_arl, n, alpha, beta, evals = 1e3) {
    #Define a cost function returning descrepancy between target and specified.
    target_p <- 1 / target_arl
    cost <- function(k) {
        chart$k <- k
        lims <- get_lims(chart, n) * n
        p <- bb.cdf(floor(lims[1]), n, alpha, beta) + 
            (1-bb.cdf(ceiling(lims[2]), n, alpha, beta))
        return((p - target_p)^2)
    }
    xs <- seq(0, n, length.out=evals)
    ys <- sapply(xs, cost)
    optim_ind <- which.min(ys)
    chart$k <- optimize(cost, c(xs[optim_ind-1], xs[optim_ind+1]))$minimum
    #TODO: Use optimize to search a small interval around this point.
    return(chart)

}
