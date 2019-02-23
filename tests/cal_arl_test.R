#!/usr/bin/Rscript
#  test_cal_arl.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.24.2018

## Test functions that calibrate in control ARL.

source('lib.R')
source('charts/re_beta_chart.R')
source('charts/x_chart.R')
source('charts/laney_chart.R')
source('charts/p_chart.R')
library(progress)

### Generate data from a beta-binomial model
n.mu <- 1000#Mean number of observations per sample (Drawn form a poisson dist)
alpha <- 10#Parameters for random effect distribution (a beta dist)
beta <- 10
max_periods <- 100#After how many periods do we stop the simulation if the method still hasn't signaled?
iters <- 1e3#How many times do we run the experiment?

##Prepare the charts
charts <- list()

#Beta Random Effects
charts$beta.re <- re_beta_chart(alpha = alpha, beta = beta)

# X chart
#TODO: The sd evaluation here is approximate because n is changing but we assume it to be fixed on its mean
mu <- bb.mean(alpha, beta, n.mu) / n.mu
sig <- sqrt(bb.var(alpha, beta, n.mu)) / n.mu
charts$x <- x_chart(mu = mu,
                    sig = sig)

#Laney chart
params <- bb_mm('laney', alpha, beta, n.mu)
charts$laney <- laney_chart(params$rho, params$sig_p, params$sig_z)

# P chart
charts$p <- p_chart(alpha / (alpha + beta))

# Calibrate the charts to have the desired in control ARL
target_arl <- 100
charts <- lapply(charts, function(chart) cal_arl(chart, target_arl, n.mu, alpha, beta))

# Prepare Storage
run.lens <- matrix(NA, nrow = iters, ncol = length(charts))

# Do the actual simulation
pb <- progress_bar$new(total = iters)
for (iter in 1:iters) {
    i <- 0
    pb$tick()
    while (i < max_periods) {
        i <- i + 1
        #Generate a sample
        n <- rpois(1, n.mu)
        rho <- rbeta(1, alpha, beta)
        x <- rbinom(1, n, rho)

        #Check each chart to see if it has signaled
        lims <- lapply(charts, function(chart) get_lims(chart, n))

        #Make sure all of our limits are attainable
        if (sum(sapply(lims, function(lim) lim[1] == 0 || lim[2] == 1)) > 0) {
            stop("A limit is either 0 or 1")
        }

        for (l in 1:length(lims)) {
            lim <- lims[[l]]
            #if this chart hasn't signalled yet, check to see if it has this time.
            if (is.na(run.lens[iter, l])) {
                if (x / n < lim[1] || x / n > lim[2]) {
                    run.lens[iter, l] <- i
                }
            }
        }

        #See if all charts have signalled.
        if (sum(is.na(run.lens[iter,])) == 0) {
            break
        }
    }
}

## These should all be near target_arl
colMeans(run.lens)
