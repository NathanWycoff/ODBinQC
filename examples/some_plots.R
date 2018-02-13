#!/usr/bin/Rscript
#  examples/some_plots.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.13.2018

## As a sanity check, visualize a run of the "shift_in_mean.R".

source('lib.R')
source('./charts/re_beta_chart.R')
source('./charts/x_chart.R')
source('./charts/laney_chart.R')
source('./charts/p_chart.R')

## Randomly shift the mean and observe each method's performance during Phase I

K <- 20#When does the shift occur?
N.mu <- 1e3#How many points do we observe on average at each time point

#Specify a certain mean and variance...
rho_ic <- 0.1
rho_oc <- 0.2
sig <- N.mu*0.4

# and get the appropriate alpha and beta
res_ic <- bb_mm('betabinom', N = N.mu, mu = N.mu*rho_ic, sig = sig)
alpha_ic <- res_ic$alpha
beta_ic <- res_ic$beta
res_oc <- bb_mm('betabinom', N = N.mu, mu = N.mu*rho_oc, sig = sig)
alpha_oc <- res_oc$alpha
beta_oc <- res_oc$beta

##Prepare the charts
charts <- list()

#Beta Random Effects
charts$beta.re <- re_beta_chart(alpha = alpha_ic, beta = beta_ic)

# X chart
#TODO: The sd evaluation here is approximate because n is changing but we assume it to be fixed on its mean
mu <- bb.mean(alpha_ic, beta_ic, N.mu) / N.mu
sig <- sqrt(bb.var(alpha_ic, beta_ic, N.mu)) / N.mu
charts$x <- x_chart(mu = mu,
                    sig = sig)

#Laney chart
params <- bb_mm('laney', alpha_ic, beta_ic, N.mu)
charts$laney <- laney_chart(params$rho, params$sig_p, params$sig_z)

# P chart
charts$p <- p_chart(alpha_ic / (alpha_ic + beta_ic))

# Calibrate the charts to have the desired in control ARL
target_arl <- 100
charts <- lapply(charts, function(chart) cal_arl(chart, target_arl, N.mu, alpha_ic, beta_ic))

# Prepare Storage
run.lens <- matrix(NA, nrow = iters, ncol = length(charts))

# Do the actual simulation
pb <- progress_bar$new(total = iters)
for (iter in 1:iters) {
    i <- 0
    pb$tick()
    while (TRUE) {
        i <- i + 1
        #Generate a sample
        n <- rpois(1, N.mu)
        if (i < K) {
            rho <- rbeta(1, alpha_ic, beta_ic)
        } else {
            rho <- rbeta(1, alpha_oc, beta_oc)
        }
        x <- rbinom(1, n, rho)

        #Check each chart to see if it has signaled
        lims <- lapply(charts, function(chart) get_lims(chart, n))

        #Make sure all of our limits are attainable
        if (sum(sapply(lims, function(lim) lim[1] == 0 && lim[2] == 1)) > 0) {
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
