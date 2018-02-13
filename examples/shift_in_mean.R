#!/usr/bin/Rscript
#  examples/shift_in_mean.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.06.2018

options(error=dump.frames)

source('lib.R')
source('./charts/re_beta_chart.R')
source('./charts/x_chart.R')
source('./charts/laney_chart.R')
source('./charts/p_chart.R')
require(progress)
require(ggplot2)

## Randomly shift the mean and observe each method's performance during Phase I

iters <- 1e1#How many times do we run the experiment?
K <- 20#When does the shift occur?
N.mu <- 1e3#How many points do we observe on average at each time point

#Specify a certain mean and variance...
rho_ic <- 0.1
sig <- N.mu*0.05

# and get the appropriate alpha and beta
res_ic <- bb_mm('betabinom', N = N.mu, mu = N.mu*rho_ic, sig = sig)
alpha_ic <- res_ic$alpha
beta_ic <- res_ic$beta

# Go over a grid of mean differences
deltas <- seq(0, 0.2, length.out = 11)
ARL_profile <- matrix(NA, ncol = 4, nrow = length(deltas))
colnames(ARL_profile) <- c("BetaBinom", "X.Chart", "Laney.Chart", "P.Chart")
rownames(ARL_profile) <- rho_ic + deltas
for (d in 1:length(deltas)) {
    ## Determine out of sample parameters
    delta <- deltas[d]
    rho_oc <- rho_ic + delta
    res_oc <- bb_mm('betabinom', N = N.mu, mu = N.mu*rho_oc, sig = sig)
    alpha_oc <- res_oc$alpha
    beta_oc <- res_oc$beta

    ##########################################################################################
    ###Prepare the charts
    charts <- list()

    #Beta Random Effects
    charts$beta.re <- re_beta_chart(alpha = alpha_ic, beta = beta_ic)

    # X chart
    #TODO: The sd evaluation here is approximate because n is changing but we assume it to be fixed on its mean
    mu <- bb.mean(alpha_ic, beta_ic, N.mu) / N.mu
    sig_xc <- sqrt(bb.var(alpha_ic, beta_ic, N.mu)) / N.mu
    charts$x <- x_chart(mu = mu,
                        sig = sig_xc)

    #Laney chart
    params <- bb_mm('laney', alpha_ic, beta_ic, N.mu)
    charts$laney <- laney_chart(params$rho, params$sig_p, params$sig_z)

    # P chart
    charts$p <- p_chart(alpha_ic / (alpha_ic + beta_ic))

    # Calibrate the charts to have the desired in control ARL
    target_arl <- 100
    charts <- lapply(charts, function(chart) cal_arl(chart, target_arl, N.mu, alpha_ic, beta_ic))


    ##########################################################################################
    ### Do the actual simulation
    run.lens <- matrix(NA, nrow = iters, ncol = length(charts))
    pb <- progress_bar$new(total = iters)
    for (iter in 1:iters) {
        i <- 0
        pb$tick()

        #Prepare some storage
        xs <- c()
        ns <- c()
        all_lims <- list()
        while (TRUE) {
            i <- i + 1
            #Generate a sample
            n <- rpois(1, N.mu)
            ns <- c(ns, n)
            if (i < K) {
                rho <- rbeta(1, alpha_ic, beta_ic)
            } else {
                rho <- rbeta(1, alpha_oc, beta_oc)
            }

            #Draw the current count of quality characteristic
            x <- rbinom(1, n, rho)
            xs <- c(xs, x)

            #Check each chart to see if it has signaled
            lims <- lapply(charts, function(chart) get_lims(chart, n))
            all_lims[[i]] <- lims

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

        ## On the last iteration, let's make a plot of each decision boundary.
        if (iter == iters) {
            pdf(paste('./images/type_1_diff_', rho_oc - rho_ic, '_sig_', sig/N.mu, '.pdf', sep = ''))
            par(mfrow=c(2,2))
            for (chart in charts) {
                plot(chart, xs, ns)
            }
            dev.off()
        }
    }

    ARL_profile[d,] <- colMeans(run.lens)
}

#Save the output to file.
capture.output(print(ARL_profile), file = './output/type_1_arlprofile.txt')
