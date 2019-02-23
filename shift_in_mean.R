#!/usr/bin/Rscript
#  examples/shift_in_mean.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.06.2018

# This file had ought to be run in the root dir of the repo.

source('lib.R')
source('sim_lib.R')
source('./charts/re_beta_chart.R')
source('./charts/x_chart.R')
source('./charts/laney_chart.R')
source('./charts/p_chart.R')

set.seed(123)

# Phase I or phase II analysis (estimate params (phase I) or treat them as known (phase II)?)
anls_type <- 1

## Randomly shift the mean and observe each method's performance 
iters <- 1e1#How many times do we run the experiment?
start_est <- 40#How many periods until we start estimating?
time_of_shift <- 45#How many periods until the shift occurs?
N.mu <- 5e3#How many points do we observe on average at each time point
rho_ic <- 0.1#In control mmean proportion of failures

#Percent increases in variation
naive_sig_ic <- sqrt(N.mu) * sqrt(rho_ic * (1-rho_ic))
ks <- c(1.001, 5.0, 20.0)#Scales expected varation under binomial model

for (k in ks) {
    # and get the shape parameters corresponding to a desired mean and variance
    res_ic <- bb_mm('betabinom', N = N.mu, mu = N.mu*rho_ic, sig = k*naive_sig_ic)
    alpha_ic <- res_ic$alpha
    beta_ic <- res_ic$beta

    #Init our charts by giving them true params and calibrating in control ARL.
    charts <- make_charts(alpha_ic, beta_ic)

    # Initialize some storage
    deltas <- seq(0, 0.05, length.out = 11)
    ARL_profile <- matrix(NA, ncol = length(charts), nrow = length(deltas))
    ARL_uncert <- matrix(NA, ncol = length(charts), nrow = length(deltas))
    colnames(ARL_profile) <- c("BetaBinomML", "BetaBinomRMM", "X.Chart", "Laney.Chart", "P.Chart")
    rownames(ARL_profile) <- rho_ic + deltas

    for (d in 1:length(deltas)) {
        ## Determine out of control parameters
        delta <- deltas[d]
        rho_oc <- rho_ic + delta
        naive_sig_oc <- sqrt(N.mu) * sqrt(rho_oc * (1 - rho_oc))
        res_oc <- bb_mm('betabinom', N = N.mu, mu = N.mu*rho_oc, sig = k*naive_sig_oc)
        alpha_oc <- res_oc$alpha
        beta_oc <- res_oc$beta

        run.lens <- matrix(NA, nrow = iters, ncol = length(charts))
        for (iter in 1:iters) {
            i <- 0

            #Prepare some storage
            xs <- c()
            ns <- c()
            all_lims <- list()
            while (TRUE) {
                i <- i + 1

                #Draw a sample
                samp <- draw_sample(i)
                x <- samp$x
                n <- samp$n
                xs <- c(xs, x)
                ns <- c(ns, n)

                #If a Phase I analysis is desired, estimate parameters
                if (anls_type == 1 && i >= start_est) {
                    charts <- lapply(charts, function(chart) est_params(chart, xs, ns))
                }

                #Determine if any limits were triggered (only do this if we have started estimating for Phase I)
                if (anls_type == 2 || i >= start_est) {
                    run.lens[iter,] <- check_lims(x, n, chart)
                }

                #See if all charts have signalled.
                if (sum(is.na(run.lens[iter,])) == 0) {
                    break
                }
            }

            ## On the last iteration, let's make a plot of each decision boundary.
            if (iter == iters) {
                pdf(paste('./images/phase_', anls_type, '_diff_', 
                          rho_oc - rho_ic, '_sig_', k, '.pdf', sep = ''))
                par(mfrow=c(2,2))
                for (chart in charts) {
                    plot(chart, xs, ns)
                }
                dev.off()
            }
        }

        #Store the results of this particular mean difference / variance combination
        ARL_profile[d,] <- colMeans(run.lens)
        ARL_uncert[d,] <- apply(run.lens, 2, sd) / sqrt(iters)
    }

    #Save the output related to a particular variance level to file.
    capture.output(print(ARL_profile), file = paste('./images/phase_', anls_type, '_diff_', 
                          rho_oc - rho_ic, '_sig_arlprofile', k, '.txt', sep = ''))
}
