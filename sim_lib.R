#!/usr/bin/Rscript
#  sim_lib.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.20.2018

## Some functions related to the simulation study###Prepare charts for this study
make_charts <- function(alpha_ic, beta_ic) {
    charts <- list()

    #Beta Random Effects, with different estimation strategies.
    charts$beta.re.mle <- re_beta_chart(alpha = alpha_ic, beta = beta_ic, 
                                    def_inf = 'mle')
    charts$beta.re.rmm <- re_beta_chart(alpha = alpha_ic, beta = beta_ic, 
                                    def_inf = 'rmm')

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

    return(charts)
}

draw_sample <- function(i) {
    #Generate a sample
    n <- rpois(1, N.mu)
    if (i < time_of_shift) {
        rho <- rbeta(1, alpha_ic, beta_ic)
    } else {
        rho <- rbeta(1, alpha_oc, beta_oc)
    }
    #Draw the current count of quality characteristic
    x <- rbinom(1, n, rho)

    return(list(x=x, n=n))
}

check_lims <- function(x, n, chart) {
    #Check each chart to see if it has signaled
    lims <- lapply(charts, function(chart) get_lims(chart, n))
    all_lims[[i]] <- lims

    #Make sure all of our limits are attainable
    if (sum(sapply(lims, function(lim) lim[1] == 0 && lim[2] == 1)) > 0) {
        stop(paste("A limit is either 0 or 1 with d = ", d, ", k = ", k, sep =""))
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

    return(run.lens[iter,])
}


