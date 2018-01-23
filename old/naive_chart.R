#!/usr/bin/R
#  naive_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.24.2017
#Make a naive p chart

#rhos_obs is a vector of observed proportions
#n is a vector of observed sample sizes
naive_p <- function(rho_obs, n) {
    #Estimate overall mean
    mean_est <- sum(rho_obs * n) / sum(n)
    m <- length(n)

    #Estimate observation variance
    var_est <- mean_est * (1-mean_est)

    #Get the variance for each observed proportion
    prop_var <- var_est / n

    #Calc 3-sigma limits
    bounds <- 3*sqrt(prop_var)
    l_lim <- pmax(0, mean_est - bounds)
    u_lim <- pmin(1, mean_est + bounds)

    #Create our plots
    max_point <- max(c(max(u_lim), max(rho_obs)))
    min_point <- min(c(min(l_lim), min(rho_obs)))

    plot(NULL, ylim = c(min_point, max_point), xlim = c(1, m), 
         main = "Naive P chart", xlab = "Index", ylab = "Proportion")

    #Check which points are in or out of control.
    in_control <- rho_obs < u_lim & rho_obs > l_lim
    in_control_col <- sapply(in_control, function(i) ifelse(i,'black','red'))

    points(rho_obs, col = in_control_col)
    points(l_lim, type = 'l', col = 'orange')
    points(u_lim, type = 'l', col = 'orange')
}
