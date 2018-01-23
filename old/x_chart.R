#!/usr/bin/R
#  x_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.24.2017
#Make a x chart, effectively treating each proportion as a realization of
#a normal variate. 

#rhos_obs is a vector of observed proportions
#n is a vector of observed sample sizes
x_chart <- function(rho_obs, n) {
    m <- length(rho_obs)
    mean_est <- sum(rho_obs * n) / sum(n)

    #Calculate sd using range
    ranges <- abs(diff(rho_obs))
    r_bar <- mean(ranges)
    sigma_p <- r_bar / 1.128

    #Calc 3-sigma limits
    bounds <- 3*sigma_p
    l_lim <- pmax(0, mean_est -rep(bounds, m))
    u_lim <- pmin(1, mean_est + rep(bounds, m))

    #Create our plots
    max_point <- max(c(max(u_lim), max(rho_obs)))
    min_point <- min(c(min(l_lim), min(rho_obs)))

    plot(NULL, ylim = c(min_point, max_point), xlim = c(1, m), 
         main = "X-Chart for Proportions", xlab = "Index", ylab = "Proportion")

    #Check which points are in or out of control.
    in_control <- rho_obs <= u_lim & rho_obs >= l_lim
    in_control_col <- sapply(in_control, function(i) ifelse(i,'black','red'))

    points(rho_obs, col = in_control_col)
    points(l_lim, type = 'l', col = 'orange')
    points(u_lim, type = 'l', col = 'orange')
}
