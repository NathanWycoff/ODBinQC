#!/usr/bin/R
#  p_chart_comparison.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.24.2017
#Bring in some plotting functions
source("naive_chart.R")
source("x_chart.R")
source("laney_chart.R")
source("gre_chart.R")

#Generate data with intragroup variation (statisticians call this a "random effect").
mu <- 1#Overall logodds mean
sigma <- 0.2#Random Effects variance
n_mu <- 10000#Mean number of observations per sample
m <- 10#How many days do we sample?

rho_obs <- c()
ns <- c()
for (i in 1:m) {
    #Draw a sample size and mean for this run
    n <- rpois(1, n_mu)
    ns <- c(ns, n)
    rho <- plogis(mu + rnorm(1,0,sigma))

    #Draw the sample for this run
    rho_obs <- c(rho_obs, rbinom(1,n,rho) / n)
}

par(mfrow=c(2,2))
naive_p(rho_obs, ns)
x_chart(rho_obs, ns)
laney_chart(rho_obs, ns)
glmm_chart(rho_obs, ns)
