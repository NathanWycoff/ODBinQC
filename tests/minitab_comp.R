#!/usr/bin/Rscript
#  tests/minitab_comp.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.30.2018

## Compare the results of our plots to minitab's output for a given dataset.
source('lib.R')
source('./charts/re_beta_chart.R')
source('./charts/x_chart.R')
source('./charts/laney_chart.R')
source('./charts/p_chart.R')

## Generate some data
seed <- 123
set.seed(seed)
alpha <- 20
beta <- 10

m <- 100
n.mu <- 1000
N <- rpois(m, n.mu)
rhos <- rbeta(m, alpha, beta)
X <- rbinom(m, N, rhos)

## Write them to a csv
write.csv(as.data.frame(cbind(X, N)), './data/minitab_test_data.csv')

## Create some charts
charts <- list()

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

## Estimate their parameters
charts <- lapply(charts, function(chart) est_params(chart, X, N))

## Get their limits
lapply(N, function(n) lapply(charts, function(chart) get_lims(chart, n)))

## Plot them!
par(mfrow=c(2,2))
plot(charts$x, X, N)#IMR Chart looks good!
plot(charts$laney, X, N)
plot(charts$p, X, N)#P chart looks good!
