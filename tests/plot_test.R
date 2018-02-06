#!/usr/bin/Rscript
#  tests/plot_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.30.2018

## See what our plots look like!
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
n.mu <- 10
N <- rpois(m, n.mu-1)^2 + 1
rhos <- rbeta(m, alpha, beta)
X <- rbinom(m, N, rhos)

## Create some charts
charts <- list()

# The beta RE chart
charts$beta.re <- re_beta_chart(alpha = alpha, beta = beta)

# X chart
#TODO: The sd evaluation here is approximate because n is changing but we assume it to be fixed on its mean
mu <- bb.mean(alpha, beta, mean(N)) / mean(N)
sig <- sqrt(bb.var(alpha, beta, mean(N))) / mean(N)
charts$x <- x_chart(mu = mu,
                    sig = sig)

#Laney chart
params <- bb_mm('laney', alpha, beta, mean(N))
charts$laney <- laney_chart(params$rho, params$sig_p, params$sig_z)

# P chart
charts$p <- p_chart(alpha / (alpha + beta))

## Estimate their parameters
charts <- lapply(charts, function(chart) est_params(chart, X, N))

## Plot them!
par(mfrow=c(2,2))
plot(charts$beta.re, X, N)
plot(charts$x, X, N)
plot(charts$laney, X, N)
plot(charts$p, X, N)
