#!/usr/bin/Rscript
#  tests/x_chart_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 06.28.2019

## Test functionality of some x_chart functions
source('./charts/x_chart.R')
source('lib.R')

n <- 100
alpha <- 20
beta <- 40
mu <- alpha / (alpha + beta)
sig <- sqrt(mu * (1-mu)) / n
chart <- x_chart(mu, sig)

target_arl <- 100
chart <- cal_arl(chart, target_arl, n, alpha, beta)

## Test that the ARL is right
iters <- 1e6
rhos <- rbeta(iters, alpha, beta)
X <- rbinom(iters, n, rhos)
rho_hats <- X / n
lims <- get_lims(chart)
mean(rho_hats <= lims[1] | rho_hats >= lims[2])
