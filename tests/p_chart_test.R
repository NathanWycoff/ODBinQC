#!/usr/bin/Rscript
#  tests/p_chart_test.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.30.2018

## Test functionality of some p_chart functions
source('./charts/p_chart.R')
source('lib.R')

## Test that we can recover rho, the probability of success, even under random effects
n <- 10000
N.mu <- 100
N <- rpois(n,N.mu-1)+1
alpha <- 2
beta <- 4
true_rho_exp <- alpha / (alpha + beta)
rhos <- rbeta(n,alpha,beta)
X <- rbinom(n,N,rhos)

chart <- p_chart()
chart <- est_params(chart, X, N)

print(chart$rho)
print(true_rho_exp)
