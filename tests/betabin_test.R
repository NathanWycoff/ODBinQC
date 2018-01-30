#!/usr/bin/Rscript
#  test_betabin_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.23.2018

## This script puts some beta-binomial functions to the test
source('./charts/re_beta_chart.R')
source('./lib.R')

alpha <- 10
beta <- 2
n <- 1000
ns <- 100

## Compare the two methods against one another, using the hypergeometric func
## vs just summing the pmf
ks <- round(seq(0,n, length.out = ns))
system.time(as <- sapply(ks, function(k) bb.cdf(k,n,alpha,beta)))

## Compare the quantile func to observed data
p <- 0.9
n <- 100
q <- bb.quantile(p, n, alpha, beta)

iters <- 10000
X <- sapply(1:iters, function(i) {
                rho <- rbeta(1,alpha,beta)
                return(rbinom(1,n,rho))
})

cat('observed')
quantile(X, p)
cat('actual')
q

## See if we can recover the appropriate beta binomial parameters to match a specifeid mean and variance.
alpha <- 2000
beta <- 100
N <- 20
mu <- bb.mean(alpha, beta, N)
sig <- sqrt(bb.var(alpha, beta, N))

bb_mm('betabinom', mu = mu, sig = sig, N = N)
