#!/usr/bin/Rscript
#  test_betabin_funcs.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.23.2018

## This script puts some beta-binomial functions to the test
library(bbmle)
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

## See if we can estimate alpha and beta given some data.
alpha <- 10
beta <- 2
N <- 100
m <- 300
rho <- rbeta(m, alpha, beta)
X <- rbinom(m, N, rho)
chart <- re_beta_chart()

# Using robust method of moments
est_params(chart, X, N)

# Using maximum likelihood
est_params(chart, X, N, type = 'mle')

## Why can we get good values for alpha and beta with sd's on proportions that don't make any sense

#Specify a mean and variance
rho <- 0.1
N <- 1e2
sig_rho_naive <- sqrt(rho * (1 - rho) / N)

mu_x <- N * rho
sig_x <- N * sig_rho_naive

k <- 1
sig_x <- k * sqrt(N) * sqrt(rho*(1-rho))

sig_x <- sqrt(N) * sqrt(0.6)


###
#Get the alpha and beta
res <- bb_mm('betabinom', N = N, mu = mu_x, sig = sig_x)
print(res)

#Sample from it
m <- 3e4
rho <- rbeta(m, res$alpha, res$beta)
X <- rbinom(m, N, rho)

#Find out what's happenin
mean(X)
sd(X)
mean(X / N)
sd(X / N)
