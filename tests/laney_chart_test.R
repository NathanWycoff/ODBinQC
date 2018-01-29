#!/usr/bin/Rscript
#  test_laney_ests.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.25.2018

## In order to evaluate ARL properties of the laney chart with known parameters,
## we need to understand how beta-binomial quantities translate to the population
## quantities mentioned in Laney's paper, namely sigma_z and sigma_p.
## This script verifies analytical derivations of these quantities through sims.

source('./charts/laney_chart.R')
source('lib.R')

## Double check that under a model with no random effects, we eventually estimate
## that the z-scale variation is indeed 1.
chart <- laney_chart()
m <- 1000
N <- rpois(m,10)
X <- rbinom(m, N, 0.5)

est_params(chart, X, N)

## Check that, as alpha, beta \to \infty with \alpha / (\alpha + \beta) fixed, we have that sig_z \to 1.
alpha <- 2e10
beta <- 1e10

bb_mm('laney', alpha, beta, 3)


## Check that for large sample sizes, the sample estimates agree with the translated
## parameters.
## Interesting note: the moving range estimate for variance is only good if p is no
## close to 0 or 1, as this skews the distribution of Z scores.
chart <- laney_chart()
alpha <- 20
beta <- 10
m <- 1e4
n.mu <- 1e4
N <- 1+rpois(m,n.mu-1)
rhos <- rbeta(m, alpha, beta)
X <- rbinom(m, N, rhos)

#Translated Params
bb_mm('laney', alpha, beta, n.mu)

#Estiamtes
est_params(chart, X, N)
