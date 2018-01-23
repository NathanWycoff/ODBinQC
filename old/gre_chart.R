#!/usr/bin/Rscript
#  gre_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 12.26.2017

## A generalized Random Effects chart for attribute data.

require(lme4)
require(MCMCglmm)

#Generate data with intragroup variation (statisticians call this a "random effect").
mu <- 1#Overall logodds mean
sigma <- 0.2#Random Effects variance
n_mu <- 1000#Mean number of observations per sample
m <- 10#How many days do we sample?

rho_obs <- c()
xs <- c()
ns <- c()
for (i in 1:m) {
    #Draw a sample size and mean for this run
    n <- rpois(1, n_mu)
    ns <- c(ns, n)
    rho <- plogis(mu + rnorm(1,0,sigma))

    #Draw the sample for this run
    x_obs <- rbinom(1,n,rho) 
    rho_obs <- c(rho_obs, x_obs / n)
    xs <- c(xs, x_obs)
}

ys <- unlist(sapply(1:m, function(i) c(rep(1,xs[i]), rep(0, ns[i]-xs[i]))))
day <- rep(1:m, ns)

#Fit the model using asymptotics
fit <- glmer(ys ~ (1 | day), family = 'binomial')
summary(fit)

predict(fit, newdata = c(), interval = 'predict')

#Fit the model using MCMC
dummy <- 1
df <- as.data.frame(cbind(ys, day, dummy))
fit <- MCMCglmm(fixed = ys ~ dummy + 0, random=~day, family = 'categorical', data = df)

#Do prediction
newdata = 
