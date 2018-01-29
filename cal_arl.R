#!/usr/bin/Rscript
#  sustained_shift.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.22.2018

## Calculate ARL for a couple of different charts while in control to calibrate chart parameters.

source('lib.R')
source('./charts/re_beta_chart.R')
source('./charts/x_chart.R')
require(progress)

## Calibrate the individual x chart to get desired ARL
#Set some params
n.mu <- 1000#Mean number of observations per sample (Drawn form a poisson dist)
alpha <- 10#Parameters for random effect distribution (a beta dist)
beta <- 2

mu <- bb.mean(alpha, beta, n.mu) / n.mu
sig <- sqrt(bb.var(alpha, beta, n.mu)) / n.mu
xc <- x_chart(mu = mu, sig = sig)

target.arl <- 20
target.p <- 1 / target.arl

cost <- function(k) {
    lims <- xc$get.lims(n.mu, k = k) * n.mu
    p <- bb.cdf(floor(lims[1]), n.mu, alpha, beta) + 
        (1-bb.cdf(ceiling(lims[2]), n.mu, alpha, beta))
    return((p - target.p)^2)
}

x.chart.optim.k <- optimize(cost, c(0, 10))$minimum

### Generate data from a beta-binomial model for charts which require simulation
n.mu <- 1000#Mean number of observations per sample (Drawn form a poisson dist)
alpha <- 10#Parameters for random effect distribution (a beta dist)
beta <- 2
max.periods <- 100#After how many periods do we stop the simulation if the method still hasn't signaled?
iters <- 1e2#How many times do we run the experiment?

#Prepare the charts
charts <- list()
charts$beta.re <- re.beta.chart(alpha = alpha, beta = beta)
#TODO: The sd evaluation here is approximate because n is changing but we assume it to be fixed on its mean
mu <- bb.mean(alpha, beta, n.mu) / n.mu
sig <- sqrt(bb.var(alpha, beta, n.mu)) / n.mu
charts$x <- x.chart(mu = mu,
                    sig = sig)

run.lens <- matrix(NA, nrow = iters, ncol = length(charts))

pb <- progress_bar$new(total = iters)
for (iter in 1:iters) {
    i <- 0
    pb$tick()
    while (TRUE) {
        i <- i + 1
        #Generate a sample
        n <- n.mu#rpois(1, n.mu)
        rho <- rbeta(1, alpha, beta)
        x <- rbinom(1, n, rho)

        #Check each chart to see if it has signaled
        #lims <- lapply(charts, function(chart) chart$get.lims(n))
        re.lims <- charts$beta.re$get.lims(n)
        x.lims <- charts$x$get.lims(n.mu, x.chart.optim.k)
        lims <- list(re.lims, x.lims)

        for (l in 1:length(lims)) {
            lim <- lims[[l]]
            #if this chart hasn't signalled yet, check to see if it has this time.
            if (is.na(run.lens[iter, l])) {
                if (x / n < lim[1] || x / n > lim[2]) {
                    run.lens[iter, l] <- i
                }
            }
        }

        #See if all charts have signalled.
        if (sum(is.na(run.lens[iter,])) == 0) {
            break
        }
    }
}

colMeans(run.lens)
colMeans(1/run.lens)

### Using Normal random effect with logit link
n.mu <- 100#Mean number of observations per sample (Drawn form a poisson dist)
eta.ic <- 0.5#In control mean log-odds of success
eta.oc <- 0.5#Out of control mean log-odds of success
sig <- 0.9#In control standard deviation
max.periods <- 100#After how many periods do we stop the simulation if the method still hasn't signaled?

#Prepare the chart
chart <- re.logit.chart(eta.ic, sig)

i <- 0
while (i < max.periods) {
    i <- i + 1
    n <- rpois(1, n.mu)
    rho <- plogis(eta.ic + rnorm(1, 0, sig))
    x <- rbinom(1, n, rho)
    lims <- chart$get.lims(x, n)
    if (x < lims[1] || x > lims[2]) {
        print("Oh shit cousin")
    }
}
