#!/usr/bin/Rscript
#  sustained_shift.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.22.2018

## A simulation comparing methods for attributes with inherant batch-to-batch variation 
## with a sustained shift in the mean occuring a random amount of time after start.
## ARL is calculated for each.

source('re_beta_chart.R')
require(progress)

### Using beta random effect for beta-binomial posterior predictive.
n.mu <- 1000#Mean number of observations per sample (Drawn form a poisson dist)
alpha <- 10#Parameters for random effect distribution (a beta dist)
beta <- 2
max.periods <- 100#After how many periods do we stop the simulation if the method still hasn't signaled?
iters <- 1e4#How many times do we run the experiment?

#Prepare the chart
chart <- re.beta.chart(alpha, beta)

run.lens <- rep(NA, iters)

pb <- progress_bar$new(total = iters)
for (iter in 1:iters) {
    i <- 0
    pb$tick()
    while (TRUE) {
        i <- i + 1
        n <- rpois(1, n.mu)
        rho <- rbeta(1, alpha, beta)
        x <- rbinom(1, n, rho)
        lims <- chart$get.lims(n)
        if (x < lims[1] || x > lims[2]) {
            run.lens[iter] <- i
            break
        }
    }
}

mean(run.lens)
mean(1/run.lens)

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


