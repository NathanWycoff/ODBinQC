#!/usr/bin/Rscript
#  glmm_chart.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.24.2017

#Make a naive p chart

#For inference on model parameters
require(rjags)

#Estimate the overall mean by simply taking the total average
#Estimate the variance using our mean estimate
#rhos_obs is a vector of observed proportions
#n is a vector of observed sample sizes
#mu_0 is the prior mean for eta, the overall mean logodds of success (norm dist).
#sigma2_0 is the prior variance for eta
#theta_0 is the prior mean for sigma2, the random effect variance (exp dist).
glmm_chart <- function(rho_obs, ns, ...) {
    m <- length(rho_obs) 

    #Do inference on the model params
    result <- model_inference(rho_obs, ns, ...)

    #Create quantiles via monte carlo
    intervals <- sapply(ns, function(n) int_rho(n, result$eta, result$sigma_eta))

    #Create the limits
    l_lim <- intervals[1,]
    u_lim <- intervals[2,]

    #Create our plots
    max_point <- max(c(max(u_lim), max(rho_obs)))
    min_point <- min(c(min(l_lim), min(rho_obs)))

    plot(NULL, ylim = c(min_point, max_point), xlim = c(1, m), 
         main = "GLMM Plot", xlab = "Index", ylab = "Proportion")

    #Check which points are in or out of control.
    in_control <- rho_obs <= u_lim & rho_obs >= l_lim
    in_control_col <- sapply(in_control, function(i) ifelse(i,'black','red'))

    points(rho_obs, col = in_control_col)
    points(l_lim, type = 'l', col = 'orange')
    points(u_lim, type = 'l', col = 'orange')
}

model_inference <- function(rho_obs, ns, mu_0 = 0, sigma2_0 = 1, theta_0 = 1) {
    #Store our data in a form jags likes
    data <- list(n_obs = rho_obs*ns, ns = ns, N = length(ns))
    hyper <- list(mu_0 = mu_0, sigma2_0 = sigma2_0, theta_0 = theta_0)

    #Specify our model
    model <- "
    model {
        #Likelihood
        for (i in 1:N) {
            #Generate the random effects
            etas[i] ~ dnorm(0, phi_eta)

            #Generate the data
            n_obs[i] ~ dbinom(phi(eta + etas[i]), ns[i])
        }

        #R vibes better with standard deviations than with precisions.
        sigma_eta <- sqrt(1/phi_eta)

        #Prior
        eta ~ dnorm(mu_0, sigma2_0)
        phi_eta ~ dexp(theta_0)
    }
    "

    #Run the actual sampler!
    sampler <- jags.model(textConnection(model), 
                          data=c(data, hyper))
    update(sampler, n.iter=100)
    output=coda.samples(model=sampler, variable.names=c("eta", "sigma_eta"),
                        n.iter=1000)
    post.mean <- colMeans(output[[1]])
    return(list(eta = post.mean[['eta']], sigma_eta = post.mean[['sigma_eta']]))
}

#Sample from the generative model
generative_model <- function(m, eta, sigma_eta, n_mu) {
    ns <- rpois(m, n_mu)
    etas <- rnorm(m, 0, sigma_eta)
    rho_obs <- rbinom(m, ns, pnorm(eta + etas)) / ns
    return(list(rho_obs = rho_obs, ns = ns))
}

#Get the interval for a rho, with eta and sigma_eta fixed
int_rho <- function(n, eta, sigma_eta, iters = 1000) {
    x <- rbinom(iters, n, pnorm(eta + rnorm(iters, 0, sigma_eta))) / n

    return(quantile(x, c(0.025, 0.975)))
}

#Test out the sampler a bit
#m <- generative_model(100,1,1,10)
#result <- model_inference(m$rho_obs, m$ns)
