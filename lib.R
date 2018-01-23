#!/usr/bin/Rscript
#  lib.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 01.22.2018

## Functions useful for multiple files in this repository.

#' Generate overdispered binomial data
#' 
#' Generate binomial data with overdispersion through random Effects using either a Beta-Binomial model or GLMM-like random effects model with normal random effects.
#' @param re.dist A scalar character vector indicating which random effects distribution to use for overdispersion. Currently, 'beta' and 'normal' are supported. A logistic link is used with the normal data.
