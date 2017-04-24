#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Helpers to simulate data according to an NMF model.
## author: kriss1@stanford.edu

## ---- simulation-helpers ----
#' Merge Default NMF options
#'
#' @param opts [list] A partially filled list of options, to fill in with
#'   defaults.
#' @return opts [list] The version of opts with defaults filled in.
#' @export
merge_nmf_opts <- function(opts = list()) {
  default_opts <- list(
    "K" = 2,
    "N" = 100,
    "P" = 75,
    "a" = 1,
    "b" = 1,
    "c" = 1,
    "d" = 1,
    "zero_inf_prob" = 0
  )
  modifyList(default_opts, opts)
}

#' Simulate NMF Data (optionally zero inflated)
#'
#' To facilitate simulation across many parameters, it's useful to have a single
#' function to generate all the quantities of interest.
#'
#' @param opts [list] A list containing parameters for simulation. Any options
#'   that are not specified will be filled in with defaults, according to
#'   merge_nmf_defaults().
#' @return A list with the latent thetas, betas, mask, and observed Y.
#' @export
nmf_sim <- function(opts) {
  opts <- merge_nmf_opts(opts)

  ## scores
  theta <- matrix(
    rgamma(opts$N * opts$K, rate = opts$a, shape = opts$b),
    opts$N, opts$K
  )

  ## factors
  beta <- matrix(
    rgamma(opts$P * opts$K, rate = opts$c, shape = opts$d),
    opts$P, opts$K
  )

  sim <- sim_from_params(theta, beta, opts$zero_inf_prob)
  list(
    "theta" = theta,
    "beta" = beta,
    "mask" = sim$mask,
    "y" = sim$y
  )
}

#' Simulate (ZI)-NMF data from theta, beta, and inflation parameters
#'
#' @param theta [matrix] An N x K matrix of scores
#' @param beta [matrix] A V x K matrix of factors
#' @param zero_inf_prob [numeric] The amount of zero inflation; i.e., the
#'   proportion in [0, 1] that will be set to zero.
#' @return y [matrix] The N x V matrix of observed (potentially zero-inflated
#'   counts).
sim_from_params <- function(theta, beta, zero_inf_prob = 0) {
  N <- nrow(theta)
  P <- nrow(beta)

  ## observations
  y <- matrix(
    rpois(N * P, theta %*% t(beta)),
    N, P
  )

  ## set some proportion to zero
  mask <- matrix(
    sample(
      c(0, 1),
      N * P,
      replace = TRUE,
      prob = c(1 - zero_inf_prob, zero_inf_prob)),
    N, P
  )
  y[mask == 1] <- 0

  list("mask" = mask, "y" = y)
}
