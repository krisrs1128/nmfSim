#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Helpers for fitting NMF / Z-NMF models using STAN.
##
## author: kriss1@stanford.edu

## ---- modeling-helpers ----
#' Merge Default Modeling Parameters
#'
#' This lets us run models with partially specified options.
#'
#' @param opts [list] A partially filled list of options, to fill in with
#'   defaults.
#' @return opts [list] The version of opts with defaults filled in.
#' @export
merge_model_opts <- function(opts = list()) {
  default_opts <- list(
    "inference" = "gibbs",
    "method" = file.path(.libPaths(), "nmfSim", "extdata", "nmf_gamma_poisson.stan"),
    "K" = 2
  )
  modifyList(default_opts, opts)
}

#' Fit a generic NMF Stan Model
#'
#' This wraps vb() and stan() in the STAN package to let us run either approach
#' using a single command.
#'
#' @param y [matrix] The data on which to fit the NMF model.
#' @param model_opts [list] A partially filled list of model fitting options.
#'   Unspecified options will be passed into merge_model_opts().
#' @param prior_opts [list] A list of prior information, required by the NMF
#'   fitting STAN code.
#' @return result [stan object] The fitted stan object.
#' @importFrom rstan stan stan_model vb extract cpp_object_initializer
#' @export
fit_model <- function(y, model_opts = list(), prior_opts = list()) {
  stan_data <- list(
    "N" = nrow(y),
    "P" = ncol(y),
    "y" = y,
    "K" = model_opts$K,
    "zero_inf_prob" = model_opts$zero_inf_prob
  )
  stan_data <- c(stan_data, prior_opts)

  if (grepl("zero", model_opts$method)) {
    stan_data$zero_inf_prob <- NULL
  }

  if (model_opts$inference == "gibbs") {
    result <- stan(file = model_opts$method, data = stan_data, chain = 1)
  } else if (model_opts$inference == "vb") {
    f <- stan_model(model_opts$method)
    result <- vb(f, stan_data)
  } else {
    stop("model_opts$inference is not recognized")
  }

  extract(result)
}
