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
    "K" = 2,
    "zero_inf_prob" = 0
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
#' @return result [stan object] The fitted stan object.
#' @importFrom rstan stan stan_model vb extract cpp_object_initializer
#' @export
fit_model <- function(y, model_opts = list()) {
  model_opts <- merge_model_opts(model_opts)
  stan_data <- list(
    "N" = nrow(y),
    "P" = ncol(y),
    "y" = y,
    "K" = model_opts$K,
    "zero_inf_prob" = model_opts$zero_inf_prob
  )

  clear_tmp()
  if (model_opts$inference == "gibbs") {
    result <- rstan::extract(stan(
                       file = model_opts$method,
                       data = stan_data,
                       chain = 1,
                       warmup = 1000,
                       iter = 1500
                     ))
  } else if (model_opts$inference == "vb") {
    f <- stan_model(model_opts$method)
    result <- rstan::extract(vb(f, stan_data, output_samples = 500))
  } else if (model_opts$inference == "bootstrap") {
    result <- bootstrap_vb(model_opts$method, data = stan_data)
  } else {
    stop("model_opts$inference is not recognized")
  }

  result
}

#' Wrapper to get posterior means in NMF
#'
#' @param samples The output of a call to vb() using the nmf code
#' @return result [list] A list with the posterior means as matrices, with
#'   factors (k) as columns
#' @export
nmf_posterior_means <- function(samples) {
  list(
    "theta_hat" = t(ldaSim::posterior_mean(samples$theta, c("k", "i"))),
    "beta_hat" = t(ldaSim::posterior_mean(samples$beta, c("k", "v"))),
    "prior_hat" = colMeans(samples$prior_params)
  )
}

#' Clear all files in tempdir()
#'
#' Hack used to avoid crashing https://groups.google.com/forum/#!category-topic/stan-users/general/Bg8HhyB6qG4
clear_tmp <- function() {
  tmp_files <- list.files(tempdir(), full.names = TRUE)
  file.remove(tmp_files)
}

#' Fit Parameteric Bootstrap for Variational Bayes
#'
#' @param method [character] The path to stan file containing model fitting
#'   code.
#' @param data [list] A list of data to input in the data{} field of the stan
#'   file.
#' @param B [integer] The number of bootstrap replicates to compute.
#' @return result [list] A list containing beta and theta fields, which is an
#'   array similar to what is output by stan(), except columns are now bootstrap
#'   replicates instead of sampling iterations
#' @export
bootstrap_vb <- function(method, data, B = 500) {
  ## First, make a VB fit, to use as the estimated parameters in the parametric
  ## bootstrap
  f <- stan_model(method)
  vb_fit <- vb(f, data, check_data = FALSE, adapt_engaged = FALSE, eta = 0.1, output_samples = B)
  samples <- extract(vb_fit)
  means0 <- nmf_posterior_means(samples)

  theta_boot <- array(0, c(B, data$N, data$K))
  beta_boot <- array(0, c(B, data$P, data$K))
  prior_boot <- matrix(0, B, 4)

  for (b in seq_len(B)) {
    cat(sprintf("Bootstrap iteration %s\n", b))

    ## Simulate according to fitted parameters
    cur_data <- data
    cur_data$y <- sim_from_params(
      means0$theta_hat,
      means0$beta_hat,
      data$zero_inf_prob
    )$y

    ## Fit another VB iteration
    cur_fit <- try(
      vb(f, cur_data, check_data = FALSE, adapt_engaged = FALSE, eta = 0.1, output_samples = B)
    )
    if (class(cur_fit) != "try-error") {
      cur_means <- nmf_posterior_means(extract(cur_fit))
      theta_boot[b,,] <- t(cur_means$theta_hat)
      beta_boot[b,,] <- t(cur_means$beta_hat)
      prior_boot[b,] <- cur_means$prior_hat
    }
  }

  list("theta" = theta_boot, "beta" = beta_boot, "prior" = prior_boot)
}
