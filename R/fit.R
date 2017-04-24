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
  } else if (model_opts$inference == "bootstrap") {
    result <- bootstrap_vb(model_opts$method, data = stan_data)
  } else {
    stop("model_opts$inference is not recognized")
  }

  rstan::extract(result)
}

nmf_posterior_means <- function(samples) {
  list(
    "theta_hat" = t(ldaSim::posterior_mean(samples$theta, c("k", "i"))),
    "beta_hat" = t(ldaSim::posterior_mean(samples$beta, c("k", "v")))
  )
}

bootstrap_vb <- function(method, data, B = 3) {

  ## First, make a VB fit, to use as the estimated parameters in the parametric
  ## bootstrap
  f <- stan_model(method)
  vb_fit <- vb(f, data)
  samples <- extract(vb_fit)

  tmp_theta <- tempfile()
  tmp_beta <- tempfile()

  for (b in seq_len(B)) {

    if (b %% 10 == 0) {
      cat(sprintf("Bootstrap iteration %s\n", b))
    }

    cur_data <- data
    cur_data$y <- sim_from_params(theta_hat, beta_hat, data$zero_inf_prob)
    cur_fit <- vb(f, cur_data)

    cur_means <- nmf_posterior_means(extract(cur_fit))

    theta_data <- melt(cur_means$theta_hat, varnames = c("i", "k", "theta"))
    theta_data$iteration = b
    beta_data <- melt(cur_means$beta_hat, varnames = c("v", "k", "beta"))
    beta_data$iteration = b

    write.table(
      theta_data,
      tmp_theta,
      append = TRUE,
      row.names = FALSE,
      col.names = b == 1
    )

    write.table(
      beta_data,
      tmp_theta,
      append = TRUE,
      row.names = FALSE,
      col.names = b == 1
    )
  }

  list(
    "theta" = readr::read_csv(tmp_theta),
    "beta" = readr::read_csv(tmp_beta)
  )
}
