#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Helpers for writing configuration files specifying experimental setup.
##
## author: kriss1@stanford.edu

#' Write experiment configurations files
#'
#' This generates the JSON file on which all the experiments will be based. It
#' just creates a JSON will the full factorial combinations of input options.
#'
#' @param factors [list] A list of one-dimensional factors on which we compute
#'   the full cross-product.
#' @param path [string] The path to which to save the configurations JSON.
#' @return NULL. side-effects: Writes the configurations JSON to path.
#' @examples
#' sim_factors <- list(
#'   "N" = c(50, 100, 200),
#'   "P" = c(75, 125),
#'   "zero_inf_prob" = c(0, 0.2, 0.5, 0.8)
#' )
#' model_factors <- list(
#'   "inference" = c("gibbs", "vb"),
#'   "method" = c("zinf_nmf", "nmf")
#' )
#' #write_configs(sim_factors, model_factors)
#' @importFrom jsonlite toJSON
#' @export
write_configs <- function(config_df,
                          n_batches = 50,
                          completed_fits = c(),
                          config_path = "config.json",
                          base_id = "fit",
                          output_dir = "./") {
  config_df$batch <- rep(seq_len(n_batches), length.out = nrow(config_df))

  config <- vector(length = nrow(config_df), mode = "list")
  sim_ix <- colnames(config_df) %in% c("N", "P", "prior_params", "zero_inf_prob")
  model_ix <- colnames(config_df) %in% c("inference", "method")

  ## reshape into a form appropriate for the config json
  for (i in seq_len(nrow(config_df))) {
    config[[i]]$sim_opts <- as.list(config_df[i, sim_ix]) %>%
      merge_nmf_opts()

    config[[i]]$model_opts <- as.list(config_df[i, model_ix]) %>%
      merge_model_opts()

    config[[i]]$output_dir <- output_dir
    config[[i]]$id <- sprintf("%s-%d", base_id, i)
    config[[i]]$batch <- config_df[i, "batch"]
  }

  remove_ids <- sapply(config, function(x) x$id) %in% gsub(".rda", "", completed_fits)
  cat(toJSON(config[!remove_ids], auto_unbox = TRUE), file = config_path)
}
