#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
## Helper functions for plotting the output of an NMF simulation experiment.

#' Get Score Means
#'
#' This groups the main scores matrix by sample index, and average posterior
#' samples.
#'
#' @param scores [data.frame] The data.frame giving the scores for true and
#'   simulated data.
#' @param grouping_vars [character vector] The names of columns on which to
#'   group. Everything else will be averaged over.
#' @return scores [data.frame] The original data.frame with the value_1 and
#'   value_2 columns averaged depending on grouping variables.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by_ summarise
#' @export
score_means <- function(scores, grouping_vars) {
  scores %>%
    group_by_(.dots = grouping_vars) %>%
    summarise(mean_1 = mean(value_1), mean_2 = mean(value_2), truth_1 = truth_1[1], truth_2 = truth_2[1])
}

#' Plot Contours and Associated Coordinate
#'
#' This is a wrapper of ggcontours in the ggscaffold package that shows the
#' coordinate associated with point clouds.
#'
#' @param plot_data [data.frame] The data frame that contains the sampled
#'   coordinates for the contours, along with the true positions (which will be
#'   indicated with text).
#' @param plot_opts [list] A list specifying plotting appearance. Pretty much
#'   the same as the usual ggcontours plot_opts, except an extra option for plot
#'   limits.
#' @return [list] A list containing two ggplot2 objects, for the faceted and
#'   unfaceted scores plots.
#' @importFrom ggscaffold ggcontours
#' @importFrom ggplot2 geom_text aes_string scale_x_sqrt scale_y_sqrt theme
#'   facet_wrap
#' @importFrom grid unit
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export
scores_contours <- function(plot_data, plot_opts) {
  p1 <- ggcontours(plot_data, plot_opts) +
    geom_text(
      data = score_means(plot_data, c(plot_opts$group, plot_opts$facet_terms)),
      aes_string(x = "mean_1", y = "mean_2", label = plot_opts$group),
      col = plot_opts$mean_col,
      size = plot_opts$text_size
    ) +
    geom_text(
      data = plot_data %>% filter(iteration == 1),
      aes_string(x = "truth_1", y = "truth_2", label = plot_opts$group),
      size = plot_opts$text_size
    ) +
    scale_x_sqrt(limits = plot_opts$x_lim, expand = c(0, 0)) +
    scale_y_sqrt(limits = plot_opts$y_lim, expand = c(0, 0)) +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  p2 <- p1 +
    facet_wrap(formula(paste0("~", plot_opts$group))) +
    theme(plot.margin = unit(c(0,0,0,0), "mm"))
  list("grouped" = p1, "coordinates" = p2)
}

#' Errors histogram
#'
#' Plot the histograms of errors associated with the scatterplots from the NMF fits.
#'
#' @param plot_data [data.frame] The data used to plot the error between truth
#'   vs. estimate across all dimensions. See the output of
#'   melt_reshaped_samples().
#' @param facet_terms [character vector] The columns on which to facet_grid the
#'   plot.
#' @param n_bins [int] The number of bins in each histogram panel. Defaults to 75.
#' @param alpha [numeric] The alpha transparency for the different factors.
#' @param colors [character vector] The colors to use for each factor.
#' @return hist_plot [ggplot] The ggplot object showing error histograms across
#'   factors and simulation configurations.
#' @importFrom ggplot2 ggplot geom_histogram aes facet_grid scale_y_continuous
#'   scale_fill_manual theme
#' @importFrom ggscaffold min_theme
#' @importFrom scales pretty_breaks
#' @export
error_histograms <- function(plot_data,
                             facet_terms = NULL,
                             n_bins = 75,
                             alpha = 0.7,
                             colors = c("#d95f02", "#7570b3")) {
  ggplot(plot_data) +
    geom_histogram(
      aes(x = sqrt(estimate) - sqrt(truth), fill = dimension, y = ..density..),
      position = "identity", alpha = alpha, bins = n_bins
    ) +
    facet_grid(formula(paste(facet_terms, collapse = "~"))) +
    scale_y_continuous(breaks = pretty_breaks(3)) +
    scale_fill_manual(values = colors) +
    min_theme() +
    theme(
      legend.position = "bottom"
    )
}
