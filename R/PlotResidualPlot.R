#' Create residual plot
#'
#' Create a figure that plots the residuals of the data. These residuals are the
#' difference between expected normal quantiles and observed quantiles.
#'
#' @param centre_width A numeric value between 0.0 and 1.0 that describes the
#'   width of the centre of the data. Can be NULL.
#' @param show_original Show residuals for original, untransformed, data in
#'   addition to transformed data.
#' @param use_alpha Use transparency for points in case the data contains many
#'   instances.
#' @param use_absolute_deviation Plot absolute deviation instead of residuals.
#' @param ggtheme `ggplot2` theme to use for the plot. If not provided,
#'   `ggplot2::theme_light` is used.
#'
#' @inheritParams assess_transformation
#'
#' @return A `ggplot2` plot object for a Q-Q plot.
#' @export
#'
#' @examples
#' x <- exp(stats::rnorm(1000))
#' transformer <- find_transformation_parameters(
#'   x = x,
#'   method = "box_cox"
#' )
#'
#' if (rlang::is_installed("ggplot2")) {
#'   plot_residual_plot(
#'     x = x,
#'     transformer = transformer
#'   )
#'
#'   # Plot only central 80% of the data.
#'   plot_residual_plot(
#'     x = x,
#'     transformer = transformer,
#'     centre_width = 0.80,
#'     show_original = FALSE
#'   )
#' }
plot_residual_plot <- function(
    x,
    transformer,
    centre_width = NULL,
    show_original = TRUE,
    use_alpha = TRUE,
    use_absolute_deviation = TRUE,
    ggtheme = NULL) {

  # Prevent CRAN NOTE due to non-standard use of variables by data.table.
  type <- residual <- z_expected <- NULL

  # Check that packages are present.
  require_package(c("ggplot2", "rlang"), purpose = "to create a residual plot")

  # Perform checks on x.
  .check_data(x)

  # Check the transformer.
  .check_transformer(transformer)

  # Create none-transformer to compare QQ-plot with original data.
  transformer_none <- methods::new("transformationNone")
  transformer_none <- .set_transformation_parameters(
    object = transformer_none,
    x = x)

  # Compute residuals for transformed data.
  residual_data_transformed <- get_residuals(
    x = x,
    transformer = transformer)
  residual_data_transformed[, "type" := "transformed"]

  # Compute residuals for original data.
  residual_data_original <- get_residuals(
    x = x,
    transformer = transformer_none)
  residual_data_original[, "type" := "original"]

  # Combine residual data for plotting.
  data <- data.table::rbindlist(
    list(residual_data_transformed, residual_data_original),
    use.names = TRUE)

  # Check ggtheme.
  ggtheme <- .check_ggtheme(ggtheme)

  # Compute alpha values for the points based on the number of instances.
  n_instances <- nrow(residual_data_transformed)

  if (n_instances < 100L || !use_alpha) {
    point_alpha <- 1.0

  } else if (n_instances > 1E5) {
    point_alpha <- 0.01

  } else {
    # Empirical decay to increasingly smaller alpha.
    point_alpha <- 1.0 - 0.99 * (log10(10000) - 2)^0.1 / 3^0.1
  }

  # Encode type.
  data$type <- factor(data$type, levels = c("original", "transformed"))

  # Keep only transformed data.
  if (!show_original) data <- data[type == "transformed", ]

  # Assign full weight to central elements. We do this to ensure that poor fits
  # don't get down-weighted in the centre, and in that way produce
  # log-likelihood values that are too optimistic.
  if (!is.null(centre_width)) {
    if (centre_width < 0.0 || centre_width > 1.0) {
      stop(paste0("centre_width should be NULL or between 0.0 and 1.0. Found: ", centre_width))
    }

    data <- data[abs(z_expected) <= stats::qnorm(0.50 + centre_width / 2), ]
    if (nrow(data) == 0) {
      stop(paste0("centre_width may be too small: no residuals were selected."))
    }
  }

  # Use absolute values, if required.
  if (use_absolute_deviation) data[, "residual" := abs(residual)]

  # Start plotting.
  p <- ggplot2::ggplot(data = data)
  p <- p + ggtheme

  if (show_original) {

    # Add points.
    p <- p + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$z_expected,
        y = .data$residual,
        colour = .data$type,
        alpha = point_alpha))

    # Add colour scale.
    p <- p + ggplot2::scale_colour_manual(
      name = "data",
      values = c("#4e79a7", "#f28e2b"),
      breaks = c("original", "transformed"))

  } else {

    # Add points.
    p <- p + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$z_expected,
        y = .data$residual,
        alpha = point_alpha))
  }

  # Remove guide for alpha.
  p <- p + ggplot2::guides(alpha = "none")

  # Set x-limits.
  if (is.null(centre_width)) {
    x_limits <- c(-3.0, 3.0)

  } else {
    x_limits <- c(stats::qnorm(0.50 - centre_width / 2), stats::qnorm(0.50 + centre_width / 2))
  }

  # Set y-limits.
  # Check that 0 is contained in the residual plot.
  y_limits <- c(NA_real_, NA_real_)
  if (all(data$residual > 0.0)) {
    y_limits <- c(0.0, NA_real_)

  } else if (all(data$residual < 0.0)) {
    y_limits <- c(NA_real_, 0.0)
  }

  # Set limit.
  p <- p + ggplot2::coord_cartesian(
    xlim = x_limits,
    ylim = y_limits)

  # Set labels.
  p <- p + ggplot2::xlab("Expected theoretical quantiles")
  p <- p + ggplot2::ylab(ifelse(use_absolute_deviation, "Absolute deviation", "Deviation"))

  return(p)
}
