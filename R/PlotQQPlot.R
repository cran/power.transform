#' Create Q-Q plot
#'
#' Create a figure that plots the expected, theoretical normal quantiles
#' (z-scores) against the observed normal quantiles (z-scores) of the data.
#'
#' @param show_original Show quantiles for original, untransformed, data in
#'   addition to transformed data.
#' @param show_identity Show identity line that indicates equivalence between
#'   expected and observed quantiles.
#' @param use_alpha Use transparency for points in case the data contains many
#'   instances.
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
#'   method = "box_cox")
#'
#' if (rlang::is_installed("ggplot2")) {
#'   plot_qq_plot(
#'     x = x,
#'     transformer = transformer
#'   )
#' }
#'
plot_qq_plot <- function(
    x,
    transformer,
    show_original = TRUE,
    show_identity = TRUE,
    use_alpha = TRUE,
    ggtheme = NULL) {

  # Prevent CRAN NOTE due to non-standard use of variables by data.table.
  type <- NULL

  # Check that packages are present.
  require_package(c("ggplot2", "rlang"), purpose = "to create a QQ-plot")

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

  # Start plotting.
  p <- ggplot2::ggplot(data = data)
  p <- p + ggtheme

  if (show_original) {

    # Add points.
    p <- p + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$z_expected,
        y = .data$z_observed,
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
        y = .data$z_observed,
        alpha = point_alpha))
  }

  # Add identity-line
  if (show_identity) {
    p <- p + ggplot2::geom_abline(
      intercept = 0.0,
      slope = 1.0)
  }

  # Remove guide for alpha.
  p <- p + ggplot2::guides(alpha = "none")

  # Set limit.
  p <- p + ggplot2::coord_cartesian(xlim = c(-3.0, 3.0))

  # Set labels.
  p <- p + ggplot2::xlab("Expected normal quantiles")
  p <- p + ggplot2::ylab("Observed normal quantiles")

  return(p)
}
