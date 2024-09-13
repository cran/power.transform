#' @include TransformationObjects.R
NULL


#' Assess normality of transformed data
#'
#' Not all data allows for a reasonable transformation to normality using power
#' transformation. For example, uniformly distributed data or multi-modal data
#' cannot be transformed to normality. This function computes a p-value for an
#' empirical goodness of fit test for central normality. A distribution is
#' centrally normal if the central 80% of the data are approximately normally
#' distributed. The null-hypothesis is that the transformed distribution is
#' centrally normal.
#'
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`.
#' @param verbose Sets verbosity of the fubction.
#' @param ... Unused arguments.
#'
#' @inheritParams power_transform
#'
#' @return p-value for empirical goodness of fit test.
#' @export
#'
#' @examples
#' x <- exp(stats::rnorm(1000))
#' transformer <- find_transformation_parameters(
#'   x = x,
#'   method = "box_cox")
#'
#' assess_transformation(
#'   x = x,
#'   transformer = transformer)
assess_transformation <- function(
    x,
    transformer,
    verbose = TRUE,
    ...) {

  # Prevent CRAN NOTE due to non-standard use of variables by data.table.
  p_observed <- alpha <- NULL

  # Compute fit data.
  residual_data <- get_residuals(
    x = x,
    transformer = transformer)

  # The test uses a central portion kappa = 0.80.
  residual_data <- residual_data[p_observed >= 0.10 & p_observed <= 0.90]

  test_statistic_value <- mean(abs(residual_data$residual))

  if (test_statistic_value < min(gof_lookup_table$test_statistic)) return(1.0)
  if (test_statistic_value > gof_lookup_table[alpha == 0.0001]$test_statistic) {

    if (verbose) message("p-value is smaller than 10^-4")

    return(0.0)
  }

  p_value <- stats::spline(
    x = gof_lookup_table$test_statistic,
    y = gof_lookup_table$alpha,
    method = "hyman",
    xout = test_statistic_value
  )$y

  return(p_value)
}



#' Compute residuals of transformation to normality
#'
#' @inheritParams assess_transformation
#'
#' @return A `data.table` containing the expected (according to a normal
#'   distribution) and observed z-scores, and their difference as residuals.
#' @export
#'
#' @examples
#' x <- exp(stats::rnorm(1000))
#' transformer <- find_transformation_parameters(
#'   x = x,
#'   method = "box_cox")
#'
#' residual_data <- get_residuals(
#'   x = x,
#'   transformer = transformer)
get_residuals <- function(
    x,
    transformer,
    ...) {

  # Perform checks on x.
  .check_data(x)

  # Check the transformer.
  .check_transformer(transformer)

  return(.get_fit_data(
    object = transformer,
    x = x))
}



#### .get_fit_data (generic) ---------------------------------------------------
setGeneric(
  ".get_fit_data",
  function(object, ...) standardGeneric(".get_fit_data"))



#### .get_fit_data (general) ---------------------------------------------------
setMethod(
  ".get_fit_data",
  signature(object = "transformationPowerTransform"),
  function(
    object,
    x,
    ...) {

    # Sort x, if required.
    if (is.unsorted(x)) x <- sort(x)

    # Perform transformation.
    y <- ..transform(
      object = object,
      x = x)

    # Compute the expected z-score.
    z_expected <- compute_expected_z(x = x)

    # Compute M-estimates for locality and scale
    robust_estimates <- huber_estimate(y, tol = 1E-3)

    # Avoid division by 0.0.
    if (robust_estimates$sigma == 0.0) robust_estimates$sigma <- 1.0

    # Compute the observed z-score.
    z_observed <- (y - robust_estimates$mu) / robust_estimates$sigma

    # Compute residuals.
    residual <- z_observed - z_expected

    return(data.table::data.table(
      "z_expected" = z_expected,
      "z_observed" = z_observed,
      "residual" = residual,
      "p_observed" = (seq_along(x) - 1 / 3) / (length(x) + 1 / 3)))
  }
)
