#' @include ParameterEstimators.R
#' @include TransformationObjects.R
NULL


# estimatorRaymaekersRobust definition -----------------------------------------
setClass(
  "estimatorRaymaekersRobust",
  contains = "estimatorGeneric")



# ..get_available_weighting_functions (general) --------------------------------
setMethod(
  "..get_available_weighting_functions",
  signature(transformer = "transformationPowerTransform", estimator = "estimatorRaymaekersRobust"),
  function(transformer, estimator, ...) {
    # This method does not allow for weighting functions, but uses its own
    # internal procedure.
    return("none")
  }
)



# ..get_default_weighting_function (general) -----------------------------------
setMethod(
  "..get_default_weighting_function",
  signature(transformer = "transformationPowerTransform", estimator = "estimatorRaymaekersRobust"),
  function(transformer, estimator, ...) {
    return("none")
  }
)



# ..optimise_transformation_parameters (Raymaekers robust) ---------------------
setMethod(
  ".optimise_transformation_parameters",
  signature(object = "estimatorRaymaekersRobust"),
  function(
    object,
    transformer,
    optimiser,
    optimisation_parameters,
    x,
    ...) {

    # Check that we do not inadvertently pass problems that do not require
    # optimisation to the optimiser.
    if (is.null(optimisation_parameters)) stop(
      "DEV: optimisation_parameters cannot be empty.")
    if (!is(transformer, "transformationPowerTransform")) stop(
      "DEV: transformer should be a valid power transformation object.")
    if (!is.numeric(x)) stop("DEV: x should be numeric.")

    # Perform method-specific checks on proper usage.
    if (!setequal(optimisation_parameters$parameter_type, "lambda")) stop(
      "DEV: Raymaekers and Rousseeuws robust method can only be used for optimising lambda.")
    if (!transformer@robust) stop(
      paste0(
        "DEV: Raymaekers and Rousseuws robust method is intended for robust optimisation, ",
        "but the robust attribute is FALSE."))

    # Sort x if necessary.
    if (is.unsorted(x)) x <- sort(x)

    # Compute z-values according to the inverse cumulative density function.
    z_expected <- compute_expected_z(x = x)

    # Step 1: Compute initial estimate for lambda, based on residuals with
    # rectified residual errors.
    lambda <- suppressWarnings(
      stats::optimise(
        .raymaekers_robust_rectified_optimisation,
        interval = c(optimisation_parameters$lower, optimisation_parameters$upper),
        x = x,
        z = z_expected,
        transformer = transformer,
        maximum = FALSE))

    lambda <- ifelse(
      is.finite(lambda$objective),
      lambda$minimum,
      1.0)

    # Step 2: Compute lambda from reweighted maximum likelihood.
    lambda <- .raymaekers_robust_reweighted_optimisation(
      lambda = lambda,
      x = x,
      transformer = transformer,
      lambda_range = c(optimisation_parameters$lower, optimisation_parameters$upper))

    # Step 3: Compute lambda from reweighted maximum likelihood again.
    lambda <- .raymaekers_robust_reweighted_optimisation(
      lambda = lambda,
      x = x,
      transformer = transformer,
      lambda_range = c(optimisation_parameters$lower, optimisation_parameters$upper))

    if (!is.finite(lambda)) lambda <- NA_real_

    return(list("lambda" = lambda, "value" = NA_real_))
  }
)



.raymaekers_robust_rectified_optimisation <- function(lambda, x, z, transformer) {

  # Update lambda parameter of the transformer locally.
  transformer@lambda <- lambda

  # Compute values after power transformation with the appropriate lambda value.
  y <- ..raymaekers_robust_rectified_transform(
    transformer = transformer,
    x = x)

  # Compute M-estimates for locality and scale
  robust_estimates <- huber_estimate(y, tol = 1E-3)

  # Check problematic values.
  if (!is.finite(robust_estimates$sigma)) return(NA_real_)
  if (robust_estimates$sigma <= .Machine$double.eps) return(NA_real_)

  # Compute residuals.
  residual <- (y - robust_estimates$mu) / robust_estimates$sigma - z

  # Compute Tukey bisquare function to truncate weights of outlier residuals.
  truncated_weights <- numeric(length(residual)) + 1.0
  valid_residuals <- which(abs(residual) <= 0.5)

  if (length(valid_residuals) > 0) {
    truncated_weights[valid_residuals] <- 1.0 - (1.0 - (residual[valid_residuals] / 0.5)^2)^3
  }

  return(sum(truncated_weights))
}



..raymaekers_robust_rectified_transform <- function(transformer,  x) {

  # The rectified transform replaces part of the transformed values by a first
  # order (linear) approximation. Linear approximation of a function f(x) at
  # point a is defined as y = f(a) + (x - a) * f'(a), with f'(a) being the
  # derivative of f(x=a). Here function f is a Box-Cox or Yeo-Johnson
  # transformation, and point a is the first or third quartile, depending on
  # lambda.

  # Find first and third quartiles.
  cut_off <- stats::quantile(x, probs = c(0.25, 0.75), names = FALSE)

  y <- numeric(length(x))

  # Perform rectified transformation
  if (transformer@lambda == 1.0) {
    # For lambda equal to 1, the mapping is linear, and no elements are
    # out-of-range and require rectification.
    out_of_range <- logical(length(x))

  } else if (transformer@lambda > 1.0) {
    # Select the cut-off value, i.e. the first quartile.
    cut_off <- cut_off[1]

    # Elements that have value below Cl (1st quartile) are rectified.
    out_of_range <- x < cut_off

  } else {
    # Lambda < 1.0.

    # Select the cut-off value, i.e. the third quartile.
    cut_off <- cut_off[2]

    # Elements that have value above Cu (3rd quartile) are rectified.
    out_of_range <- x > cut_off
  }

  if (any(out_of_range)) {
    # Linear approximation to out-of-range elements.
    y[out_of_range] <- ..transform(object = transformer, x = cut_off) +
      (x[out_of_range] - cut_off) * ..first_derivative(object = transformer, x = cut_off)
  }

  # Map elements that do not require rectification using the normal Box-Cox
  # transformation.
  if (any(!out_of_range)) {
    y[!out_of_range] <- ..transform(object = transformer, x = x[!out_of_range])
  }

  return(y)

}



.raymaekers_robust_reweighted_optimisation <- function(lambda, x, transformer, lambda_range) {

  # Update lambda parameter of the transformer locally.
  transformer@lambda <- lambda

  # Compute transformed values based on lambda.
  y <- ..transform(
    object = transformer,
    x = x)

  # Compute M-estimates for locality and scale
  robust_estimates <- huber_estimate(y, tol = 1E-3)

  # Check problematic values.
  if (!is.finite(robust_estimates$sigma)) return(NA_real_)
  if (robust_estimates$sigma <= .Machine$double.eps) return(NA_real_)

  # Compute weights.
  weights <- as.numeric(abs(y - robust_estimates$mu) / robust_estimates$sigma <= stats::qnorm(0.99))
  if (sum(weights) == 0.0) return(NA_real_)

  # Compute new optimal lambda.
  lambda <- suppressWarnings(
    stats::optimise(
      ..raymaekers_robust_log_likelihood,
      interval = lambda_range,
      x = x,
      w = weights,
      transformer = transformer,
      maximum = TRUE))

  lambda <- ifelse(
    is.finite(lambda$objective),
    lambda$maximum,
    NA_real_)

  return(lambda)
}



..raymaekers_robust_log_likelihood <- function(
    lambda,
    transformer,
    x,
    w) {

  # Update lambda parameter of the transformer locally.
  transformer@lambda <- lambda

  # Transform x under the provided lambda.
  y <- ..transform(
    object = transformer,
    x = x)

  if (any(!is.finite(y))) return(NA_real_)

  # Compute log-likelihood
  llf <- .log_likelihood(
    transformer = transformer,
    x = x,
    y = y,
    w = w)

  return(llf)
}



..estimators_raymaekers_robust <- function() {
  return("raymaekers_robust")
}
