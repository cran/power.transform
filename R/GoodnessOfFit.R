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
#' This function is a wrapper around `ecn.test`.
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

  # Perform test.
  h <- ecn.test(x = x, transformer = transformer)

  if (is.na(h$p_value)) {
    message("p-value could not be determined.")
  } else if (h$p_value > 0.0001 && verbose) {
    message("p-value is smaller than 0.0001")
  }

  return(h$p_value)
}








#' Empirical central normality test
#'
#' Assesses central normality of input data using an empirical test.
#'
#' @param x vector of input data, of at least length 5.
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`. Optional, if present residuals are
#'   determined from x after transformation.
#'
#' @return list with mean absolute error (`tau`), critical value (at
#'   significance level = 0.95) of the test statistic (`tau_critical`) and
#'   p-value (`p_value`) for the empirical central normality test.
#' @export
ecn.test <- function(x, transformer = NULL) {
  return(cn.test(x = x, transformer = transformer, robust = TRUE))
}




#' Central normality test
#'
#' Assesses central normality of input data.
#'
#' @param x vector of input data, of at least length 5.
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`. Optional, if present residuals are
#'   determined from x after transformation.
#' @param robust Determines the use of the strict test for central normality
#'   (`FALSE`) or the more robust version where test statistics where determined
#'   from normally distributed data with 10% outliers (`TRUE`).
#'
#' @return list with mean absolute error (`tau`), critical value (at
#'   significance level = 0.95) of the test statistic (`tau_critical`) and
#'   p-value (`p_value`) for (empirical) central normality test.
#' @export
cn.test <- function(x, transformer = NULL, robust = FALSE) {

  # Prevent CRAN NOTE due to non-standard use of variables by data.table.
  p_observed <- alpha <- NULL

  # To use get_residuals we need to have a transformer. We just obfuscate this
  # by creating a dud transformer.
  if (is.null(transformer)) {
    transformer <- find_transformation_parameters(x = x, method = "none")
  }

  # Compute residual data.
  residual_data <- get_residuals(
    x = x,
    transformer = transformer
  )

  # Select look-up table.
  if (!robust) {
    lookup_table <- central_normality_lookup_table

  } else {
    lookup_table <- empirical_central_normality_lookup_table
  }

  # The test uses a central portion kappa = 0.80.
  residual_data <- residual_data[p_observed >= 0.1 & p_observed <= 0.9]

  # Compute mean residual error for the central portion.
  test_statistic_value <- mean(abs(residual_data$residual))

  # Lookup critical tau at alpha = 0.95.
  test_statistic_critical <- .interpolate_2d(
    list("n" = length(x), "alpha" = 0.05),
    data = lookup_table
  )

  # Lookup p-value.
  p_value <- .interpolate_alpha(
    n_lookup = length(x),
    tau_lookup = test_statistic_value,
    data = lookup_table
  )

  if (length(x) < 5L) {
    p_value <- NA_real_
  }

  return(list(
    "tau" = test_statistic_value,
    "tau_critical" = test_statistic_critical,
    "p_value" = p_value
  ))
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



.interpolate_2d <- function(x_interp, data, limited = TRUE) {
  # Bi-linear interpolation. data: data.table with three columns that is
  # organised as a 3D grid. x_interp: named list with two elements of length 1.
  # limited: whether values outside the grid should be limited to the extremes
  # of the grid.

  # Get names corresponding to x and y variables.
  x_name <- names(x_interp)[1L]
  y_name <- names(x_interp)[2L]
  z_name <- setdiff(colnames(data), c(x_name, y_name))

  # Get interpolation position.
  x_0 <- x_interp[[x_name]][1L]
  y_0 <- x_interp[[y_name]][1L]

  # Find nearest points to interpolation position. Using native data.table
  # routines is faster.
  x <- unique(data[, mget(x_name)][order(get(x_name))])[[x_name]]
  y <- unique(data[, mget(y_name)][order(get(y_name))])[[y_name]]

  # Lower and upper values of x.
  x_l <- utils::tail(x[x < x_0], n = 1L)
  x_u <- utils::head(x[x >= x_0], n = 1L)

  if (length(x_l) == 0L) {
    if (!limited && x_u != x_0) return(NA_real_)
    x_l <- x_u
  }
  if (length(x_u) == 0L) {
    if (!limited) return(NA_real_)
    x_u <- x_l
  }

  # Lower and upper values of y.
  y_l <- utils::tail(y[y < y_0], n = 1L)
  y_u <- utils::head(y[y >= y_0], n = 1L)

  if (length(y_l) == 0L) {
    if (!limited && y_u != y_0) return(NA_real_)
    y_l <- y_u
  }
  if (length(y_u) == 0L) {
    if (!limited) return(NA_real_)
    y_u <- y_l
  }

  # Find values.
  z_ll <- data[get(x_name) == x_l & get(y_name) == y_l, get(z_name)][1L]
  z_ul <- data[get(x_name) == x_u & get(y_name) == y_l, get(z_name)][1L]
  z_lu <- data[get(x_name) == x_l & get(y_name) == y_u, get(z_name)][1L]
  z_uu <- data[get(x_name) == x_u & get(y_name) == y_u, get(z_name)][1L]

  if (x_l == x_u) {
    z_0l <- z_ll
    z_0u <- z_uu

  } else {
    z_0l <- (x_u - x_0) / (x_u - x_l) * z_ll + (x_0 - x_l) / (x_u - x_l) * z_ul
    z_0u <- (x_u - x_0) / (x_u - x_l) * z_lu + (x_0 - x_l) / (x_u - x_l) * z_uu
  }

  if (y_l == y_u) {
    z_00 <- z_0l

  } else {
    z_00 <- (y_u - y_0) / (y_u - y_l) * z_0l + (y_0 - y_l) / (y_u - y_l) * z_0u
  }

  return(z_00)
}



.interpolate_alpha <- function(n_lookup, tau_lookup, data) {
  # Specific interpolation for alpha values based on n and tau. Since tau is not
  # on a grid (it's the surface), we need to interpolate tau (in the
  # lookup-table) for n and every alpha in the lookup-table first. Then we need
  # to interpolate alpha for the provided tau.

  # Prevent CRAN NOTE due to non-standard use of variables by data.table.
  n <- tau <- NULL

  # Interpolate tau for n_lookup.
  data <- data[
    ,
    list(
      "tau" = stats::spline(
        x = n,
        y = tau,
        method = "fmm",
        xout = n_lookup
      )$y
    ),
    by = c("alpha")
  ]

  # Interpolate alpha for tau_lookup
  alpha <- stats::spline(
    x = data$tau,
    y = data$alpha,
    method = "fmm",
    xout = tau_lookup
  )$y

  # Check if lookup is outside the domain of tau.
  if(tau_lookup > max(data$tau)) alpha <- 0.0
  if(tau_lookup < min(data$tau)) alpha <- 1.0

  # Limit alpha to (0, 1) range.
  if (alpha < 0.0) alpha <- 0.0
  if (alpha > 1.0) alpha <- 1.0

  return(alpha)
}
