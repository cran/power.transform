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
#' This function is a wrapper around `ecn_test`.
#'
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`.
#' @param kappa Central portion of the distribution.
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
#'   method = "box_cox"
#' )
#'
#' assess_transformation(
#'   x = x,
#'   transformer = transformer
#' )
assess_transformation <- function(
    x,
    transformer,
    kappa = 0.80,
    verbose = TRUE,
    ...
) {

  # Perform test. Note that assess_transformation is interested only on the
  # p-value of the test, hence we try to capture any errors, and handle them
  # gracefully.
  h <- tryCatch(
    ecn_test(
      x = x,
      transformer = transformer,
      kappa = kappa,
      ...
    ),
    error = identity
  )

  if (is(h, "simpleError")) h <- list("p_value" = NA_real_)

  if (is.na(h$p_value)) {
    message("p-value could not be determined.")
  } else if (h$p_value > 0.0001 && verbose) {
    message("p-value is smaller than 0.0001")
  }

  return(h$p_value)
}




#' Empirical central normality test
#'
#' Assesses central normality of input data using an empirical test. The test
#' has the null hypothesis that the input data were sampled from a central
#' normal distribution.
#'
#' @param x vector of input data, of at least length 5.
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`. Optional, if present residuals are
#'   determined from x after transformation.
#' @param tau test statistic value. Usually computed from `x`, but if can be
#'   provided instead of `x`. `tau` cannot be negative.
#' @param n number of samples. Usually the length of `x`, but can be provided
#'   directly with `tau` and `kappa`. `n` must be 5 or greater.
#' @param kappa central portion of the distribution. For `kappa = 1.0` the test
#'   is equivalent to a normality test, such as the Shapiro-Wilk test.
#' @return list with the test statistic (`tau`), the critical value of the test
#'   statistic (at significance level = 0.95) (`tau_critical`) and p-value
#'   (`p_value`) for the empirical central normality test.
#' @export
ecn_test <- function(
    x = NULL,
    transformer = NULL,
    tau = NULL,
    n = NULL,
    kappa = 0.8
) {
  # Prevent CRAN NOTE due to non-standard use of variables by data.table.
  p_observed <- alpha <- alpha_filter <- NULL

  # Checks on input data.
  if (!is.numeric(x) && !is.null(x)) {
    rlang::abort(paste0(
      "x is expected to be a numeric vector of length 5 or longer, or NULL. Found: ",
      paste_s(class(x))
    ))

  }

  if (is.numeric(x)) {
    if (length(x) < 5L) {
      rlang::abort(paste0(
        "x is expected to be a numeric vector of length 5 or longer. Found: ",
        length(x)
      ))
    }

  }

  if (is.null(x) && is.null(tau)) {
    rlang::abort("Either x or tau must be provided. Both were not set.")
  }

  # Checks on tau.
  if (!is.numeric(tau) && !is.null(tau)) {
    rlang::abort(paste0(
      "tau is expected to be a single, positive number, or NULL. Found: ",
      paste_s(class(tau))
    ))
  }

  if (is.numeric(tau)) {
    if (length(tau) != 1L) {
      rlang::abort(paste0(
        "tau is expected to be a single, positive number. Found: ",
        length(tau), " numbers."
      ))
    }

    if (tau < 0.0) {
      rlang::abort(paste0(
        "tau is expected to be a single, positive number. Found: ",
        tau
      ))
    }
  }

  # Checks on n.
  if (!is.numeric(n) && !is.null(n)) {
    rlang::abort(paste0(
      "n is expected to be a single number greater than or equal to 5, or NULL. Found: ",
      paste_s(class(n))
    ))
  }

  if (is.numeric(n)) {
    if (length(n) != 1L) {
      rlang::abort(paste0(
        "n is expected to be a single number greater than or equal to 5. Found: ",
        length(n), " numbers."
      ))
    }

    if (n > 5.0) {
      rlang::abort(paste0(
        "n is expected to be a single number greater than or equal to 5. Found: ",
        n
      ))
    }
  }

  # Checks on kappa.
  if (!is.numeric(kappa)) {
    rlang::abort(paste0(
      "kappa is expected to be a single number between 0.5 and 1.0. Found: ",
      paste_s(class(kappa))
    ))
  }

  if (length(kappa) != 1L) {
    rlang::abort(paste0(
      "kappa is expected to be a single number between 0.5 and 1.0. Found: ",
      length(kappa), " numbers."
    ))
  }

  if (kappa < 0.5 || kappa > 1.0) {
    rlang::abort(paste0(
      "kappa is expected to be a single number between 0.5 and 1.0. Found: ",
      kappa
    ))
  }


  # We need to compute tau_lookup, if it is not provided through tau.
  tau_lookup <- tau
  n_lookup <- n
  kappa_lookup <- kappa
  if (is.null(tau_lookup)) {
    # To use get_residuals we need to have a transformer. We just obfuscate this
    # by creating a dud transformer.
    if (is.null(transformer)) {
      transformer <- find_transformation_parameters(x = x, method = "none")
    }

    # Compute residual data.
    residual_data <- get_residuals(
      x = x,
      transformer = transformer,
      kappa = kappa_lookup
    )

    # The test uses a central portion kappa = 0.80.
    lower_bound <- (1.0 - kappa_lookup) / 2.0
    upper_bound <- 1.0 - (1.0 - kappa_lookup) / 2.0
    residual_data <- residual_data[p_observed >= lower_bound & p_observed <= upper_bound]

    # Compute mean residual error for the central portion.
    tau_lookup <- mean(abs(residual_data$residual))

    # Set n_lookup.
    n_lookup <- length(x)
  }

  # Additional check on n_lookup. If tau was provided by the user, this value
  # is required.
  if (is.null(n_lookup)) {
    rlang::abort(paste0(
      "n is expected to be a single number greater than or equal to 5. Found: ",
      length(n), " numbers."
    ))
  }

  # Import from package data.
  data <- data.table::copy(ecn_lookup_table)

  # Step 1: filter table to include only nearest values for n, and kappa.
  n_close <- levels(data$n)[c(
    max(which(as.numeric(levels(data$n)) <= n_lookup)),
    min(which(as.numeric(levels(data$n)) >= n_lookup))
  )]
  kappa_close <- levels(data$kappa)[c(
    max(which(as.numeric(levels(data$kappa)) <= kappa_lookup)),
    min(which(as.numeric(levels(data$kappa)) >= kappa_lookup))
  )]

  # Select data.
  new_data <- data[n %in% n_close & kappa %in% kappa_close]

  ..alpha_filter <- function(tau, tau_lookup) {
    x <- logical(length(tau))

    # Find those elements that are closest to tau_lookup.
    y <- tau - tau_lookup
    x[which(y == min(y[y >= 0.0]))] <- TRUE
    x[which(y == max(y[y <= 0.0]))] <- TRUE

    return(x)
  }

  # Determine which alpha values are of interest given the value of tau, and
  # filter again to avoid superfluous interpolation.
  new_data[, "alpha_filter" := ..alpha_filter(tau, tau_lookup), by = c("n", "kappa")]
  selected_alpha <- new_data[alpha_filter == TRUE]$alpha
  selected_alpha <- c(min(selected_alpha), max(selected_alpha))

  # Filter alpha values, but always include 0.05 because we need it for the
  # critical test statistic.
  new_data <- new_data[data.table::between(alpha, selected_alpha[1L], selected_alpha[2L]) | alpha == "0.05"]

  # Convert n, and kappa to numeric values. They were originally
  # stored as factors to facilitate lookup and compress the data.
  new_data$n <- as.numeric(levels(new_data$n)[as.integer(new_data$n)])
  new_data$kappa <- as.numeric(levels(new_data$kappa)[as.integer(new_data$kappa)])

  # For each alpha, interpolate tau given n, and kappa.
  tau_interp <- sapply(
    split(new_data, new_data$alpha, drop = TRUE),
    function(data, n_lookup, tau_lookup) {
      return(.interp_2d(
        f = data$tau,
        x = data$n,
        y = data$kappa,
        x_c = n_lookup,
        y_c = kappa_lookup
      ))
    },
    n_lookup = n_lookup,
    tau_lookup = tau_lookup
  )
  tau_data <- data.table::data.table(
    tau = tau_interp,
    alpha = factor(names(tau_interp), levels = levels(new_data$alpha))
  )

  # In tau_data, interpolate over tau to find the corresponding alpha value.
  # Interpolate alpha for tau_lookup
  alpha_out <- stats::spline(
    x = tau_data$tau,
    y = as.numeric(levels(tau_data$alpha)[as.integer(tau_data$alpha)]),
    method = "fmm",
    xout = tau_lookup
  )$y

  # Check if lookup is outside the domain of tau.
  if(tau_lookup > max(tau_data$tau)) alpha_out <- 0.0
  if(tau_lookup < min(tau_data$tau)) alpha_out <- 1.0

  # Limit alpha to (0, 1) range.
  if (alpha_out < 0.0) alpha_out <- 0.0
  if (alpha_out > 1.0) alpha_out <- 1.0

  # Determine critical tau.
  tau_critical <- tau_data[alpha == "0.05"]$tau

  return(list(
    "tau" = tau_lookup,
    "tau_critical" = tau_critical,
    "p_value" = alpha_out
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
    kappa = 0.80,
    ...
) {

  # Perform checks on x.
  .check_data(x)

  # Check the transformer.
  .check_transformer(transformer)

  return(.get_fit_data(
    object = transformer,
    x = x,
    kappa = kappa
  ))
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
    kappa = 0.8,
    ...
  ) {

    # Sort x, if required.
    if (is.unsorted(x)) x <- sort(x)

    # Perform transformation.
    y <- ..transform(
      object = object,
      x = x)

    # Compute the expected z-score.
    z_expected <- compute_expected_z(x = x)

    # Set k for the robust Huber's M-estimates based on kappa.
    cutoff_k <- stats::qnorm(0.5 + kappa / 2.0)
    if (is.infinite(cutoff_k)) cutoff_k <- 10.0

    # Compute M-estimates for locality and scale.
    robust_estimates <- huber_estimate(y, k = cutoff_k, tol = 1E-3)

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
      "p_observed" = (seq_along(x) - 1.0 / 3.0) / (length(x) + 1.0 / 3.0)
    ))
  }
)



.interp_2d <- function(f, x, y, x_c, y_c) {
  # Helper function for 2d interpolation.
  # f, x, y are numeric vectors of the same length, with f the values of the
  # lattice points, and x and y their coordinates. x_c, y_c are the coordinates
  # of the interpolation point.

  ..get_coord <- function(x, x_c) {
    x_0 <- min(x)
    x_1 <- max(x)
    if (x_1 == x_0) {
      x_d <- x_0
    } else {
      x_d <- (x_c - x_0) / (x_1 - x_0)
    }

    return(x_d)
  }

  ..get_value <- function(f, x, y) {
    f_out <- numeric(4L)
    f[1L] <- f[which(x == min(x) & y == min(y))]  # 00
    f[2L] <- f[which(x == max(x) & y == min(y))]  # 10
    f[3L] <- f[which(x == min(x) & y == max(y))]  # 01
    f[4L] <- f[which(x == max(x) & y == max(y))]  # 11

    return(f)
  }

  # Find values at each lattice point, in an organised manner.
  f <- ..get_value(f, x, y)

  # Find coordinates to interpolate at.
  x_d <- ..get_coord(x, x_c)
  y_d <- ..get_coord(y, y_c)

  f_out <-
    f[1L] * (1.0 - x_d) * (1.0 - y_d) +
    f[2L] * x_d         * (1.0 - y_d) +
    f[3L] * (1.0 - x_d) * y_d         +
    f[4L] * x_d         * y_d

  return(f_out)
}



.interp_3d <- function(f, x, y, z, x_c, y_c, z_c) {
  # Helper function for 3d interpolation.
  # f, x, y, z are numeric vectors of the same length, with f the values of
  # the lattice points, and x, y, and z their coordinates. x_c, y_c, z_c are
  # the coordinates of the interpolation point.

  ..get_coord <- function(x, x_c) {
    x_0 <- min(x)
    x_1 <- max(x)
    if (x_1 == x_0) {
      x_d <- x_0
    } else {
      x_d <- (x_c - x_0) / (x_1 - x_0)
    }

    return(x_d)
  }

  ..get_value <- function(f, x, y, z) {
    f_out <- numeric(8L)
    f[1L] <- f[which(x == min(x) & y == min(y) & z == min(z))]  # 000
    f[2L] <- f[which(x == max(x) & y == min(y) & z == min(z))]  # 100
    f[3L] <- f[which(x == min(x) & y == max(y) & z == min(z))]  # 010
    f[4L] <- f[which(x == max(x) & y == max(y) & z == min(z))]  # 110
    f[5L] <- f[which(x == min(x) & y == min(y) & z == max(z))]  # 001
    f[6L] <- f[which(x == max(x) & y == min(y) & z == max(z))]  # 101
    f[7L] <- f[which(x == min(x) & y == max(y) & z == max(z))]  # 011
    f[8L] <- f[which(x == max(x) & y == max(y) & z == max(z))]  # 111

    return(f)
  }

  # Find values at each lattice point, in an organised manner.
  f <- ..get_value(f, x, y, z)

  # Find coordinates to interpolate at.
  x_d <- ..get_coord(x, x_c)
  y_d <- ..get_coord(y, y_c)
  z_d <- ..get_coord(z, z_c)

  f_out <-
    f[1L] * (1.0 - x_d) * (1.0 - y_d) * (1.0 - z_d) +
    f[2L] * x_d         * (1.0 - y_d) * (1.0 - z_d) +
    f[3L] * (1.0 - x_d) * y_d         * (1.0 - z_d) +
    f[4L] * x_d         * y_d         * (1.0 - z_d) +
    f[5L] * (1.0 - x_d) * (1.0 - y_d) * z_d +
    f[6L] * x_d         * (1.0 - y_d) * z_d +
    f[7L] * (1.0 - x_d) * y_d         * z_d +
    f[8L] * x_d         * y_d         * z_d

  return(f_out)
}
