
#' @title Huber M-estimate
#'
#' @description Iteratively computes M-estimates for location and scale. These
#'   are robust estimates of the mean and standard deviation of the data.
#'
#' @param x Vector of numeric values for which the location and scale should be
#'   estimated.
#' @param k Numeric value > 0 that the determines the value beyond which the
#'   signal is winsorized.
#' @param tol Tolerance for the iterative procedure.
#'
#' @return list with location estimate `"mu"` and scale estimate `"sigma"`.
#' @export
huber_estimate <- function(x, k = 1.28, tol = 1E-4) {
  # k=1.28 is based on Wilcox RR. Introduction to Robust Estimation and
  # Hypothesis Testing. Academic Press; 2011.

  # Filter missing values
  x <- x[is.finite(x)]

  if (length(x) == 0) return(list("mu" = NA_real_, "sigma" = NA_real_))

  # Initial estimates for the estimate mu and scale sigma
  mu_0 <- stats::median(x)
  sigma_0 <- stats::mad(x)

  # Check that sigma_0 is not 0.
  if (sigma_0 == 0.0) sigma_0 <- stats::sd(x)

  # Check that there is an initial estimate for scale.
  if (sigma_0 == 0.0) return(list("mu" = mu_0, "sigma" = 0.0))

  # Determine theta, i.e. the probability corresponding to k, according to a
  # Gaussian distribution.
  theta <- 2.0 * stats::pnorm(k) - 1.0
  beta <- theta + k^2 * (1 - theta) - 2.0 * k * stats::dnorm(k)

  for (ii in seq_len(50)) {
    # Apply Huber's rho, that basically winsorizes values |x| >= k * sigma_0
    xx <- pmin(pmax(mu_0 - k * sigma_0, x), mu_0 + k * sigma_0)

    # Compute updated mean and scale estimate.
    mu_1 <- sum(xx) / length(xx)
    sigma_1 <- sqrt(sum((xx - mu_1)^2) / (beta * (length(xx) - 1)))

    # Check that mu_1 and sigma_1 are finite.
    if (!is.finite(mu_1) || !is.finite(sigma_1)) return(list("mu" = NA_real_, "sigma" = NA_real_))

    # Check convergence.
    if (abs(mu_0 - mu_1) < tol * mu_0 && abs(sigma_0 - sigma_1) < tol * sigma_0) break

    # Values for next iteration
    mu_0 <- mu_1
    sigma_0 <- sigma_1
  }

  return(list("mu" = mu_0, "sigma" = sigma_0))
}



compute_expected_z <- function(x) {

  # Check if x is sorted, and sort otherwise.
  if (is.unsorted(x)) stop(paste0("DEV: x is expected to be sorted in ascending order."))

  # Compute expected quantile.
  p <- (seq_along(x) - 1 / 3) / (length(x) + 1 / 3)

  # Set up a data.table.
  data <- data.table::data.table(
    "x" = x,
    "p" = p)

  # Average quantile for when x has multiple values. Though this necessitates
  # using the data.table package, this is by far the fastest implementation.
  data[, "p_group" := mean(p), by = "x"]

  # Compute z-scores.
  z <- stats::qnorm(p = data$p_group)

  return(z)
}




is_package_installed <- function(name) {
  # Try to obtain the package version. This perhaps the cleanest way to check
  # whether a package exists. require and requireNameSpace attach and load
  # packages, which is not required here. The find.package documentation
  # actively discourages its use to identify whether a package is installed.
  #
  # Originally from the familiar R package, under the EUPL license.

  if (length(name) == 0) return(TRUE)

  installed_version <- rlang::env_cache(
    power_transform_global_env,
    paste0(name, "_version"),
    default = tryCatch(
      utils::packageVersion(name),
      error = identity
    )
  )

  return(!inherits(installed_version, "error"))
}



require_package <- function(x, purpose = NULL, call = rlang::caller_env(), ...) {

  # Check whether packages are installed, without loading the
  # packages.
  package_loaded <- sapply(x, is_package_installed)

  # Skip further analysis if all packages are present.
  if (all(package_loaded)) return(invisible(TRUE))

  # Find all packages that are missing.
  x <- x[!package_loaded]

  # Select unique packages.
  x <- unique(x)

  # Write error message.
  message_str <- paste0(
    "The following package",
    ifelse(length(x) > 1, "s have", " has"),
    " to be installed",
    ifelse(is.null(purpose), ": ", paste0(" ", purpose, ": ")),
    paste_s(x), ".")

  rlang::abort(message_str, call = call)
}



paste_s <- function(...) {
  # Function to collapse a series of strings into a summation in the form
  # "element_1, element_2, ..., and element_n".
  #
  # Originally from the familiar R package, under the EUPL license.
  dots <- c(...)

  if (length(dots) > 2) {
    # For more than 2 elements, split into an initial and final section.
    initial_string <- paste0(utils::head(dots, n = length(dots) - 2L), collapse = ", ")
    final_string <- paste0(utils::tail(dots, n = 2L), collapse = " and ")

    return(paste0(c(initial_string, final_string), collapse = ", "))

  } else if (length(dots) == 2) {
    # For exactly 2 elements, combine with "and".
    return(paste0(dots, collapse = " and "))

  } else {
    # For only one element, return as is.
    return(paste0(dots))
  }
}


#' Random Values from the Asymmetric Generalised Normal Distribution
#'
#' Draws random values from an asymmetric generalised normal distribution.
#'
#' @param n number of instances
#' @param location central location of the distribution
#' @param scale scale of the distribution. Must be strictly positive: `scale >
#'   0.0`
#' @param alpha value between 0.0 and 1.0 that determines the skewness of the
#'   distribution. `alpha > 0.5` creates a distribution with a negative skew
#'   (left-skewed), i.e. the left tail of the distribution is elongated, and the
#'   bulk of the distribution is located to the right. `alpha < 0.5` creates a
#'   distribution with a positive skew (right-skewed), i.e. the right tail of
#'   the distribution is elongated, and the bulk of the distribution is located
#'   to the left. For `alpha = 0.0`, the distribution does not have a skew.
#' @param beta Strictly positive value (`beta > 0.0`) that determines the
#'   overall shape of the generalised normal distribution. For `beta = 1`, an
#'   asymmetric Laplace distribution is used. `beta = 2` draws values according
#'   to an asymmetric normal distribution. For large `beta` the distribution
#'   will approximate the uniform distribution.
#'
#' @details Random values drawn according to an asymmetric generalised normal
#'   distribution. Here the asymmetric generalised normal distribution is a
#'   symmetric general normal distribution, that is made asymmetric using the
#'   procedure described by Gijbels et al. To generate random values we use the
#'   quantile function of the symmetric generalised normal distribution that was
#'   derived by M. Griffin.
#'
#'   The default parameter values produce values as if drawn from the standard
#'   normal distribution with \eqn{\sigma = \sqrt{2}}, that is, the standard
#'   deviation is not \eqn{\sqrt{2}} instead of \eqn{1}.
#'
#' @return One or more numeric values drawn from the asymmetric generalised
#'   normal distribution.
#' @export
#'
#' @references
#' 1. Gijbels I, Karim R, Verhasselt A. Quantile Estimation in a Generalized
# Asymmetric Distributional Setting. Stochastic Models, Statistics and Their
# Applications. Springer International Publishing; 2019. pp. 13–40
#'
#' 1. Griffin M (2018). gnorm: Generalized Normal/Exponential Power Distribution.
# R package version 1.0.0
#'
#' @examples
#' # Draw values from a standard normal distribution.
#' x <- power.transform::ragn(n = 10000, scale = 1/sqrt(2))
#' hist(x, 50)
#'
#' # Draw values from a left-skewed normal distribution.
#' x <- power.transform::ragn(n = 10000, scale = 1/sqrt(2), alpha = 0.8)
#' hist(x, 50)
#'
#' # Draw values from a right-skewed normal distribution.
#' x <- power.transform::ragn(n = 10000, scale = 1/sqrt(2), alpha = 0.2)
#' hist(x, 50)
#'
#' # Draw values from a standard laplace distribution.
#' x <- power.transform::ragn(n = 10000, scale = 1/sqrt(2), beta = 1.0)
#' hist(x, 50)
ragn <- function(n, location = 0, scale = 1, alpha = 0.5, beta = 2) {
  # Random values drawn according to an asymmetric generalised normal
  # distribution. Here the asymmetric generalised normal distribution is a
  # symmetric general normal distribution, that is made asymmetric using the
  # procedure described by Gijbels et al. To generate random values we require
  # the quantile function of the symmetric generalised normal distribution,
  # which is was derived by M. Griffin.

  # alpha needs to be in (0, 1).
  if (alpha < 0.0 || alpha > 1.0) stop("alpha should be between 0.0 and 1.0.")

  # Beta cannot be zero or negative.
  if (beta <= 0.0) stop("beta must be strictly positive.")

  # Scale cannot be zero or negative.
  if (scale <= 0.0) stop("scale must be strictly positive.")

  # Draw n values from an uniform distribution. Avoid extreme tails for
  # numerical reasons.
  p <- stats::runif(
    n = n,
    min = sqrt(.Machine$double.eps),
    max = 1 - sqrt(.Machine$double.eps))

  # Avoid extreme alpha. There are singularities at 0 and 1.
  if (alpha < sqrt(.Machine$double.eps)) {
    alpha <- sqrt(.Machine$double.eps)

  } else if (alpha > 1.0 - sqrt(.Machine$double.eps)) {
    alpha <- 1.0 - sqrt(.Machine$double.eps)
  }

  # Identify values of p which are less and greater than alpha.
  p_less <- which(p <= alpha)
  p_greater <- which(p > alpha)

  # Compute p for the quantile function of the symmetric generalised normal
  # distribution.
  p_sgnd <- numeric(n)
  p_sgnd[p_less] <- p[p_less] / (2.0 * alpha)
  p_sgnd[p_greater] <- (1.0 + p[p_greater] - 2.0 * alpha) / (2.0 * (1.0 - alpha))

  # Compute values according to the quantile function of the symmetric
  # generalised normal distribution.
  x_sgnd <- sign(p_sgnd - 0.5) * stats::qgamma(
    p = 2.0 * abs(p_sgnd - 0.5),
    shape = 1.0 / beta,
    scale = 1.0)^(1.0 / beta)

  # Populate array of values x.
  x <- numeric(n)
  x[p_less] <- location + scale / (1.0 - alpha) * x_sgnd[p_less]
  x[p_greater] <- location + scale / alpha * x_sgnd[p_greater]

  return(x)
}
