#' @include ParameterEstimators.R
NULL



# estimatorMaximumLikelihoodEstimation definition ------------------------------
setClass(
  "estimatorMaximumLikelihoodEstimation",
  contains = "estimatorGeneric",
  slots = list("method" = "character"),
  prototype = list("method" = "mle"))



# .compute_objective (MLE) -----------------------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorMaximumLikelihoodEstimation"),
  function(object, transformer, x, ...) {

    # Make sure that x is sorted.
    if (is.unsorted(x)) x <- sort(x)

    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if (!all(is.finite(w))) return(NA_real_)

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

    # Normalise by weights. This is to prevent the optimiser from cheating by
    # finding fitting parameters where most of the weights will be (close to)
    # zero.
    llf <- llf / sum(w)

    # Log-likelihood should be maximised, hence return -llf because optimisers
    # use minimisation.
    return(-llf)
  }
)



# ..get_default_optimiser_control (MLE) ----------------------------------------
setMethod(
  "..get_default_optimiser_control",
  signature(object = "estimatorMaximumLikelihoodEstimation"),
  function(object, optimiser, ...) {

    if (optimiser %in% c("direct", "direct-l", "mlsl")) {
      parameters <- list("xtol_rel" = 1E-3, "maxeval" = 300)

    } else if (optimiser %in% c("bobyqa")) {
      parameters <- list("xtol_rel" = 1E-6)

    } else {
      parameters <- list("xtol_rel" = 1E-3, "ftol_abs" = 1E-5)
    }

    return(parameters)
  }
)



.log_likelihood <- function(transformer, x, y, w) {

  # Compute the sum of the weights.
  sum_w <- sum(w)
  if (sum_w == 0.0) return(NA_real_)

  # Compute the weighted estimates of the mean mu and variance sigma squared
  # for y.
  mu_hat <- sum(w * y) / sum_w
  sigma_hat_squared <- sum(w * (y - mu_hat)^2) / sum_w

  # Log-likelihood cannot be estimated if sigma is NaN, or equals 0.0.
  if (!is.finite(sigma_hat_squared)) return(NA_real_)
  if (sigma_hat_squared <= .Machine$double.eps) return(NA_real_)

  # Compute the log likelihood under the assumption that the transformed
  # variable y follows the normal distribution.
  llf <- ..log_likelihood(
    object = transformer,
    x = x,
    y = y,
    w = w,
    mu_hat = mu_hat,
    sigma_hat_squared = sigma_hat_squared
  )

  return(llf)
}



..estimators_mle <- function() {
  return(c("mle"))
}
