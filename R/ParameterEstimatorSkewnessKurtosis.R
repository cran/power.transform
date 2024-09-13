#' @include ParameterEstimators.R
NULL



# estimatorSkewnessKurtosis definition -----------------------------------------
setClass(
  "estimatorSkewnessKurtosis",
  contains = "estimatorGeneric")

# estimatorJarqueBera definition -----------------------------------------------
setClass(
  "estimatorJarqueBera",
  contains = "estimatorSkewnessKurtosis",
  slots = list("method" = "character"),
  prototype = list("method" = "jarque_bera"))

# estimatorDAgostino definition ------------------------------------------------
setClass(
  "estimatorDAgostino",
  contains = "estimatorSkewnessKurtosis",
  slots = list("method" = "character"),
  prototype = list("method" = "dagostino"))



# .compute_objective (general SK) ---------------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorSkewnessKurtosis"),
  function(object, transformer, x, ...) {

    # Make sure that x is sorted.
    if (is.unsorted(x)) x <- sort(x)

    # Transform feature values
    y <- ..transform(
      object = transformer,
      x = x)

    if (any(!is.finite(y))) return(NA_real_)

    # Compute weights.
    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if (!all(is.finite(w))) return(NA_real_)

    sum_w <- sum(w)
    if (sum_w == 0) return(NA_real_)

    # Compute (weighted) mean and variance, because this is used for skewness
    # and kurtosis.
    mu <- sum(w * y) / sum_w
    sigma_squared <- sum(w * (y - mu)^2) / sum_w

    # Check problematic values.
    if (!is.finite(sigma_squared)) return(NA_real_)
    if (sigma_squared <= .Machine$double.eps) return(NA_real_)

    # Compute (weighted) skewness and kurtosis.
    skewness <- (1.0 / sigma_squared^(3 / 2)) * (sum(w * (y - mu)^3)) / sum_w
    kurtosis <- (1.0 / sigma_squared^2) * (sum(w * (y - mu)^4)) / sum_w

    return(list(
      "skewness" = skewness,
      "kurtosis" = kurtosis,
      "n" = sum_w))
  }
)



# .compute_objective (Jarque-Bera) ---------------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorJarqueBera"),
  function(object, transformer, x, ...) {
    # Based on the Jarque-Bera test statistic.

    # Compute skewness and kurtosis first.
    moments <- callNextMethod()

    if (any(is.na(moments))) return(NA_real_)

    # Compute the interesting part of the test statistic.
    t <- moments$skewness^2 + 0.25 * (moments$kurtosis - 3.0)^2

    # Statistic should be minimised for better fits.
    return(t)
  }
)



# .compute_objective (D'Agostino) ----------------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorDAgostino"),
  function(object, transformer, x, ...) {
    # Based on the D'Agostino's K-squared test statistic.

    # Compute skewness and kurtosis first.
    moments <- callNextMethod()

    if (any(is.na(moments))) return(NA_real_)
    n <- moments$n

    # Division by 0 is possible for n = 2 and n = 3, and the kurtosis-based
    # statistic doesn't work well for n less than 20.
    if (n < 20) return(NA_real_)

    # Skewness test statistic following D'Agostino and Pearson (1973).
    s <- moments$skewness

    beta_1 <- s * sqrt(((n + 1.0) * (n + 3.0)) / (6.0 * (n - 2.0)))
    beta_2 <- 3.0 * (n^2 + 27.0 * n - 70.0) * (n + 1.0) * (n + 3.0) / ((n - 2.0) * (n + 5.0) * (n + 7.0) * (n + 9.0))
    alpha <- sqrt(2.0 / (sqrt(2 * beta_2 - 2) - 2.0 ))
    delta <- 1 / sqrt(log(sqrt(-1 + sqrt(2 * beta_2 - 2))))

    z_s <- delta * log(beta_1 / alpha + sqrt(beta_1^2 / alpha^2 + 1))

    if (!is.finite(z_s)) return(NA_real_)

    # Kurtosis test statistic following Anscombe and Glynn (1983).
    k <- moments$kurtosis

    beta_1 <- 3.0 * (n - 1.0) / (n + 1.0)
    beta_2 <- 24.0 * n * (n - 2.0) * (n - 3.0) / ((n + 1.0)^2 * (n + 3.0) * (n + 5.0))
    beta_3 <- sqrt(6.0 * (n + 3.0) * (n + 5.0) / (n * (n - 2.0) * (n - 3.0))) * 6.0 * (n^2 - 5.0 * n + 2.0) / ((n + 7.0) * (n + 9.0))
    alpha_1 <- 6.0 + 8.0 / beta_3 * (2.0 / beta_3 + sqrt(1 + 4.0 / beta_3^2))
    alpha_2 <- (k - beta_1) / sqrt(beta_2)

    z_k <- sqrt(9.0 * alpha_1 / 2.0) * (1.0 - 2.0 / (9.0 * alpha_1) - ((1.0 - 2.0 / alpha_1) / (1.0 + alpha_2 * sqrt(2.0 / (alpha_1 - 4.0))))^(1 / 3))

    if (!is.finite(z_k)) return(NA_real_)

    # Compute test statistic.
    t <- z_s^2 + z_k^2

    # Statistic should be minimised for better fits.
    return(t)
  }
)



# ..get_default_optimiser_control (Jarque-Bera) --------------------------------
setMethod(
  "..get_default_optimiser_control",
  signature(object = "estimatorJarqueBera"),
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



# ..get_default_optimiser_control (D'Agostino)) --------------------------------
setMethod(
  "..get_default_optimiser_control",
  signature(object = "estimatorDAgostino"),
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



..estimators_jarque_bera <- function() {
  return(c("jarque_bera"))
}

..estimators_dagostino <- function() {
  return(c("dagostino"))
}
