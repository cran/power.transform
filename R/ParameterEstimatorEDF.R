#' @include ParameterEstimators.R
NULL



# estimatorEmpiricalDistributionFunction definition ----------------------------
setClass(
  "estimatorEmpiricalDistributionFunction",
  contains = "estimatorGeneric")

# estimatorAndersonDarling definition ------------------------------------------
setClass(
  "estimatorAndersonDarling",
  contains = "estimatorEmpiricalDistributionFunction",
  slots = list("method" = "character"),
  prototype = list("method" = "anderson_darling"))

# estimatorCramervonMises definition -------------------------------------------
setClass(
  "estimatorCramervonMises",
  contains = "estimatorEmpiricalDistributionFunction",
  slots = list("method" = "character"),
  prototype = list("method" = "cramer_von_mises"))

# estimatorKolmogorovSmirnov definition -----------------------------------------
setClass(
  "estimatorKolmogorovSmirnov",
  contains = "estimatorEmpiricalDistributionFunction",
  slots = list("method" = "character"),
  prototype = list("method" = "kolmogorov_smirnov"))



# .compute_objective (general EDF) ---------------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorEmpiricalDistributionFunction"),
  function(object, transformer, x, ...) {

    # Prevent NOTE due to non-standard evaluation in data.table.
    p <- NULL

    # Make sure that x is sorted.
    if (is.unsorted(x)) x <- sort(x)

    # Transform feature values
    y <- ..transform(
      object = transformer,
      x = x
    )

    if (any(!is.finite(y))) return(NA_real_)

    # Compute (merged) empirical probabilities. Average empirical probabilities
    # for when x has multiple values. Though this necessitates using the
    # data.table package, this is by far the fastest implementation.
    p_empirical <- (seq_along(x) - 1 / 3) / (length(x) + 1 / 3)

    # Set up a data.table.
    data <- data.table::data.table("x" = x, "p" = p_empirical)
    data[, "p_group" := mean(p), by = c("x")]
    p_empirical <- data$p_group

    # Determine mu and sigma for putative normal distribution from the data.
    # This is necessary to compute the probabilities expected by theory.
    if (transformer@robust) {
      # Compute M-estimates for locality and scale
      estimates <- huber_estimate(y, tol = 1E-3)

    } else {
      estimates <- list("mu" = mean(y), "sigma" = stats::sd(y))
    }

    # Check problematic values.
    if (!is.finite(estimates$sigma)) return(NA_real_)
    if (estimates$sigma <= .Machine$double.eps) return(NA_real_)

    # Compute expected probabilities according to the cumulative density
    # function of the normal distribution parametrised by mu and sigma.
    p_expected <- stats::pnorm(
      q = y,
      mean = estimates$mu,
      sd = estimates$sigma)

    return(list(
      "p_empirical" = p_empirical,
      "p_expected" = p_expected))
  }
)



# .compute_objective (Anderson-Darling) ----------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorAndersonDarling"),
  function(object, transformer, x, ...) {
    # Based on the Anderson-Darling test statistic.

    # Get general EDF data first.
    p <- callNextMethod()

    if (!is.list(p)) return(NA_real_)

    # Prevent numerical issues due to very small or very large expected values.
    p_expected <- p$p_expected
    p_expected[p_expected <= 1E-4] <- 1E-4
    p_expected[p_expected >= 1.0 - 1E-4] <- 1.0 - 1E-4

    # Get weights from the weighting function..
    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if (!all(is.finite(w))) return(NA_real_)

    sum_w <- sum(w)
    if (sum_w == 0.0) return(NA_real_)

    # Get weights for Anderson-Darling.
    w_ad <- 1.0 / (p_expected * (1.0 - p_expected))

    # Compute (weighted) statistic. Multiplication with 1E2 prevents values from
    # becoming vanishingly small.
    t <- sum(w * w_ad * (1E2 * (p$p_empirical - p$p_expected))^2)

    # Normalise by weights. This is to prevent the optimiser from cheating by
    # finding fitting parameters where most of the weights will be (close to)
    # zero.
    t <- t / sum_w

    # Statistic should be minimised for better fits.
    return(t)
  }
)



# .compute_objective (Cramér-von Mises) ----------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorCramervonMises"),
  function(object, transformer, x, ...) {
    # Based on the Cramér-von Mises test statistic.

    # Get general EDF data first.
    p <- callNextMethod()

    if (!is.list(p)) return(NA_real_)

    # Get weights from the weighting function..
    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if (!all(is.finite(w))) return(NA_real_)

    sum_w <- sum(w)
    if (sum_w == 0.0) return(NA_real_)

    # Compute (weighted) statistic. Multiplication with 1E2 is not really
    # required, but does not really hurt either.
    t <- sum(w * (1E2 * (p$p_empirical - p$p_expected))^2)

    # Normalise by weights. This is to prevent the optimiser from cheating by
    # finding fitting parameters where most of the weights will be (close to)
    # zero.
    t <- t / sum_w

    # Statistic should be minimised for better fits.
    return(t)
  }
)



# .compute_objective (Kolmogorov-Smirnov) --------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorKolmogorovSmirnov"),
  function(object, transformer, x, ...) {
    # Based on the Kolmogorov-Smirnov test statistic.

    # Get general EDF data first.
    p <- callNextMethod()

    if (!is.list(p)) return(NA_real_)

    # Get weights from the weighting function..
    w <- .get_weights(
      object = object,
      transformer = transformer,
      x = x)

    if (!all(is.finite(w))) return(NA_real_)
    if (all(w == 0.0)) return(NA_real_)

    # Compute absolute error between empirical and theoretic probabilities.
    t <- 1E2 * abs(p$p_empirical - p$p_expected)

    # Find maximum value. The weights are only used for finding the maximum
    # value, and unlike other EDF-based tests we do not normalise by the sum of
    # weights.
    t <- t[which.max(w * t)]

    # Statistic should be minimised for better fits.
    return(t)
  }
)



# ..get_default_optimiser_control (Anderson-Darling) ---------------------------
setMethod(
  "..get_default_optimiser_control",
  signature(object = "estimatorAndersonDarling"),
  function(object, optimiser, ...) {

    if (optimiser %in% c("direct", "direct-l", "mlsl")) {
      parameters <- list("xtol_rel" = 1E-3, "maxeval" = 500)

    } else if (optimiser %in% c("bobyqa")) {
      parameters <- list("xtol_rel" = 1E-6)

    } else {
      parameters <- list("xtol_rel" = 1E-3, "ftol_abs" = 1E-3)
    }

    return(parameters)
  }
)



# ..get_default_optimiser_control (Cramér-von Mises) ---------------------------
setMethod(
  "..get_default_optimiser_control",
  signature(object = "estimatorCramervonMises"),
  function(object, optimiser, ...) {

    if (optimiser %in% c("direct", "direct-l", "mlsl")) {
      parameters <- list("xtol_rel" = 1E-3, "maxeval" = 500)

    } else if (optimiser %in% c("bobyqa")) {
      parameters <- list("xtol_rel" = 1E-6)

    } else {
      parameters <- list("xtol_rel" = 1E-3, "ftol_abs" = 1E-3)
    }

    return(parameters)
  }
)



# ..get_default_optimiser_control (Kolmogorov-Smirnov) -------------------------
setMethod(
  "..get_default_optimiser_control",
  signature(object = "estimatorKolmogorovSmirnov"),
  function(object, optimiser, ...) {

    if (optimiser %in% c("direct", "direct-l", "mlsl")) {
      parameters <- list("xtol_rel" = 1E-3, "maxeval" = 500)

    } else if (optimiser %in% c("bobyqa")) {
      parameters <- list("xtol_rel" = 1E-6)

    } else {
      parameters <- list("xtol_rel" = 1E-3, "ftol_abs" = 1E-3)
    }

    return(parameters)
  }
)



..estimators_anderson_darling <- function() {
  return(c("anderson_darling"))
}

..estimators_cramer_von_mises <- function() {
  return(c("cramer_von_mises"))
}

..estimators_kolmogorov_smirnov <- function() {
  return(c("kolmogorov_smirnov"))
}
