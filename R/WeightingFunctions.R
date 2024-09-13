.set_weighting_function <- function(
    transformer,
    estimator,
    weighting_function,
    weighting_function_parameters) {

  # Find those weighting functions that are allowed for the particular
  # combination of transformer and estimator.
  available_weighting_functions <- ..get_available_weighting_functions(
    transformer = transformer,
    estimator = estimator)

  if (is.null(weighting_function)) weighting_function <- ..get_default_weighting_function(
    transformer = transformer,
    estimator = estimator)

  if (!weighting_function %in% available_weighting_functions) {
    stop(paste0(
      "The desired weighting function ", weighting_function,
      " is not available."))
  }

  # Create weighting function object.
  if (weighting_function == "none") {
    weighting_object <- methods::new("weightingMethodNone")

  } else if (weighting_function == "empirical_probability_step") {
    weighting_object <- methods::new("weightingMethodEmpiricalProbabilityStep")

  } else if (weighting_function == "empirical_probability_triangle") {
    weighting_object <- methods::new("weightingMethodEmpiricalProbabilityTriangle")

  } else if (weighting_function == "empirical_probability_cosine") {
    weighting_object <- methods::new("weightingMethodEmpiricalProbabilityCosine")

  } else if (weighting_function == "transformed_step") {
    weighting_object <- methods::new("weightingMethodTransformedStep")

  } else if (weighting_function == "transformed_triangle") {
    weighting_object <- methods::new("weightingMethodTransformedTriangle")

  } else if (weighting_function == "transformed_cosine") {
    weighting_object <- methods::new("weightingMethodTransformedCosine")

  } else if (weighting_function == "residual_step") {
    weighting_object <- methods::new("weightingMethodResidualStep")

  } else if (weighting_function == "residual_triangle") {
    weighting_object <- methods::new("weightingMethodResidualTriangle")

  } else if (weighting_function == "residual_cosine") {
    weighting_object <- methods::new("weightingMethodResidualCosine")

  } else {
    stop(paste0("DEV: the intended weighting method cannot be set: ", weighting_function))
  }

  # Set parameters of the weighting functions.
  default_parameters <- ..get_default_weighting_parameters(
    object = weighting_object,
    transformer = transformer,
    estimator = estimator)

  .check_weighting_function_parameters(
    x = weighting_function_parameters,
    default_parameters = default_parameters)

  for (parameter in names(default_parameters)){
    if (parameter %in% names(weighting_function_parameters)) {
      slot(weighting_object, parameter) <- weighting_function_parameters[[parameter]]

    } else {
      slot(weighting_object, parameter) <- default_parameters[[parameter]]
    }
  }

  return(weighting_object)
}


setClass("weightingSourceGeneric")
setClass("weightingSourceNone", contains = "weightingSourceGeneric")
setClass("weightingSourceEmpiricalProbability", contains = "weightingSourceGeneric")
setClass("weightingSourceTransformed", contains = "weightingSourceGeneric")
setClass("weightingSourceResidual", contains = "weightingSourceGeneric")

setClass("weightingFunctionGeneric")
setClass("weightingFunctionNone", contains = "weightingFunctionGeneric")
setClass(
  "weightingFunctionStep",
  contains = "weightingFunctionGeneric",
  slots = list("k1" = "numeric"),
  prototype = list("k1" = NA_real_))
setClass(
  "weightingFunctionTriangle",
  contains = "weightingFunctionGeneric",
  slots = list("k1" = "numeric", "k2" = "numeric"),
  prototype = list("k1" = NA_real_, "k2" = NA_real_))
setClass(
  "weightingFunctionCosine",
  contains = "weightingFunctionGeneric",
  slots = list("k1" = "numeric", "k2" = "numeric"),
  prototype = list("k1" = NA_real_, "k2" = NA_real_))

# Weighting method classes are a combination of weighting source and function.
setClass(
  "weightingMethodNone",
  contains = c("weightingSourceNone", "weightingFunctionNone"))
setClass(
  "weightingMethodEmpiricalProbabilityStep",
  contains = c("weightingSourceEmpiricalProbability", "weightingFunctionStep"))
setClass(
  "weightingMethodEmpiricalProbabilityTriangle",
  contains = c("weightingSourceEmpiricalProbability", "weightingFunctionTriangle"))
setClass(
  "weightingMethodEmpiricalProbabilityCosine",
  contains = c("weightingSourceEmpiricalProbability", "weightingFunctionCosine"))
setClass(
  "weightingMethodTransformedStep",
  contains = c("weightingSourceTransformed", "weightingFunctionStep"))
setClass(
  "weightingMethodTransformedTriangle",
  contains = c("weightingSourceTransformed", "weightingFunctionTriangle"))
setClass(
  "weightingMethodTransformedCosine",
  contains = c("weightingSourceTransformed", "weightingFunctionCosine"))
setClass(
  "weightingMethodResidualStep",
  contains = c("weightingSourceResidual", "weightingFunctionStep"))
setClass(
  "weightingMethodResidualTriangle",
  contains = c("weightingSourceResidual", "weightingFunctionTriangle"))
setClass(
  "weightingMethodResidualCosine",
  contains = c("weightingSourceResidual", "weightingFunctionCosine"))



# .get_weights (general weighting) ---------------------------------------------
setMethod(
  ".get_weights",
  signature(object = "weightingFunctionGeneric"),
  function(object, transformer, x, ...) {
    # Weights are set by first processing the data, and then applying the
    # weight function to the processed input data.

    y <- .compute_weight_input(
      object = object,
      transformer = transformer,
      x = x)

    w <- .apply_weight_function(
      object = object,
      x = y)

    return(w)
  }
)



# .compute_weight_input (generic) ----------------------------------------------
setGeneric(
  ".compute_weight_input",
  function(object, ...) standardGeneric(".compute_weight_input"))



# .compute_weight_input (none) -------------------------------------------------
setMethod(
  ".compute_weight_input",
  signature(object = "weightingSourceNone"),
  function(object, x, ...) {
    return(x)
  }
)



# .compute_weight_input (empirical probability) --------------------------------
setMethod(
  ".compute_weight_input",
  signature(object = "weightingSourceEmpiricalProbability"),
  function(object, x, ...) {
    # Check if x is sorted.
    if (is.unsorted(x)) stop(paste0("DEV: x is expected to be sorted in ascending order."))

    # Compute empirical probabilities of each point
    p <- (seq_along(x) - 1 / 3) / (length(x) + 1 / 3)

    # Centralise and map to [-1, 1] range.
    p <- 2.0 * (p - 0.5)

    return(p)
  }
)



# .compute_weight_input (transformed) ------------------------------------------
setMethod(
  ".compute_weight_input",
  signature(object = "weightingSourceTransformed"),
  function(object, transformer, x, ...) {
    # Find transformed feature values.
    y <- ..transform(
      object = transformer,
      x = x)

    # Approximate Huber's M-estimates for locality and scale.
    robust_estimates <- huber_estimate(y, tol = 1E-3)

    # Check problematic values.
    if (!is.finite(robust_estimates$sigma)) return(NA_real_)
    if (robust_estimates$sigma <= .Machine$double.eps) return(NA_real_)

    return((y - robust_estimates$mu) / robust_estimates$sigma)
  }
)



# .compute_weight_input (residual) ---------------------------------------------
setMethod(
  ".compute_weight_input",
  signature(object = "weightingSourceResidual"),
  function(object, transformer, x, ...) {
    # Find transformed feature values.
    y <- ..transform(
      object = transformer,
      x = x)

    # Compute the expected z-score.
    z_expected <- compute_expected_z(x = x)

    # Approximate Huber's M-estimates for locality and scale.
    robust_estimates <- huber_estimate(y, tol = 1E-3)

    # Check problematic values.
    if (!is.finite(robust_estimates$sigma)) return(NA_real_)
    if (robust_estimates$sigma <= .Machine$double.eps) return(NA_real_)

    # Compute the observed z-score.
    z_observed <- (y - robust_estimates$mu) / robust_estimates$sigma

    # Compute residuals.
    residual <- z_observed - z_expected

    return(residual)
  }
)



# .apply_weight_function (generic) ---------------------------------------------
setGeneric(
  ".apply_weight_function",
  function(object, ...) standardGeneric(".apply_weight_function"))



# .apply_weight_function (none) ------------------------------------------------
setMethod(
  ".apply_weight_function",
  signature(object = "weightingFunctionNone"),
  function(object, x, ...) {
    # All weights are 1.0.
    w <- rep_len(1.0, length(x))

    return(w)
  }
)



# .apply_weight_function (step) ------------------------------------------------
setMethod(
  ".apply_weight_function",
  signature(object = "weightingFunctionStep"),
  function(object, x, ...) {
    # Set weights using a step window.

    # k1 should be 0 or greater.
    if (object@k1 < 0.0) return(NA_real_)

    # Initialise weights.
    w <- rep_len(1.0, length(x))

    # Set step window.
    w[abs(x) > object@k1] <- 0.0

    return(w)
  }
)



# .apply_weight_function (triangular) ------------------------------------------
setMethod(
  ".apply_weight_function",
  signature(object = "weightingFunctionTriangle"),
  function(object, x, ...) {
    # Set weights using a triangular window.

    # k1 should be 0 or greater, and k2 cannot be smaller than k1.
    if (object@k2 < object@k1) return(NA_real_)
    if (object@k1 < 0.0) return(NA_real_)

    # Initialise weights.
    w <- rep_len(1.0, length(x))

    # Set weights of elements between k1 and k2.
    if (object@k1 != object@k2) {
      lobe_elements <- which(abs(x) >= object@k1 & abs(x) <= object@k2)
      w[lobe_elements] <- 1.0 - (abs(x[lobe_elements]) - object@k1) / (object@k2 - object@k1)
    }

    # Set weights of elements greater than k2.
    w[abs(x) > object@k2] <- 0.0

    return(w)
  }
)



# .apply_weight_function (tapered cosine) --------------------------------------
setMethod(
  ".apply_weight_function",
  signature(object = "weightingFunctionCosine"),
  function(object, x, ...) {
    # Set weights using a tapered cosine window.

    # k1 should be 0 or greater, and k2 cannot be smaller than k1.
    if (object@k2 < object@k1) return(NA_real_)
    if (object@k1 < 0.0) return(NA_real_)

    # Initialise weights.
    w <- rep_len(1.0, length(x))

    # Set weights of elements between k1 and k2.
    if (object@k1 != object@k2) {
      lobe_elements <- which(abs(x) >= object@k1 & abs(x) <= object@k2)
      w[lobe_elements] <- 0.5 + 0.5 * cos((abs(x[lobe_elements]) - object@k1) / (object@k2 - object@k1) * pi)
    }

    # Set weights of elements greater than k2.
    w[abs(x) > object@k2] <- 0.0

    return(w)
  }
)



# ..get_default_weighting_parameters (generic) ---------------------------------
setGeneric(
  "..get_default_weighting_parameters",
  function(object, ...) standardGeneric("..get_default_weighting_parameters"))



..weight_function_external <- function(
    x,
    weight_function,
    ...) {
  # Used for externally showing weight functions. Not used internally by the
  # package, but a convenience function for the manuscript.

  if (weight_function == "step") {
    object <- methods::new("weightingFunctionStep", ...)

  } else if (weight_function == "triangle") {
    object <- methods::new("weightingFunctionTriangle", ...)

  } else if (weight_function == "cosine") {
    object <- methods::new("weightingFunctionCosine", ...)

  } else {
    stop("weight_function was not recognised.")
  }

  return(.apply_weight_function(
    object = object,
    x = x))
}



..weighting_functions_all <- function() {
  return(c(
    "none",
    "empirical_probability_step",
    "empirical_probability_triangle",
    "empirical_probability_cosine",
    "transformed_step",
    "transformed_triangle",
    "transformed_cosine",
    "residual_step",
    "residual_triangle",
    "residual_cosine"))
}
