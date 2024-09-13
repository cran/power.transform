#' @include WeightingFunctions.R

# ..get_default_weighting_parameters (none) ------------------------------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodNone"),
  function(object, transformer, estimator, ...) {

    # No weights required.
    return(NULL)
  }
)



# ..get_default_weighting_parameters (emp prob, step, general) -----------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodEmpiricalProbabilityStep"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 0.90)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "empirical_probability-step" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (emp prob, triangle, general) -------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodEmpiricalProbabilityTriangle"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 0.85, "k2" = 0.95)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "empirical_probability-triangle" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
        default_values$k2 <- optimal_values$k2
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (emp prob, cosine, general) ---------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodEmpiricalProbabilityCosine"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 0.85, "k2" = 0.95)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "empirical_probability-cosine" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
        default_values$k2 <- optimal_values$k2
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (transformed, step, general) --------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodTransformedStep"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 1.96)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "transformed-step" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (transformed, triangle, general) ----------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodTransformedTriangle"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 0.50, "k2" = 8.00)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "transformed-triangle" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
        default_values$k2 <- optimal_values$k2
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (transformed, cosine, general) ------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodTransformedCosine"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 0.50, "k2" = 8.00)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "transformed-cosine" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
        default_values$k2 <- optimal_values$k2
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (residual, step, general) -----------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodResidualStep"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 2.00)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "residual-step" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (residual, triangle, general) -------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodResidualTriangle"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 0.50, "k2" = 8.00)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "residual-triangle" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
        default_values$k2 <- optimal_values$k2
      }
    }

    return(default_values)
  }
)



# ..get_default_weighting_parameters (residual, cosine, general) ---------------
setMethod(
  "..get_default_weighting_parameters",
  signature(object = "weightingMethodResidualCosine"),
  function(object, transformer, estimator, ...) {

    # Prevent NOTE due to non-standard evaluation.
    name <- method <- estimation_method <- NULL

    default_values <- list("k1" = 0.50, "k2" = 8.00)

    # Check for known values that where obtained for the manuscript.
    if (..requires_shift_scale_optimisation(transformer)) {
      optimal_values <- two_sided_function_parameters[
        name == "residual-cosine" & method == transformer@method & estimation_method == estimator@method
      ]

      if (nrow(optimal_values) == 1) {
        default_values$k1 <- optimal_values$k1
        default_values$k2 <- optimal_values$k2
      }
    }

    return(default_values)
  }
)
