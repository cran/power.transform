#' @include TransformationObjects.R
NULL

# transformationBoxCox definition ----------------------------------------------

#' Box-Cox transformation object
#'
#' This class is used for Box-Cox transformations.
#'
#' @slot method Main transformation method, i.e. `"box_cox"`.
#' @slot robust Indicates whether a robust version of the Box-Cox transformation
#'   is used to set transformation parameters. The value depends on the `robust`
#'   argument of the `find_transformation_parameters` function.
#' @slot lambda Numeric lambda parameter for the Box-Cox transformation.
#' @slot shift Numeric shift parameter for the Box-Cox transformation. The value
#'   depends on the data used for setting transformation parameters. If all data
#'   are strictly positive, `shift` has a value of `0.0`. When negative or zero
#'   values are present, data are shifted to be strictly positive. If
#'   `invariant=TRUE` in the `find_transformation_parameters` function,
#'   `lambda`, `shift` and `scale` parameters are optimised simultaneously.
#' @slot scale Numeric scale parameter for the Box-Cox transformation. If
#'   `invariant=TRUE` in the `find_transformation_parameters` function,
#'   `lambda`, `shift` and `scale` parameters are optimised simultaneously.
#'   Otherwise, the `scale` parameter has a value of `1.0`.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @seealso [find_transformation_parameters]
#' @export
#' @rdname transformation_box_cox

setClass(
  "transformationBoxCox",
  contains = "transformationPowerTransform",
  slots = list(
    "method" = "character",
    "robust" = "logical",
    "lambda" = "numeric",
    "shift" = "numeric",
    "scale" = "numeric"
  ),
  prototype = list(
    "method" = "box_cox",
    "robust" = FALSE,
    "lambda" = 1.0,
    "shift" = 0.0,
    "scale" = 1.0
  )
)



# transformationBoxCoxInvariant definition -------------------------------------

#' @rdname transformation_box_cox
#' @export
setClass(
  "transformationBoxCoxInvariant",
  contains = "transformationBoxCox")



# .set_transformation_parameters (Box-Cox) -------------------------------------
setMethod(
  ".set_transformation_parameters",
  signature(object = "transformationBoxCox"),
  function(
    object,
    x,
    lambda,
    estimation_method = "mle",
    weighting_function = NULL,
    weighting_function_parameters = NULL,
    optimiser = NULL,
    backup_use_default = TRUE,
    ...
  ) {

    # Set lambda range. If lambda is NULL, set a very wide range.
    if (is.null(lambda)) lambda <- ..get_default_lambda_range(object = object)

    # Pre-scale and shift x. For conventional Box-Cox transformations, this does
    # nothing. For scale and shift-invariant Box-Cox transformations, the data
    # are standardised, and corresponding scale and shift parameters are
    # exported.
    init_transform_data <- ..standardise_data(
      object = object,
      x = x
    )
    x <- init_transform_data$x
    init_shift <- init_transform_data$shift
    init_scale <- init_transform_data$scale

    # Set minimum shift, in case any negative values are present.
    object <- ..set_minimum_shift(
      object = object,
      x = x
    )

    # Set lambda, in case a fixed lambda is provided.
    object <- ..set_lambda(
      object = object,
      lambda = lambda
    )

    # Get optimisation parameters to initialise and configure the optimiser.
    global_optimisation_parameters <- ..optimisation_parameters(
      object = object,
      x = x,
      lambda = lambda,
      estimation_method = estimation_method
    )

    # Skip optimisation if there is nothing to optimise.
    if (is.null(global_optimisation_parameters)) {
      # Add package version.
      object <- .set_version(object = object)

      object@complete <- TRUE

      return(object)
    }

    # Initialise the estimator.
    estimator <- .set_estimator(
      transformer = object,
      estimation_method = estimation_method,
      weighting_function = weighting_function,
      weighting_function_parameters = weighting_function_parameters
    )

    # Set optimiser.
    if (is.null(optimiser)) {
      optimiser <- ..get_default_optimiser(
        object = estimator,
        transformer = object
      )
    }

    # Optimise transformation parameters in a global search.
    optimised_parameters <- .optimise_transformation_parameters(
      object = estimator,
      transformer = object,
      x = x,
      optimiser = optimiser,
      optimisation_parameters = global_optimisation_parameters,
      ...
    )

    # Update parameters for local search.
    local_optimisation_parameters <- ..local_optimisation_parameters(
      object = object,
      global_parameters = global_optimisation_parameters,
      selected_parameters = optimised_parameters
    )

    if (!is.null(local_optimisation_parameters)) {
      # Optimise transformation parameters in a local search.
      local_optimised_parameters <- .optimise_transformation_parameters(
        object = estimator,
        transformer = object,
        x = x,
        optimiser = optimiser,
        optimisation_parameters = local_optimisation_parameters,
        ...
      )

      # Only update optimal parameters obtained during global search if those
      # obtained during local search are better.
      if (!is.null(local_optimised_parameters$value) && !is.null(optimised_parameters$value)) {
        if (is.finite(local_optimised_parameters$value) && is.finite(optimised_parameters$value)) {
          if (local_optimised_parameters$value < optimised_parameters$value) {
            optimised_parameters <- local_optimised_parameters
          }
        }
      } else if (!is.null(local_optimised_parameters$value)) {
        optimised_parameters <- local_optimised_parameters
      }
    }

    if (!is.finite(optimised_parameters$lambda) && estimation_method != "mle") {
      # Fallback option in case a non-MLE method fails.
      estimator <- .set_estimator(
        transformer = object,
        estimation_method = "mle",
        weighting_function = weighting_function,
        weighting_function_parameters = weighting_function_parameters
      )

      # Optimise transformation parameters.
      optimised_parameters <- .optimise_transformation_parameters(
        object = estimator,
        transformer = object,
        x = x,
        optimiser = optimiser,
        optimisation_parameters = global_optimisation_parameters,
        ...
      )
    }

    # Update shift, scale and lambda values with the optimised parameters.
    if (!is.null(optimised_parameters$shift)) {
      if (is.finite(optimised_parameters$shift)) {
        object@shift <- optimised_parameters$shift * init_scale + init_shift

      } else if (!backup_use_default) {
        object@shift <- optimised_parameters$shift * init_scale + init_shift

      } else {
        object@shift <- 0.0
      }
    }

    if (!is.null(optimised_parameters$scale)) {
      if (is.finite(optimised_parameters$scale)) {
        object@scale <- optimised_parameters$scale * init_scale

      } else if (!backup_use_default) {
        object@scale <- optimised_parameters$scale * init_scale

      } else {
        object@scale <- 1.0
      }
    }

    if (!is.null(optimised_parameters$lambda)) {
      if (is.finite(optimised_parameters$lambda)) {
        object@lambda <- optimised_parameters$lambda

      } else if (!backup_use_default) {
        object@lambda <- optimised_parameters$lambda

      } else {
        object@lambda <- 1.0
      }
    }

    # Add package version.
    object <- .set_version(object = object)

    object@complete <- TRUE

    return(object)
  }
)



# .transform (Box-Cox) ---------------------------------------------------------
setMethod(
  ".transform",
  signature(object = "transformationBoxCox"),
  function(object, x, oob_action = "na", ...) {

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Check oob_action
    .check_oob_action(oob_action)

    # Find non-positive and non-finite instances.
    na_entries <- which(x <= object@shift | !is.finite(x))

    if (length(na_entries) > 0 && oob_action == "na") {
      rlang::warn(
        message = paste0(
          "Box-cox power transforms are only defined for strictly positive values. ",
          length(na_entries), " zero, negative, NA or infinite values were found."),
        class = "power_transform_transform_invalid_values"
      )

      y[na_entries] <- NA_real_

    } else if (length(na_entries) > 0 && oob_action == "valid") {
      lower_bound <- ..get_value_bounds(object = object)[1]
      invalid_values <- is.finite(x) & x < lower_bound

      if (any(invalid_values)) {
        rlang::warn(
          message = paste0(
            "Box-cox power transforms are only defined for strictly positive values. ",
            length(na_entries), " zero or negative values were found."),
          class = "power_transform_transform_invalid_values"
        )

        y[invalid_values] <- lower_bound
      }
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if (length(valid_entries) > 0) {
      y[valid_entries] <- ..transform(
        object = object,
        x = x[valid_entries])
    }

    return(y)
  }
)



# ..transform (Box-Cox) --------------------------------------------------------
setMethod(
  "..transform",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    # Subtract shift and divide by scale.
    x <- (x - object@shift) / object@scale

    if (object@lambda == 0) {
      y <- log(x)

    } else {
      y <- (x^object@lambda - 1) / object@lambda
    }

    return(y)
  }
)


# .revert_transformation (Box-Cox) ---------------------------------------------
setMethod(
  ".revert_transform",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    # Revert transformation.
    if (object@lambda == 0) {
      y <- exp(x)

    } else {
      y <- (x * object@lambda + 1)^(1 / object@lambda)
    }

    # Apply scale and shift.
    y <- y * object@scale + object@shift

    return(y)
  })



# ..requires_shift_scale_optimisation (Box-Cox (invariant)) --------------------
setMethod(
  "..requires_shift_scale_optimisation",
  signature(object = "transformationBoxCoxInvariant"),
  function(object, ...) {
    return(TRUE)
  }
)



# ..get_default_shift_range (Box-Cox) ------------------------------------------
setMethod(
  "..get_default_shift_range",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {
    # Find the value that brings the entire distribution to 0. We need to
    # increment slightly to avoid x containing 0s.
    max_value <- min(x, na.rm = TRUE)
    max_value_offset <- 1.0
    min_value_offset <- 1.0

    # Find the typical, non-zero, distance between values. NA values should be
    # removed.
    dx <- unique(diff(sort(x, na.last = NA)))
    dx <- dx[dx > 0.0]

    # It shouldn't happen that all values are the same - but better check it.
    if (length(dx) > 0) max_value_offset <- stats::median(dx)

    # The increment should not grow too much.
    if (max_value_offset > 0.5) max_value_offset <- 0.5

    # Update min_value_offset.
    min_value_offset <- stats::median(x - max_value, na.rm = TRUE)
    if (min_value_offset < 1.0) min_value_offset <- 1.0

    # Set minimum shift value.
    min_value <- max_value - min_value_offset

    # Set maximum shift value.
    max_value <- max_value - max_value_offset

    # Set shift range. Occasionally, the values not be sorted.
    shift_range <- sort(c(min_value, min_value, max_value))

    return(shift_range)
  }
)



# ..get_default_scale_range (Box-Cox) ------------------------------------------
setMethod(
  "..get_default_scale_range",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    scale <- stats::IQR(x)
    if (scale == 0.0) scale <- stats::sd(x)
    if (!is.finite(scale)) scale <- 1.0
    if (scale == 0.0) scale <- 1.0

    # Limit scaling to half the scale up to twice the scale.
    scale_range <- c(scale / 2.0, scale, scale * 2.0)

    return(scale_range)
  }
)



# ..get_default_lambda_range (Box-Cox) -----------------------------------------
setMethod(
  "..get_default_lambda_range",
  signature(object = "transformationBoxCox"),
  function(object, ...) {
    return(c(-100.0, 100.0))
  }
)



# ..standardise_data (Box-Cox (invariant)) -------------------------------------
setMethod(
  "..standardise_data",
  signature(object = "transformationBoxCoxInvariant"),
  function(object, x, ...) {

    shift <- stats::median(x)
    if (!is.finite(shift)) shift <- 0.0

    scale <- stats::IQR(x)
    if (scale == 0.0) scale <- stats::sd(x)
    if (!is.finite(scale)) scale <- 1.0
    if (scale == 0.0) scale <- 1.0

    # Set shift so that the minimum value after scaling is 1.
    shift <- (min(x) / scale - 1.0) * scale

    # Shift and scale x.
    x <- (x - shift) / scale

    return(list(
      "x" = x,
      "shift" = shift,
      "scale" = scale
    ))
  }
)



# ..set_minimum_shift (Box-Cox) ------------------------------------------------
setMethod(
  "..set_minimum_shift",
  signature(object = "transformationBoxCox"),
  function(object, x, ...) {

    # Shift is necessary to avoid transformation with negative or zero values.
    if (all(x > 0.0)) return(object)

    # Warn to notify that zero or negative values are present. Avoid warning if
    # shifts are optimised.
    if (!is(object, "transformationBoxCoxInvariant")) {
      rlang::warn(
        message = paste0(
          "Box-cox power transforms are only defined for strictly positive values. ",
          "One or more zero or negative values were found in \"x\". ",
          "The values are shifted to induce strictly positive values."),
        class = "power_transform_transform_invalid_values"
      )
    }

    # Find shift range.
    x_range <- ..get_default_shift_range(
      object = object,
      x = x)

    # Set shift to its maximum value.
    object@shift <- max(x_range)

    return(object)
  }
)



# ..get_value_bounds (Box-Cox) -------------------------------------------------
setMethod(
  "..get_value_bounds",
  signature(object = "transformationBoxCox"),
  function(object, ...) {
    # By default, the minimum valid value is indicated by shift attribute. All
    # values that are equal to or smaller than the shift attribute, are invalid.
    minimum_shift <- object@shift

    # Increase the minimum shift by a small amount to create the valid lower
    # bound.
    return(c(minimum_shift + 1E-8, Inf))
  }
)



# ..set_lambda (Box-Cox) ------------------------------------------------------
setMethod(
  "..set_lambda",
  signature(object = "transformationBoxCox"),
  function(object, lambda, ...) {

    # Only update lambda if it is not a range.
    if (length(lambda) != 1) return(object)

    # Update lambda.
    object@lambda <- lambda

    return(object)
  }
)



# ..optimisation_parameters (Box-Cox) ------------------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationBoxCox"),
  function(object, lambda, ...) {

    if (length(lambda) == 1) return(NULL)

    return(
      list(
        "initial" = mean(lambda),
        "lower" = min(lambda),
        "upper" = max(lambda),
        "parameter_type" = "lambda"))
  }
)



# ..optimisation_parameters (Box-Cox (invariant)) ------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationBoxCoxInvariant"),
  function(
    object,
    x,
    lambda,
    estimation_method,
    ...
  ) {

    # Set up x-range.
    x_range <- ..get_default_shift_range(
      object = object,
      x = x
    )

    # Set up scale.
    s_range <- ..get_default_scale_range(
      object = object,
      x = x
    )

    if (estimation_method %in% c("currently_no_method")) {
      if (length(lambda) == 1) {
        return(
          list(
            "lower" = c(x_range[1]),
            "initial" = c(x_range[2]),
            "upper" = c(x_range[3]),
            "parameter_type" = c("shift")
          )
        )

      } else {
        return(
          list(
            "lower" = c(x_range[1], min(lambda)),
            "initial" = c(x_range[2], mean(lambda)),
            "upper" = c(x_range[3], max(lambda)),
            "parameter_type" = c("shift", "lambda")
          )
        )
      }
    } else {
      if (length(lambda) == 1) {
        return(
          list(
            "lower" = c(x_range[1], s_range[1]),
            "initial" = c(x_range[2], s_range[2]),
            "upper" = c(x_range[3], s_range[3]),
            "parameter_type" = c("shift", "scale")
          )
        )

      } else {
        return(
          list(
            "lower" = c(x_range[1], s_range[1], min(lambda)),
            "initial" = c(x_range[2], s_range[2], mean(lambda)),
            "upper" = c(x_range[3], s_range[3], max(lambda)),
            "parameter_type" = c("shift", "scale", "lambda")
          )
        )
      }
    }
  }
)



# ..log_likelihood (Box-Cox) ---------------------------------------------------
setMethod(
  "..log_likelihood",
  signature(object = "transformationBoxCox"),
  function(
    object,
    x,
    y,
    w,
    mu_hat,
    sigma_hat_squared,
    ...
  ) {
    # Ensure that data are standardised in the same manner as for
    # transformation.
    x <- (x - object@shift) / object@scale

    # Compute the log likelihood under the assumption that the transformed
    # variable y follows the normal distribution. See manuscript appendix for
    # derivation.
    llf <- -sum(w) / 2.0 * log(2 * pi * sigma_hat_squared) +
           -1.0 / (2.0 * sigma_hat_squared) * sum(w * (y - mu_hat)^2) +
           -sum(w) * log(object@scale) +
           (object@lambda - 1.0) * sum(w * log(x))

    return(llf)
  }
)



# ..first_derivative (Box-Cox) -------------------------------------------------
setMethod(
  "..first_derivative",
  signature(object = "transformationBoxCox"),
  function(object, x) {
    return(x^(object@lambda - 1.0))
  }
)



# ..get_available_estimators (Box-Cox) -----------------------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationBoxCox"),
  function(object, ...) {

    available_estimators <- ..estimators_all()

    # Only allow Raymaekers and Rousseeuw's method for robust optimisation.
    if (!object@robust) {
      available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())
    }

    return(available_estimators)
  }
)


# ..get_available_estimators (Box-Cox (invariant)) -----------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationBoxCoxInvariant"),
  function(object, ...) {

    available_estimators <- ..estimators_all()

    # Raymaekers and Rousseeuw's method for robust optimisation is not suited
    # for simultaneous estimation of shift, scale and lambda parameters.
    available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())

    return(available_estimators)
  }
)



# show (Box-Cox) ---------------------------------------------------------------
setMethod(
  "show",
  signature("transformationBoxCox"),
  function(object) {

    if (object@shift != 0.0 && object@scale != 1.0) {
      invariant_str <- "shift- and scale-invariant "
    } else if (object@shift != 0.0) {
      invariant_str <- "shift-invariant "
    } else if (object@scale != 1.0) {
      invariant_str <- "scale-invariant "
    } else {
      invariant_str <- ""
    }

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      invariant_str,
      "Box-Cox transformation object"
    )

    if (object@complete) {
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")

      if (object@shift != 0.0) cat(paste0("  shift: ", object@shift, "\n"))
      if (object@scale != 1.0) cat(paste0("  scale: ", object@scale, "\n"))

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)
