#' @include TransformationObjects.R
NULL

# transformationYeoJohnson definition ------------------------------------------

#' Yeo-Johnson transformation object
#'
#' This class is used for Yeo-Johnson transformations.
#'
#' @slot method Main transformation method, i.e. `"yeo_johnson"`.
#' @slot robust Indicates whether a robust version of the Yeo-Johnson
#'   transformation is used to set transformation parameters. The value depends
#'   on the `robust` argument of the `find_transformation_parameters` function.
#' @slot lambda Numeric lambda parameter for the Yeo-Johnson transformation.
#' @slot shift Numeric shift parameter for the Yeo-Johnson transformation. If
#'   `invariant=TRUE` in the `find_transformation_parameters` function,
#'   `lambda`, `shift` and `scale` parameters are optimised simultaneously.
#'   Otherwise, the `shift` parameter has a value of `0.0`.
#' @slot scale Numeric scale parameter for the Yeo-Johnson transformation. If
#'   `invariant=TRUE` in the `find_transformation_parameters` function,
#'   `lambda`, `shift` and `scale` parameters are optimised simultaneously.
#'   Otherwise, the `scale` parameter has a value of `1.0`.
#' @slot complete Indicates whether transformation parameters were set.
#'
#' @seealso [find_transformation_parameters]
#' @export
#' @rdname transformation_yeo_johnson

setClass(
  "transformationYeoJohnson",
  contains = "transformationPowerTransform",
  slots = list(
    "method" = "character",
    "robust" = "logical",
    "lambda" = "numeric",
    "shift" = "numeric",
    "scale" = "numeric"
  ),
  prototype = list(
    "method" = "yeo_johnson",
    "robust" = FALSE,
    "lambda" = 1.0,
    "shift" = 0.0,
    "scale" = 1.0
  )
)



# transformationYeoJohnsonInvariant definition ---------------------------------

#' @rdname transformation_yeo_johnson
#' @export
setClass(
  "transformationYeoJohnsonInvariant",
  contains = "transformationYeoJohnson")



# .set_transformation_parameters (Yeo-Johnson) ---------------------------------
setMethod(
  ".set_transformation_parameters",
  signature(object = "transformationYeoJohnson"),
  function(
    object,
    x,
    lambda,
    estimation_method = "mle",
    weighting_function = NULL,
    weighting_function_parameters = NULL,
    optimiser = NULL,
    backup_use_default = TRUE,
    ...) {

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

    if (!is.finite(optimised_parameters$lambda) & estimation_method != "mle") {
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



# .transform (Yeo-Johnson) -----------------------------------------------------
setMethod(
  ".transform",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {

    # Set up vector to copy new values to.
    y <- numeric(length(x))

    # Find non-finite instances.
    na_entries <- which(!is.finite(x))
    if (length(na_entries) > 0) {
      rlang::warn(
        message = paste0(
          "Yeo-Johnson transforms are only defined for finite values. ",
          length(na_entries), " NA or infinite values were found."),
        class = "power_transform_transform_invalid_values"
      )

      y[na_entries] <- NA_real_
    }

    # Find remaining valid instances.
    valid_entries <- setdiff(seq_along(x), na_entries)

    # Perform transformation.
    if (length(valid_entries) > 0) y[valid_entries] <- ..transform(
      object = object,
      x = x[valid_entries])

    return(y)
  }
)



# ..transform (Yeo-Johnson) ----------------------------------------------------
setMethod(
  "..transform",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {

    # Subtract shift and divide by scale.
    x <- (x - object@shift) / object@scale

    # Determine positive and negative elements of the input vector
    pos_index <- x >= 0
    neg_index <- x < 0

    # Initialise y.
    y <- rep(NA_real_, length(x))

    if (any(pos_index)) {
      if (object@lambda == 0.0) {
        y[pos_index] <- log1p(x[pos_index])

      } else {
        y[pos_index] <- ((x[pos_index] + 1)^object@lambda - 1) / object@lambda
      }
    }

    if (any(neg_index)) {
      if (object@lambda == 2.0) {
        y[neg_index] <- -log1p(-x[neg_index])

      } else {
        y[neg_index] <- -((-x[neg_index] + 1)^(2 - object@lambda) - 1) / (2 - object@lambda)
      }
    }

    return(y)
  }
)



# .revert_transform (Yeo-Johnson) ----------------------------------------------
setMethod(
  ".revert_transform",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {

    # Copy output
    y <- x

    # Determine positive and negative elements of the input vector
    pos_index <- x >= 0 & is.finite(x)
    neg_index <- x < 0 & is.finite(x)

    if (any(pos_index)) {
      if (object@lambda != 0) {
        y[pos_index] <- ((x[pos_index] * object@lambda + 1)^(1 / object@lambda) - 1)

      } else {
        y[pos_index] <- exp(x[pos_index]) - 1
      }
    }

    if (any(neg_index)) {
      if (object@lambda != 2) {
        y[neg_index] <- 1 - (x[neg_index] * (object@lambda - 2) + 1)^(1 / (2 - object@lambda))

      } else {
        y[neg_index] <- 1 - exp(-x[neg_index])
      }
    }

    # Apply scale and shift.
    y <- y * object@scale + object@shift

    return(y)
  })



# ..requires_shift_scale_optimisation (Yeo-Johnson (invariant)) ----------------
setMethod(
  "..requires_shift_scale_optimisation",
  signature(object = "transformationYeoJohnsonInvariant"),
  function(object, ...) {
    return(TRUE)
  }
)



# ..get_default_shift_range (Yeo-Johnson) --------------------------------------
setMethod(
  "..get_default_shift_range",
  signature(object = "transformationYeoJohnson"),
  function(object, x, ...) {

    # Compute skewness. This prevent local optimisers from becoming stuck in a
    # local optimum at the edge of the search grid.
    mu <- sum(x) / length(x)
    sigma_squared <- sum((x - mu)^2) / length(x)

    # Compute (weighted) skewness and kurtosis.
    skewness <- (1.0 / sigma_squared^(3 / 2)) * (sum((x - mu)^3)) / length(x)
    if (is.na(skewness)) skewness <- 0.0

    if (skewness < 0.0) {
      x_initial <- stats::quantile(x, 0.95)

    } else {
      x_initial <- stats::quantile(x, 0.05)
    }

    # Set shift range. Allow for moving the values entirely negative or
    # positive.
    shift_range <- c(min(x) - 1.0, x_initial, max(x) + 1.0)

    return(shift_range)
  }
)



# ..get_default_scale_range (Yeo-Johnson) --------------------------------------
setMethod(
  "..get_default_scale_range",
  signature(object = "transformationYeoJohnson"),
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



# ..get_default_lambda_range (Yeo-Johnson) -------------------------------------
setMethod(
  "..get_default_lambda_range",
  signature(object = "transformationYeoJohnson"),
  function(object, ...) {
    return(c(-100.0, 100.0))
  }
)


# ..standardise_data (Yeo-Johnson (invariant)) ---------------------------------
setMethod(
  "..standardise_data",
  signature(object = "transformationYeoJohnsonInvariant"),
  function(object, x, ...) {

    shift <- stats::median(x)
    if (!is.finite(shift)) shift <- 0.0

    scale <- stats::IQR(x)
    if (scale == 0.0) scale <- stats::sd(x)
    if (!is.finite(scale)) scale <- 1.0
    if (scale == 0.0) scale <- 1.0

    # if (all(x > 0.0)) {
    #   # Set shift so that the minimum value after scaling is 1, if all values
    #   # are positive.
    #   shift <- (min(x) / scale - 1.0) * scale
    # }

    # Shift and scale x.
    x <- (x - shift) / scale

    return(list(
      "x" = x,
      "shift" = shift,
      "scale" = scale
    ))
  }
)


# ..set_lambda (Yeo-Johnson) ---------------------------------------------------
setMethod(
  "..set_lambda",
  signature(object = "transformationYeoJohnson"),
  function(object, lambda, ...) {

    # Only update lambda if it is not a range.
    if (length(lambda) != 1) return(object)

    # Update lambda.
    object@lambda <- lambda

    return(object)
  }
)



# ..optimisation_parameters (Yeo-Johnson) --------------------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnson"),
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



# ..optimisation_parameters (Yeo-Johnson (invariant)) --------------------------
setMethod(
  "..optimisation_parameters",
  signature(object = "transformationYeoJohnsonInvariant"),
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
    x_range <- unname(x_range)

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



# ..log_likelihood (Yeo-Johnson) -----------------------------------------------
setMethod(
  "..log_likelihood",
  signature(object = "transformationYeoJohnson"),
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
           (object@lambda - 1.0) * sum(w * sign(x) * log1p(abs(x)))


    return(llf)
  }
)



# ..first_derivative (Yeo-Johnson) ---------------------------------------------
setMethod(
  "..first_derivative",
  signature(object = "transformationYeoJohnson"),
  function(object, x) {
    return((1 + abs(x))^(sign(x) * (object@lambda - 1)))
  }
)



# ..get_available_estimators (Yeo-Johnson) -------------------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationYeoJohnson"),
  function(object, ...) {

    available_estimators <- ..estimators_all()

    # Only allow Raymaekers and Rousseeuw's method for robust optimisation.
    if (!object@robust) {
      available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())
    }

    return(available_estimators)
  }
)



# ..get_available_estimators (Yeo-Johnson (invariant)) -------------------------
setMethod(
  "..get_available_estimators",
  signature(object = "transformationYeoJohnsonInvariant"),
  function(object, ...) {

    available_estimators <- ..estimators_all()

    # Raymaekers and Rousseeuw's method for robust optimisation is not suited
    # for simultaneous estimation of shift and lambda parameters.
    available_estimators <- setdiff(available_estimators, ..estimators_raymaekers_robust())

    return(available_estimators)
  }
)



# show (Yeo-Johnson) -----------------------------------------------------------
setMethod(
  "show",
  signature(object = "transformationYeoJohnson"),
  function(object) {

    str <- paste0(
      "A ", ifelse(object@robust, "robust ", ""),
      ifelse(object@shift != 0.0, "shift-invariant ", ""),
      "Yeo-Johnson transformation object"
    )

    if (object@complete) {
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")

      if (object@shift != 0.0) cat(paste0("  shift: ", object@shift, "\n"))

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)



# show (Yeo-Johnson (invariant)) -----------------------------------------------
setMethod(
  "show",
  signature(object = "transformationYeoJohnsonInvariant"),
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
      "Yeo-Johnson transformation object"
    )

    if (object@complete) {
      cat(paste0(str, ".\n"))
      cat("  lambda: ", object@lambda, "\n")
      cat("  shift: ", object@shift, "\n")
      cat("  scale: ", object@scale, "\n")

    } else {
      cat(paste0(str, " with unset transformation parameters.\n"))
    }
  }
)
