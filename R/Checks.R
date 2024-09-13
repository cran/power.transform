.check_data <- function(x) {
  # Perform checks on x.
  if (length(x) == 0) {
    rlang::abort("x does not contain any values.")
  }

  if (is.factor(x)) {
    rlang::abort("x is categorical, and power transformations are not applicable.")
  }

  if (!is.numeric(x)) {
    rlang::abort("x does not contain numeric values.")
  }

  if (all(!is.finite(x))) {
    rlang::abort("x only contains NA or inf values.")
  }

  return(invisible(TRUE))
}



.check_transformer <- function(x) {

  if (!is(x, "transformationPowerTransform")) {
    rlang::abort(paste0(
      "The transformer object does not have the expected class. ",
      "Expected: transformationPowerTransform (or subclass). ",
      "Found: ", class(x)[1]))
  }

  if (!x@complete) {
    rlang::abort(paste0(
      "The transformer object did not have all fitting parameters set."))
  }

  return(invisible(TRUE))
}



.check_lambda_range <- function(x) {
  # This checks the lambda argument for find_parameters.

  # NULL is a valid value.
  if (is.null(x)) return(invisible(TRUE))

  # Otherwise, must be length 2, numeric, finite and sorted.
  if (!length(x) %in% c(1L, 2L)) {
    rlang::abort(paste0("lambda should consist of 1 or 2 numeric values. ", length(x), " values were found."))
  }

  if (!is.numeric(x)) {
    rlang::abort("lambda should consist of 1 or 2 numeric values. The values are not numeric.")
  }

  if (any(!is.finite(x))) {
    rlang::abort("lambda should consist of 1 or 2 numeric values. One or both values are not finite.")
  }

  if (length(x) == 2) {
    if (diff(x) == 0.0) {
      rlang::abort("If lambda consists of 2 numeric values, these can not be identical.")
    }

    if (is.unsorted(x)) {
      rlang::abort("If lambda consists of 2 numeric values, the values should be ordered by increasing value.")
    }
  }

  return(invisible(TRUE))
}



.check_lambda_value <- function(x) {
  # This checks the lambda argument for the mutator.
  if (!length(x) == 1) {
    rlang::abort(paste0("lambda should be a single, finite, numeric value. ", length(x), " values were found."))
  }

  if (!is.numeric(x)) {
    rlang::abort("lambda should be a single, finite, numeric value. Found: ", paste_s(class(x)))
  }

  if (!is.finite(x)) {
    rlang::abort(paste0(
      "lambda should be a single, finite, numeric value. ",
      "Found: a numeric value that is not finite (", x, ")"))
  }

  return(invisible(TRUE))
}



.check_shift_value <- function(x) {
  # This checks the shift argument for the mutator.
  if (!length(x) == 1) {
    rlang::abort(paste0("shift should be a single, finite, numeric value. ", length(x), " values were found."))
  }

  if (!is.numeric(x)) {
    rlang::abort("shift should be a single, finite, numeric value. Found: ", paste_s(class(x)))
  }

  if (!is.finite(x)) {
    rlang::abort(paste0(
      "shift should be a single, finite, numeric value. ",
      "Found: a numeric value that is not finite (", x, ")"))
  }

  return(invisible(TRUE))
}



.check_scale_value <- function(x) {
  # This checks the scale argument for the mutator.
  if (!length(x) == 1) {
    rlang::abort(paste0("scale should be a single, finite, numeric value. ", length(x), " values were found."))
  }

  if (!is.numeric(x)) {
    rlang::abort("scale should be a single, finite, numeric value. Found: ", paste_s(class(x)))
  }

  if (!is.finite(x)) {
    rlang::abort(paste0(
      "scale should be a single, finite, numeric value. ",
      "Found: a numeric value that is not finite (", x, ")"))
  }

  return(invisible(TRUE))
}



.check_oob_action <- function(x) {
  if (length(x) != 1) {
    rlang::abort(paste0(
      "One of the following should be provided as the oob_action argument: \"na\", or \"valid\".",
      "Found: ", length(x), " arguments."))
  }

  if (!any(x %in% c("na", "valid"))) {
    rlang::abort(paste0("One of the following should be provided as the oob_action argument: \"na\", or \"valid\".",
    "Found: ", x))
  }

  return(invisible(TRUE))
}



.check_gof_test_p_value <- function(x, descriptor) {
  # NULL is a valid value.
  if (is.null(x)) return(invisible(TRUE))

  if (length(x) != 1) {
    rlang::abort(paste0(
      "The ", descriptor, " should consist of a single, numeric value. ",
      length(x), " values were provided."))
  }

  if (!is.numeric(x)) {
    rlang::abort(paste0(
      "The ", descriptor, " should consist of a single, numeric value. ",
      "The provided value is not numeric: ", paste_s(class(x))))
  }

  if (x > 1.0 || x < 0.0) {
    rlang::abort(paste0(
      "The ", descriptor, " should be a value between 0.0 and 1.0. Found: ", x))
  }

  return(invisible(TRUE))
}



.check_weighting_function_parameters <- function(x, default_parameters) {

  # Skip if there is nothing to check.
  if (length(x) == 0) return(invisible(TRUE))

  # Throw an error if the parameters are not named.
  if (is.null(names(x))) {
    rlang::abort(paste0(
      "Names are currently missing from weighting_function_parameters. ",
      "Parameters should be named, e.g. weighting_function_parameters = list(\"k1\" = 0.5). "))
  }

  unknown_parameters <- setdiff(names(x), names(default_parameters))
  # Throw a rlang::warn if all provided parameters could not be matched.
  if (length(unknown_parameters) > 0) {
    if (length(default_parameters) > 0) {
      rlang::warn(paste0(
        "One or more weighting function parameters could not be matched: ",
        paste_s(unknown_parameters), " ",
        "The following parameters can be set: ",
        paste_s(names(default_parameters))
      ))

    } else {
      rlang::warn(paste0(
        "One or more weighting function parameters could not be matched: ",
        paste_s(unknown_parameters), " ",
        "The weighting function does not require any parameters."
      ))
    }
  }

  return(invisible(TRUE))
}
