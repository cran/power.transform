
#' Set transformation parameters
#'
#' `find_transformation_parameters` is used to find optimal parameters for
#' univariate transformation to normality.
#'
#' @param x A vector with numeric values.
#' @param method One of the following methods for power transformation:
#'
#'   * `box_cox`: Transformation using the Box-Cox transformation (Box and Cox,
#'   1964). The Box-Cox transformation requires that all data are strictly
#'   positive. Features that contain zero or negative values cannot be
#'   transformed using this transformation. In their work, Box and Cox define a
#'   shifted variant. We use this variant to shift values to a strictly positive
#'   range, when negative values are present. The Box-Cox transformation relies
#'   on a single parameter lambda, which is estimated through maximisation of
#'   the log-likelihood function corresponding to a normal distribution.
#'
#'   * `yeo_johnson`:Transformation using the Yeo-Johnson
#'   transformation (Yeo and Johnson, 2000). Unlike the Box-Cox transformation,
#'   the Yeo-Johnson transformation allows for negative and positive values.
#'   Like the Box-Cox transformation, this transformation relies on a single
#'   parameter lambda, which is estimated through maximisation of the
#'   log-likelihood function corresponding to a normal distribution.
#'
#'   * `none`: A fall-back method that will not transform values.
#'
#' @param robust Flag for using a robust version of Box-Cox or Yeo-Johnson
#'   transformation, as defined by Raymaekers and Rousseeuw (2021). This version
#'   is less sensitive in the presence outliers.
#' @param invariant Flag for using a version of Box-Cox or Yeo-Johnson
#'   transformation that simultaneously optimises location and scale in addition
#'   to the lambda parameter.
#' @param lambda Single lambda value, or range of lambda values that should be
#'   considered. Default: c(4.0, 6.0). Can be `NULL` to force optimisation
#'   without a constraint in lambda values.
#' @param empirical_gof_normality_p_value Significance value for the empirical
#'   goodness-of-fit test for central normality. The p-value is computed through
#'   the `assess_transformation` function. By setting this parameter to a
#'   numeric value other than `NULL`, the transformation will be rejected when
#'   the p-value of the test is below the significance value.
#' @param ... Unused parameters.
#'
#' @return A transformer object that can be used to transform values.
#' @export
#'
#' @references 1. Yeo, I. & Johnson, R. A. A new family of power transformations
#'   to improve normality or symmetry. Biometrika 87, 954–959 (2000).
#'
#'   1. Box, G. E. P. & Cox, D. R. An analysis of transformations. J. R. Stat.
#'   Soc. Series B Stat. Methodol. 26, 211–252 (1964).
#'
#'   1. Raymaekers, J., Rousseeuw,  P. J. Transforming variables to central
#'   normality. Mach Learn. (2021).
#'
#' @examples
#' x <- exp(stats::rnorm(1000))
#' transformer <- find_transformation_parameters(
#'   x = x,
#'   method = "box_cox")
find_transformation_parameters <- function(
    x,
    method = "yeo_johnson",
    robust = TRUE,
    invariant = TRUE,
    lambda = c(-4.0, 6.0),
    empirical_gof_normality_p_value = NULL,
    ...) {

  # Check transformation methods.
  if (!method %in% c("box_cox", "yeo_johnson", "none")) {
    stop(paste0(
      "The method argument should be one of \"box_cox\", \"yeo_johnson\" or \"none\". ",
      "Found: ", method))
  }

  if (method != "none") {
    # Perform checks on x.
    .check_data(x)

    # Remove NA or inf values.
    x <- x[is.finite(x)]

    # Check number of unique values.
    n_unique_values <- length(unique(x))
    if (n_unique_values <= 3 && method != "none") {
      rlang::warn(
        message = paste0(
          "x contains ", n_unique_values, " unique values, and ",
          "power transformation is not performed. power.transform requires ",
          "four or more unique values."),
        class = "power_transform_no_transform"
      )

      method <- "none"
    }

    if (n_unique_values > 3 && n_unique_values <= 10) {
      rlang::warn(
        message = paste0(
          "x contains ", n_unique_values, " unique values. Is the feature a categorical feature? ",
          "power.transform will treat the feature as numeric since it has at least four unique values."),
        class = "power_transform_few_unique_values"
      )
    }

    # Sort x.
    x <- sort(x)
  }

  # Create transformation objects.
  if (method == "none") {
    object <- methods::new("transformationNone")

  } else if (method == "box_cox") {
    object <- methods::new(
      "transformationBoxCox",
      robust = robust
    )

    if (invariant) {
      object <- methods::new(
        "transformationBoxCoxInvariant",
        object
      )
    }

    # Check lambda.
    .check_lambda_range(lambda)

  } else if (method == "yeo_johnson") {
    object <- methods::new(
      "transformationYeoJohnson",
      robust = robust
    )

    if (invariant) {
      object <- methods::new(
        "transformationYeoJohnsonInvariant",
        object
      )
    }

    # Check lambda.
    .check_lambda_range(lambda)

  } else {
    stop(
      paste0("Encountered an unknown transformation method: ", method)
    )
  }

  # Set transformation parameters.
  object <- .set_transformation_parameters(
    object = object,
    x = x,
    lambda = lambda,
    ...
  )

  # Reject transformations in case the goodness-of-fit test p-value is below the
  # required significance value.
  if (!is(object, "transformationNone") && !is.null(empirical_gof_normality_p_value)) {
    .check_gof_test_p_value(
      x = empirical_gof_normality_p_value,
      descriptor = "empirical goodness-of-fit test significance level (empirical_gof_normality_p_value)")

    gof_test_p <- assess_transformation(
      x = x,
      transformer = object,
      verbose = FALSE)

    if (is.na(gof_test_p)) {
      rlang::warn(
        message = paste0(
          "The p-value of the transformed data (", gof_test_p, ") could not be determined, likely due to an internal error. ",
          "The transformation is rejected, and data are kept as is."
        ),
        class = "power_transform_no_transform"
      )

      object <- methods::new("transformationNone")

      object <- .set_transformation_parameters(
        object = object,
        x = x,
        lambda = lambda,
        ...
      )

    } else if (gof_test_p < empirical_gof_normality_p_value) {
      rlang::warn(
        message = paste0(
          "The p-value of the transformed data (", gof_test_p, ") is below the required ",
          "significance level (", empirical_gof_normality_p_value, "). ",
          "The transformation is rejected, and data are kept as is."),
        class = "power_transform_no_transform")

      object <- methods::new("transformationNone")

      object <- .set_transformation_parameters(
        object = object,
        x = x,
        lambda = lambda,
        ...)
    }
  }

  return(object)
}



#' Transform values
#'
#' `power_transform` transforms numeric values to normality.
#'
#' @param x A vector with numeric values that should be transformed to
#'   normality.
#' @param transformer A transformer object created using
#'   `find_transformation_parameters`. If `NULL`, a transformer is generated
#'   internally.
#' @param oob_action Action that should be taken when out-of-bounds values are
#'   encountered in `x`. This can for example be 0 or negative values for
#'   Box-Cox transformations.
#'
#'   * `na` (default): replaces out-of-bounds values by NA values.
#'
#'   * `valid`: replaces out-of-bounds values by the closest valid boundary
#'   values.
#'
#'   This argument has no effect for Yeo-Johnson transformations.
#'
#' @param ... Parameters passed on to `find_transformation_parameters`.
#'
#' @inheritDotParams find_transformation_parameters
#'
#' @return A vector of transformed values of `x`.
#' @export
#'
#' @seealso [find_transformation_parameters]
#' @examples
#' x <- exp(stats::rnorm(1000))
#' y <- power_transform(
#'   x = x,
#'   method = "box_cox")
power_transform <- function(
    x,
    transformer = NULL,
    oob_action = "na",
    ...) {

  # Create a transformer.
  if (is.null(transformer)) {
    transformer <- do.call(
      find_transformation_parameters,
      c(list("x" = x), list(...)))
  }

  # Check the transformer.
  .check_transformer(transformer)

  # Check that x is numeric.
  if (!is(transformer, "transformationNone")) {
    if (!is.numeric(x)) {
      stop("x does not contain numeric values.")
    }
    if (!is.finite(get_lambda(transformer))) {
      stop("Transformer was not fitted correctly.")
    }
  }

  # Transform data using the transformer.
  y <- .transform(
    object = transformer,
    x = x,
    oob_action = oob_action)

  return(y)
}



#' Revert transformation
#'
#' `revert_power_transform` reverts the transformation of numeric values to
#' normality.
#'
#' @param y A vector with numeric values that was previously transformed to
#'   normality.
#' @param transformer A transformer object created using
#'   `find_transformation_parameters` that was used to transform the values to
#'   normality previously. Cannot be `NULL`.
#'
#' @return A vector of values.
#' @export
#'
#' @examples
#' x0 <- exp(stats::rnorm(1000))
#'
#' transformer <- find_transformation_parameters(
#'   x = x0,
#'   method = "box_cox")
#'
#' y <- power_transform(
#'   x = x0,
#'   transformer = transformer)
#'
#' x1 <- revert_power_transform(
#'   y = y,
#'   transformer = transformer)
revert_power_transform <- function(
    y,
    transformer) {

  if (missing(transformer)) {
    stop(paste0(
      "A transformer object is required to revert the transformation."
    ))
  }

  # Check the transformer.
  .check_transformer(transformer)

  # Revert transformation.
  x <- .revert_transform(
    object = transformer,
    x = y)

  return(x)
}
