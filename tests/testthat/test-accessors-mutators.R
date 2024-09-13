
for (method in c("box_cox", "yeo_johnson", "none")) {
  testthat::test_that(
    "Accessing and changing lambda parameters functions correctly.", {

      # Generate data.
      x <- stats::rnorm(n = 1000L)

      # Create the transformer.
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        estimation_method = "mle"
      )

      if (method == "none") {
        # Check that lambda parameter cannot be set.
        testthat::expect_warning(
          transformer_new <- power.transform::set_lambda(
            object = transformer,
            lambda = 2.0),
          class = "power_transform_no_attribute"
        )

        # Expect that the transformer is not changed.
        testthat::expect_equal(transformer, transformer_new)

        # Check that lambda parameter cannot be read.
        testthat::expect_warning(
          lambda_value <- power.transform::get_lambda(
            object = transformer_new),
          class = "power_transform_no_attribute"
        )

        testthat::expect_equal(lambda_value, NA_real_)

      } else {
        transformer_new <- power.transform::set_lambda(
          object = transformer,
          lambda = 2.0
        )

        lambda_value <- power.transform::get_lambda(
          object = transformer_new
        )

        testthat::expect_equal(lambda_value, 2.0)
      }
    }
  )

  testthat::test_that(
    "Accessing and changing shift parameters functions correctly.", {

      # Generate data.
      x <- stats::rnorm(n = 1000L)

      # Create the transformer. Note that the data are not normally distributed.
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        estimation_method = "mle"
      )

      if (method == "none") {
        # Check that shift parameter cannot be set.
        testthat::expect_warning(
          transformer_new <- power.transform::set_shift(
            object = transformer,
            shift = 3.0),
          class = "power_transform_no_attribute"
        )

        # Expect that the transformer is not changed.
        testthat::expect_equal(transformer, transformer_new)

        # Check that shift parameter cannot be read.
        testthat::expect_warning(
          shift_value <- power.transform::get_shift(
            object = transformer_new),
          class = "power_transform_no_attribute"
        )

        testthat::expect_equal(shift_value, NA_real_)

      } else {
        transformer_new <- power.transform::set_shift(
          object = transformer,
          shift = 3.0
        )

        shift_value <- power.transform::get_shift(
          object = transformer_new
        )

        testthat::expect_equal(shift_value, 3.0)
      }
    }
  )

  testthat::test_that(
    "Accessing and changing scale parameters functions correctly.", {

      # Generate data.
      x <- stats::rnorm(n = 1000L)

      # Create the transformer. Note that the data are not normally distributed.
      transformer <- power.transform::find_transformation_parameters(
        x = x,
        method = method,
        estimation_method = "mle"
      )

      if (method == "none") {
        # Check that shift parameter cannot be set.
        testthat::expect_warning(
          transformer_new <- power.transform::set_scale(
            object = transformer,
            scale = 3.0),
          class = "power_transform_no_attribute"
        )

        # Expect that the transformer is not changed.
        testthat::expect_equal(transformer, transformer_new)

        # Check that shift parameter cannot be read.
        testthat::expect_warning(
          scale_value <- power.transform::get_scale(
            object = transformer_new),
          class = "power_transform_no_attribute"
        )

        testthat::expect_equal(scale_value, NA_real_)

      } else {
        transformer_new <- power.transform::set_scale(
          object = transformer,
          scale = 3.0
        )

        scale_value <- power.transform::get_scale(
          object = transformer_new
        )

        testthat::expect_equal(scale_value, 3.0)
      }
    }
  )
}

