parameter_list <- list()
ii <- 1
for (method in c("box_cox", "yeo_johnson", "none")) {
  for (robust in c(FALSE, TRUE)) {
    for (invariant in c(FALSE, TRUE)) {
      parameter_list[[ii]] <- list(
        "method" = method,
        "robust" = robust,
        "invariant" = invariant,
        "estimation_method" = "mle"
      )

      ii <- ii + 1
    }
  }
}

# Iterate over different sample numbers.
for (n in c(100, 1000, 10000)) {
  # Set seed.
  set.seed(19L)

  # Draw samples according to a normal distribution.
  x_normal <- stats::rnorm(n)

  # Make exponentially distributed.
  x_exponential <- exp(x_normal)

  # Draw bimodal with large separation
  x_bimodal_large <- c(
    stats::rnorm(n / 2, mean = 0.0),
    stats::rnorm(n / 2, mean = 6.0)
  )

  # Draw bimodal with intermediate separation
  x_bimodal_intermediate <- c(
    stats::rnorm(n / 2, mean = 0.0),
    stats::rnorm(n / 2, mean = 3.0)
  )

  # Draw bimodal with small separation
  x_bimodal_small <- c(
    stats::rnorm(n / 2, mean = 0.0),
    stats::rnorm(n / 2, mean = 1.0)
  )

  # Draw uniform.
  x_uniform <- stats::runif(n, min = 1E-5)

  # Draw Gamma (left skewed)
  x_gamma_left <- stats::rgamma(n, shape = 2.0, scale = 2.0)

  # Draw Gamma (monotonic decreasing)
  x_gamma_monotonic <- stats::rgamma(n, shape = 1.0, scale = 2.0)

  # Iterate over all parameter sets.
  for (ii in seq_along(parameter_list)) {
    testthat::test_that(
      paste0(
        "Assessing transformation goodness generates the correct results. ",
        "(", ii,
        "; n: ", n,
        "; method: ", parameter_list[[ii]]$method,
        "; robust: ", parameter_list[[ii]]$robust,
        "; invariant: ", parameter_list[[ii]]$invariant, ")"
      ),
      {
        # Normally distributed data --------------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_normal),
            parameter_list[[ii]]
          )
        ))

        # Compute the p-value.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_normal,
          transformer = transformer
        ))

        if (n == 100) {
          testthat::expect_gt(p_value, 0.10)

        } else {
          testthat::expect_gt(p_value, 0.50)
        }

        # Exponentially distributed data ---------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_exponential),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_exponential,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none") {
          testthat::expect_lt(p_value, 0.01)
        } else {
          testthat::expect_gt(p_value, 0.05)
        }


        # Left-skewed gamma distribution ---------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_gamma_left),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_gamma_left,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none") {
          testthat::expect_lt(p_value, 0.20)
        } else {
          testthat::expect_gt(p_value, 0.05)
        }

        # Monotonically decreasing gamma distribution --------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_gamma_monotonic),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_gamma_monotonic,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none") {
          testthat::expect_lt(p_value, 0.01)
        } else {
          testthat::expect_gt(p_value, 0.05)
        }

        # Bi-modal distribution (large separation) -----------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_bimodal_large),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_bimodal_large,
          transformer = transformer
        ))

        testthat::expect_lt(p_value, 0.01)

        # Bi-modal distribution (intermediate separation) ----------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_bimodal_intermediate),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_bimodal_intermediate,
          transformer = transformer
        ))

        testthat::expect_lt(p_value, 0.50)

        # Bi-modal distribution (small separation) -----------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_bimodal_small),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_bimodal_small,
          transformer = transformer
        ))

        testthat::expect_gt(p_value, 0.10)

        # Uniform distribution -------------------------------------------------

        # Create the transformer.
        transformer <- suppressWarnings(do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_uniform),
            parameter_list[[ii]]
          )
        ))

        # Compute the fraction of instances which exceed the threshold.
        p_value <- suppressMessages(power.transform::assess_transformation(
          x = x_uniform,
          transformer = transformer
        ))

        if (parameter_list[[ii]]$method == "none" && n == 100) {
          testthat::expect_lt(p_value, 0.05)
        } else {
          testthat::expect_gt(p_value, 0.01)
        }
      }
    )
  }
}


# Assess empirical goodness of fit test for rejecting non-normal distributions.
testthat::test_that(
  "Non-normal transformations are correctly rejected.", {

    # Generate data.
    x_bimodal_large <- c(
      stats::rnorm(500L, mean = 0.0),
      stats::rnorm(500L, mean = 6.0)
    )

    # Create the transformer. Note that the data are not normally distributed.
    transformer_default <- power.transform::find_transformation_parameters(
      x = x_bimodal_large,
      method = "yeo_johnson",
      estimation_method = "mle"
    )

    testthat::expect_s4_class(transformer_default, "transformationYeoJohnsonInvariant")

    testthat::expect_warning(
      transformer_rejected <- power.transform::find_transformation_parameters(
        x = x_bimodal_large,
        method = "yeo_johnson",
        estimation_method = "mle",
        empirical_gof_normality_p_value = 0.05),
      class = "power_transform_no_transform"
    )

    testthat::expect_s4_class(transformer_rejected, "transformationNone")
  }
)


