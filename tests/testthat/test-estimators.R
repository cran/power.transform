parameter_list <- list()
ii <- 1
for(method in c("box_cox", "yeo_johnson")){

  # Set estimation methods.
  estimation_methods <- setdiff(power.transform:::..estimators_all(), power.transform:::..estimators_raymaekers_robust())
  for(estimation_method in estimation_methods){

    parameter_list[[ii]] <- list(
      "method" =method,
      "robust" = FALSE,
      "invariant" = TRUE,
      "estimation_method" = estimation_method
    )

    ii <- ii + 1
  }
}



# Set seed.
set.seed(19L)

# Draw 1000 normally
x <- stats::rnorm(1000)

# Iterate over all parameter sets.
for(ii in seq_along(parameter_list)){

  # All-positive values -----------------------------------------------------
  x_positive <- exp(x)

  testthat::test_that(
    paste0(
      "Transforming all-positive values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_list[[ii]]$method,
      "; estimation_method: ", parameter_list[[ii]]$estimation_method, ")"
    ),
    {
      # Create the transformer.
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args = c(
          list("x" = x_positive),
          parameter_list[[ii]]
        )
      )

      if(parameter_list[[ii]]$method == "box_cox"){
        testthat::expect_true(abs(transformer@lambda) < 0.2)
        testthat::expect_true(abs(transformer@shift) < 0.1)

      } else if(parameter_list[[ii]]$method == "yeo_johnson"){
        testthat::expect_true(abs(transformer@lambda + 0.2) < 0.2)
      }
    }
  )

  # Some negative values -------------------------------------------------------
  x_part_negative <- exp(x) - 1

  testthat::test_that(
    paste0(
      "Transforming partially negative values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_list[[ii]]$method,
      "; estimation_method: ", parameter_list[[ii]]$estimation_method, ")"
    ),
    {
      # Create the transformer
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args = c(
          list("x" = x_part_negative),
          parameter_list[[ii]]
        )
      )

      # Check lambda values.
      if(parameter_list[[ii]]$method == "box_cox"){
        testthat::expect_true(abs(transformer@lambda) < 0.2)
        testthat::expect_true(abs(transformer@shift + 1.0) < 0.1)

      } else if(parameter_list[[ii]]$method == "yeo_johnson"){
        testthat::expect_true(abs(transformer@lambda + 0.2) < 0.2)
      }
    }
  )


  # All-negative values --------------------------------------------------------
  x_negative <- exp(x) - exp(max(x)) - 1E-8

  testthat::test_that(
    paste0(
      "Transforming completely negative values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_list[[ii]]$method,
      "; estimation_method: ", parameter_list[[ii]]$estimation_method, ")"
    ),
    {
      # Create the transformer
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args = c(
          list("x" = x_negative),
          parameter_list[[ii]]
        )
      )

      if(parameter_list[[ii]]$method == "box_cox"){
        testthat::expect_true(abs(transformer@lambda) < 0.2)

      } else if(parameter_list[[ii]]$method == "yeo_johnson"){
        testthat::expect_true(abs(transformer@lambda + 0.2) < 0.2)
      }
    }
  )


  # Some NA values -------------------------------------------------------------
  x_some_na <- exp(x)
  x_some_na[c(1,2)] <- NA_real_

  testthat::test_that(
    paste0(
      "Transforming all-positive values, with some NA values, generates the correct results. ",
      "(", ii,
      "; method: ", parameter_list[[ii]]$method,
      "; estimation_method: ", parameter_list[[ii]]$estimation_method, ")"
    ),
    {
      # Create the transformer.
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args = c(
          list("x" = x_some_na),
          parameter_list[[ii]]
        )
      )

      # Check lambda values.
      if(parameter_list[[ii]]$method == "box_cox"){
        testthat::expect_true(abs(transformer@lambda) < 0.2)
        testthat::expect_true(abs(transformer@shift) < 0.1)

      } else if(parameter_list[[ii]]$method == "yeo_johnson"){
        testthat::expect_true(abs(transformer@lambda + 0.2) < 0.2)
      }
    }
  )


  # Some infinite values -------------------------------------------------------
  x_some_inf <- exp(x)
  x_some_inf[c(1,2)] <- Inf

  testthat::test_that(
    paste0(
      "Transforming all-positive values, with some Inf values, generates the correct results. ",
      "(", ii,
      "; method: ", parameter_list[[ii]]$method,
      "; estimation_method: ", parameter_list[[ii]]$estimation_method, ")"
    ),
    {
      # Create the transformer.
      transformer <- do.call(
        power.transform::find_transformation_parameters,
        args = c(
          list("x" = x_some_inf),
          parameter_list[[ii]]
        )
      )

      # Check lambda values.
      if(parameter_list[[ii]]$method == "box_cox"){
        testthat::expect_true(abs(transformer@lambda) < 0.2)
        testthat::expect_true(abs(transformer@shift) < 0.1)

      } else if(parameter_list[[ii]]$method == "yeo_johnson"){
        testthat::expect_true(abs(transformer@lambda + 0.2) < 0.2)
      }
    }
  )

  # Few unique values (< 10) ---------------------------------------------------
  x_few_unique <- ceiling(x - min(x) + 1E-8)

  testthat::test_that(
    paste0(
      "Transforming vector with fewer than 10 unique values generates the correct results. ",
      "(", ii,
      "; method: ", parameter_list[[ii]]$method,
      "; estimation_method: ", parameter_list[[ii]]$estimation_method, ")"
    ),
    {
      # Creating the transformer should throw a warning, but otherwise function
      # normally.
      testthat::expect_warning(
        transformer <- do.call(
          power.transform::find_transformation_parameters,
          args = c(
            list("x" = x_few_unique),
            parameter_list[[ii]]
          )
        ),
        class = "power_transform_few_unique_values"
      )

      # Check lambda values.
      if(parameter_list[[ii]]$method == "box_cox"){
        testthat::expect_true(abs(transformer@lambda - 1.0) < 0.2)

      } else if(parameter_list[[ii]]$method == "yeo_johnson"){
        testthat::expect_true(abs(transformer@lambda -1.0) < 0.2)
      }
    }
  )
}
