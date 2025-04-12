..estimator_wrapper <- function(
    par,
    transformer,
    estimator,
    x,
    parameter_type,
    verbose = FALSE,
    ...) {

  if ("lambda" %in% parameter_type) {
    transformer@lambda <- par[which(parameter_type == "lambda")]
  }
  if ("shift" %in% parameter_type) {
    transformer@shift <- par[which(parameter_type == "shift")]
  }
  if ("scale" %in% parameter_type) {
    transformer@scale <- par[which(parameter_type == "scale")]
  }

  # Compute target value.
  target <- .compute_objective(
    object = estimator,
    transformer = transformer,
    x = x,
    ...
  )

  if (verbose) {
    cat(paste0(
      "target: ", target, "; ",
      paste0(
        mapply(
          function(type, value) (paste0(type, ": ", value)),
          type = parameter_type,
          value = par),
        collapse = "; "),
      "\n"))
  }

  return(target)
}



# estimatorGeneric definition --------------------------------------------------
setClass(
  "estimatorGeneric",
  slots = list(
    "method" = "character",
    "weighting_method" = "ANY"),
  prototype = list(
    "method" = "none",
    "weighting_method" = NULL))



# .compute_objective (generic) -------------------------------------------------
setGeneric(
  ".compute_objective",
  function(object, ...) standardGeneric(".compute_objective"))



# .compute_objective (general) -------------------------------------------------
setMethod(
  ".compute_objective",
  signature(object = "estimatorGeneric"),
  function(object, ...) {

    # Estimators should be specified for each estimator - an error is thrown to
    # ensure that this general method is never used.
    rlang::abort(paste0(
      "DEV: missing .compute_objective method for the ",
      paste_s(class(object)), " class."
    ))
  }
)



# .optimise_transformation_parameters (generic) --------------------------------
setGeneric(
  ".optimise_transformation_parameters",
  function(object, ...) standardGeneric(".optimise_transformation_parameters"))



# ..optimise_transformation_parameters (general) -------------------------------
setMethod(
  ".optimise_transformation_parameters",
  signature(object = "estimatorGeneric"),
  function(
    object,
    transformer,
    optimiser = NULL,
    optimisation_parameters,
    x,
    optimiser_control = NULL,
    verbose = FALSE,
    ...) {

    # Check that we do not inadvertently pass problems that do not require
    # optimisation to the optimiser.
    if (is.null(optimisation_parameters)) {
      rlang::abort("DEV: optimisation_parameters cannot be empty.")
    }
    if (!is(transformer, "transformationPowerTransform")) {
      rlang::abort("DEV: transformer should be a valid power transformation object.")
    }
    if (!is.numeric(x)) {
      rlang::abort("DEV: x should be numeric.")
    }

    # Check fall-back option.
    if (
      !is_package_installed("nloptr") &&
      optimiser %in% c("direct", "direct-l", "mlsl", "bobyqa", "subplex", "nelder-mead")
    ) {
      rlang::warn(
        message = paste0(
          "The nloptr package is required to optimise power transformation parameters using the ",
          optimiser, " algoritm. stats::optim is used as a fallback option."),
        class = "power_transform_package_missing",
        .frequency = "once",
        .frequency_id = "power_transform_package_missing_nloptr"
      )

      optimiser <- "optim-nelder-mead"
    }

    # Set default optimiser values.
    if (is.null(optimiser_control)) optimiser_control <- ..get_default_optimiser_control(
      object = object,
      optimiser = optimiser)

    if (optimiser == "direct-l") {
      # DIRECT-L algorithm
      #
      # D. R. Jones, C. D. Perttunen, and B. E. Stuckmann, “Lipschitzian
      # optimization without the lipschitz constant,” J. Optimization Theory and
      # Applications, vol. 79, p. 157 (1993).
      #
      # J. M. Gablonsky and C. T. Kelley, “A locally-biased form of the DIRECT
      # algorithm," J. Global Optimization, vol. 21 (1), p. 27-37 (2001).
      results <- tryCatch(
        nloptr::directL(
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type,
          nl.info = verbose,
          verbose = verbose
        ),
        error = identity
      )

    } else if (optimiser == "direct") {
      # DIRECT algorithm
      #
      # D. R. Jones, C. D. Perttunen, and B. E. Stuckmann, “Lipschitzian
      # optimization without the lipschitz constant,” J. Optimization Theory and
      # Applications, vol. 79, p. 157 (1993).
      results <- tryCatch(
        nloptr::direct(
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type,
          nl.info = verbose,
          verbose = verbose
        ),
        error = identity
      )

    } else if (optimiser == "mlsl") {
      # Multi-Level Single-Linkage algorithm (MLSL)
      #
      # A. H. G. Rinnooy Kan and G. T. Timmer, "Stochastic global optimization
      # methods," Mathematical Programming, vol. 39, p. 27-78 (1987).
      #
      # Sergei Kucherenko and Yury Sytsko, "Application of deterministic
      # low-discrepancy sequences in global optimization," Computational
      # Optimization and Applications, vol. 30, p. 297-318 (2005).
      results <- tryCatch(
        nloptr::mlsl(
          x0 = optimisation_parameters$initial,
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type,
          nl.info = verbose,
          verbose = verbose
        ),
        error = identity
      )

    } else if (optimiser == "bobyqa") {
      # BOBYQA algorithm
      #
      # M. J. D. Powell, "The BOBYQA algorithm for bound constrained
      # optimization without derivatives," Department of Applied Mathematics and
      # Theoretical Physics, Cambridge England, technical report NA2009/06
      # (2009)
      results <- tryCatch(
        nloptr::bobyqa(
          x0 = optimisation_parameters$initial,
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type,
          nl.info = verbose,
          verbose = verbose),
        error = identity
      )

    } else if (optimiser == "subplex") {
      # SUBPLEX algorithm
      #
      # T. Rowan, “Functional Stability Analysis of Numerical
      # Algorithms”, Ph.D. thesis, Department of Computer Sciences, University
      # of Texas at Austin, 1990.
      results <- tryCatch(
        nloptr::sbplx(
          x0 = optimisation_parameters$initial,
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type,
          nl.info = verbose,
          verbose = verbose),
        error = identity
      )

    } else if (optimiser == "nelder-mead") {
      # Nelder-Mead simplex algorithm
      #
      # J. A. Nelder and R. Mead, “A simplex method for function minimization,”
      # The Computer Journal 7, p. 308-313 (1965).
      #
      # M. J. Box, “A new method of constrained optimization and a comparison
      # with other methods,” Computer J. 8 (1), 42-52 (1965).
      results <- tryCatch(
        nloptr::neldermead(
          x0 = optimisation_parameters$initial,
          fn = ..estimator_wrapper,
          lower = optimisation_parameters$lower,
          upper = optimisation_parameters$upper,
          control = optimiser_control,
          transformer = transformer,
          estimator = object,
          x = x,
          parameter_type = optimisation_parameters$parameter_type,
          nl.info = verbose,
          verbose = verbose
        ),
        error = identity
      )

    } else if (optimiser == "optim-nelder-mead") {
      # Fall-back optimiser in case nloptr is not available. The Nelder-Mead
      # algorithm in stats::optim does not yield results as consistent as the
      # nloptr optimisers.
      results <- stats::optim(
        par = optimisation_parameters$initial,
        fn = ..estimator_wrapper,
        gr = NULL,
        transformer = transformer,
        estimator = object,
        x = x,
        parameter_type = optimisation_parameters$parameter_type,
        verbose = verbose,
        control = list("abstol" = 1E-5, "reltol" = 1E-5)
      )

    } else {
      stop(paste0("Optimiser not recognised: ", optimiser))
    }

    # Check for known errors. If an unknown error is encountered, raise this
    # error to the user.
    if (inherits(results, "error")) {
      if (!results$message %in% c("objective in x0 returns NA")) {
        rlang::abort(results$message)
      }
    }

    # Initialise a parameter list to return to the calling process.
    parameter_list <- lapply(
      optimisation_parameters$parameter_type,
      function(x) (NA_real_)
    )
    names(parameter_list) <- optimisation_parameters$parameter_type

    if (!inherits(results, "error")) {
      # Insert optimal parameter values into the parameter list, if any.
      for (ii in seq_along(optimisation_parameters$parameter_type)) {
        parameter_value <- results$par[ii]
        if (is.finite(parameter_value)) {
          parameter_list[[optimisation_parameters$parameter_type[ii]]] <- parameter_value
        }
      }

      # Add target value
      parameter_list$value <- results$value
    }

    return(parameter_list)
  }
)



# ..get_available_weighting_functions (generic) --------------------------------
setGeneric(
  "..get_available_weighting_functions",
  function(transformer, estimator, ...) standardGeneric("..get_available_weighting_functions"))



# ..get_available_weighting_functions (general) --------------------------------
setMethod(
  "..get_available_weighting_functions",
  signature(transformer = "transformationPowerTransform", estimator = "estimatorGeneric"),
  function(transformer, estimator, ...) {

    # If transformation methods is not robust, do not use weights.
    if (!transformer@robust) return("none")

    return(..weighting_functions_all())
  }
)



# ..get_default_weighting_function (generic) -----------------------------------
setGeneric(
  "..get_default_weighting_function",
  function(transformer, estimator, ...) standardGeneric("..get_default_weighting_function"))



# ..get_default_weighting_function (general) -----------------------------------
setMethod(
  "..get_default_weighting_function",
  signature(transformer = "transformationPowerTransform", estimator = "estimatorGeneric"),
  function(transformer, estimator, ...) {

    if (transformer@robust) {
      if (transformer@method == "box_cox") return("empirical_probability_step")

      return("empirical_probability_cosine")

    } else {
      return("none")
    }
  }
)



# .get_weights (generic) -------------------------------------------------------
setGeneric(
  ".get_weights",
  function(object, ...) standardGeneric(".get_weights"))



# .get_weights (general) -------------------------------------------------------
setMethod(
  ".get_weights",
  signature(object = "estimatorGeneric"),
  function(object, transformer, x, ...) {
    return(.get_weights(
      object = object@weighting_method,
      transformer = transformer,
      x = x))
  }
)



# ..get_default_optimiser (generic) --------------------------------------------
setGeneric(
  "..get_default_optimiser",
  function(object, ...) standardGeneric("..get_default_optimiser"))



# ..get_default_optimiser (general) --------------------------------------------
setMethod(
  "..get_default_optimiser",
  signature(object = "estimatorGeneric"),
  function(object, transformer, ...) {

    if (..requires_shift_scale_optimisation(transformer)) {
      return("bobyqa")

    } else {
      return("bobyqa")
    }
  }
)



# ..get_default_optimiser_control (generic) ------------------------------------
setGeneric(
  "..get_default_optimiser_control",
  function(object, ...) standardGeneric("..get_default_optimiser_control"))



.set_estimator <- function(
  transformer,
  estimation_method,
  weighting_function,
  weighting_function_parameters) {

  # Check estimation method.
  available_estimators <- ..get_available_estimators(object = transformer)
  if (!estimation_method %in% available_estimators) {
    stop(paste0(
      "The desired estimator ", estimation_method, " is not available for ",
      ifelse(transformer@robust, "robust ", ""),
      "optimisation for ", paste_s(class(transformer)), " objects."))
  }

  # Create estimator instance.
  if (estimation_method %in% ..estimators_mle()) {
    estimator <- methods::new("estimatorMaximumLikelihoodEstimation")

  } else if (estimation_method %in% ..estimators_raymaekers_robust()) {
    estimator <- methods::new("estimatorRaymaekersRobust")

  } else if (estimation_method %in% ..estimators_anderson_darling()) {
    estimator <- methods::new("estimatorAndersonDarling")

  } else if (estimation_method %in% ..estimators_cramer_von_mises()) {
    estimator <- methods::new("estimatorCramervonMises")

  } else if (estimation_method %in% ..estimators_kolmogorov_smirnov()) {
    estimator <- methods::new("estimatorKolmogorovSmirnov")

  } else if (estimation_method %in% ..estimators_jarque_bera()) {
    estimator <- methods::new("estimatorJarqueBera")

  } else if (estimation_method %in% ..estimators_dagostino()) {
    estimator <- methods::new("estimatorDAgostino")

  } else {
    stop(paste0("DEV: unknown estimation_method: ", estimation_method))
  }

  estimator@weighting_method <- .set_weighting_function(
    transformer = transformer,
    estimator = estimator,
    weighting_function = weighting_function,
    weighting_function_parameters = weighting_function_parameters)

  return(estimator)
}



..estimators_all <- function() {
  # Create a list with all estimators implemented in power.transform. Update
  # when more estimators are added.

  return(c(
    ..estimators_mle(),
    ..estimators_raymaekers_robust(),
    ..estimators_anderson_darling(),
    ..estimators_cramer_von_mises(),
    # ..estimators_kolmogorov_smirnov(), <-- The optimiser has trouble with this one for larger sample sizes.
    ..estimators_jarque_bera(),
    ..estimators_dagostino()
  ))
}
