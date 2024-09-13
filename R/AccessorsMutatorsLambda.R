#' @include TransformationObjects.R
#' @include TransformationBoxCox.R
#' @include TransformationYeoJohnson.R

# set_lambda (generic) ---------------------------------------------------------

#' @title Set lambda value
#'
#' @description Set the lambda value of a transformer object.
#'
#' @param object Transformer object
#' @param lambda Lambda value
#' @param ... Unused arguments
#'
#' @return Transformer object with updated lambda value.
#' @export
#'
#' @rdname lambda-mutator-method
setGeneric(
  "set_lambda",
  function(object, lambda, ...) standardGeneric("set_lambda"))


#' @rdname lambda-mutator-method
setMethod(
  "set_lambda",
  signature(object = "transformationPowerTransform"),
  function(object, lambda, ...) {
    rlang::warn(
      message = "The current transformer object does not have a lambda value that can be set.",
      class = "power_transform_no_attribute"
    )

    return(object)
  }
)


#' @rdname lambda-mutator-method
setMethod(
  "set_lambda",
  signature(object = "transformationBoxCox"),
  function(object, lambda, ...) {

    .check_lambda_value(x = lambda)

    object <- ..set_lambda(
      object = object,
      lambda = lambda)

    return(object)
  }
)


#' @rdname lambda-mutator-method
setMethod(
  "set_lambda",
  signature(object = "transformationYeoJohnson"),
  function(object, lambda, ...) {

    .check_lambda_value(x = lambda)

    object <- ..set_lambda(
      object = object,
      lambda = lambda)

    return(object)
  }
)



# get_lambda (generic) ---------------------------------------------------------

#' @title Get lambda value
#'
#' @description Get the lambda value of a transformer object.
#'
#' @param object Transformer object
#' @param ... Unused arguments
#'
#' @return Lambda value of the transformer.
#' @export
#'
#' @rdname lambda-accessor-method
setGeneric(
  "get_lambda",
  function(object, ...) standardGeneric("get_lambda"))


#' @rdname lambda-accessor-method
setMethod(
  "get_lambda",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
    rlang::warn(
      message = "The current transformer object does not have a lambda value that can be read.",
      class = "power_transform_no_attribute"
    )

    return(NA_real_)
  }
)


#' @rdname lambda-accessor-method
setMethod(
  "get_lambda",
  signature(object = "transformationBoxCox"),
  function(object, ...) {

    return(object@lambda)
  }
)


#' @rdname lambda-accessor-method
setMethod(
  "get_lambda",
  signature(object = "transformationYeoJohnson"),
  function(object, ...) {

    return(object@lambda)
  }
)
