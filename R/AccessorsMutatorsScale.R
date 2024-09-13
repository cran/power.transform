#' @include TransformationObjects.R
#' @include TransformationBoxCox.R
#' @include TransformationYeoJohnson.R

# set_scale (generic) ----------------------------------------------------------

#' @title Set scale value
#'
#' @description Set the scale value of a transformer object.
#'
#' @param object Transformer object
#' @param scale scale value
#' @param ... Unused arguments
#'
#' @return Transformer object with updated scale value.
#' @export
#'
#' @rdname scale-mutator-method
setGeneric(
  "set_scale",
  function(object, scale, ...) standardGeneric("set_scale"))


#' @rdname scale-mutator-method
setMethod(
  "set_scale",
  signature(object = "transformationPowerTransform"),
  function(object, scale, ...) {
    rlang::warn(
      message = "The current transformer object does not have a scale value that can be set.",
      class = "power_transform_no_attribute"
    )

    return(object)
  }
)


#' @rdname scale-mutator-method
setMethod(
  "set_scale",
  signature(object = "transformationBoxCox"),
  function(object, scale, ...) {

    .check_scale_value(x = scale)

    object@scale <- scale

    return(object)
  }
)


#' @rdname scale-mutator-method
setMethod(
  "set_scale",
  signature(object = "transformationYeoJohnson"),
  function(object, scale, ...) {

    .check_scale_value(x = scale)

    object@scale <- scale

    return(object)
  }
)



# get_scale (generic) ---------------------------------------------------------

#' @title Get scale value
#'
#' @description Get the scale value of a transformer object.
#'
#' @param object Transformer object
#' @param ... Unused arguments
#'
#' @return scale value of the transformer.
#' @export
#'
#' @rdname scale-accessor-method
setGeneric(
  "get_scale",
  function(object, ...) standardGeneric("get_scale"))


#' @rdname scale-accessor-method
setMethod(
  "get_scale",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
    rlang::warn(
      message = "The current transformer object does not have a scale value that can be read.",
      class = "power_transform_no_attribute"
    )

    return(NA_real_)
  }
)


#' @rdname scale-accessor-method
setMethod(
  "get_scale",
  signature(object = "transformationBoxCox"),
  function(object, ...) {

    return(object@scale)
  }
)


#' @rdname scale-accessor-method
setMethod(
  "get_scale",
  signature(object = "transformationYeoJohnson"),
  function(object, ...) {

    return(object@scale)
  }
)
