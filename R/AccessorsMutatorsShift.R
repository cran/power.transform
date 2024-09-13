#' @include TransformationObjects.R
#' @include TransformationBoxCox.R
#' @include TransformationYeoJohnson.R

# set_shift (generic) ----------------------------------------------------------

#' @title Set shift value
#'
#' @description Set the shift value of a transformer object.
#'
#' @param object Transformer object
#' @param shift Shift value
#' @param ... Unused arguments
#'
#' @return Transformer object with updated shift value.
#' @export
#'
#' @rdname shift-mutator-method
setGeneric(
  "set_shift",
  function(object, shift, ...) standardGeneric("set_shift"))


#' @rdname shift-mutator-method
setMethod(
  "set_shift",
  signature(object = "transformationPowerTransform"),
  function(object, shift, ...) {
    rlang::warn(
      message = "The current transformer object does not have a shift value that can be set.",
      class = "power_transform_no_attribute"
    )

    return(object)
  }
)


#' @rdname shift-mutator-method
setMethod(
  "set_shift",
  signature(object = "transformationBoxCox"),
  function(object, shift, ...) {

    .check_shift_value(x = shift)

    object@shift <- shift

    return(object)
  }
)


#' @rdname shift-mutator-method
setMethod(
  "set_shift",
  signature(object = "transformationYeoJohnson"),
  function(object, shift, ...) {

    .check_shift_value(x = shift)

    object@shift <- shift

    return(object)
  }
)



# get_shift (generic) ---------------------------------------------------------

#' @title Get shift value
#'
#' @description Get the shift value of a transformer object.
#'
#' @param object Transformer object
#' @param ... Unused arguments
#'
#' @return shift value of the transformer.
#' @export
#'
#' @rdname shift-accessor-method
setGeneric(
  "get_shift",
  function(object, ...) standardGeneric("get_shift"))


#' @rdname shift-accessor-method
setMethod(
  "get_shift",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
    rlang::warn(
      message = "The current transformer object does not have a shift value that can be read.",
      class = "power_transform_no_attribute"
    )

    return(NA_real_)
  }
)


#' @rdname shift-accessor-method
setMethod(
  "get_shift",
  signature(object = "transformationBoxCox"),
  function(object, ...) {

    return(object@shift)
  }
)


#' @rdname shift-accessor-method
setMethod(
  "get_shift",
  signature(object = "transformationYeoJohnson"),
  function(object, ...) {

    return(object@shift)
  }
)
