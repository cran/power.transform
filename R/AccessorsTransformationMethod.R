# get_transformation_method (generic) ---------------------------------------------------------

#' @title Get transformation method
#'
#' @description Get the transformation method of a transformer object.
#'
#' @param object Transformer object
#' @param ... Unused arguments
#'
#' @return Transformation method
#' @export
#'
#' @rdname transformation-method-accessor-method
setGeneric(
  "get_transformation_method",
  function(object, ...) standardGeneric("get_transformation_method"))


#' @rdname transformation-method-accessor-method
setMethod(
  "get_transformation_method",
  signature(object = "transformationPowerTransform"),
  function(object, ...) {
    return(object@method)
  }
)
