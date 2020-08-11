#' Applies method values() to different output classes
#'
#' @description Enables the internal function \code{values}, used in the \code{mr_allmethods} function.
#' @docType methods
#' @name values
#' @param object Object (could be an object of class "WeightedMedian", "Egger", or "IVW").
#'
#' @keywords internal
NULL

#' @rdname values
setMethod("values",
          "WeightedMedian",
          function(object){
            return(c(object@Estimate,
                     object@StdError,
                     object@CILower,
                     object@CIUpper,
                     object@Pvalue
                     ))
          }
)

#--------------------------------------------------------------------------------------------

#' @rdname values
setMethod("values",
          "IVW",
          function(object){
            return(c(object@Estimate,
                     object@StdError,
                     object@CILower,
                     object@CIUpper,
                     object@Pvalue
            ))
          }
)

#--------------------------------------------------------------------------------------------

#' @rdname values
setMethod("values",
          "Egger",
          function(object){
            return(rbind(c(object@Estimate,
                     object@StdError.Est,
                     object@CILower.Est,
                     object@CIUpper.Est,
                     object@Pvalue.Est),
                     c(object@Intercept,
                       object@StdError.Int,
                       object@CILower.Int,
                       object@CIUpper.Int,
                       object@Pvalue.Int)))
          }
)
