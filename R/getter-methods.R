#' Applies method $ to different classes
#'
#' Enables slots of objects in this package to be referenced easily.
#' @docType methods
#' @name getter
#' @param x Object.
#' @param name Name of slot. 
#'
#' @keywords internal
NULL

#' @rdname getter
setMethod("$",
          "WeightedMedian",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "IVW",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "Egger",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MRAll",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MRInput",
          function(x, name) slot(x, name))