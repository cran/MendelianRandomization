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
          "MRMBE",
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

#' @rdname getter
setMethod("$",
          "MaxLik",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MVIVW",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MVEgger",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MRMVInput",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MRHetPen",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MRConMix",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MVMedian",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MVLasso",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MRLasso",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "DIVW",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MRcML",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "PIVW",
          function(x, name) slot(x, name))

#' @rdname getter
setMethod("$",
          "MVMRcML",
          function(x, name) slot(x, name))

