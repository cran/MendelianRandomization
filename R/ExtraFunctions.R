# Confidence Interval Calculator

#' Calculate confidence intervals using the normal distribution
#'
#' @description Internal function for calculating confidence intervals using the normal distribution.
#'
#' @param type "l" for lower, "u" for upper.
#' @param mean Causal estimate.
#' @param se Standard error of estimate.
#' @param alpha Significance level.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Numeric value of confidence interval limit.
#'
#' @examples ci_normal(type = "l", mean = 0, se = 1, alpha = 0.05)
#'
#' @export

ci_normal <- function(type, mean, se, alpha){
  x <- 1 - alpha/2

  if(type == "l") return(mean - qnorm(x)*se)
  else if (type == "u") return(mean + qnorm(x)*se)
}

#' Calculate confidence intervals using the t-distribution
#'
#' @description Internal function for calculating confidence intervals using the t-distribution.
#'
#' @param type "l" for lower, "u" for upper.
#' @param mean Causal estimate.
#' @param se Standard error of estimate.
#' @param df Number of degrees of freedom.
#' @param alpha Significance level.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Numeric value of confidence interval limit.
#'
#' @examples ci_t(type = "l", mean = 0, se = 1, df = 8, alpha = 0.05)
#'
#' @export

ci_t <- function(type, mean, se, df, alpha){
  x <- 1 - alpha/2

  if(type == "l") return(mean - qt(x, df = df)*se)
  else if (type == "u") return(mean + qt(x, df = df)*se)
}

# Rounding

#' Produce nicely rounded numbers
#'
#' @description Internal function for putting extra zeros on numbers (if needed).
#'
#' @param number A number to be rounded.
#' @param places Number of decimal places.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Returns a nicely formatted number. For example, 0.3 to three decimal places would be rendered as 0.300.
#'
#' @examples decimals(number = 0.3, places = 3)
#'
#' @export


decimals <- function(number, places) format(round(number, places), nsmall = places)

# Name of variable
# var.name <- function(var) deparse(substitute(var))

#' Capitalize a word
#'
#' @description Internal function for capitalizing a word.
#'
#' @param x A word.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Returns a capitalized word.
#'
#' @examples simpleCap(x = "weighted")
#'
#' @export

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}
