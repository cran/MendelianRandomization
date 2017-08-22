# Required Functions
#' Weighted median function
#'
#' @description Internal function for calculating weighted median estimate (or simple median estimate if weights are all equal).
#'
#' @param theta Ratio causal estimates for each genetic variant considered individually.
#' @param weights Weights.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Causal estimate.
#'
#' @examples weighted.median(theta = chdlodds/ldlc, weights = abs(chdloddsse/ldlc))
#'
#' @export


weighted.median <- function(theta, weights){

  theta.sorted <- theta[order(theta)]
  weights.sorted <- weights[order(theta)]
  cumsum <- cumsum(weights.sorted) - 0.5*weights.sorted
  cumsum <- cumsum/sum(weights.sorted)
  k <- length(which(cumsum < 0.5))

  ratio <- (0.5 - cumsum[k])/(cumsum[k+1] - cumsum[k])
  weighted.estimate <- theta.sorted[k] + (theta.sorted[k+1] - theta.sorted[k])*ratio

  return(weighted.estimate)
}

#' Weighted median standard error function
#'
#' @description Internal function for calculating standard error of weighted median estimate (or simple median estimator if weights are all equal) using bootstrapping. The number of iterations and initial value of the random seed can also be set.
#'
#' @param Bx A numeric vector of beta-coefficient values for genetic associations with the exposure.
#' @param Bxse The standard errors associated with the beta-coefficients in \code{Bx}.
#' @param By A numeric vector of beta-coefficient values for genetic associations with the outcome.
#' @param Byse The standard errors associated with the beta-coefficients in \code{By}.
#' @param weights Weights.
#' @param iter The number of bootstrap samples to generate when calculating the standard error.
#' @param seed The random seed to use when generating the bootstrap samples (for reproducibility). If set to \code{NA}, the random seed will not be set (for example, if the function is used as part of a larger simulation).
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Causal estimate.
#'
#' @examples weighted.median.boot.se(Bx = ldlc, By = chdlodds, Bxse = ldlcse, Byse = chdloddsse,
#' weights = chdloddsse, iter = 100, seed = 314)
#'
#' @export

weighted.median.boot.se <- function(Bx, By, Bxse, Byse, weights, iter, seed){

  med <- NULL
  snps <- length(By)

if (!is.na(seed)) { set.seed(seed) }

  for(i in 1:iter){
    Bx.boot <- rnorm(snps, mean = Bx, sd = Bxse)
    By.boot <- rnorm(snps, mean = By, sd = Byse)

    theta.boot <- By.boot/Bx.boot
    med[i] <- weighted.median(theta.boot, weights)
  }
  return(sd(med))
}


#' @docType methods
#' @rdname mr_median

setMethod("mr_median",
          signature = c(object = "MRInput"),
          function(object, weighting = "weighted", 
                   distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265){

            if(length(object@betaX) < 3){
              cat("Method requires data on >2 variants.")
              return()
            } else {

            Bx = object@betaX
            By = object@betaY
            Theta = object@betaY/object@betaX
            Bxse = object@betaXse
            Byse = object@betaYse

            Simple = rep(1/length(Bx), length(Bx))
            Weighted = (Bx/Byse)^2

            if(weighting %in% c("simple", "weighted", "penalized") & distribution %in% c("normal", "t-dist")){

              if(weighting == "simple"){
                thetaWM <- weighted.median(Theta, Simple)
                seBoot <- weighted.median.boot.se(Bx, By, Bxse, Byse, Simple, iterations, seed)
                if(distribution == "normal"){
                  ciLower <- ci_normal("l", thetaWM, seBoot, alpha)
                  ciUpper <- ci_normal("u", thetaWM, seBoot, alpha)
                  pval <- 2*pnorm(-abs(thetaWM/seBoot))
                } else if (distribution == "t-dist"){
                  ciLower <- ci_t("l", thetaWM, seBoot, length(Bx) - 1, alpha)
                  ciUpper <- ci_t("u", thetaWM, seBoot, length(Bx) - 1, alpha)
                  pval <- 2*pt(-abs(thetaWM/seBoot), df=length(Bx) - 1)
                }

              } else if (weighting == "weighted"){
                thetaWM <- weighted.median(Theta, Weighted)
                seBoot <- weighted.median.boot.se(Bx, By, Bxse, Byse, Weighted, iterations, seed)
                if(distribution == "normal"){
                  ciLower <- ci_normal("l", thetaWM, seBoot, alpha)
                  ciUpper <- ci_normal("u", thetaWM, seBoot, alpha)
                  pval <- 2*pnorm(-abs(thetaWM/seBoot))
                } else if (distribution == "t-dist"){
                  ciLower <- ci_t("l", thetaWM, seBoot, length(Bx) - 1, alpha)
                  ciUpper <- ci_t("u", thetaWM, seBoot, length(Bx) - 1, alpha)
                  pval <- 2*pt(-abs(thetaWM/seBoot), df=length(Bx) - 1)
                }

              } else if (weighting == "penalized"){
                penalty <- pchisq(Weighted*(Theta - weighted.median(Theta, Weighted))^2, df = 1, lower.tail = FALSE)
                pen.weights <- Weighted*pmin(1, penalty*20)

                thetaWM <- weighted.median(Theta, pen.weights)
                seBoot <- weighted.median.boot.se(Bx, By, Bxse, Byse, pen.weights, iterations, seed)
                if(distribution == "normal"){
                  ciLower <- ci_normal("l", thetaWM, seBoot, alpha)
                  ciUpper <- ci_normal("u", thetaWM, seBoot, alpha)
                  pval <- 2*pnorm(-abs(thetaWM/seBoot))
                } else if (distribution == "t-dist"){
                  ciLower <- ci_t("l", thetaWM, seBoot, length(Bx) - 1, alpha)
                  ciUpper <- ci_t("u", thetaWM, seBoot, length(Bx) - 1, alpha)
                  pval <- 2*pt(-abs(thetaWM/seBoot), df=length(Bx) - 1)
                }
              }

              return(new("WeightedMedian",
                         Type = weighting,
                         Exposure = object@exposure,
                         Outcome = object@outcome,

                         Estimate = thetaWM,
                         StdError = seBoot,
                         CILower = ciLower,
                         CIUpper = ciUpper,
                         Alpha = alpha,

                         Pvalue = as.numeric(pval),
                         SNPs = length(By))
              )

            } else {
              cat("Weighting must be one of : simple, weighted, penalized. \n")
              cat("Distribution must be one of : normal, t-dist." )
            }

          }}
)