# Calculating the inverse-variance weighted estimate (and relevant statistics) using the input data

# Extra Function

#' Calculates p-values for penalization of weights
#'
#' @description Internal function for calculating penalized weights in conjunction with \code{r.weights}.
#'
#' These weights are used in either the \code{mr_ivw} or \code{mr_egger} functions when \code{penalized = TRUE}, or in the \code{mr_median} function when \code{method = "penalized"}.
#'
#' @param .Bx Beta-coefficient for genetic association with the risk factor.
#' @param .Bxse Standard error of genetic association with the risk factor.
#' @param .By Beta-coefficient for genetic association with the outcome.
#' @param .Byse Standard error of genetic association with the outcome.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return P-values corresponding to heterogeneity test for each genetic variant in turn (based on a chi-squared-1 distribution).
#'
#' @examples penalised.weights(ldlc, ldlcse, chdlodds, chdloddsse)
#'
#' @export

penalised.weights <- function(.Bx, .Bxse, .By, .Byse) {
  theta <- median(.By/.Bx)
  pen.weights <- pchisq((.Bx^2/.Byse^2)*(.By/.Bx - theta)^2, df = 1, lower.tail = F)

  return(pen.weights)
}

#' Calculates p-values for penalization of weights with second-order weights
#'
#' @description Internal function for calculating penalized weights in conjunction with \code{r.weights} with \code{weights} equal to \code{"delta"}.
#'
#' These weights are used in either the \code{mr_ivw} or \code{mr_egger} functions when \code{penalized = TRUE}, or in the \code{mr_median} function when \code{method = "penalized"}.
#'
#' @param .Bx Beta-coefficient for genetic association with the risk factor.
#' @param .Bxse Standard error of genetic association with the risk factor.
#' @param .By Beta-coefficient for genetic association with the outcome.
#' @param .Byse Standard error of genetic association with the outcome.
#' @param .psi The correlation between the genetic associations with the exposure and the association with the outcome for each variant resulting from sample overlap.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return P-values corresponding to heterogeneity test for each genetic variant in turn (based on a chi-squared-1 distribution).
#'
#' @examples penalised.weights.delta(ldlc, ldlcse, chdlodds, chdloddsse, 0)
#'
#' @export

penalised.weights.delta <- function(.Bx, .Bxse, .By, .Byse, .psi) {
#  theta <- sum(.By/.Bx*(.Byse^2/.Bx^2+.By^2*.Bxse^2/.Bx^4-2*.psi*.By*.Bxse*.Byse/.Bx^3)^-1)/sum((.Byse^2/.Bx^2+.By^2*.Bxse^2/.Bx^4-2*.psi*.By*.Bxse*.Byse/.Bx^3)^-1)
  theta <- median(.By/.Bx)
  pen.weights <- pchisq(((.Byse^2/.Bx^2+.By^2*.Bxse^2/.Bx^4-2*.psi*.By*.Bxse*.Byse/.Bx^3)^-1)*(.By/.Bx - theta)^2, df = 1, lower.tail = F)

  return(pen.weights)
}


#' Calculates p-values for penalization of weights
#'
#' @description Internal function for calculating penalized weights in conjunction with \code{penalised.weights}.
#'
#' These weights are used in either the \code{mr_ivw} or \code{mr_egger} functions when \code{penalized = TRUE}.
#'
#' @param .Byse Standard error of genetic association with the outcome.
#' @param .pen.weights Factors for penalizing weights.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Penalized weights.
#'
#' @examples r.weights(chdloddsse, penalised.weights(ldlc, ldlcse, chdlodds, chdloddsse))
#'
#' @export

r.weights <- function(.Byse, .pen.weights){
  return(.Byse^(-2)*pmin(1, .pen.weights*100))
}

#' Calculates p-values for penalization of weights with second-order weights
#'
#' @description Internal function for calculating penalized weights in conjunction with \code{penalised.weights} with \code{weights} equal to \code{"delta"}.
#'
#' These weights are used in either the \code{mr_ivw} or \code{mr_egger} functions when \code{penalized = TRUE}.
#'
#' @param .Bx Beta-coefficient for genetic association with the risk factor.
#' @param .Bxse Standard error of genetic association with the risk factor.
#' @param .By Beta-coefficient for genetic association with the outcome.
#' @param .Byse Standard error of genetic association with the outcome.
#' @param .psi The correlation between the genetic associations with the exposure and the association with the outcome for each variant resulting from sample overlap.
#' @param .pen.weights Factors for penalizing weights.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Penalized weights.
#'
#' @examples r.weights.delta(ldlc, ldlcse, chdlodds, chdloddsse, 0,
#'   penalised.weights.delta(ldlc, ldlcse, chdlodds, chdloddsse, 0))
#'
#' @export

r.weights.delta <- function(.Bx, .Bxse, .By, .Byse, .psi, .pen.weights){
  return((.Byse^2 + .By^2*.Bxse^2/.Bx^2-2*.psi*.By*.Bxse*.Byse/.Bx)^-1*pmin(1, .pen.weights*20))
}


#' @docType methods
#' @rdname mr_ivw

setMethod("mr_ivw",
          "MRInput",
          function(object,
                   model = "default",
                   robust = FALSE,
                   penalized = FALSE,
                   weights = "simple",
                   psi = 0,
                   correl = FALSE,
                   distribution = "normal",
                   alpha = 0.05, ...){

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            rho = object@correlation

            nsnps <- length(Bx)

            if(model == "default"){
              if(nsnps < 4) {
                model <- "fixed"
              } else {
                model <- "random"
              }
            }

            if(!is.na(sum(rho))) { correl = TRUE }

              if(model %in% c("random", "fixed") & distribution %in% c("normal", "t-dist") & weights %in% c("simple", "delta")){


            if(correl == T){

              if(is.na(sum(rho))){

                cat("Correlation matrix not given.")

              } else {

                  omega <- Byse%o%Byse*rho

                  thetaIVW <- as.numeric(solve(t(Bx)%*%solve(omega)%*%Bx)*t(Bx)%*%solve(omega)%*%By)
                    rse <- By - thetaIVW*Bx

                    if(model == "random") {
                      thetaIVWse <- sqrt(solve(t(Bx)%*%solve(omega)%*%Bx))*max(sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-1)),1)
                    } else if (model == "fixed"){
                      thetaIVWse <- sqrt(solve(t(Bx)%*%solve(omega)%*%Bx))
                    }
                  }

                  robust <- FALSE
                  penalized <- FALSE
                  correlation <- TRUE
                  weights <- "simple"

                  if(distribution == "normal"){
                    ciLower <- ci_normal("l", thetaIVW, thetaIVWse, alpha)
                    ciUpper <- ci_normal("u", thetaIVW, thetaIVWse, alpha)
                  } else if (distribution == "t-dist"){
                    ciLower <- ci_t("l", thetaIVW, thetaIVWse, nsnps - 1, alpha)
                    ciUpper <- ci_t("u", thetaIVW, thetaIVWse, nsnps - 1, alpha)
                  }

                  rse.corr = sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-1))
                  heter.stat <- (length(Bx) - 1)*(rse.corr^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = length(Bx)-1, lower.tail = F)
  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=length(Bx)-1) }

 fstat = sum(((Bx/Bxse)%*%solve(chol(rho)))^2)/nsnps

                  return(new("IVW",
                             Model = model,
                             Exposure = object@exposure,
                             Outcome = object@outcome,

                             Robust = robust,
                             Penalized = penalized,
                             Correlation = object@correlation,

                             Estimate = as.numeric(thetaIVW),
                             StdError = as.numeric(thetaIVWse),
                             CILower =  as.numeric(ciLower),
                             CIUpper = as.numeric(ciUpper),

                             SNPs = nsnps,
                             Pvalue = as.numeric(pvalue),

                             Alpha = alpha,

                             RSE = as.numeric(rse.corr),
                             Heter.Stat = c(heter.stat,
                             pvalue.heter.stat),
                             Fstat = fstat))
            } else {

            if(nsnps == 1){

              thetaIVW <- By/Bx
if (weights == "simple") {    thetaIVWse <- abs(Byse/Bx)  }
if (weights == "delta")  {    thetaIVWse <- sqrt(Byse^2/Bx^2+By^2*Bxse^2/Bx^4-2*psi*By*Bxse*Byse/Bx^3) }
              rse <- 1
              pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse))

              robust <- FALSE
              penalized <- FALSE
              correlation <- FALSE

              if(distribution == "normal"){
                ciLower <- ci_normal("l", thetaIVW, thetaIVWse, alpha)
                ciUpper <- ci_normal("u", thetaIVW, thetaIVWse, alpha)
              } else if (distribution == "t-dist"){
                ciLower <- ci_t("l", thetaIVW, thetaIVWse, length(Bx) - 1, alpha)
                ciUpper <- ci_t("u", thetaIVW, thetaIVWse, length(Bx) - 1, alpha)
              }
    fstat = (Bx/Bxse)^2

              return(new("IVW",
                         Model = model,
                         Exposure = object@exposure,
                         Outcome = object@outcome,

                         Robust = robust,
                         Penalized = penalized,
                         Correlation = object@correlation,

                         Estimate = thetaIVW,
                         StdError = thetaIVWse,
                         CILower =  ciLower,
                         CIUpper = ciUpper,

                         SNPs = length(By),
                         Pvalue = pvalue,

                         Alpha = alpha,

                         RSE = rse,
                         Heter.Stat = c(NaN, NaN),
                         Fstat = fstat))

            } else if (nsnps > 1) {
                if(robust == TRUE){
                  if(penalized == TRUE){
                    # method : robust and penalized
if (weights == "simple") {  pen.weights <- penalised.weights(Bx, Bxse, By, Byse)
                    penalised.robust.summary <- summary(lmrob(By ~ Bx - 1, weights = r.weights(Byse, pen.weights), k.max = 500, maxit.scale = 500, ...))  }
if (weights == "delta") {  pen.weights <- penalised.weights.delta(Bx, Bxse, By, Byse, psi)
                    penalised.robust.summary <- summary(lmrob(By ~ Bx - 1, weights = r.weights.delta(Bx, Bxse, By, Byse, psi, pen.weights), k.max = 500, maxit.scale = 500, ...)) }


                    thetaIVW <- penalised.robust.summary$coef[1]

                    if(model == "random") thetaIVWse <- penalised.robust.summary$coef[1,2]/min(penalised.robust.summary$sigma, 1)
                    else thetaIVWse <- penalised.robust.summary$coef[1,2]/penalised.robust.summary$sigma

  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=length(Bx)-1) }
                    rse <- penalised.robust.summary$sigma

                    heter.stat <- NaN
                    pvalue.heter.stat <- NaN

                  } else {
                    # method : robust (not penalised)
if (weights == "simple") {
    robust.summary <- summary(lmrob(By ~ Bx - 1, weights = Byse^(-2),  k.max = 500, maxit.scale = 500, ...)) }
if (weights == "delta") {
    robust.summary <- summary(lmrob(By ~ Bx - 1, weights = (Byse^2 + By^2*Bxse^2/Bx^2-2*psi*By*Bxse*Byse/Bx)^-1,  k.max = 500, maxit.scale = 500, ...)) }
                    thetaIVW <- robust.summary$coef[1]

                    if(model == "random") thetaIVWse <- robust.summary$coef[1,2] / min(robust.summary$sigma, 1)
                    else thetaIVWse <- robust.summary$coef[1,2] / robust.summary$sigma

  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=length(Bx)-1) }
                    rse <- robust.summary$sigma

                    heter.stat <- (length(Bx) - 1)*(rse^2)
                    pvalue.heter.stat <- pchisq(heter.stat, df = length(Bx)-1, lower.tail = F)

                  }

                } else if (penalized == TRUE & robust == FALSE){
                  # method : penalized (not robust)
if (weights == "simple") {  pen.weights <- penalised.weights(Bx, Bxse, By, Byse)
                  penalised.summary <- summary(lm(By ~ Bx - 1, weights = r.weights(Byse, pen.weights), ...)) }
if (weights == "delta") {  pen.weights <- penalised.weights.delta(Bx, Bxse, By, Byse, psi) 
                  penalised.summary <- summary(lm(By ~ Bx - 1, weights = r.weights.delta(Bx, Bxse, By, Byse, psi, pen.weights), ...)) }


                  thetaIVW <- penalised.summary$coef[1]

                  if(model == "random") thetaIVWse <- penalised.summary$coef[1,2]/min(penalised.summary$sigma, 1)
                  else thetaIVWse <- penalised.summary$coef[1,2]/penalised.summary$sigma

  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=length(Bx)-1) }
                  rse <- penalised.summary$sigma

                  heter.stat <- NaN
                  pvalue.heter.stat <- NaN

                } else {
                  # method : not robust or penalised
if (weights == "simple") { summary <- summary(lm(By ~ Bx - 1, weights = Byse^(-2), ...))}
if (weights == "delta") { summary <- summary(lm(By ~ Bx - 1, weights =(Byse^2 + By^2*Bxse^2/Bx^2-2*psi*By*Bxse*Byse/Bx)^-1, ...))}

                  thetaIVW <- summary$coef[1]

                  if(model == "random") thetaIVWse <- summary$coef[1,2]/min(summary$sigma,1)
                  else thetaIVWse <- summary$coef[1,2]/summary$sigma

  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=length(Bx)-1) }
                  rse <- summary$sigma

                  heter.stat <- (length(Bx) - 1)*(rse^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = length(Bx)-1, lower.tail = F)

                }

                if(distribution == "normal"){
                  ciLower <- ci_normal("l", thetaIVW, thetaIVWse, alpha)
                  ciUpper <- ci_normal("u", thetaIVW, thetaIVWse, alpha)
                } else if (distribution == "t-dist"){
                  ciLower <- ci_t("l", thetaIVW, thetaIVWse, length(Bx) - 1, alpha)
                  ciUpper <- ci_t("u", thetaIVW, thetaIVWse, length(Bx) - 1, alpha)
                }

  fstat = sum((Bx/Bxse)^2)/nsnps

                return(new("IVW",
                           Model = model,
                           Exposure = object@exposure,
                           Outcome = object@outcome,

                           Robust = robust,
                           Penalized = penalized,
                           Correlation = object@correlation,

                           Estimate = thetaIVW,
                           StdError = thetaIVWse,
                           CILower =  ciLower,
                           CIUpper = ciUpper,

                           SNPs = length(By),
                           Pvalue = pvalue,

                           Alpha = alpha,

                           RSE = rse,
                           Heter.Stat = c(heter.stat,
                           pvalue.heter.stat),
                           Fstat = fstat))

            } else {

              cat("No SNPs detected.")

            } } }
 else {
                cat("Model type must be one of: default, random, fixed. \n")
                cat("Distribution must be one of : normal, t-dist. \n")
                cat("Weights must be one of : simple, delta. \n")
                cat("See documentation for details. \n")
              }
          }
)