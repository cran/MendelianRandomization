# Calculating the MR-Egger estimate (and relevant statistics) using the input data

#' Calculates confidence intervals for the MR-Egger method
#'
#' @description Internal function for calculating confidence intervals for the MR-Egger method.
#'
#' @param type "l" for lower, "u" for upper.
#' @param dist "normal" for normal distribution, "t-dist" for t-distribution.
#' @param mean Causal estimate.
#' @param se Standard error of estimate.
#' @param df Degrees of freedom (for t-distribution).
#' @param .rse Residual standard error.
#' @param .alpha Significance level.
#'
#' @details The slight complication of this function is that when the estimate of the residual standard error (RSE) is less than one, the t-distribution confidence interval is calculated as either the confidence interval using a normal deviate and setting the RSE to 1, or a t-distribution deviate and using the estimated RSE. The wider of the two intervals is reported. This ensures that under-dispersion is not doubly penalized, while also making sure that the estimate is no more precise than that from a fixed-effect analysis.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Numeric value of confidence interval limit.
#'
#' @examples egger.bounds(type = "l", dist = "normal", .theta = 0, .thetase = 1, df = 0, .rse = 1, .alpha = 0.05)
#'
#' @export


egger.bounds <- function(type, dist, .theta, .thetase, df = 0, .rse, .alpha){

  x <- 1 - .alpha/2

  if(type == "u"){

    if (dist == "normal") return(.theta + qnorm(x)*.thetase)
    else if (dist == "t-dist") return(ifelse(.rse < 1,
                                           max(.theta + qnorm(x)*.thetase, .theta + qt(x, df = df)*.thetase*.rse),
                                           .theta + qt(x, df = df)*.thetase))

  } else if (type == "l") {

    if (dist == "normal") return(.theta - qnorm(x)*.thetase)
    else if (dist == "t-dist") return(ifelse(.rse < 1,
                                           min(.theta - qnorm(x)*.thetase, .theta - qt(x, df = df)*.thetase*.rse),
                                           .theta - qt(x, df = df)*.thetase))

  } else { return(NA) }
}

# Calculating the estimate (and relevant statistics) using all the methods specified
#' @docType methods
#' @rdname mr_egger

setMethod("mr_egger",
          "MRInput",
          function(object,
                   robust = FALSE,
                   penalized = FALSE,
                   correl = FALSE,
                   distribution = "normal",
                   alpha = 0.05, ...){

            if(length(object@betaX) < 3){
              cat("Method requires data on >2 variants.")
              return()
            } else if(distribution %in% c("normal", "t-dist")){

            # Ensure that the betaX estimate has positive values
            By = sign(object@betaX)*object@betaY
            rho = object@correlation
            Bx = abs(object@betaX)
            Bxse = object@betaXse
            Byse = object@betaYse


            nsnps <- length(Bx)

            if(!is.na(sum(rho))) { correl = TRUE }

            if(correl == TRUE){

              if(is.na(sum(rho))){

                cat("Correlation matrix not given.")

              } else {

            rho = object@correlation*(sign(object@betaX)%o%sign(object@betaX))

                omega <- Byse%o%Byse*rho
                theta.vals <- solve(t(cbind(rep(1, nsnps), Bx))%*%solve(omega)%*%cbind(rep(1, nsnps), Bx))%*%t(cbind(rep(1, nsnps), Bx))%*%solve(omega)%*%By
                thetaE <- theta.vals[2]
                thetaInter <- theta.vals[1]
                # first component is intercept term, second component is slope term (causal estimate)
                rse <- By - thetaInter - thetaE*Bx
                  rse.corr = as.numeric(sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-2)))

                  sigma <- solve(t(cbind(rep(1, nsnps), Bx))%*%solve(omega)%*%cbind(rep(1, nsnps), Bx))*max(sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-2)),1)
                  thetaEse <- sqrt(sigma[2,2])/min(1, rse.corr)
                  thetaInterse <- sqrt(sigma[1,1])/min(1, rse.corr)



                ciUpper <- egger.bounds("u", distribution, thetaE, thetaEse, nsnps - 2, rse.corr, alpha)
                ciLower <- egger.bounds("l", distribution, thetaE, thetaEse, nsnps - 2, rse.corr, alpha)
                ciUpperInter <- egger.bounds("u", distribution, thetaInter, thetaInterse, nsnps - 2, rse.corr, alpha)
                ciLowerInter <- egger.bounds("l", distribution, thetaInter, thetaInterse, nsnps - 2, rse.corr, alpha)

                heter.stat <- sum(((1/Byse)*(By - thetaInter - thetaE*Bx))^2)
                pvalue.heter.stat <- pchisq(heter.stat, df = nsnps - 2, lower.tail = F)

if (distribution=="normal") { pleiotropic.pvalue <- 2*(1-pnorm(abs(thetaInter/thetaInterse))) }
if (distribution=="t-dist") { 
                pleiotropic.pvalue <- ifelse(rse.corr < 1,
                                        max(2*(1-pnorm(abs(thetaInter/thetaInterse))),
                                            2*(1-pt(abs(thetaInter/thetaInterse/rse.corr),
                                                    df = length(Bx)-2))),

                                        2*(1-pt(abs(thetaInter/thetaInterse),
                                                df = length(Bx)-2)))
 }


if (distribution=="normal") { causal.pvalue <- 2*(1-pnorm(abs(thetaE/thetaEse))) }
if (distribution=="t-dist") { 
                causal.pvalue <- ifelse(rse.corr < 1,
                                        max(2*(1-pnorm(abs(thetaE/thetaEse))),
                                            2*(1-pt(abs(thetaE/thetaEse/rse.corr),
                                                    df = length(Bx)-2))),

                                        2*(1-pt(abs(thetaE/thetaEse),
                                                df = length(Bx)-2)))
 }


                  heter.stat <- (length(Bx) - 2)*(rse.corr^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = length(Bx)-2, lower.tail = F)

                return(new("Egger",
                           Model = "random",
                           Exposure = object@exposure,
                           Outcome = object@outcome,

                           Robust = robust,
                           Penalized = penalized,
                           Correlation = rho,

                           Estimate = thetaE,
                           StdError.Est = thetaEse,
                           CILower.Est = ciLower,
                           CIUpper.Est = ciUpper,
                           Pvalue.Est = causal.pvalue,

                           Intercept = thetaInter,
                           StdError.Int = thetaInterse,
                           CILower.Int = ciLowerInter,
                           CIUpper.Int = ciUpperInter,
                           Pvalue.Int = pleiotropic.pvalue,

                           Pleio.pval = pleiotropic.pvalue,
                           Causal.pval = causal.pvalue,

                           Alpha = alpha,

                           SNPs = nsnps,
                           RSE = as.numeric(rse.corr),
                           Heter.Stat = c(heter.stat, pvalue.heter.stat),
                           I.sq = NaN ))
              }

            } else {

                if(robust == TRUE){
                  if(penalized == TRUE){
                    # method for random : robust and penalized
                    summary.intercept <- summary(lm(By~Bx, weights = 1/(Byse^2), ...))
                    intercept <- summary.intercept$coef[1,1]
                    interceptSE <- summary.intercept$coef[1,2]

                    pvalue.intercept <- 2*(1-pt(abs(intercept/interceptSE), df = nsnps))
                    slope <- summary.intercept$coef[2,1]
                    slopeSE <- summary.intercept$coef[2,2]

                    pen.weights <- pchisq((1/Byse^2)*(By - intercept - slope*Bx)^2, df = 1, lower.tail = F)
                    r.weights <- Byse^(-2)*pmin(1, pen.weights*100)


                    penalised.summary <- summary(lmrob(By ~ Bx, weights = r.weights,  k.max = 500, maxit.scale = 500, ...))
                    rse <- penalised.summary$sigma

                    thetaE <- penalised.summary$coef[2,1]
                    thetaEse <- penalised.summary$coef[2,2]/min(penalised.summary$sigma, 1)
                    ciUpper <- egger.bounds("u", distribution, thetaE, thetaEse, nsnps - 2, rse, alpha)
                    ciLower <- egger.bounds("l", distribution, thetaE, thetaEse, nsnps - 2, rse, alpha)

                    thetaInter <- penalised.summary$coef[1,1]
                    thetaInterse <- penalised.summary$coef[1,2]/min(penalised.summary$sigma, 1)
                    ciUpperInter <- egger.bounds("u", distribution, thetaInter, thetaInterse, nsnps - 2, rse, alpha)
                    ciLowerInter <- egger.bounds("l", distribution, thetaInter, thetaInterse, nsnps - 2, rse, alpha)

                    heter.stat <- NaN
                    pvalue.heter.stat <- NaN

                  } else {
                    # method for random : robust (not penalised)
                    robust.summary <- summary(lmrob(By ~ Bx, weights = (Byse)^(-2),  k.max = 500, maxit.scale = 500, ...))
                    pvalueInt <- robust.summary$coef[1,4]
                    pvalue <- robust.summary$coef[2,4]
                    rse <- robust.summary$sigma

                    thetaE <- robust.summary$coef[2,1]
                    thetaEse <- robust.summary$coef[2,2]/min(robust.summary$sigma,1)
                    ciUpper <- egger.bounds("u", distribution, thetaE, thetaEse, length(Bx) - 2, rse, alpha)
                    ciLower <- egger.bounds("l", distribution, thetaE, thetaEse, length(Bx) - 2, rse, alpha)

                    thetaInter <- robust.summary$coef[1,1]
                    thetaInterse <- robust.summary$coef[1,2]/min(robust.summary$sigma, 1)
                    ciUpperInter <- egger.bounds("u", distribution, thetaInter, thetaInterse, length(Bx) - 2, rse, alpha)
                    ciLowerInter <- egger.bounds("l", distribution, thetaInter, thetaInterse, length(Bx) - 2, rse, alpha)

                    heter.stat <- sum(((1/Byse)*(By - thetaInter - thetaE*Bx))^2)
                    pvalue.heter.stat <- pchisq(heter.stat, df = length(Bx) - 2, lower.tail = F)

                  }

                } else if (penalized == TRUE & robust == FALSE){
                  # method for random : penalized (not robust)
                  summary.intercept <- summary(lm(By~Bx, weights = 1/(Byse^2), ...))
                  intercept <- summary.intercept$coef[1,1]
                  interceptSE <- summary.intercept$coef[1,2]
                  pvalue.intercept <- 2*(1-pt(abs(intercept/interceptSE), df = length(Bx)))
                  slope <- summary.intercept$coef[2,1]
                  slopeSE <- summary.intercept$coef[2,2]

                  pen.weights <- pchisq((1/Byse^2)*(By - intercept - slope*Bx)^2, df = 1, lower.tail = F)
                  r.weights <- Byse^(-2)*pmin(1, pen.weights*100)


                  penalised.summary <- summary(lm(By ~ Bx, weights = r.weights, ...))
                  pvalueInt <- penalised.summary$coef[1,4]
                  pvalue <- penalised.summary$coef[2,4]
                  rse <- penalised.summary$sigma

                  thetaE <- penalised.summary$coef[2,1]
                  thetaEse <- penalised.summary$coef[2,2]/min(penalised.summary$sigma, 1)
                  ciUpper <- egger.bounds("u", distribution, thetaE, thetaEse, length(Bx) - 2, rse, alpha)
                  ciLower <- egger.bounds("l", distribution, thetaE, thetaEse, length(Bx) - 2, rse, alpha)

                  thetaInter <- penalised.summary$coef[1,1]
                  thetaInterse <- penalised.summary$coef[1,2]/min(penalised.summary$sigma, 1)
                  ciUpperInter <- egger.bounds("u", distribution, thetaInter, thetaInterse, length(Bx) - 2, rse, alpha)
                  ciLowerInter <- egger.bounds("l", distribution, thetaInter, thetaInterse, length(Bx) - 2, rse, alpha)

                  heter.stat <- NaN
                  pvalue.heter.stat <- NaN

                } else {
                  # method for random : not robust or penalised
                  summary <- summary(lm(By ~ Bx, weights = (Byse)^(-2), ...))
                  pvalueInt <- summary$coef[1,4]
                  pvalue <- summary$coef[2,4]
                  rse <- summary$sigma

                  thetaE <- summary$coef[2,1]
                  thetaEse <- summary$coef[2,2]/min(summary$sigma,1)
                  ciUpper <- egger.bounds("u", distribution, thetaE, thetaEse, length(Bx) - 2, rse, alpha)
                  ciLower <- egger.bounds("l", distribution, thetaE, thetaEse, length(Bx) - 2, rse, alpha)

                  thetaInter <- summary$coef[1,1]
                  thetaInterse <- summary$coef[1,2]/min(summary$sigma,1)
                  ciUpperInter <- egger.bounds("u", distribution, thetaInter, thetaInterse, length(Bx) - 2, rse, alpha)
                  ciLowerInter <- egger.bounds("l", distribution, thetaInter, thetaInterse, length(Bx) - 2, rse, alpha)

                  heter.stat <- sum(((1/Byse)*(By - thetaInter - thetaE*Bx))^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = length(Bx) - 2, lower.tail = F)

                }

if (distribution=="normal") { pleiotropic.pvalue <- 2*(1-pnorm(abs(thetaInter/thetaInterse))) }
if (distribution=="t-dist") { 
                pleiotropic.pvalue <- ifelse(rse < 1,
                                        max(2*(1-pnorm(abs(thetaInter/thetaInterse))),
                                            2*(1-pt(abs(thetaInter/thetaInterse/rse),
                                                    df = length(Bx)-2))),

                                        2*(1-pt(abs(thetaInter/thetaInterse),
                                                df = length(Bx)-2)))
 }

if (distribution=="normal") { causal.pvalue <- 2*(1-pnorm(abs(thetaE/thetaEse))) }
if (distribution=="t-dist") { 
                causal.pvalue <- ifelse(rse < 1,
                                        max(2*(1-pnorm(abs(thetaE/thetaEse))),
                                            2*(1-pt(abs(thetaE/thetaEse/rse),
                                                    df = length(Bx)-2))),

                                        2*(1-pt(abs(thetaE/thetaEse),
                                                df = length(Bx)-2)))
 }


Q = sum((Bxse/Byse)^-2*(Bx/Byse-weighted.mean(Bx/Byse, w=(Bxse/Byse)^-2))^2)
Isq = max(0, (Q-(length(Bx)-1))/Q)

                return(new("Egger",
                           Model = "random",
                           Exposure = object@exposure,
                           Outcome = object@outcome,

                           Robust = robust,
                           Penalized = penalized,
                           Correlation = matrix(),

                           Estimate = thetaE,
                           StdError.Est = thetaEse,
                           CILower.Est = ciLower,
                           CIUpper.Est = ciUpper,
                           Pvalue.Est = causal.pvalue,

                           Intercept = thetaInter,
                           StdError.Int = thetaInterse,
                           CILower.Int = ciLowerInter,
                           CIUpper.Int = ciUpperInter,
                           Pvalue.Int = pleiotropic.pvalue,

                           Pleio.pval = pleiotropic.pvalue,
                           Causal.pval = causal.pvalue,

                           Alpha = alpha,

                           SNPs = length(object@betaY),
                           RSE = rse,
                           Heter.Stat = c(heter.stat, pvalue.heter.stat),
                           I.sq = as.numeric(Isq) ))
              }
 } else {      cat("Distribution must be one of : normal, t-dist. \n")
                cat("See documentation for details. \n")
              }
 }
)