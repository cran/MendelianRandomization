#' Calculates log-likelihood with uncorrelated variants in two-sample setting (no correlation from sample overlap)
#'
#' @description Internal function for calculating log-likelihood for for maximum-likelihood method.
#'
#' @param param Parameters in the model.
#' @param x Genetic associations with the exposure.
#' @param sigmax Standard errors of associations with the exposure.
#' @param y Genetic associations with the outcome.
#' @param sigmay Standard errors of genetic associations with the outcome.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Log-likelihood.
#'
#' @export

loglikelihood <- function(param, x, sigmax, y, sigmay) {                    # log-likelihood function
  return(1/2*sum((x-param[1:length(x)])^2/sigmax^2)+1/2*
   sum((y-param[length(x)+1]*param[1:length(x)])^2/sigmay^2)) }

#' Calculates log-likelihood with correlated variants in two-sample setting (no correlation from sample overlap)
#'
#' @description Internal function for calculating log-likelihood for for maximum-likelihood method.
#'
#' @param param Parameters in the model.
#' @param x Genetic associations with the exposure.
#' @param Taux Precision matrix for associations with the exposure.
#' @param y Genetic associations with the outcome.
#' @param Tauy Precision matrix for associations with the outcome.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Log-likelihood.
#'
#' @export

loglikelihoodcorrel <- function(param, x, Taux, y, Tauy) {               # log-likelihood function
   return(1/2*t(x-param[1:length(x)])%*%Taux%*%(x-param[1:length(x)])+1/2*
     t(y-param[length(x)+1]*param[1:length(x)])%*%Tauy%*%
      (y-param[length(x)+1]*param[1:length(x)])) }

#' Calculates log-likelihood with correlation from sample overlap
#'
#' @description Internal function for calculating log-likelihood for for maximum-likelihood method.
#'
#' @param param Parameters in the model.
#' @param x Genetic associations with the exposure.
#' @param y Genetic associations with the outcome.
#' @param Tauxy Precision matrix for associations with the exposure and outcome.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Log-likelihood.
#'
#' @export

loglikelihoodrhocorrel <- function(param, x, y, Tauxy) {               # log-likelihood function
   return(1/2*t(c(x-param[1:length(x)], y-param[length(x)+1]*param[1:length(x)]))%*%Tauxy%*%
        c(x-param[1:length(x)], y-param[length(x)+1]*param[1:length(x)])) }


#' @docType methods
#' @rdname mr_maxlik

setMethod("mr_maxlik",
          "MRInput",
          function(object,
                   model = "default", 
                   correl = FALSE,
                   psi = 0,
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

              if(model %in% c("random", "fixed") & distribution %in% c("normal", "t-dist")){

if (psi == 0) {
 if (correl == FALSE) {

opt <- optim(c(Bx, sum(Bx*By/Byse^2)/sum(Bx^2/Byse^2)),
            loglikelihood, x=Bx, sigmax=Bxse, y=By, sigmay=Byse, hessian=TRUE, control = list(maxit=25000), ...)
               # optimization command

    }
 else if (correl == TRUE) {
         Sigmax = Bxse%o%Bxse*rho
         Sigmay = Byse%o%Byse*rho
         Taux.in = solve(Sigmax); Tauy.in = solve(Sigmay)

opt <- optim(c(Bx, sum(Bx*By/Byse^2)/sum(Bx^2/Byse^2)),
            loglikelihoodcorrel, x=Bx, Taux=Taux.in, y=By, Tauy=Tauy.in, hessian=TRUE, control = list(maxit=25000), ...)
  }
 }
else if (psi != 0) {
 if (correl == FALSE) {  rho = diag(rep(1, nsnps))  }
         Sigmax  = Bxse%o%Bxse*rho
         Sigmay  = Byse%o%Byse*rho
         Sigmaxy = Bxse%o%Byse*rho*psi

         Sigmaxyall = rbind(cbind(Sigmax, Sigmaxy), cbind(Sigmaxy, Sigmay))

Tauxy.in = solve(Sigmaxyall)

opt <- optim(c(Bx, sum(Bx*By/Byse^2)/sum(Bx^2/Byse^2)),
            loglikelihoodrhocorrel, x=Bx, y=By, Tauxy=Tauxy.in, hessian=TRUE, control = list(maxit=25000), ...)
         }
               # optimization command

thetaML   <- opt$par[nsnps+1]
if (model == "fixed") {
thetaMLse <- sqrt(solve(opt$hessian)[nsnps+1,nsnps+1]) }
if (model == "random") {
thetaMLse <- sqrt(solve(opt$hessian)[nsnps+1,nsnps+1] * max(2*opt$value/(nsnps-1), 1)) }

rse = sqrt(2*opt$value/(nsnps-1))

                  if(distribution == "normal"){
                    ciLower <- ci_normal("l", thetaML, thetaMLse, alpha)
                    ciUpper <- ci_normal("u", thetaML, thetaMLse, alpha)
                  } else if (distribution == "t-dist"){
                    ciLower <- ci_t("l", thetaML, thetaMLse, nsnps - 1, alpha)
                    ciUpper <- ci_t("u", thetaML, thetaMLse, nsnps - 1, alpha)
                  }

                  heter.stat <- 2*opt$value
                  pvalue.heter.stat <- pchisq(2*opt$value, df=length(Bx)-1, lower.tail=FALSE)
  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaML/thetaMLse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaML/thetaMLse), df=length(Bx)-1) }

                  return(new("MaxLik",
                             Model = model,
                             Exposure = object@exposure,
                             Outcome = object@outcome,

                             Correlation = object@correlation,
                             Psi = as.numeric(psi),

                             Estimate = as.numeric(thetaML),
                             StdError = as.numeric(thetaMLse),
                             CILower =  as.numeric(ciLower),
                             CIUpper = as.numeric(ciUpper),

                             SNPs = nsnps,
                             Pvalue = as.numeric(pvalue),

                             Alpha = alpha,

                             RSE = as.numeric(rse),
                             Heter.Stat = c(heter.stat,
                             pvalue.heter.stat)))
  }
 else {
                cat("Model type must be one of: default, random, fixed. \n")
                cat("Distribution must be one of : normal, t-dist. \n")
                cat("See documentation for details. \n")
              }

})



