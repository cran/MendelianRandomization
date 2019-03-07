#' @docType methods
#' @rdname mr_mvivw

setMethod("mr_mvivw",
          "MRMVInput",
          function(object,
                   model = "default",
                   correl = FALSE,
                   distribution = "normal",
                   alpha = 0.05, ...){

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            rho = object@correlation

            nsnps <- dim(Bx)[1]

            if(model == "default"){
              if(nsnps < 4) {
                model <- "fixed"
              } else {
                model <- "random"
              }
            }

            if(!is.na(sum(rho))) { correl = TRUE }

              if(model %in% c("random", "fixed") & distribution %in% c("normal", "t-dist")){


            if(correl == TRUE){

              if(is.na(sum(rho))){

                cat("Correlation matrix not given.")

              } else {

                  omega <- Byse%o%Byse*rho

                  thetaIVW <- solve(t(Bx)%*%solve(omega)%*%Bx)%*%t(Bx)%*%solve(omega)%*%By
                    rse <- By - Bx%*%thetaIVW

                    if(model == "random") {
                      thetaIVWse <- sqrt(diag(solve(t(Bx)%*%solve(omega)%*%Bx)))*max(sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-dim(Bx)[2])),1)
                    } else if (model == "fixed"){
                      thetaIVWse <- sqrt(diag(solve(t(Bx)%*%solve(omega)%*%Bx)))
                    }

                  correlation <- TRUE

                  if(distribution == "normal"){
                    ciLower <- ci_normal("l", thetaIVW, thetaIVWse, alpha)
                    ciUpper <- ci_normal("u", thetaIVW, thetaIVWse, alpha)
                  } else if (distribution == "t-dist"){
                    ciLower <- ci_t("l", thetaIVW, thetaIVWse, nsnps - dim(Bx)[2], alpha)
                    ciUpper <- ci_t("u", thetaIVW, thetaIVWse, nsnps - dim(Bx)[2], alpha)
                  }

                  rse.corr = sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-dim(Bx)[2]))
                  heter.stat <- (dim(Bx)[1] - dim(Bx)[2])*(rse.corr^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = dim(Bx)[1]-dim(Bx)[2], lower.tail = F)
  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=nsnps-dim(Bx)[2]) }

                  return(new("MVIVW",
                             Model = model,
                             Exposure = object@exposure,
                             Outcome = object@outcome,

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
                             pvalue.heter.stat)))
            } } else {
                  # not correlated
   summary <- summary(lm(By ~ Bx - 1, weights = Byse^(-2), ...))

                  thetaIVW <- summary$coef[,1]

                  if(model == "random") thetaIVWse <- summary$coef[,2]/min(summary$sigma,1)
                  else thetaIVWse <- summary$coef[,2]/summary$sigma

  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaIVW/thetaIVWse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaIVW/thetaIVWse), df=dim(Bx)[1]-dim(Bx)[2]) }
                  rse <- summary$sigma

                  heter.stat <- (nsnps - dim(Bx)[2])*(rse^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = nsnps-dim(Bx)[2], lower.tail = F)

                if(distribution == "normal"){
                  ciLower <- ci_normal("l", thetaIVW, thetaIVWse, alpha)
                  ciUpper <- ci_normal("u", thetaIVW, thetaIVWse, alpha)
                } else if (distribution == "t-dist"){
                  ciLower <- ci_t("l", thetaIVW, thetaIVWse, dim(Bx)[1]-dim(Bx)[2], alpha)
                  ciUpper <- ci_t("u", thetaIVW, thetaIVWse, dim(Bx)[1]-dim(Bx)[2], alpha)
                }

                return(new("MVIVW",
                           Model = model,
                           Exposure = object@exposure,
                           Outcome = object@outcome,

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
                           pvalue.heter.stat)))

            } } 
 else {
                cat("Model type must be one of: default, random, fixed. \n")
                cat("Distribution must be one of : normal, t-dist. \n")
                cat("See documentation for details. \n")
              }
          }
)
