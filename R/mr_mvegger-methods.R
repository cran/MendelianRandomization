#' @docType methods
#' @rdname mr_mvegger

setMethod("mr_mvegger",
          "MRMVInput",
          function(object,
                   orientate = 1,
                   correl = FALSE,
                   distribution = "normal",
                   alpha = 0.05){

if (orientate %in% 1:dim(object@betaX)[2]) { orientAte = orientate} else {orientAte = 1}

            orient = sign(object@betaX)[,orientAte]
            Bx = object@betaX*orient
            By = object@betaY*orient
            Bxse = object@betaXse
            Byse = object@betaYse
            rho = object@correlation

            nsnps <- dim(Bx)[1]


            if(!is.na(sum(rho))) { correl = TRUE }

              if(distribution %in% c("normal", "t-dist")){


            if(correl == TRUE){

              if(is.na(sum(rho))){

                cat("Correlation matrix not given.")

              } else {

            rho = object@correlation*(orient%o%orient)
                  omega <- Byse%o%Byse*rho

                  thetaEgger <- solve(t(cbind(Bx, rep(1, nsnps)))%*%solve(omega)%*%cbind(Bx, rep(1, nsnps)))%*%t(cbind(Bx, rep(1, nsnps)))%*%solve(omega)%*%By
                    rse <- By - cbind(Bx, rep(1, nsnps))%*%thetaEgger
                  thetaEggerse <- sqrt(diag(solve(t(cbind(Bx, rep(1, nsnps)))%*%solve(omega)%*%cbind(Bx, rep(1, nsnps)))))*max(sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-dim(Bx)[2]-1)),1)

                  correlation <- TRUE

                  if(distribution == "normal"){
                    ciLower <- ci_normal("l", thetaEgger, thetaEggerse, alpha)
                    ciUpper <- ci_normal("u", thetaEgger, thetaEggerse, alpha)
                  } else if (distribution == "t-dist"){
                    ciLower <- ci_t("l", thetaEgger, thetaEggerse, nsnps - dim(Bx)[2]-1, alpha)
                    ciUpper <- ci_t("u", thetaEgger, thetaEggerse, nsnps - dim(Bx)[2]-1, alpha)
                  }

                  rse.corr = sqrt(t(rse)%*%solve(omega)%*%rse/(nsnps-dim(Bx)[2]-1))
                  heter.stat <- (dim(Bx)[1] - dim(Bx)[2]-1)*(rse.corr^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = dim(Bx)[1]-dim(Bx)[2]-1, lower.tail = F)
  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaEgger/thetaEggerse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaEgger/thetaEggerse), df=nsnps-dim(Bx)[2]-1) }

                  return(new("MVEgger",
                             Model = "random",
                             Orientate = orientAte,
                             Exposure = object@exposure,
                             Outcome = object@outcome,

                             Correlation = object@correlation,

                             Estimate = as.numeric(thetaEgger)[1:dim(Bx)[2]],
                             StdError.Est = as.numeric(thetaEggerse)[1:dim(Bx)[2]],
                             CILower.Est =  as.numeric(ciLower)[1:dim(Bx)[2]],
                             CIUpper.Est = as.numeric(ciUpper)[1:dim(Bx)[2]],
                             Pvalue.Est = as.numeric(pvalue)[1:dim(Bx)[2]],

                             Intercept = as.numeric(thetaEgger)[dim(Bx)[2]+1],
                             StdError.Int = as.numeric(thetaEggerse)[dim(Bx)[2]+1],
                             CILower.Int =  as.numeric(ciLower)[dim(Bx)[2]+1],
                             CIUpper.Int = as.numeric(ciUpper)[dim(Bx)[2]+1],
                             Pvalue.Int = as.numeric(pvalue)[dim(Bx)[2]+1],

                             Alpha = alpha,
                             SNPs = nsnps,

                             RSE = as.numeric(rse.corr),
                             Heter.Stat = c(heter.stat,
                             pvalue.heter.stat)))
            } } else {
                  # not correlated
   summary <- summary(lm(By ~ Bx, weights = Byse^(-2)))

                  thetaEgger <- summary$coef[,1]
     thetaEggerse <- summary$coef[,2]/min(summary$sigma,1)

  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaEgger/thetaEggerse)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaEgger/thetaEggerse), df=dim(Bx)[1]-dim(Bx)[2]) }
                  rse <- summary$sigma

                  heter.stat <- (nsnps - dim(Bx)[2]-1)*(rse^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = nsnps-dim(Bx)[2]-1, lower.tail = F)

                if(distribution == "normal"){
                  ciLower <- ci_normal("l", thetaEgger, thetaEggerse, alpha)
                  ciUpper <- ci_normal("u", thetaEgger, thetaEggerse, alpha)
                } else if (distribution == "t-dist"){
                  ciLower <- ci_t("l", thetaEgger, thetaEggerse, dim(Bx)[1]-dim(Bx)[2]-1, alpha)
                  ciUpper <- ci_t("u", thetaEgger, thetaEggerse, dim(Bx)[1]-dim(Bx)[2]-1, alpha)
                }

                return(new("MVEgger",
                           Model = "random",
                           Orientate = orientAte,
                           Exposure = object@exposure,
                           Outcome = object@outcome,

                           Correlation = object@correlation,

                             Estimate = as.numeric(thetaEgger)[2:(dim(Bx)[2]+1)],
                             StdError.Est = as.numeric(thetaEggerse)[2:(dim(Bx)[2]+1)],
                             CILower.Est =  as.numeric(ciLower)[2:(dim(Bx)[2]+1)],
                             CIUpper.Est = as.numeric(ciUpper)[2:(dim(Bx)[2]+1)],
                             Pvalue.Est = as.numeric(pvalue)[2:(dim(Bx)[2]+1)],

                             Intercept = as.numeric(thetaEgger)[1],
                             StdError.Int = as.numeric(thetaEggerse)[1],
                             CILower.Int =  as.numeric(ciLower)[1],
                             CIUpper.Int = as.numeric(ciUpper)[1],
                             Pvalue.Int = as.numeric(pvalue)[1],

                           Alpha = alpha,
                           SNPs = nsnps,

                           RSE = rse,
                           Heter.Stat = c(heter.stat,
                           pvalue.heter.stat)))

            } } 
 else {

                cat("Distribution must be one of : normal, t-dist. \n")
                cat("See documentation for details. \n")
              }
          }
)
