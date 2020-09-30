#' @docType methods
#' @rdname mr_mvmedian

setMethod(mr_mvmedian,
          "MRMVInput",
          function(object,
                   distribution = "normal",
                   alpha = 0.05, 
                   iterations = 10000,
                   seed = 314159265){
            
if( exists(".Random.seed") ) {
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
 }


if (!is.na(seed)) { set.seed(seed) }


            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            
            nsnps <- dim(Bx)[1]

              if(distribution %in% c("normal", "t-dist")){


            qr_mod = rq(By ~ Bx - 1, weights = Byse^-2)

              #Boot
              if (iterations > 0){
                est = sapply(1:iterations, function(i){
                  p = length(By)
                  k = dim(Bx)[2]
                  Bxboot = sapply(1:k, function(j){rnorm(p, Bx[, j], Bxse[, j])})
                  Byboot = rnorm(p, By, Byse)
                  rq(Byboot ~ Bxboot - 1, weights = Byse^-2)$coefficients
                })
                se = apply(est, 1, sd)
              }
              else {se = NA}
              thetaMed = qr_mod$coefficients
              thetaMedse = se
              
              if(distribution == "normal"){
                ciLower <- ci_normal("l", thetaMed, thetaMedse, alpha)
                ciUpper <- ci_normal("u", thetaMed, thetaMedse, alpha)
              } else if (distribution == "t-dist"){
                ciLower <- ci_t("l", thetaMed, thetaMedse, nsnps - dim(Bx)[2], alpha)
                ciUpper <- ci_t("u", thetaMed, thetaMedse, nsnps - dim(Bx)[2], alpha)
              }
              
              if (distribution == "normal") { pvalue <- 2*pnorm(-abs(thetaMed/thetaMedse)) }
              if (distribution == "t-dist") { pvalue <- 2*pt(-abs(thetaMed/thetaMedse), df=nsnps-dim(Bx)[2]) }
              return(new("MVMedian",
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         
                         Estimate = as.numeric(thetaMed),
                         StdError = as.numeric(thetaMedse),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),

                         Alpha = alpha,
                         Pvalue = as.numeric(pvalue),
                         SNPs = nsnps))
 }
 else {

                cat("Distribution must be one of : normal, t-dist. \n")
                cat("See documentation for details. \n")
              }
          }
)


