# Required Functions
#' Mode-based estimate (Hartwig) estimation function
#'
#' @description Internal function for calculating mode-based estimate.
#'
#' @param BetaIV.in Ratio causal estimates for each genetic variant considered individually.
#' @param seBetaIV.in Standard errors of ratio causal estimates.
#' @param phi.in Bandwidth multiplication factor.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Mode-based estimate.
#'
#' @examples mbe_est(BetaIV.in = chdlodds/ldlc, seBetaIV.in = abs(chdloddsse/ldlc), phi.in = 1)
#'
#' @export

mbe_est <- function(BetaIV.in, seBetaIV.in, phi.in) { 
    s <- 0.9*(min(sd(BetaIV.in), mad(BetaIV.in)))/length(BetaIV.in)^(1/5)
    
    #Standardised weights
    s.weights <- seBetaIV.in^-2/sum(seBetaIV.in^-2)
   
    #Define the actual bandwidth
    h <- s*phi.in
    #Compute the smoothed empirical density function
    densityIV <- density(BetaIV.in, weights=s.weights, bw=h)
    #Extract the point with the highest density as the point estimate 
    beta_est <- densityIV$x[densityIV$y==max(densityIV$y)]
    return(beta_est) }

#' Mode-based estimate (Hartwig) bootstrap function
#'
#' @description Internal function for calculating standard error of mode-based estimate.
#'
#' @param BetaIV.in Ratio causal estimates for each genetic variant considered individually.
#' @param seBetaIV.in Standard errors of ratio causal estimates.
#' @param iterations.in Number of bootstrap iterations.
#' @param phi.in Bandwidth multiplication factor.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Bootstrapped mode-based estimates.
#'
#' @examples mbe_boot(BetaIV.in = chdlodds/ldlc, seBetaIV.in = abs(chdloddsse/ldlc),
#'    weighting.in = "simple", iterations.in = 100)
#'
#' @export

mbe_boot <- function(BetaIV.in, seBetaIV.in, weighting.in, iterations.in, phi.in) {
    beta.boot <- NULL
    
    for(i in 1:iterations.in) {
      BetaIV.boot      <- rnorm(length(BetaIV.in), mean=BetaIV.in, sd=seBetaIV.in)
  if (weighting.in=="weighted")   { beta.boot[i] <- mbe_est(BetaIV.in=BetaIV.boot, seBetaIV.in=seBetaIV.in, phi.in=phi.in) }
  if (weighting.in=="unweighted") { beta.boot[i] <- mbe_est(BetaIV.in=BetaIV.boot, seBetaIV.in=rep(1, length(BetaIV.in)), phi.in=phi.in) }
      }

    return(beta.boot)
  }



#' @docType methods
#' @rdname mr_mbe

setMethod("mr_mbe",
          "MRInput",
          function(object,
                   weighting = "weighted",
                   stderror = "delta",
                   phi = 1,
                   seed = 314159265,
                   iterations = 10000,
                   distribution = "normal",
                   alpha = 0.05){

if( exists(".Random.seed") ) {
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
 }

if (!is.na(seed)) { set.seed(seed) }

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse

            nsnps <- length(Bx)

       if(weighting %in% c("weighted", "unweighted") & stderror %in% c("simple", "delta")){

  BetaIV   <- By/Bx    
  #SEs of ratio estimates
  if (stderror=="simple") { seBetaIV <- Byse/abs(Bx) }
  if (stderror=="delta")  { seBetaIV <- sqrt(Byse^2/abs(Bx)^2 + By^2*Bxse^2/Bx^4) }

  if (weighting=="weighted")   { betaMBE <- mbe_est(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, phi.in = phi) }
  if (weighting=="unweighted") { betaMBE <- mbe_est(BetaIV.in=BetaIV, seBetaIV.in=rep(1, length(BetaIV)), phi.in = phi) }

  betaMBE.boot <- mbe_boot(BetaIV.in=BetaIV, seBetaIV.in=seBetaIV, weighting.in = weighting, iterations.in = iterations, phi.in = phi)
  seMBE <- mad(betaMBE.boot)
  
  if (distribution == "normal") { pvalue <- 2*pnorm(-abs(betaMBE/seMBE)) }
  if (distribution == "t-dist") { pvalue <- 2*pt(-abs(betaMBE/seMBE), df=length(Bx)-1) }


                if(distribution == "normal"){
                  ciLower <- ci_normal("l", betaMBE, seMBE, alpha)
                  ciUpper <- ci_normal("u", betaMBE, seMBE, alpha)
                } else if (distribution == "t-dist"){
                  ciLower <- ci_t("l", betaMBE, seMBE, length(Bx) - 1, alpha)
                  ciUpper <- ci_t("u", betaMBE, seMBE, length(Bx) - 1, alpha)
                }


                  return(new("MRMBE",
                             Exposure = object@exposure,
                             Outcome = object@outcome,
                             Weighting = as.character(weighting),
                             StdErr = as.character(stderror),
                             Phi = as.numeric(phi),

                             Estimate = as.numeric(betaMBE),
                             StdError = as.numeric(seMBE),
                             CILower =  as.numeric(ciLower),
                             CIUpper = as.numeric(ciUpper),

                             SNPs = nsnps,
                             Pvalue = as.numeric(pvalue),
                             Alpha = alpha))
            } else {
                cat("Weighting must be one of: weighted, unweighted. \n")
                cat("Stderror must be one of : simple, delta. \n")
                cat("See documentation for details. \n")
              }
          }
)
