#' @docType methods
#' @rdname mr_conmix

setMethod("mr_conmix",
          "MRInput",
          function(object,
                   psi    = 0,
                   CIMin  = -1,
                   CIMax  = 1,
                   CIStep = 0.001,
                   alpha = 0.05){

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse

            nsnps = length(Bx)

            ratio = By/Bx; ratio.se = abs(Byse/Bx);
   if (psi < 0 | psi == 0) {   psi = 1.5*sd(ratio)  }
            theta = seq(from = CIMin, to = CIMax, by = CIStep)
            iters = length(theta)
lik=NULL
 for (j1 in 1:iters) {
  lik.inc = exp(-(theta[j1]-ratio)^2/2/ratio.se^2)/sqrt(2*pi*ratio.se^2)
  lik.exc = exp(-ratio^2/2/(psi^2+ratio.se^2))/(sqrt(2*pi*(psi^2+ratio.se^2)))
  valid = (lik.inc>lik.exc)*1
  lik[j1] = prod(c(lik.inc[valid==1], lik.exc[valid==0]))
  if (which.max(lik)==length(lik)) { valid.best = valid }
 }
  phi = ifelse(sum(valid.best)<1.5, 1,
    max(sqrt(sum(((ratio[valid.best==1]-weighted.mean(ratio[valid.best==1],
           ratio.se[valid.best==1]^-2))^2*
           ratio.se[valid.best==1]^-2))/(sum(valid.best)-1)), 1))
 loglik = log(lik)

whichin = which(2*loglik>(2*max(loglik)-qchisq(1-alpha, df=1)*phi^2))
   # provides an index of estimate values in the 95% confidence interval
            betaConMix = CIMin+CIStep*(which.max(loglik)-1)
   # modal estimate
            CIRange    = CIMin+CIStep*(whichin-1);
            CILower <- c(min(CIRange), CIRange[which(diff(CIRange)>1.01*CIStep)+1])
            CIUpper <- c(CIRange[which(diff(CIRange)>1.01*CIStep)], max(CIRange))


                  return(new("MRConMix",
                             Exposure = object@exposure,
                             Outcome = object@outcome,
                             Psi = as.numeric(psi),

                             Estimate = as.numeric(betaConMix),
                             CIRange  = as.numeric(CIRange),
                             CILower  = as.numeric(CILower),
                             CIUpper  = as.numeric(CIUpper),

                             CIMin    = as.numeric(CIMin),
                             CIMax    = as.numeric(CIMax),
                             CIStep   = as.numeric(CIStep),

                             SNPs = nsnps,
                             Alpha = alpha))
   }
)
