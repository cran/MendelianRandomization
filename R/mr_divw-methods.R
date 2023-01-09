# Calculating the debiased inverse-variance weighted estimate (and relevant statistics) using the input data

#' @docType methods
#' @rdname mr_divw

setMethod("mr_divw",
          signature = c(object = "MRInput"),
          function(object, over.dispersion = TRUE, alpha = 0.05, diagnostics=FALSE){

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            SE.ratio<-object@betaXse/object@betaYse

            beta_dIVW<-sum(By*Bx/Byse^2)/sum((Bx^2-Bxse^2)/Byse^2)
            mu<-Bx/Bxse
            condition<-(mean(mu^2)-1)*sqrt(length(Bx))
            tau.square<-ifelse(over.dispersion==FALSE,0,
                               max(0,sum(((By-beta_dIVW*Bx)^2-Byse^2-beta_dIVW^2*Bxse^2)/Byse^2)/sum(Byse^(-2))))
            V1<-sum((SE.ratio^2*mu^2+beta_dIVW^2*SE.ratio^4*(mu^2+1)+tau.square*SE.ratio^2/Byse^2*mu^2))
            V2<-sum((SE.ratio^2*(mu^2-1)))
            se_dIVW<-sqrt(V1/V2^2)
            c_alpha<-qnorm(alpha/2,lower.tail = FALSE)
            ciLower<-beta_dIVW-c_alpha*se_dIVW
            ciUpper<-beta_dIVW+c_alpha*se_dIVW
            pval <- 2*pnorm(-abs(beta_dIVW/se_dIVW))

            if(diagnostics){
              t<-(By-beta_dIVW*Bx)/sqrt(Byse^2+tau.square+beta_dIVW^2*Bxse^2)
              qqnorm(t)
              qqline(t)
            }

            return(new("DIVW",
                       Over.dispersion = over.dispersion,
                       Exposure = object@exposure,
                       Outcome = object@outcome,

                       Estimate = beta_dIVW,
                       StdError = se_dIVW,
                       CILower = ciLower,
                       CIUpper = ciUpper,
                       Alpha = alpha,

                       Pvalue = as.numeric(pval),
                       SNPs = length(By),
                       Condition = condition)
            )

          }
)
