# Required Functions
#' Generate bootstrap samples for the bootstrapping Fieller's confidence interval of the penalized inverse-variance weighted (pIVW) method
#'
#' @description Internal function of the penalized inverse-variance weighted (pIVW) method, which generates bootstrap samples for the bootstrapping Fieller's confidence interval.
#'
#' @param object An \code{MRInput} object.
#' @param beta_hat The causal effect estimate.
#' @param tau2 The estimated variance of the horizontal pleiotropy.
#' @param lambda The penalty parameter in the pIVW estimator. By default, \code{lambda=1}.
#' @param n_boot The sample size of the bootstrap samples. By default, \code{n_boot=1000}.
#' @param seed_boot The seed for random sampling in the bootstrap method. By default, \code{seed_boot=1}.
#'
#' @keywords internal
#'
#' @return \item{z_b}{A vector containing the bootstrap samples for the bootstrapping Fieller's confidence interval.}
#'
#' @export

BF_dist = function(object,beta_hat=0,tau2=0,lambda=1,n_boot=1000,seed_boot=1){

            if( exists(".Random.seed") ) {
              old <- .Random.seed
              on.exit( { .Random.seed <<- old } )
            }


  i = 0
  set.seed(seed_boot)
  z_b = numeric()

  while(i < n_boot){
    Bx = object@betaX
    Byse  = object@betaYse
    Bxse  = object@betaXse
    alpha = rnorm(length(Bx),0, sqrt(tau2))
    By_hat = rnorm(length(Bx),Bx*beta_hat+alpha, Byse)
    Bx_hat = rnorm(length(Bx),Bx, Bxse)

    v1 = sum(Byse^(-4)*(By_hat^2*Bx_hat^2-(By_hat^2-Byse^2-tau2)*(Bx_hat^2-Bxse^2)))
    v2 = sum(Byse^(-4)*(4*(Bx_hat^2-Bxse^2)*Bxse^2 + 2*Bxse^4))
    v12= sum(2*Byse^(-4)*By_hat*Bx_hat*Bxse^2)
    mu1= sum(Byse^(-2)*By_hat*Bx_hat)
    mu2= sum(Byse^(-2)*(Bx_hat^2-Bxse^2))
    mu2p = mu2/2+sign(mu2)*sqrt(mu2^2/4+lambda*v2)
    mu1p = mu1+(v12/v2)*(mu2p-mu2)
    w = mu2p/(2*mu2p-mu2)
    v1p = v1
    v2p = w^2*v2
    v12p = w*v12

    temp = (mu1p-beta_hat*mu2p)^2/(v1p-2*beta_hat*v12p+beta_hat^2*v2p)
    if(temp<0){
      next
    }else{
      z_b = c(temp,z_b)
      i = i+1
    }
  }
  z_b = sort(z_b)
  return(z_b)
}



# Calculating the penalized inverse-variance weighted estimate (and relevant statistics) using the input data

#' @docType methods
#' @rdname mr_pivw

setMethod("mr_pivw",
          signature = c(object = "MRInput"),
          function(object, lambda=1, over.dispersion=TRUE, delta=0, sel.pval=NULL, Boot.Fieller=TRUE, alpha=0.05){

            Bx = object@betaX
            Bxse  =object@betaXse
            p = length(Bx)

            if(lambda<0){
              cat("\'lambda\' cannot be smaller than zero.","\n")
              return()
            }

            if(delta<0){
              cat("\'delta\' cannot be smaller than zero.","\n")
              return()
            }

            if(alpha<=0){
              alpha = 0.05
              cat("\'alpha\' provided is less than or equal to zero. \'alpha\ is set to be 0.05.","\n")
            }

            if(delta>0){
              if(is.null(sel.pval)){
                delta=0
                cat("\'sel.pval\' is not provided. \'delta\' is set to be zero and no IV selection conducted.","\n")
              }else if(length(sel.pval)!=p){
                delta=0
                cat("The length of \'sel.pval\' doesn't match the number of snps in the MRInput object. \'delta\' is set to be zero and no IV selection conducted.","\n")
              }
            }



            if(delta==0){
              kappa = mean(Bx^2/Bxse^2)-1
              eta = kappa*sqrt(p)
            }else{
              sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
              if(sum(sel)==0){
                cat("No snps is kept in the anaysis after IV selection. Please try a smaller \'delta\' value.","\n")
                return()
              }
              kappa = mean(Bx[sel]^2/Bxse[sel]^2)-1
              eta = kappa*sqrt(sum(sel))
              sel.z = qnorm(sel.pval/2,lower.tail = FALSE)
              q = pnorm(sel.z-delta,lower.tail = TRUE)+pnorm(-sel.z-delta,lower.tail = TRUE)
              psi2 = sum(((Bx/Bxse)^4-6*(Bx/Bxse)^2+3)*q*(1-q))/sum(sel)
              eta = eta/max(1,sqrt(psi2))
            }

            if(delta!=0 & over.dispersion!=TRUE){
              sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
              object@betaY = object@betaY[sel]
              object@betaX = object@betaX[sel]
              object@betaYse = object@betaYse[sel]
              object@betaXse = object@betaXse[sel]
            }

            By = object@betaY
            Bx = object@betaX
            Byse  = object@betaYse
            Bxse  = object@betaXse

            v2 = sum(Byse^(-4)*(4*(Bx^2-Bxse^2)*Bxse^2 + 2*Bxse^4))
            v12 = sum(2*Byse^(-4)*By*Bx*Bxse^2)
            mu1 = sum(Byse^(-2)*By*Bx)
            mu2 = sum(Byse^(-2)*(Bx^2-Bxse^2))
            mu2p = mu2/2 + sign(mu2)*sqrt(mu2^2/4+lambda*v2)
            mu1p = mu1 + (v12/v2)*(mu2p-mu2)
            beta_pIVW = mu1p/mu2p

            if(over.dispersion==TRUE){
              tau2 = sum(((By-beta_pIVW*Bx)^2 - Byse^2 -beta_pIVW^2*Bxse^2) * Byse^(-2))/sum(Byse^(-2))
              if (tau2 < 0){tau2 = 0}
            }else{
              tau2 = 0
            }

            if(delta!=0 & over.dispersion==TRUE){
              sel = sel.pval<2*pnorm(delta,lower.tail = FALSE)
              object@betaY = object@betaY[sel]
              object@betaX = object@betaX[sel]
              object@betaYse = object@betaYse[sel]
              object@betaXse = object@betaXse[sel]

              By = object@betaY
              Bx = object@betaX
              Byse = object@betaYse
              Bxse = object@betaXse

              v2 = sum(Byse^(-4)*(4*(Bx^2-Bxse^2)*Bxse^2 + 2*Bxse^4))
              v12 = sum(2*Byse^(-4)*By*Bx*Bxse^2)
              mu1 = sum(Byse^(-2)*By*Bx)
              mu2 = sum(Byse^(-2)*(Bx^2-Bxse^2))
              mu2p = mu2/2 + sign(mu2)*sqrt(mu2^2/4+lambda*v2)
              mu1p = mu1 + (v12/v2)*(mu2p-mu2)
              beta_pIVW = mu1p/mu2p
            }
            se_pIVW = 1/abs(mu2p)*sqrt(sum((Bx^2/Byse^2)*(1+tau2*Byse^(-2))
                                           + beta_pIVW^2*(Bxse^2/Byse^4)*(Bx^2+Bxse^2)))


            if(Boot.Fieller==TRUE){

              z_b = BF_dist(object,beta_pIVW,tau2,lambda)
              v1p = sum(Byse^(-4)*(By^2*Bx^2-(By^2-Byse^2-tau2)*(Bx^2-Bxse^2)))
              w = mu2p/(2*mu2p-mu2)
              v2p = w^2*v2
              v12p = w*v12
              z0 = mu1p^2/v1p
              pval = sum(z0<z_b)/length(z_b)

              qt = z_b[round(length(z_b)*(1-alpha))]
              A = mu2p^2-qt*v2p
              B = 2*(qt*v12p - mu1p*mu2p)
              C = mu1p^2 -qt*v1p
              D = B^2-4*A*C

              if(A>0){
                r1 = (-B-sqrt(D))/2/A
                r2 = (-B+sqrt(D))/2/A
                CI = cbind(r1,r2)
                CI_ll = r1
                CI_ul = r2
              }else if (D>0){
                r1 = (-B-sqrt(D))/2/A
                r2 = (-B+sqrt(D))/2/A
                CI1 = c(-Inf,r1)
                CI2 = c(r2,Inf)
                CI = rbind(CI1,CI2)
                CI_ll = c(-Inf,r2)
                CI_ul = c(r1,Inf)
              }else{
                CI = cbind(-Inf,Inf)
                CI_ll = -Inf
                CI_ul = Inf
              }

            }else{
              pval = 2*pnorm(abs(beta_pIVW),0,se_pIVW,lower.tail = FALSE)
              CI_ll = beta_pIVW+qnorm(alpha/2)*se_pIVW
              CI_ul = beta_pIVW-qnorm(alpha/2)*se_pIVW
              CI = cbind(CI_ll,CI_ul)
            }


            return(new("PIVW",
                       Over.dispersion = over.dispersion,
                       Boot.Fieller = Boot.Fieller,
                       Lambda = lambda,
                       Delta = delta,
                       Exposure = object@exposure,
                       Outcome = object@outcome,

                       Estimate = beta_pIVW,
                       StdError = se_pIVW,
                       CILower = CI_ll,
                       CIUpper = CI_ul,
                       Alpha = alpha,

                       Pvalue = as.numeric(pval),
                       Tau2 = tau2,
                       SNPs = length(By),
                       Condition = eta)
            )




          }
)

