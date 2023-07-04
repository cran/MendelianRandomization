#' @docType methods
#' @rdname mr_mvcML

setMethod("mr_mvcML",
          "MRMVInput",
          function(object,
                   n,
                   DP = TRUE,
                   rho_mat = diag(ncol(object@betaX)+1),
                   K_vec = 0:(ceiling(nrow(object@betaX)/2)),
                   random_start = 0,
                   num_pert = 200,
                   min_theta_range = -0.5,
                   max_theta_range = 0.5,
                   maxit = 100,
                   alpha = 0.05, 
                   seed = 314159265){
            
if( exists(".Random.seed") ) {
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
 }


if (!is.na(seed)) { set.seed(seed) }


            Bx = object@betaX
            By = as.matrix(object@betaY)
            Bxse = object@betaXse
            Byse = object@betaYse
            
            nsnps <- dim(Bx)[1]
            Sig_inv_l = invcov_mvmr(Bxse, Byse, rho_mat)

            # run MVMRcML without data perturbation, i.e. only with BIC selection #
            if(!DP){
              res = MVmr_cML(b_exp = Bx, b_out = By, se_bx = Bxse, Sig_inv_l = Sig_inv_l,  
              n = n, K_vec = K_vec, random_start = random_start + 1, min_theta_range = min_theta_range, max_theta_range = max_theta_range, 
              maxit = maxit, thres = 1e-4)
              theta = res$BIC_theta
              se_theta = MVcML_SdTheta(b_exp = Bx, b_out = By, Sig_inv_l = Sig_inv_l,
                theta = theta, zero_ind = setdiff(1:nsnps,res$BIC_invalid))
              pvalue = 2*pnorm(-abs(theta/se_theta))
              ciLower = theta - qnorm(1-alpha/2)*se_theta
              ciUpper = theta + qnorm(1-alpha/2)*se_theta
              if(res$Converge == 1){
                stop("MVMRcML-BIC failed to converge, please try using multiple random starting points (random_start) or increasing number of iterations (maxit).")
              }
              return(new("MVMRcML",
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         DP = FALSE,
                         Estimate = as.numeric(theta),
                         StdError = as.numeric(se_theta),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),
                         Alpha = alpha,
                         Pvalue = as.numeric(pvalue),
                         BIC_invalid = as.numeric(res$BIC_invalid),
                         K_hat = as.numeric(res$Khat),
                         SNPs = nsnps))

            }

            ## run MVMRcML with data perturbation 
            res = MVmr_cML_DP(b_exp = Bx, b_out = By, se_bx = Bxse, Sig_inv_l = Sig_inv_l,
              n = n, K_vec = K_vec, random_start = random_start + 1, num_pert = num_pert, min_theta_range = min_theta_range, max_theta_range = max_theta_range,
              maxit = maxit, thres = 1e-4)
            theta = res$BIC_DP_theta
            se_theta = res$BIC_DP_se
            pvalue = 2*pnorm(-abs(theta/se_theta))
            ciLower = theta - qnorm(1-alpha/2)*se_theta
            ciUpper = theta + qnorm(1-alpha/2)*se_theta
            if(res$eff_DP_B < 100){
              warning("The number of data perturbation is too small, results may not be stable. num_pert >= 100 is recommended.")
            }
            if(res$eff_DP_B < num_pert/2){
              warning("Less than half of the perturbations converged, results may not be stable. You may try increasing the number of random starting points (random_start).")
            }

              return(new("MVMRcML",
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         DP = TRUE,
                         Estimate = as.numeric(theta),
                         StdError = as.numeric(se_theta),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),
                         Alpha = alpha,
                         Pvalue = as.numeric(pvalue),
                         BIC_invalid = as.numeric(res$BIC_invalid),
                         K_hat = as.numeric(res$DP_ninvalid),
                         eff_DP_B = res$eff_DP_B,
                         SNPs = nsnps))
            
      } 
)


#' Profile likelihood of valid IVs
#'
#' calculate the profile likelihood of valid IVs up to a constant
#'
#' @keywords internal
#'
pl <- function(x,b_exp_v,b_out_v,Sig_inv_v){
  k = ncol(b_exp_v)
  m_valid = length(b_out_v)
  pll = 0
  for(i in 1:m_valid){
    W = Sig_inv_v[[i]]
    b_exp_i = b_exp_v[i,]
    b_out_i = b_out_v[i]
    beta = c(b_exp_i,b_out_i)
    B = W[-(k+1),] %*% beta + c(W[(k+1),] %*% beta) * x
    A = W[-(k+1),-(k+1)] + W[-(k+1),k+1] %*% t(x) + x %*% t(W[k+1,-(k+1)]) + W[k+1,k+1] * x %*% t(x)
    bhat_xi = solve(A) %*% B
    b_i = c(bhat_xi,t(bhat_xi) %*% x)
    pll = pll  -1/2 * t(b_i) %*% W %*% b_i + t(beta) %*% W %*% b_i
  }
  return(-pll)
}

#' Generate the list of inverse of covariance matrices used in \code{MVMR-cML-DP}
#'
#' @param se_bx A m*L matrix of standard errors of SNP-exposure association
#' @param se_by A vector of standard errors of SNP-outcome association
#' @param rho_mat The correlation matrix among the L exposures and the outcome
#' @return A list of inverse of covariance matrices with respect to each genetic variant, retaining the ordering in \code{se_bx}
#'
#' @export
#'
invcov_mvmr <- function(se_bx, se_by, rho_mat){
  Sig_l = Sig_inv_l = list()
  m = nrow(se_bx)
  for(i in 1:m){
    Sig = rho_mat*crossprod(t(c(se_bx[i,],se_by[i])))
    Sig_l[[i]] = Sig
    Sig_inv_l[[i]] = solve(Sig)
  }
  return(Sig_inv_l)
}

#' Standard error estimate for MVMR-cML-BIC
#'
#' This is based on the profile likelihood of the set of valid IVs, which is not robust to uncertainty in model selection.
#'
#' @param b_exp A matrix of SNP effects on the exposure variable.
#' @param b_out A vector of SNP effects on the outcome variable.
#' @param Sig_inv_l A list of inverse of covariance matrix.
#' @param theta A vector of final estimates of causal effect of each exposure by MVMR-cML-BIC obtained from \code{MVmr_cML_DP}.
#' @param zero_ind A vector of the index of valid IVs.
#' @param r_vec A vector of estimated horizontal pleiotropic effects.
#'
#'
#' @import numDeriv
#' @return A vector
#'
#' @export
#'
MVcML_SdTheta <- function(b_exp,b_out,Sig_inv_l,theta,zero_ind,r_vec=NULL){
  if(!is.null(r_vec)){
      zero_ind = which(r_vec==0)
    }
  H = hessian(pl,x=theta,b_exp_v = b_exp[zero_ind,,drop=FALSE], b_out_v=b_out[zero_ind],Sig_inv_v =Sig_inv_l[zero_ind])
  se_theta = sqrt(diag(solve(H)))
  return(se_theta)
}


