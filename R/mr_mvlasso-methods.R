#' @docType methods
#' @rdname mr_mvlasso
#' 
setMethod(mr_mvlasso,
  "MRMVInput",
  function(object,
           orientate = 1,
           distribution = "normal",
           alpha = 0.05,
           lambda = numeric(0)){
    
    K = dim(object@betaX)[2]
    if (orientate %in% 1:K) { orientAte = orientate} else {orientAte = 1}
    
    orient = sign(object@betaX)[,orientAte]
    Bx = object@betaX*orient
    By = object@betaY*orient
    Bxse = object@betaXse
    Byse = object@betaYse
    nsnps = dim(Bx)[1]
    
    if(distribution %in% c("normal", "t-dist")){
    
      S = diag(Byse^-2)
      b = S^(1/2) %*% Bx
      Pb = b %*% solve(t(b) %*% b, t(b))
      xlas = (diag(nsnps) - Pb) %*% S^(1/2)
      ylas = (diag(nsnps) - Pb) %*% S^(1/2) %*% By
      
      if (length(lambda) != 0) {
        las_fit = glmnet(xlas, ylas, intercept = FALSE, lambda = lambda)
        las_mod = list(fit = las_fit$beta[, 1], lambda = lambda)
      } else {
        las_fit = glmnet(xlas, ylas, intercept = FALSE)
        lamseq = sort(las_fit$lambda)
        lamlen = length(lamseq)
        rse = sapply(1:lamlen, function(j){
          av = which(las_fit$beta[, (lamlen - j + 1)] == 0)
          mod = lm(S[av, av]^(1/2) %*% By[av] ~ S[av, av]^(1/2) %*% Bx[av, ] - 1)
          c(sqrt(t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)), length(av))
        })
        rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen-1)]
        het = which(rse[1, 2:lamlen] > 1 & rse_inc > ((qchisq(0.95, 1) / rse[2, 2:lamlen])))
        if (length(het) == 0){
          lam_pos = lamlen
        } else {
          lam_pos = min(het)
        }
        num_valid = rse[2, ]
        min_lam_pos = min(which(num_valid > K))
        if (lam_pos < min_lam_pos){lam_pos = min_lam_pos}
        las_mod = list(fit = las_fit$beta[, (lamlen - lam_pos + 1)], lambda = lamseq[lam_pos])
      }
      a = las_mod$fit
      e = By - a
      est = solve(t(Bx) %*% S %*% Bx, t(Bx) %*% S %*% e)
      v = which(a == 0)
      if (length(v) > K){
        post_mod = summary(lm(By[v] ~ Bx[v, ] - 1, weights = Byse[v]^-2))
        post_est = post_mod$coef[, 1]
        post_se = post_mod$coef[, 2] / min(post_mod$sigma, 1)
      } else if (length(v) == K) {
        post_est = rep(NA, K)
        post_se = rep(NA, K)
        cat("Specified value of lambda results in the same number of valid instruments as risk factors. Post-lasso method cannot be performed.")
      } else {
        post_est = rep(NA, K)
        post_se = rep(NA, K)
        cat("Specified value of lambda results in fewer valid instruments than risk factors. Post-lasso method cannot be performed.")
      }
      if(distribution == "normal"){
        ciLower <- ci_normal("l", post_est, post_se, alpha)
        ciUpper <- ci_normal("u", post_est, post_se, alpha)
      } else if (distribution == "t-dist"){
        ciLower <- ci_t("l", post_est, post_se, length(v)-K, alpha)
        ciUpper <- ci_t("u", post_est, post_se, length(v)-K, alpha)
      }
      
      if (distribution == "normal") { pvalue = 2*pnorm(-abs(post_est/post_se)) }
      if (distribution == "t-dist") { pvalue = 2*pt(-abs(post_est/post_se), df=length(v)-K) }
    
      return(new("MVLasso",
                 Orientate = orientate,
                 Exposure = object@exposure,
                 Outcome = object@outcome,
            
                 Estimate = as.numeric(post_est),
                 StdError = as.numeric(post_se),
                 CILower =  as.numeric(ciLower),
                 CIUpper = as.numeric(ciUpper),
                 
                 Alpha = alpha,
                 Pvalue = as.numeric(pvalue),
                 SNPs = nsnps,
                 
                 RegEstimate = as.numeric(est),
                 RegIntercept = as.numeric(a),
                 Valid = length(v),
                 ValidSNPs= as.character(object$snps[v]),
                 Lambda = las_mod$lambda))
    } else {
      cat("Distribution must be one of : normal, t-dist. \n")
      cat("See documentation for details. \n")
    }
  }
)