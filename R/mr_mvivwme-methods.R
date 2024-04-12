# Extra Function

#' Sampling from multivariate normal distribution
#'
#' @description Sampling from multivariate normal distribution
#'
#' @param n Number of draws.
#' @param mu Mean vector.
#' @param Sigma Covariance matrix.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Random vector of length n.
#'
#' @examples mv_norm(1, 0.2, matrix(0.1))
#' 
#' @export

mv_norm = function(n, mu, Sigma){
  d = dim(Sigma)[1]
  if (length(mu) != d){
    stop('mu and Sigma must be the same dimension.')
  }
  A = chol(Sigma)
  Z = sapply(1:d, function(j){rnorm(n, 0, 1)})
  M = matrix(rep(mu, n), nrow = n, byrow = TRUE)
  X = drop(M) + Z %*% A
}


#' @docType methods
#' @rdname mr_mvivwme

setMethod("mr_mvivwme",
          "MRMVInput",
          function(object,
                   model = "default",
                   correl = FALSE,
                   correl.x = NULL,
                   distribution = "normal",
                   alpha = 0.05,
				   max_iter = 100,
				   no_ini = 1,
				   seed = 20201201,
				   ...){
				   
			if( exists(".Random.seed") ) {
			 old <- .Random.seed
			 on.exit( { .Random.seed <<- old } )
			}

			if (!is.na(seed)) { set.seed(seed) }



            bxhat = object@betaX
            byhat = object@betaY
            sebx = object@betaXse
            seby = object@betaYse
            rho = object@correlation # LD matrix
			
			p = length(byhat) # nsnps
			K = dim(bxhat)[2]
			S = diag(seby^-2)


            # if genetic variant correlation matrix is not specified, assume they are uncorrelated
			
             if(!is.na(sum(rho))){
              S <- diag(1/seby) %*% rho %*% diag(1/seby)
             } else {
			  if(correl){
			   cat("Genetic variant correlation matrix not given.")
			  }
			 }
			
			if(is.null(correl.x)){
				SigX = lapply(1:p, function(j){diag(sebx[j, ]^2, length(sebx[j, ]))}) 
			} else {
			SigX = lapply(1:p, function(j){
				S1 = diag(sebx[j, ])
				S1 %*% correl.x %*% S1
			})
			}
			
			
            
            if(model == "default"){
              if(p < 4) {
                model <- "fixed"
              } else {
                model <- "random"
              }
            }

			  
			if(!(model %in% c("random", "fixed") & distribution %in% c("normal", "t-dist"))){  
                cat("Model type must be one of: default, random, fixed. \n")
                cat("Distribution must be one of : normal, t-dist. \n")
                cat("See documentation for details. \n")
				return(NULL)
				}
  
			  
			  
  l = matrix(nrow = max_iter, ncol = no_ini)
  thest = matrix(nrow = K, ncol = no_ini)
  for (k in 1:no_ini){
    bxtilde = t(sapply(1:p, function(j){mv_norm(1, bxhat[j, ], SigX[[j]])}))
    for (i in 1:100){
      thest[, k] = solve(t(bxtilde) %*% S %*% bxtilde, t(bxtilde) %*% S %*% byhat)
      l[i, k] = -0.5 * sum(sapply(1:p, function(j){
        (byhat[j] - t(bxhat[j, ]) %*% thest[, k])^2 / (seby[j]^2 + t(thest[, k]) %*% SigX[[j]] %*% thest[, k])
      }))
      for (j in 1:p){
        bxtilde[j, ] = t(solve(thest[, k] %*% t(thest[, k]) / seby[j]^2 + solve(SigX[[j]]), byhat[j] * thest[, k] / seby[j]^2 + solve(SigX[[j]], bxhat[j, ])))
      }
      if (i > 1){
        if (abs(l[i] - l[(i-1)]) < 1e-4) {break}
      }
    }
  }
  k0 = which.max(apply(as.matrix(l[is.na(l[, 1]) == F, ]), 2, max))
  th = thest[, k0]
  
  v = sapply(1:p, function(j){seby[j]^2 + t(th) %*% SigX[[j]] %*% th})
  e = sapply(1:p, function(j){byhat[j] - bxhat[j, ] %*% th})
  t = sapply(1:p, function(j){e[j] / sqrt(v[j])})
  
  dt = sapply(1:p, function(j){(-v[j] * bxhat[j, ] - e[j] * SigX[[j]] %*% th) / v[j]^(3/2)})
  B = dt %*% t(dt)
  
  dt2 = vector(length = p, mode = "list")
  for (j in 1:p){
    dt2[[j]] = matrix(nrow = K, ncol = K)
    S = SigX[[j]] %*% th
    for (k in 1:K){
      for (l in 1:K){
        dt2[[j]][k, l] = v[j]^(-3/2) * (-2 * S[l] * bxhat[j, k] + S[k] * bxhat[j, l] - e[j] * SigX[[j]][k, l]) +
          3 * v[j]^(-5/2) *(v[j] * bxhat[j, k] + e[j] * S[k]) * S[l]
      }
    }
  }
  a = Reduce('+', lapply(1:p, function(j){c(t[j]) * dt2[[j]]}))
  A = (dt %*% t(dt) + a)
  
  Var = solve(A, B) %*% t(solve(A))

  thetaIVW <- th
  
					resids <- byhat - bxhat %*% th
					rmse <- as.numeric(sqrt(t(resids) %*% diag(seby^-2) %*% resids / (p - K)))

                    if(model == "random") {
                      thetaIVWse <- sqrt(diag(Var))*max(rmse, 1)
                    } else if (model == "fixed"){
                      thetaIVWse <- sqrt(diag(Var))
                    }

                  if(distribution == "normal"){
                    ciLower <- ci_normal("l", thetaIVW, thetaIVWse, alpha)
                    ciUpper <- ci_normal("u", thetaIVW, thetaIVWse, alpha)
                  } else if (distribution == "t-dist"){
                    ciLower <- ci_t("l", thetaIVW, thetaIVWse, p - K, alpha)
                    ciUpper <- ci_t("u", thetaIVW, thetaIVWse, p - K, alpha)
                  }

                  heter.stat <- (p - K)*(rmse^2)
                  pvalue.heter.stat <- pchisq(heter.stat, df = p - K, lower.tail = F)
				  
  if (distribution == "normal") { pvalue <- 2*pnorm(abs(thetaIVW/thetaIVWse), lower.tail = FALSE) }
  if (distribution == "t-dist") { pvalue <- 2*pt(abs(thetaIVW/thetaIVWse), df = p - K, lower.tail = FALSE) }
        
                  return(new("MVIVWME",
                             Model = model,
                             Exposure = object@exposure,
                             Outcome = object@outcome,
							 
							 Correlation = object@correlation,

                             Estimate = as.numeric(thetaIVW),
                             StdError = as.numeric(thetaIVWse),
                             CILower =  as.numeric(ciLower),
                             CIUpper = as.numeric(ciUpper),

                             SNPs = p,
                             Pvalue = as.numeric(pvalue),

                             Alpha = alpha,

							 RSE = rmse,
                             Heter.Stat = c(heter.stat,
                             pvalue.heter.stat)
							 ))

}
)
