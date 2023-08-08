# Extra Function

#' Calculates conditional F-statistic for each risk factor using summarized data
#'
#' @description Internal function for calculating conditional F-statistics for the output \code{condFstat}.
#'
#' @param .Bx Beta-coefficient for genetic associations with the risk factors.
#' @param .Bxse Standard error of genetic associations with the risk factors.
#' @param .nx sample sizes used to compute genetic associations with the risk factors.
#' @param .ld genetic variant correlation matrix
#' @param .cor.x Risk factor correlation matrix
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Conditional F-statistic for each exposure
#'
#' @examples condFstat(.Bx = cbind(ldlc, hdlc, trig), .Bxse = cbind(ldlcse, hdlcse, trigse),
#'  .nx = rep(17723,3), .ld=diag(length(ldlc)), .cor.x=diag(3))
#' 
#' @export

condFstat <- function(.Bx, .Bxse, .nx, .ld, .cor.x){
  bx=.Bx; sx=.Bxse; nx=.nx; ld=.ld; cor.x=.cor.x; J=nrow(bx); K=ncol(bx)
  
  # exposures quantities of interest
  ax <- matrix(NA,nrow=J,ncol=K); Ax <- list()
  for (k in 1:K){ax[,k] <- 1/((nx[k]*sx[,k]^2)+bx[,k]^2)}
  for (k in 1:K){Ax[[k]] <- (sqrt(ax[,k])%*%t(sqrt(ax[,k])))*ld}
  Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)
  sqrt.Ax <- function(k){
    evec <- eigen(Ax[[k]])$vectors; eval <- eigen(Ax[[k]])$values
    return((evec%*%diag(sqrt(eval))%*%t(evec)))
  }
  sqrt.Ax <- lapply(1:K,sqrt.Ax)
  
  SigX <- function(k,l){
    solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))*(cor.x[k,l]-(as.numeric(t(Bx[,k])%*%solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))%*%Bx[,l])))*(1/(sqrt(nx[k]-J+1)*sqrt(nx[l]-J+1)))
  }
  
  SigX2 <- function(m){
    SigX2a <- list()
    for (m1 in 1:K){SigX2a[[m1]] <- SigX(m1,m)}
    SigX2a <- do.call(rbind, SigX2a)
    return(SigX2a)
  }
  
  SigX3 <- list()
  for (m1 in 1:K){SigX3[[m1]] <- SigX2(m1)}
  SigX3 <- do.call(cbind, SigX3)
  # check: SigX3[(2*J+1):(3*J),(1*J+1):(2*J)] == SigX(3,2)
  SigX <- SigX3; rm(SigX2, SigX3)
  gamX_est <- function(k){as.vector(solve(Ax[[k]])%*%Bx[,k])}
  gamX_est <- sapply(1:K,gamX_est)
  
  # compute conditional F-statistics
  if(K>2){
  condF <- function(j){
    SigXX <- SigX[c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))])),c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))]))]
    g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]%*%tet))}
    Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
    tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
    Om.nr <- function(tet){as.matrix(cbind(diag(J),kronecker(t(-tet),diag(J)))%*%SigXX%*%t(cbind(diag(J),kronecker(t(-tet),diag(J)))))}
    Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
    G <- -gamX_est[,-j]
    DQ.nr <- function(tet){2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet))}
    condF <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(J-K+1)
    return(condF)
  }
  condF <- sapply(1:K,condF)
  }
  
  if(K==2){
    condF <- function(j){
      SigXX <- SigX[c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))])),c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))]))]
      g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]*tet))}
      Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
      tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
      Om.nr <- function(tet){as.matrix(cbind(diag(J),kronecker(t(-tet),diag(J)))%*%SigXX%*%t(cbind(diag(J),kronecker(t(-tet),diag(J)))))}
      Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
      G <- -gamX_est[,-j]
      DQ.nr <- function(tet){as.numeric(2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet)))}
      condF <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(J-K+1)
      return(condF)
    }
    condF <- sapply(1:K,condF)
  }
  return(condF)
}

#' @docType methods
#' @rdname mr_mvivw

setMethod("mr_mvivw",
          "MRMVInput",
          function(object,
                   model = "default",
                   robust = FALSE,
                   correl = FALSE,
                   correl.x = NULL,
                   nx = NA,
                   distribution = "normal",
                   alpha = 0.05, ...){

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            rho = object@correlation

            nsnps <- dim(Bx)[1]
            
            if(sum(is.na(nx))>0){nx <- rep(NA,dim(Bx)[2]); 
            } else { nx=nx }
            
            if(length(nx)==1){nx <- rep(nx,dim(Bx)[2])}

            # if risk factor correlation matrix is not specified, assume they are uncorrelated
            if(missing(correl.x)){cor.x <- diag(dim(Bx)[2])
            } else {cor.x <- correl.x}
            
            # if genetic variant correlation matrix is not specified, assume they are uncorrelated
            if(!is.na(sum(rho))){ld <- rho
            } else {ld <- diag(dim(Bx)[1])}
            
            # compute conditional F-statistics if nx is supplied
            if(sum(is.na(nx)) > 0){condF = as.numeric(rep(NA, dim(Bx)[2]))}
            if(sum(is.na(nx)) == 0){
      condF = condFstat(Bx, Bxse, nx, ld, cor.x)
            if(sum(condF<0)>0) { condF = as.numeric(rep(NA, dim(Bx)[2])); cat("Conditional F statistics did not converge to positive values - should the sample sizes be larger?\n") } }
            
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

                cat("Genetic variant correlation matrix not given.")

              } else {

                  robust <- FALSE
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

                             Robust = robust,
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
                             pvalue.heter.stat),
                             CondFstat = condF))

            } } else {
                  # not correlated

if (robust == TRUE) {
   summary <- summary(lmrob(By ~ Bx - 1, weights = Byse^(-2), k.max = 500, maxit.scale = 500, ...)) }

else {

   summary <- summary(lm(By ~ Bx - 1, weights = Byse^(-2), ...))
     }

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

                           Robust = robust,
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
                           pvalue.heter.stat),
                           CondFstat = condF))
                  
            } } 
 else {
                cat("Model type must be one of: default, random, fixed. \n")
                cat("Distribution must be one of : normal, t-dist. \n")
                cat("See documentation for details. \n")
              }
          }
)
