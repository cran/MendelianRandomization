#' @docType methods
#' @rdname mr_mvgmm

setMethod("mr_mvgmm",
          "MRMVInput",
          function(object, nx, ny, cor.x=NULL, robust=TRUE, alpha=0.05, ...){
            
            bx = object@betaX
            by = object@betaY
            sx = object@betaXse
            sy = object@betaYse
            if(sum(is.na(object@correlation)) != 0){ld = diag(nrow(bx))}
            if(sum(is.na(object@correlation)) == 0){ld = object@correlation}
            p <- nrow(bx); K <- ncol(bx)
            
            if(missing(cor.x)){cor.x <- diag(K); mis.cor.x <- FALSE
            } else {cor.x <- cor.x; mis.cor.x <- TRUE}
            
            if(length(nx)==1){nx <- rep(nx,K)}

            # outcome quantities of interest
            ay <- 1/((ny*sy^2)+by^2)
            Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
            By <- ay*by
            SigY <- solve(Ay)*(1-as.numeric(t(By)%*%solve(Ay)%*%By)); SigY <- SigY*(1/(ny-p+1))
            gamY_est <- as.vector(solve(Ay)%*%By)
            
            # exposures quantities of interest
            ax <- matrix(NA,nrow=p,ncol=K); Ax <- list()
            for (k in 1:K){ax[,k] <- 1/((nx[k]*sx[,k]^2)+bx[,k]^2)}
            for (k in 1:K){Ax[[k]] <- (sqrt(ax[,k])%*%t(sqrt(ax[,k])))*ld}
            Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)
            sqrt.Ax <- function(k){
              evec <- eigen(Ax[[k]])$vectors; eval <- eigen(Ax[[k]])$values
              return((evec%*%diag(sqrt(eval))%*%t(evec)))
            }
            sqrt.Ax <- lapply(1:K,sqrt.Ax)
            
            SigX <- function(k,l){
              solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))*(cor.x[k,l]-(as.numeric(t(Bx[,k])%*%solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))%*%Bx[,l])))*(1/(sqrt(nx[k]-p+1)*sqrt(nx[l]-p+1)))
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
            # check: SigX3[(2*p+1):(3*p),(1*p+1):(2*p)] == SigX(3,2)
            SigX <- SigX3; rm(SigX2, SigX3)
            gamX_est <- function(k){as.vector(solve(Ax[[k]])%*%Bx[,k])}
            gamX_est <- sapply(1:K,gamX_est)
            
            # conditional F-test
            if(K>2){
              condF <- function(j){
                SigXX <- SigX[c(((((j-1)*p)+1):(j*p)),((1:(p*K))[-((((j-1)*p)+1):(j*p))])),c(((((j-1)*p)+1):(j*p)),((1:(p*K))[-((((j-1)*p)+1):(j*p))]))]
                g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]%*%tet))}
                Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
                tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
                Om.nr <- function(tet){as.matrix(cbind(diag(p),kronecker(t(-tet),diag(p)))%*%SigXX%*%t(cbind(diag(p),kronecker(t(-tet),diag(p)))))}
                Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
                G <- -gamX_est[,-j]
                DQ.nr <- function(tet){2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet))}
                condF <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(p-K+1)
                return(condF)
              }
              condF <- sapply(1:K,condF)
            }
            
            if(K==2){
              condF <- function(j){
                SigXX <- SigX[c(((((j-1)*p)+1):(j*p)),((1:(p*K))[-((((j-1)*p)+1):(j*p))])),c(((((j-1)*p)+1):(j*p)),((1:(p*K))[-((((j-1)*p)+1):(j*p))]))]
                g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]*tet))}
                Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
                tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
                Om.nr <- function(tet){as.matrix(cbind(diag(p),kronecker(t(-tet),diag(p)))%*%SigXX%*%t(cbind(diag(p),kronecker(t(-tet),diag(p)))))}
                Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
                G <- -gamX_est[,-j]
                DQ.nr <- function(tet){as.numeric(2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet)))}
                condF <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(p-K+1)
                return(condF)
              }
              condF <- sapply(1:K,condF)
            }
            
            
            if(robust == FALSE){
              # non-robust LIML estimate
              g <- function(tet){as.vector(gamY_est - (gamX_est%*%tet))}
              Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
              tet.gg <- nlminb(rep(0,K),objective=Q.gg)$par
              Om.nr <- function(tet){as.matrix(SigY+(kronecker(t(tet),diag(p))%*%SigX%*%t(kronecker(t(tet),diag(p)))))}
              Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
              G <- -gamX_est
              DQ.nr <- function(tet){2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet))}
              liml.nr <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$par
              var.liml.nr <- as.matrix(solve(t(G)%*%solve(Om.nr(liml.nr))%*%G))
              Q.nr <- Q.nr(liml.nr)
              ciLower <- ci_normal("l", liml.nr, sqrt(diag(var.liml.nr)), alpha)
              ciUpper <- ci_normal("u", liml.nr, sqrt(diag(var.liml.nr)), alpha)
              pvalue.heter.stat <- pchisq(Q.nr, df = p-K, lower.tail = F)
              pvalue <- 2*pnorm(-abs(liml.nr/sqrt(diag(var.liml.nr))))
              if(sum(condF<0)>0) { condF = as.numeric(rep(NA, dim(Bx)[2])); cat("Conditional F statistics did not converge to positive values - should the sample sizes be larger?\n") }
              
              return(new("MVGMM",
                         robust = FALSE,
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         Correlation = ld,
                         ExpCorrelation = mis.cor.x,
                         
                         CondFstat = as.numeric(condF),
                         Estimate = as.numeric(liml.nr),
                         StdError = as.numeric(sqrt(diag(var.liml.nr))),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),
                         
                         Pvalue = as.numeric(pvalue),
                         
                         Alpha = alpha,
                         
                         Heter.Stat = c(Q.nr,
                                        pvalue.heter.stat)))
            }
            
            if(robust == TRUE){
              # non-robust LIML estimate based on an incorrect weighting matrix (estimate should be consistent)
              g <- function(tet){as.vector(gamY_est - (gamX_est%*%tet))}
              Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
              tet.gg <- nlminb(rep(0,K),objective=Q.gg)$par
              Om.nr <- function(tet){as.matrix(SigY+(kronecker(t(tet),diag(p))%*%SigX%*%t(kronecker(t(tet),diag(p)))))}
              Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
              G <- -gamX_est
              DQ.nr <- function(tet){2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet))}
              liml.nr <- nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$par
              var.liml.nr <- as.matrix(solve(t(G)%*%solve(Om.nr(liml.nr))%*%G))
              Q.nr <- Q.nr(liml.nr)
              
              # estimating the overdispersion variance parameter
              Om <- function(tet,kappa){as.matrix(SigY+(diag(p)*kappa/ny)+(kronecker(t(tet),diag(p))%*%SigX%*%t(kronecker(t(tet),diag(p)))))}
              Q.kap <- function(kappa){as.numeric(t(g(liml.nr))%*%solve(Om(liml.nr,kappa))%*%g(liml.nr))-(p-K)}
              kap.sq <- seq(-50,50,0.1)
              Q.sq <- vector(,length=length(kap.sq))
              for (q in 1:length(kap.sq)){Q.sq[q] <- Q.kap(kap.sq[q])}
              if(sum(Q.sq>0)==0){kappa.est <- 0; default <- 1}
              if(sum(Q.sq>0)>0){kappa.est <- uniroot(Q.kap,c(kap.sq[max(which(Q.sq>0))],20), extendInt="yes", ...)$root; default <- 0}
              
              if(default == 1){cat("overdispersion parameter could not be estimated by uniroot. Default value of 0 selected, and unrobust estimates are given. \n")}
              
              # robust LIML estimate based on a plug-in overdispersion parameter estimate
              Om.r <- function(tet){as.matrix(SigY+(diag(p)*max(0,kappa.est)/ny)+(kronecker(t(tet),diag(p))%*%SigX%*%t(kronecker(t(tet),diag(p)))))}
              Q <- function(tet){as.numeric(t(g(tet))%*%solve(Om.r(tet))%*%g(tet))}
              DQ <- function(tet){2*as.matrix(t(G)%*%solve(Om.r(tet))%*%g(tet))}
              liml <- nlminb(liml.nr,objective=Q,gradient=DQ)$par
              var.liml <- as.matrix(solve(t(G)%*%solve(Om.r(liml))%*%G))
              
              ciLower <- ci_normal("l", liml, sqrt(diag(var.liml)), alpha)
              ciUpper <- ci_normal("u", liml, sqrt(diag(var.liml)), alpha)
              pvalue.heter.stat <- NA
              pvalue <- 2*pnorm(-abs(liml/sqrt(diag(var.liml))))
              if(sum(condF<0)>0) { condF = as.numeric(rep(NA, dim(Bx)[2])); cat("Conditional F statistics did not converge to positive values - should the sample sizes be larger?\n") }
              
              return(new("MVGMM",
                         robust = TRUE,
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         Correlation = ld,
                         ExpCorrelation = mis.cor.x,
                         
                         CondFstat = as.numeric(condF),
                         Estimate = as.numeric(liml),
                         StdError = as.numeric(sqrt(diag(var.liml))),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),
                         Overdispersion = kappa.est,
                         
                         Pvalue = as.numeric(pvalue),
                         
                         Alpha = alpha,
                         Heter.Stat = as.numeric(Q(liml))
              )
              )
            }
            
            else {
              cat("The robust argument must be TRUE or FALSE. \n")
              cat("See documentation for details. \n")
            }
  }
)