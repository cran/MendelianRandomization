#' @docType methods
#' @rdname mr_pcgmm

setMethod("mr_pcgmm",
          "MRInput",
          function(object, nx, ny, r = NULL, thres=0.999, robust=TRUE, alpha=0.05, ...){
            
            bx = object@betaX
            by = object@betaY
            sx = object@betaXse
            sy = object@betaYse
            ld = object@correlation
            if(sum(is.na(object@correlation)) != 0){cat("Genetic variant correlation (linkage disequilibrium) matrix is required but not provided. \n")}
            p <- length(bx)
            
            
            # estimate principal components
            Phi <- ((abs(bx)/sy)%*%t(abs(bx)/sy))*ld
            if(is.na(thres)) { thres <- 0.999 }
            if(missing(r)){r <- which(cumsum(prcomp(Phi,scale=FALSE)$sdev^2/sum((prcomp(Phi,scale=FALSE)$sdev^2)))>thres)[1]
            } else {r <- r}
            
            pc.check <- ifelse((r<2),1,0)
            if(pc.check==1){cat("More principal components needed to perform robust analysis. Consider a larger number for r or a larger threshold. \n")}
            
            lambda <- sqrt(p)*prcomp(Phi,scale=FALSE)$rotation[,1:r]
            evec <- eigen((t(lambda)%*%lambda))$vectors
            eval <- eigen((t(lambda)%*%lambda))$values
            lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
            dim(lambda) <- c(p,r)
            
            # outcome quantities of interest
            ay <- 1/((ny*sy^2)+by^2)
            Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
            By <- ay*by
            Ay.f <- t(lambda)%*%Ay%*%lambda; By.f <- as.vector(t(lambda)%*%By)
            SigY <- solve(Ay.f)*(1-as.numeric(t(By.f)%*%solve(Ay.f)%*%By.f)); SigY <- SigY*(1/(ny-r+1))
            gamY_est <- as.vector(solve(Ay.f)%*%By.f)
            
            # exposure quantities of interest
            ax <- 1/((nx*sx^2)+bx^2)
            Ax <- (sqrt(ax)%*%t(sqrt(ax)))*ld
            Bx <- ax*bx
            Ax.f <- t(lambda)%*%Ax%*%lambda; Bx.f <- as.vector(t(lambda)%*%Bx)
            SigX <- solve(Ax.f)*(1-as.numeric(t(Bx.f)%*%solve(Ax.f)%*%Bx.f)); SigX <- SigX*(1/(nx-r+1))
            gamX_est <- as.vector(solve(Ax.f)%*%Bx.f)
            
            # F-test statistic
            fstat = as.numeric(t(gamX_est)%*%solve(SigX)%*%gamX_est)/r
            
            
            if(robust == FALSE & pc.check == 0){
              # non-robust LIML estimate
              g <- function(tet){as.vector(gamY_est - (gamX_est*tet))}
              Om.nr <- function(tet){as.matrix(SigY + (SigX*(tet^2)))}
              Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
              G <- -gamX_est
              
              init.val <- seq(-1,1,0.2)
              Q.init <- vector(,length=length(init.val))
              for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q.nr, method="Brent",lower=-1e2,upper=1e2)$value}
              liml.nr <- optim(init.val[which.min(Q.init)[[1]]], Q.nr, method="Brent",lower=-1e2,upper=1e2)$par
              var.liml.nr <- as.numeric(solve(t(G)%*%solve(Om.nr(liml.nr))%*%G))
              Q.nr <- Q.nr(liml.nr)
              ciLower <- ci_normal("l", liml.nr, sqrt(var.liml.nr), alpha)
              ciUpper <- ci_normal("u", liml.nr, sqrt(var.liml.nr), alpha)
              pvalue.heter.stat <- pchisq(Q.nr, df = r-1, lower.tail = F)
              pvalue <- 2*pnorm(-abs(liml.nr/sqrt(var.liml.nr)))
              
              return(new("PCGMM",
                         robust = FALSE,
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         Correlation = object@correlation,
                         
                         Estimate = as.numeric(liml.nr),
                         StdError = as.numeric(sqrt(var.liml.nr)),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),
                         Fstat = fstat,
                         
                         PCs = r,
                         Pvalue = as.numeric(pvalue),
                         
                         Alpha = alpha,
                         
                         Heter.Stat = c(Q.nr,
                                        pvalue.heter.stat)))
            }
            
            if(robust == TRUE & pc.check == 0){
              # non-robust LIML estimate based on an incorrect weighting matrix (estimate should be consistent)
              g <- function(tet){as.vector(gamY_est - (gamX_est*tet))}
              Om.nr <- function(tet){as.matrix(SigY + (SigX*(tet^2)))}
              Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
              G <- -gamX_est
              init.val <- seq(-1,1,0.2)
              Q.init <- vector(,length=length(init.val))
              for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q.nr, method="Brent",lower=-1e2,upper=1e2)$value}
              liml.nr <- optim(init.val[which.min(Q.init)[[1]]], Q.nr, method="Brent",lower=-1e2,upper=1e2)$par
              
              # estimating the overdispersion variance parameter
              Om <- function(tet,kappa){as.matrix(SigY+(diag(r)*kappa/ny)+(SigX*(tet^2)))}
              Q.kap <- function(kappa){as.numeric(t(g(liml.nr))%*%solve(Om(liml.nr,kappa))%*%g(liml.nr))-(r-1)}
              kap.sq <- seq(-50,50,0.1)
              Q.sq <- vector(,length=length(kap.sq))
              for (q in 1:length(kap.sq)){Q.sq[q] <- Q.kap(kap.sq[q])}
              if(sum(Q.sq>0)==0){kappa.est <- 0; default <- 1}
              if(sum(Q.sq>0)>0){kappa.est <- uniroot(Q.kap,c(kap.sq[max(which(Q.sq>0))],20), extendInt="yes", ...)$root; default <- 0}
              
              if(default == 1){cat("overdispersion parameter could not be estimated by uniroot. Default value of 0 selected, and unrobust estimates are given. \n")}
              
              # robust LIML estimate based on a plug-in overdispersion parameter estimate
              Om.r <- function(tet){as.matrix(SigY+(diag(r)*max(0,kappa.est)/ny)+(SigX*(tet^2)))}
              Q <- function(tet){as.numeric(t(g(tet))%*%solve(Om.r(tet))%*%g(tet))}
              Q.init <- vector(,length=length(init.val))
              for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value}
              liml <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
              var.liml <- as.numeric(solve(t(G)%*%solve(Om.r(liml))%*%G))
              
              ciLower <- ci_normal("l", liml, sqrt(var.liml), alpha)
              ciUpper <- ci_normal("u", liml, sqrt(var.liml), alpha)
              pvalue.heter.stat <- pchisq(as.numeric(Q(liml)), df = r-1, lower.tail = F)
              pvalue <- 2*pnorm(-abs(liml/sqrt(var.liml)))
              
              return(new("PCGMM",
                         robust = TRUE,
                         Exposure = object@exposure,
                         Outcome = object@outcome,
                         Correlation = object@correlation,
                         
                         Estimate = as.numeric(liml),
                         StdError = as.numeric(sqrt(var.liml)),
                         CILower =  as.numeric(ciLower),
                         CIUpper = as.numeric(ciUpper),
                         Fstat = fstat,
                         Overdispersion = kappa.est,
                         
                         PCs = r,
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