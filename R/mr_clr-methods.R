#' @docType methods
#' @rdname mr_clr

setMethod("mr_clr",
          "MRInput",
          function(object, nx, ny, alpha=0.05,
                   CIMin  = -10,
                   CIMax  = 10,
                   CIStep = 0.01){
            
            bx = object@betaX
            by = object@betaY
            sx = object@betaXse
            sy = object@betaYse
            if(sum(is.na(object@correlation)) != 0){ld = diag(length(bx))}
            if(sum(is.na(object@correlation)) == 0){ld = object@correlation}
            p <- length(bx)
            
            # outcome quantities of interest
            ay <- 1/((ny*sy^2)+by^2)
            Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
            By <- ay*by
            Cov.y <- solve(Ay)*(1-as.numeric(t(By)%*%solve(Ay)%*%By)); Cov.y <- Cov.y*(1/(ny-p+1))
            beta.y <- as.vector(solve(Ay)%*%By)
            
            # exposure quantities of interest
            ax <- 1/((nx*sx^2)+bx^2)
            Ax <- (sqrt(ax)%*%t(sqrt(ax)))*ld
            Bx <- ax*bx
            Cov.x <- solve(Ax)*(1-as.numeric(t(Bx)%*%solve(Ax)%*%Bx)); Cov.x <- Cov.x*(1/(nx-p+1))
            beta.x <- as.vector(solve(Ax)%*%Bx)
            
            # function to compute identification-robust tests
            mrTests = function(beta.x,beta.y,Cov.x,Cov.y,beta0,alpha = 0.05, pval = FALSE){
              Cov.x_inv = solve(Cov.x); Cov.y_inv = solve(Cov.y); k = length(beta.x)
              evec.Sn <- eigen(Cov.y + beta0^2 * Cov.x)$vectors
              eval.Sn <- eigen(Cov.y + beta0^2 * Cov.x)$values
              sqrt.inv.Sn <- solve(evec.Sn%*%diag(sqrt(eval.Sn))%*%t(evec.Sn))
              evec.Rn <-  eigen(beta0^2 * Cov.y_inv + Cov.x_inv)$vectors
              eval.Rn <-  eigen(beta0^2 * Cov.y_inv + Cov.x_inv)$values
              sqrt.inv.Rn <- solve(evec.Rn%*%diag(sqrt(eval.Rn))%*%t(evec.Rn))
              Sn = sqrt.inv.Sn %*% (beta.y - beta0* beta.x)
              Rn = sqrt.inv.Rn %*% (beta0*Cov.y_inv %*% beta.y + Cov.x_inv %*% beta.x)
              Qs = c(t(Sn) %*% Sn); Qr = c(t(Rn) %*% Rn); Qsr = c(t(Sn) %*% Rn)
              mrAR = Qs; mrK = Qsr^2/Qr; mrCLR = 1/2*(Qs - Qr + ((Qs + Qr)^2 - 4 * (Qs*Qr - Qsr^2))^(1/2))
              rej1 = mrAR > qchisq(1-alpha, k); rej2 = mrK > qchisq(1 - alpha, 1)
              if(k == 1)  { rej3 = mrCLR > qchisq(1 - alpha, 1) } else {
                K = exp(lgamma(k/2)- lgamma( (k - 1)/2 )) / (pi^(1/2))
                integrand <- function(x){pchisq( (Qr + mrCLR)/(1 + Qr*x^2*(mrCLR)^(-1) ),k) * (1 - x^2)^((k - 3)/2)}
                pvalue = 1 - 2*K*integrate(integrand,lower=0,upper=1,subdivisions=5000,stop.on.error=FALSE)$value
                rej3 = pvalue < alpha }
              if(pval == TRUE) { # return p values of mrAR, mrK, and mrCLR test
                pval_AR = pchisq(mrAR, k, lower.tail = FALSE)
                pval_K = pchisq(mrK, 1, lower.tail = FALSE)
                pval_CLR = pvalue
                return(c(pval_AR, pval_K, pval_CLR))
              } else return(as.numeric(c(rej1,rej2,rej3))) # return the conclusions for mrAR, mrK, and mrCLR test; 1:reject H_0; 0: do not reject H_0
            }
            
            # compute confidence intervals by inverting identification-robust tests
            mrCI = function(beta.x, beta.y, Cov.x, Cov.y, betaset=seq(CIMin,CIMax,CIStep),alpha=0.05){
              temp = matrix(NA, nrow = 3, ncol = length(betaset)); result = list("mrAR"=NULL,"mrK" = NULL, "mrCLR" = NULL); types = c("mrAR","mrK","mrCLR")
              for(i in 1:length(betaset)) temp[,i] = mrTests(beta.x, beta.y, Cov.x, Cov.y, betaset[i], alpha = alpha)
              for(i in 1:3){
                diff = temp[i,-1] - temp[i,-length(betaset)]
                left = which(diff == -1) + 1; right = which(diff == 1)
                if(temp[i,1] == 0) left = c(1,left)
                if(temp[i,length(betaset)] == 0) right = c(right, length(betaset))
                # cat(types[i],":")
                # if(length(left) == 0) cat("empty; ") else{ for(j in 1:length(left)) cat("[", betaset[left[j]],",",betaset[right[j]],"]; ") }
                result[[i]] =  list("left" = betaset[left], "right" = betaset[right]) }
              return( result ) # return confidence intervals from mrAR, mrK and mrCLR
            }
            
            res <- mrCI(beta.x=beta.x, beta.y=beta.y, Cov.x=Cov.x, Cov.y=Cov.y, alpha=alpha)

           
            return(new("CLR",
                       Exposure = object@exposure,
                       Outcome = object@outcome,
                       Correlation = object@correlation,
                       ARlower = res$mrAR$left, ARupper = res$mrAR$right,
                       Klower = res$mrK$left, Kupper = res$mrK$right,
                       CLRlower = res$mrCLR$left, CLRupper = res$mrCLR$right,
                       CIMin    = as.numeric(CIMin),
                       CIMax    = as.numeric(CIMax),
                       CIStep   = as.numeric(CIStep),

                       Alpha = alpha))
          }
)