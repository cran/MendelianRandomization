# Required Functions
#' Prior weight function
#'
#' @description Internal function for prior weights.
#'
#' @param model.size Size of subset.
#' @param N.obs Number of genetic variants.
#' @param prob.valid.inst Probability (prior) of a genetic variant being a valid instrument.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Prior weight.
#'
#' @examples model.prior(8, 15, 0.5)
#'
#' @export

model.prior = function(model.size, N.obs, prob.valid.inst){
  pr = (prob.valid.inst^model.size)*(1-prob.valid.inst)^(N.obs-model.size)
  return(pr)
}

#' Heterogeneity-penalized weight function
#'
#' @description Internal function for calculating heterogeneity-penalized weights.
#'
#' @param prob.valid.inst Probability (prior) of a genetic variant being a valid instrument.
#' @param bx Genetic associations with risk factor.
#' @param by Genetic associations with outcome.
#' @param byse Standard errors of genetic associations with outcome.
#'
#' @keywords internal
#'
#' @details None.
#'
#' @return Heterogeneity-penalized weight.
#'
#' @examples het.weight(0.5, calcium, fastgluc, fastglucse)
#'
#' @export

het.weight = function(prob.valid.inst, bx, by, byse){
  J = length(by);
  theta.est = by/bx;
  theta.se = abs(byse/bx);
  tmp.1 = by/byse;
  tmp.2 = bx/byse;
  theta.se.sq = theta.se^2;
  log.theta.se = log(theta.se);
  est  = seest = vector("numeric", 2^J-1);
  het.weight = vector("numeric", 2^J-1);
   #
  count = 0;
 for(n in 1:J){
    perms = choose(J,n);
    inc = sparseMatrix(i=as.vector(t(replicate(n,1:perms))),
                       j=as.vector(t(getall(iterpc(J,n,c(1:J))))),
                       x=1, dims = c(perms,J));
    # sparse binary inclusion matrix
    # 1 denotes an instrument is included in the model
    # each row represents a particular model
    est.sum = inc%*%(theta.est/theta.se.sq);
    recip.var.ivw = inc%*%(1/theta.se.sq);
    est.ivw = est.sum/recip.var.ivw;
    est[(count+1):(count+perms)] = est.ivw;
  if(n>1){
     tmp = t(replicate(J, as.vector(est.ivw)));
   if(n<J){
     psi.hat = sqrt((1/(n-1))*rowSums(t(t(inc)*(tmp.1^2 - 2*tmp*(tmp.1*tmp.2) +
                    (tmp^2)*(tmp.2^2)))))
           }
      else{
        psi.hat = sqrt((1/(n-1))*sum(tmp.1^2 - 2*tmp*(tmp.1*tmp.2) +
                    (tmp^2)*(tmp.2^2)));
           }
      psi.hat[which(psi.hat<1)] = 1;
      seest[(count+1):(count+perms)] = psi.hat/sqrt(recip.var.ivw);
          }
  else if(n==1){
      seest[(count+1):(count+perms)] = inc%*%theta.se;
                }
     #
  if(n>1){
      het.exponent = rowSums(inc*t(t(t(t(inc)*theta.est) -
            as.vector(est.ivw))^2/theta.se.sq));
      het.weight[(count+1):(count+perms)] =
                 exp(-(inc%*%(log.theta.se)+0.5*het.exponent))*
                 model.prior(n,J,prob.valid.inst);
          }
    count = count+perms;
  }  # ends for loop
  newlist = list(het.weight, est, seest);
  return(newlist)
}


#' @docType methods
#' @rdname mr_hetpen

setMethod("mr_hetpen",
          "MRInput",
          function(object,
                   prior  = 0.5,
                   CIMin  = -1,
                   CIMax  = 1,
                   CIStep = 0.001,
                   alpha = 0.05){

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse

            nsnps <- length(Bx)

      if (nsnps > 30) { cat("Too many genetic variants for this version of the method to run. We are working on a solution!\n") }
    else {
      if (nsnps > 25) { cat("Method is likely to take a long time to run... Please be patient - meanwhile, we are working on a solution!\n") }

            results = het.weight(prior, Bx, By, Byse);
            het.weight.norm = results[[1]]/sum(results[[1]]);
               # normalized heterogeneity-penalized weights
            sumlik=NULL
            point = matrix(seq(CIMin, CIMax, CIStep), ncol = 1);
               #
            sumlik = vapply(point,function(i){sum(het.weight.norm*dnorm(rep(i,length(het.weight.norm)), results[[2]], results[[3]]))}, 1);
            whichin = which(2*log(sumlik)>(2*max(log(sumlik))-qchisq(1-alpha, df=1)));
   # provides an index of estimate values in the 95% confidence interval
            betaHetPen = CIMin+CIStep*(which.max(log(sumlik))-1)
   # modal estimate
            CIRange    = CIMin+CIStep*(whichin-1);
            CILower <- c(min(CIRange), CIRange[which(diff(CIRange)>1.01*CIStep)+1])
            CIUpper <- c(CIRange[which(diff(CIRange)>1.01*CIStep)], max(CIRange))


                  return(new("MRHetPen",
                             Exposure = object@exposure,
                             Outcome = object@outcome,
                             Prior = as.numeric(prior),

                             Estimate = as.numeric(betaHetPen),
                             CIRange  = as.numeric(CIRange),
                             CILower  = as.numeric(CILower),
                             CIUpper  = as.numeric(CIUpper),

                             CIMin    = as.numeric(CIMin),
                             CIMax    = as.numeric(CIMax),
                             CIStep   = as.numeric(CIStep),

                             SNPs = nsnps,
                             Alpha = alpha))
             }
        }
)
