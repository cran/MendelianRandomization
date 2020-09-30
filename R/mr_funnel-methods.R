#' @docType methods
#' @rdname mr_funnel

setMethod("mr_funnel",
          "MRInput",
          function(object,
                   CI = TRUE){
            
            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            
            r = By / Bx
            precision = abs(Bx/ Byse)
            CI_lower = (By-qnorm(0.975)*Byse)/Bx
            CI_upper = (By+qnorm(0.975)*Byse)/Bx
            
            dframe = data.frame(r, precision, CI_lower, CI_upper)
            X_label = "Causal estimate (95% CI)"
            Y_label = "Precision"
            ivw_estimate = mr_ivw(mr_input(Bx, Bxse, By, Byse))$Estimate
            
            if (CI) {
              funnel = ggplot(data= dframe, aes(x= r, y= precision, xmin=CI_lower, xmax=CI_upper)) +
                geom_point() +
                geom_linerange() +
                geom_vline(xintercept = 0) +
                geom_vline(xintercept = ivw_estimate, lty =2) +
                ylab(Y_label) + xlab(X_label) +
                theme_classic()
              
              print(funnel)
            } else {
              cat("Plotting without confidence intervals")
              funnel = ggplot(data= dframe, aes(x= r, y= precision, xmin=CI_lower, xmax=CI_upper)) +
                geom_point() +
                geom_vline(xintercept = 0) +
                geom_vline(xintercept = ivw_estimate, lty =2) +
                ylab(Y_label) + xlab(X_label) +
                theme_classic()
              
              print(funnel)
            }
            

          })