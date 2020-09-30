#' @docType methods
#' @rdname mr_loo

setMethod("mr_loo",
          "MRInput",
          function(object,
                   alpha = 0.05)
            {
            
            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            
            snps = object@snps
            #In case of named variants, below code
            if("snp" %in% snps) {
              snps <- paste("snp", 1:length(Bx), sep = "-")
            } else {
              snps <- snps
            }
            
            n = length(Bx)

            
            
            estimates = sapply(1:n, function(j){mr_ivw(mr_input(Bx[-j], Bxse[-j], By[-j], Byse[-j]), alpha = alpha)$Estimate})
            CI_lower = sapply(1:n, function(j){mr_ivw(mr_input(Bx[-j], Bxse[-j], By[-j], Byse[-j]), alpha = alpha)$CILower})
            CI_upper = sapply(1:n, function(j){mr_ivw(mr_input(Bx[-j], Bxse[-j], By[-j], Byse[-j]), alpha = alpha)$CIUpper})
            
            dframe = data.frame(snps, estimates, CI_lower, CI_upper)

            
            ivw_output = mr_ivw(mr_input(Bx,Bxse,By,Byse), alpha = alpha)
            ivw_estimate = ivw_output$Estimate
            ivw_lower = ivw_output$CILower
            ivw_upper = ivw_output$CIUpper
            
            ivw_row = data.frame("IVW estimate", ivw_estimate, ivw_lower, ivw_upper)
            names(ivw_row) = names(dframe)
            
            dframe = rbind(dframe, ivw_row)
            
            Interval_type <- paste(100*(1-alpha), "% CI)", sep = "")
            Y_label <-  paste("Leave-one-out causal estimate (", Interval_type, sep = "")
            
            dframe$snps = factor(dframe$snps, levels = rev(dframe$snps))
            
            forest <- ggplot(data = dframe, aes(y=snps, x=estimates, xmin=CI_lower, xmax = CI_upper)) +
              geom_point(shape = 15) +
              geom_point(data = ivw_row, aes(y=snps, x = estimates), shape = 23, fill = 'yellow', size = 3) +
              geom_linerange() +
              geom_vline(xintercept = 0, lty= 2) +
              ylab("Variants") + xlab(Y_label) +
              theme_classic()
            
            print(forest)
          })