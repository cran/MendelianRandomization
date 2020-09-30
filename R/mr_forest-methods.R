#' @docType methods
#' @rdname mr_forest

setMethod("mr_forest",
          signature = "MRInput",
          function(object,
                   alpha = 0.05,
                   snp_estimates = TRUE,
                   methods = "ivw",
                   ordered = FALSE){
            
            
            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse
            
            snps = object@snps
            #In case of named variants, below code
            if("snp" %in% snps) {
              snps <- paste("snp", 1:length(Bx), sep = "_")
            } else {
              snps <- snps
            }
            

            if (snp_estimates) {
              #Calculate the estimates and CIs
              estimates = By / Bx
              CI_range = qnorm(1-alpha/2)
              #Note Steve asked for / Bx but am using abs cause want CI_lower < CI_upper
              CI_lower = estimates - (CI_range * Byse)/abs(Bx)
              CI_upper = estimates + (CI_range * Byse) / abs(Bx)
              #Create the dataframe with these values
              dframe = data.frame(snps, estimates, CI_lower, CI_upper)
            } else {
              #Don't calculate the SNPs
              #Still need to create the data frame so we can append
              dframe = data.frame(matrix(ncol = 4, nrow = 0))
              names(dframe) = c("snps", "estimates", "CI_lower", "CI_upper")
            }
            
            #Calculate and implement the rows.
            
            #Define a data frame to contain the outputs of the various methods
            calculated = data.frame(matrix(ncol = 4, nrow = 0))
            names(calculated) = c("snps", "estimates", "CI_lower", "CI_upper")
            #IVW Method
            if ("ivw" %in% methods) {
              ivw_output = mr_ivw(mr_input(Bx,Bxse,By,Byse), alpha = alpha)
              ivw_estimate = ivw_output$Estimate
              ivw_CI_lower = ivw_output$CILower
              ivw_CI_upper = ivw_output$CIUpper
              
              ivw_row = data.frame("IVW estimate", ivw_estimate, ivw_CI_lower, ivw_CI_upper)
              names(ivw_row) = names(dframe)
              calculated = rbind(calculated, ivw_row)
            }
            #Median Method
            if ("median" %in% methods) {
              median_output = mr_median(mr_input(Bx, Bxse, By, Byse), weighting="simple", alpha= alpha)
              median_estimate = median_output$Estimate
              median_CI_lower = median_output$CILower
              median_CI_upper = median_output$CIUpper
              
              median_row = data.frame("Simple median estimate", median_estimate, median_CI_lower, median_CI_upper)
              names(median_row) = names(dframe)
              calculated = rbind(calculated, median_row)
            } 
            if("wmedian" %in% methods) {
              wmedian_output = mr_median(mr_input(Bx, Bxse,By, Byse), weighting="weighted", alpha = alpha)
              wmedian_estimate = wmedian_output$Estimate
              wmedian_CI_lower = wmedian_output$CILower
              wmedian_CI_upper = wmedian_output$CIUpper
              
              wmedian_row = data.frame("Weighted median estimate", wmedian_estimate, wmedian_CI_lower, wmedian_CI_upper)
              names(wmedian_row) = names(dframe)
              calculated = rbind(calculated, wmedian_row)
            } 
            if ("egger" %in% methods) {
              egger_output = mr_egger(mr_input(Bx, Bxse, By, Byse), alpha = alpha)
              egger_estimate = egger_output$Estimate
              egger_CI_lower = egger_output$CILower.Est
              egger_CI_upper = egger_output$CIUpper.Est
              
              egger_row = data.frame("MR-Egger estimate", egger_estimate, egger_CI_lower, egger_CI_upper)
              names(egger_row) = names(dframe)
              calculated = rbind(calculated, egger_row)
            } 
            if ("maxlik" %in% methods) {
              maxlik_output = mr_maxlik(mr_input(Bx, Bxse, By, Byse), alpha = alpha)
              maxlik_estimate = maxlik_output$Estimate
              maxlik_CI_lower = maxlik_output$CILower
              maxlik_CI_upper = maxlik_output$CIUpper
              
              maxlik_row = data.frame("Maximum likelihood estimate", maxlik_estimate, maxlik_CI_lower, maxlik_CI_upper)
              names(maxlik_row) = names(dframe)
              calculated = rbind(calculated, maxlik_row)
            } 
            if ("mbe" %in% methods) {
              mbe_output = mr_mbe(mr_input(Bx, Bxse, By, Byse), alpha = alpha)
              mbe_estimate = mbe_output$Estimate
              mbe_CI_lower = mbe_output$CILower
              mbe_CI_upper = mbe_output$CIUpper
              
              mbe_row = data.frame("Mode-based estimate", mbe_estimate, mbe_CI_lower, mbe_CI_upper)
              names(mbe_row) = names(dframe)
              calculated = rbind(calculated, mbe_row)
            } 
            if ("conmix" %in% methods) {
              conmix_output = mr_conmix(mr_input(Bx, Bxse, By, Byse), alpha = alpha)
              conmix_estimate = conmix_output$Estimate
              conmix_CI_lower = conmix_output$CILower
              conmix_CI_upper = conmix_output$CIUpper
              
              conmix_row = data.frame("Contamination mixture estimate", conmix_estimate, conmix_CI_lower, conmix_CI_upper)
              names(conmix_row) = names(dframe)
              calculated = rbind(calculated, conmix_row)
            }
            
            #Deal with ordered
            #This is really nasty for some reason to do with how ggplot likes to allow you to order things
            #Also need to reorder factors backwards as plotting x,y flipped
            if (ordered == TRUE){
              #This works
              #The idea is that we first order the estimates, then reverse the order
              #Then add the calculated ones to the beginning which is actually the end
              factor_order = rev(dframe$snps[order(dframe$estimates)])
              dframe$snps = factor(dframe$snps, levels = factor_order)
              dframe = rbind(dframe, calculated)
              dframe$snps = factor(dframe$snps, levels= c(calculated$snps, factor_order))
            } else{
              #This works simply on the merit that these need to be reversed
              dframe = rbind(dframe, calculated)
              dframe$snps = factor(dframe$snps, rev(dframe$snps))
            }
            
            Interval_type <- paste(100*(1-alpha), "% CI)", sep = "")
            Y_label <-  paste("Causal estimate (", Interval_type, sep = "")
            
            if (!snp_estimates) {
              X_label <- "Methods"
            } else {
              X_label <- "Variants"
            }
            
            
            forest <- ggplot(data = dframe, aes(y=snps, x=estimates, xmin=CI_lower, xmax = CI_upper)) +
              geom_point(shape = 15) +
              geom_point(data = calculated, aes(y=snps, x = estimates), shape = 18, size = 3) +
              geom_linerange() +
              geom_vline(xintercept = 0, lty= 2) +
              ylab(X_label) + xlab(Y_label) +
              theme_classic()
            
            print(forest)
            
          })