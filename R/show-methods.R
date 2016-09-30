

setMethod("show",
          "MRInput",
          function(object){

            exposure <- object@exposure
            outcome <- object@outcome
            snps <- object@snps

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse

            if("snp" %in% snps) {
              snps <- paste("snp", 1:length(Bx), sep = "_")
            } else {
              snps <- snps
            }

            betaDF <- data.frame(snps, Bx, Bxse, By, Byse)

            colnames(betaDF) <- c("SNP",
                                  paste(exposure, ".beta", sep = ""),
                                  paste(exposure, ".se", sep = ""),
                                  paste(outcome, ".beta", sep = ""),
                                  paste(outcome, ".se", sep = ""))

            print(betaDF, digits = 3)
          })

#--------------------------------------------------------------------------------------------

setMethod("show",
          "WeightedMedian",
          function(object){

            Interval_type <- paste(100*(1-object@Alpha), "% CI", sep = "")
            Statistic <- c("Method", "Estimate", "Std Error", Interval_type, "", "p-value")

            Median_type <- paste(simpleCap(object@Type), " median method", sep = "")

            Value <- c(Median_type,
                       decimals(object@Estimate,3),
                       decimals(object@StdError,3),
                       paste(decimals(object@CILower, 3), ",", sep = ""),
                       decimals(object@CIUpper,3),
                       decimals(object@Pvalue, 3))

            output.table <- data.frame(matrix(Value, nrow = 1))
            colnames(output.table) <- Statistic

            cat("\n",Median_type, "\n\n")
            cat("Number of Variants :", object@SNPs, "\n")

            cat("------------------------------------------------------------------\n")
            print(output.table, quote = F, row.names = FALSE, justify= "left")
            cat("------------------------------------------------------------------\n")
          }
)


#--------------------------------------------------------------------------------------------

setMethod("show",
          "IVW",
          function(object){

            Interval_type <- paste(100*(1-object@Alpha), "% CI", sep = "")
            Statistic <- c("Method", "Estimate", "Std Error", Interval_type, "", "p-value")

            Value <- c("IVW", decimals(c(object@Estimate, object@StdError),3),
                       paste(decimals(object@CILower, 3), ",", sep = ""), decimals(c(object@CIUpper, object@Pvalue), 3))

            output.table <- data.frame(matrix(Value, nrow = 1))
            colnames(output.table) <- Statistic
            correlation <- ifelse(sum(is.na(object@Correlation)) == 0,
                                  "correlated", "uncorrelated")
            penalized <- ifelse(object@Penalized == TRUE, "Weights of genetic variants with heterogeneous causal estimates have been penalized. ", "")
            robust <- ifelse(object@Robust == TRUE, "Robust regression used.", "")

            cat("\nInverse-variance weighted method\n")
            cat("(variants ", correlation, ", ", object@Model, "-effect model)\n\n" , sep = "")

            cat("Number of Variants :", object@SNPs, "\n")
            cat(penalized, robust, "\n", sep = "")

            cat("------------------------------------------------------------------\n")
            print(output.table, quote = F, row.names = FALSE, justify = "left")
            cat("------------------------------------------------------------------\n")

            cat("Residual standard error = ", decimals(object@RSE, 3), "\n")
            if(object@Model == "fixed") { cat("Residual standard error is set to 1 in calculation of confidence interval by fixed-effect assumption.\n") }
            if(object@RSE<1) { cat("Residual standard error is set to 1 in calculation of confidence interval when its estimate is less than 1.\n") }
            if(is.na(object@Heter.Stat[1])) {
              cat("Heterogeneity is not calculated when weights are penalized, or when there is only one variant in the analysis.")
            } else {
            cat("Heterogeneity test statistic = ", decimals(object@Heter.Stat[1],4), " on ", object@SNPs -1,
                    " degrees of freedom, (p-value = ", decimals(object@Heter.Stat[2], 4),")\n", sep = "")
            }
          }
)

#--------------------------------------------------------------------------------------------

setMethod("show",
          "MaxLik",
          function(object){

            Interval_type <- paste(100*(1-object@Alpha), "% CI", sep = "")
            Statistic <- c("Method", "Estimate", "Std Error", Interval_type, "", "p-value")

            Value <- c("MaxLik", decimals(c(object@Estimate, object@StdError),3),
                       paste(decimals(object@CILower, 3), ",", sep = ""), decimals(c(object@CIUpper, object@Pvalue), 3))

            output.table <- data.frame(matrix(Value, nrow = 1))
            colnames(output.table) <- Statistic
            correlation <- ifelse(sum(is.na(object@Correlation)) == 0,
                                  "correlated", "uncorrelated")

            cat("\nMaximum-likelihood method\n")
            cat("(variants ", correlation, ", ", object@Model, "-effect model)\n\n" , sep = "")

            cat("Number of Variants :", object@SNPs, "\n")

            cat("------------------------------------------------------------------\n")
            print(output.table, quote = F, row.names = FALSE, justify = "left")
            cat("------------------------------------------------------------------\n")

            cat("Residual standard error = ", decimals(object@RSE, 3), "\n")
            if(object@Model == "fixed") { cat("Residual standard error is set to 1 in calculation of confidence interval by fixed-effect assumption.\n") }
            if(object@RSE<1) { cat("Residual standard error is set to 1 in calculation of confidence interval when its estimate is less than 1.\n") }
            if(object@Heter.Stat[1] < 1e-16) {
              cat("Heterogeneity is not calculated when there is only one variant in the analysis.")
            } else {
            cat("Heterogeneity test statistic = ", decimals(object@Heter.Stat[1],4), " on ", object@SNPs -1,
                    " degrees of freedom, (p-value = ", decimals(object@Heter.Stat[2], 4),")\n", sep = "")
            }
          }
)

#--------------------------------------------------------------------------------------------

setMethod("show",
          "Egger",
          function(object){
            Interval_type <- paste(100*(1-object@Alpha), "% CI", sep = "")
            Statistic <- c("Method", "Estimate", "Std Error", Interval_type, "", "p-value")

            Value <- c("MR-Egger", decimals(c(object@Estimate,
                                              object@StdError.Est), 3),
                               paste(decimals(object@CILower.Est, 3), ",", sep = ""),
                                   decimals(c(object@CIUpper.Est,
                                              object@Pvalue.Est), 3),

                       "(intercept)", decimals(c(object@Intercept,
                                                 object@StdError.Int), 3),
paste(decimals(object@CILower.Int, 3), ",", sep = ""),
                                                 
                                      decimals(c(object@CIUpper.Int,
                                                 object@Pvalue.Int), 3))

            output.table <- data.frame(matrix(Value, nrow = 2, byrow = T))
            colnames(output.table) <- Statistic

            correlation <- ifelse(sum(is.na(object@Correlation)) == 0, "correlated", "uncorrelated")
            penalized <- ifelse(object@Penalized == TRUE, "Weights of genetic variants with heterogeneous causal estimates have been penalized. ", "")
            robust <- ifelse(object@Robust == TRUE, "Robust model used.", "")

            cat("\nMR-Egger method\n")
            cat("(variants ", correlation, ", ", object@Model, "-effect model)\n\n" , sep = "")

            cat("Number of Variants = ", object@SNPs, "\n")
            cat(penalized, robust, "\n", sep = "")

            cat("------------------------------------------------------------------\n")
            print(output.table, quote = F, row.names = FALSE, justify= "left")
            cat("------------------------------------------------------------------\n")

            cat("Residual Standard Error : ", decimals(object@RSE, 3), "\n")
            if(object@RSE<1) { cat("Residual standard error is set to 1 in calculation of confidence interval when its estimate is less than 1.\n") }
            if(is.na(object@Heter.Stat[1])) {
              cat("Heterogeneity not calculated when weights are penalized.\n")
            } else {
              cat("Heterogeneity test statistic = ", decimals(object@Heter.Stat[1],4), " on ", object@SNPs - 2,
                  " degrees of freedom, (p-value = ", decimals(object@Heter.Stat[2], 4),")\n", sep = "")
            }
              if(!is.nan(object@I.sq)) {
cat("I^2_GX statistic: ", decimals(object@I.sq*100, 1), "%\n", sep="") }
          }
)
#--------------------------------------------------------------------------------------------

setMethod("show",
          "MRAll",
          function(object){

            df <- slot(object, "Values")
            df[,2:6] <- decimals(df[,2:6], 3)

            space <- rep("", 6)

            if(object@Method == "all"){
              df <- rbind(df[1:3,],
                          space,
                          df[4:7,],
                          space,
                          df[8:15,])
            } else {
              df <- df
            }

            #print(object@Values, digits = 2, justify = "left")
            stargazer(df, type = "text", summary = FALSE, rownames = FALSE)
          }
)

#--------------------------------------------------------------------------------------------
