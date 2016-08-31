#' @docType methods
#' @rdname mr_plot

setMethod(f = "mr_plot",
          signature = "MRInput",
          function(object, error = TRUE, line = "ivw", orientate = FALSE){

            Bx = object@betaX
            By = object@betaY
            Bxse = object@betaXse
            Byse = object@betaYse

            signedBy = object@betaY*sign(object@betaX)
            signedBx = abs(object@betaX)


 if (orientate == TRUE) {
            By = object@betaY*sign(object@betaX)
            Bx = abs(object@betaX)
            }


            betaDF <- data.frame(Bx, By, Bxse, Byse)

            if(error == TRUE){
                plot_ly(x = Bx, y = By, showlegend = FALSE, data=betaDF, 
                        error_y = list(array = qnorm(0.975)*Byse, opacity = 0.2),
                        error_x = list(array = qnorm(0.975)*Bxse, opacity = 0.2),
                        type = "scatter",
                        mode = "markers",
                        hoverinfo = "text",
                        text = paste( "(", Bx, ",", By ,")", "<br>", "SNP = ", object@snps)) %>%
 plotly::layout(xaxis = list(title = paste("Genetic association with", object@exposure)),
                       yaxis = list(title = paste("Genetic association with", object@outcome)))

              if(line == "ivw"){
                  add_trace(y = lm(By ~ Bx - 1, weights = Byse^-2)$coef[1]*c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), x = c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), name = "ivw", showlegend = FALSE)
              } else if (line == "egger"){
                  add_trace(y = lm(signedBy ~ signedBx, weights = Byse^-2)$coef[1] + lm(signedBy ~ signedBx, weights = Byse^-2)$coef[2]*c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), x = c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), name = "egger", showlegend = FALSE) }
                else { add_trace(NULL) }
            } else if (error == FALSE) {
                plot_ly(x = Bx, y = By,  showlegend = FALSE, data = betaDF,
                        type = "scatter",
                        mode = "markers",
                        hoverinfo = "text",
                        text = paste( "(", Bx, ",", By ,")", "<br>", "SNP = ", object@snps)) %>%
                plotly::layout(xaxis = list(title = paste("Genetic association with", object@exposure)),
                       yaxis = list(title = paste("Genetic association with", object@outcome)))

              if(line == "ivw"){
                  add_trace(y = lm(By ~ Bx - 1, weights = Byse^-2)$coef[1]*c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), x = c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), name = "ivw", showlegend = FALSE)
              } else if (line == "egger"){
                  add_trace(y = lm(signedBy ~ signedBx, weights = Byse^-2)$coef[1] + lm(signedBy ~ signedBx, weights = Byse^-2)$coef[2]*c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), x = c(min(Bx-1.96*Bxse), Bx, max(Bx+1.96*Bxse)), name = "egger", showlegend = FALSE) }
               else { add_trace(NULL) }
            } else {
              cat("Note that error must be one of: TRUE, FALSE.")
            }
          })




#' @docType methods
#' @rdname mr_plot

setMethod(f = "mr_plot",
          signature = "MRAll",
          function(object){

            df <- slot(object, "Values")
            n <- nrow(df)

            new.df <- data.frame(matrix(NA, nrow = n, ncol = 5))

             if(object@Method == "median"){
               new.df <- df[,1:2]
               new.df$Intercept <- rep(0, 3)
               new.df$Color <- c("Blue", "Green", "Red")
             } else if (object@Method == "ivw"  | object@Method == "main"){
               new.df <- df[,1:2]
               new.df$Intercept <- rep(0, 4)
               new.df$Color <- c("Blue", "Green", "Purple", "Red")
             } else if (object@Method == "egger"){
               new.df <- df[c(1,3,5,7),1:2]
               new.df$Intercept <- df[c(2,4,6,8),2]
               new.df$Color <- c("Blue", "Green", "Purple", "Red")
             } else {
              new.df <- df[c(1:8, 10, 12, 14), 1:2]
              new.df$Intercept <- c(rep(0,7), df[c( 9, 11, 13, 15), 2])
            }




 if ( object@Method!="all" ) {
with(data=NULL, {        ggplot(data = NULL, aes(x = object@Data@betaX, y = object@Data@betaY)) +
              geom_point() +
              geom_hline(yintercept = 0, color = "red", alpha = 0.2) +
              geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
              geom_abline(data = new.df, aes(intercept = Intercept, slope = Estimate, color = Color), linetype = "solid",
                          show.legend = T, size = 1) +
              scale_colour_discrete(name  ="Model Type", labels = new.df$Method) +
#              scale_linetype_manual(values = rep(2:6, times = 3), name  ="Model Type", labels = new.df$Method) +
              xlab(paste("Genetic association with", object@Data@exposure)) +
              ylab(paste("Genetic association with", object@Data@outcome)) +
              #ggtitle("Comparison of all methods") +
              theme(
                plot.title = element_text(size = rel(1.5), face = "bold"),
                # Background
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
                panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
                legend.key = element_rect(fill = "white")
              )
} ) # close with
  }
 else {
 with(data=NULL, {           ggplot(data = NULL, aes(x = object@Data@betaX, y = object@Data@betaY)) +
              geom_point() +
              geom_hline(yintercept = 0, color = "red", alpha = 0.2) +
              geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
              geom_abline(data = new.df, aes(intercept = Intercept, slope = Estimate, 
    color = Method,
   linetype = Method),
                          show.legend = TRUE, size = 1) +
              scale_colour_manual(name="Model Type", labels = new.df$Method, 
              breaks = c("Simple median", "Weighted median", "Penalized weighted median", 
                         "IVW", "Penalized IVW", "Robust IVW", "Penalized Robust IVW",
                         "MR-Egger", "Penalized MR-Egger", "Robust MR-Egger", "Penalized Robust MR-Egger"),
              values = c("Red",   "Red",    "Green", "Green",  "Purple", "Purple", "Blue",   "Blue",  "Blue",  "Red",    "Green")) +
                       # IVW      MR-Egger  PenIVW   PenEgger  PenRobIVW PenRobEgg PenWMed   RobIVW   RobEgg    SimpleMed WeightedMed
              scale_linetype_manual(name="Model Type", labels = new.df$Method,
              breaks = c("Simple median", "Weighted median", "Penalized weighted median", 
                         "IVW", "Penalized IVW", "Robust IVW", "Penalized Robust IVW",
                         "MR-Egger", "Penalized MR-Egger", "Robust MR-Egger", "Penalized Robust MR-Egger"),
              values = c("solid", "dashed", "solid", "dashed", "solid",  "dashed", "dotted", "solid", "dashed", "dotted", "dotted")) +
              xlab(paste("Genetic association with", object@Data@exposure)) +
              ylab(paste("Genetic association with", object@Data@outcome)) +
              #ggtitle("Comparison of all methods") +
              theme(
                plot.title = element_text(size = rel(1.5), face = "bold"),
                # Background
                panel.background = element_rect(fill = "white"),
                panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
                panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
                legend.key = element_rect(fill = "white")
              ) }
 ) # close with
  }


          })



