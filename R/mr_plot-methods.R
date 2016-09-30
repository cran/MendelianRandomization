#' @docType methods
#' @rdname mr_plot

setMethod(f = "mr_plot",
          signature = "MRInput",
          function(object, error = TRUE, line = "ivw", 
                   orientate = FALSE, interactive = TRUE, labels = FALSE) {
            
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
            
              # Create the initial plot
              
              generic_plot <- ggplot(data = NULL, aes(x = Bx, y = By, text = object@snps)) +
                geom_point(colour = "black", alpha = 0.5, size = 2) +
                
                theme(panel.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
                      panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
                      
                      axis.text = element_text(size = 10),
                      axis.title = element_text(size = 14),
                      axis.text.x  = element_text(margin = margin(b = 7)),
                      axis.text.y  = element_text(margin = margin(l = 5))) +
                
                xlab(paste("Genetic association with", object@exposure)) +
                ylab(paste("Genetic association with", object@outcome))  + 
                geom_hline(yintercept=0) + 
                geom_vline(xintercept=0)
              
              if(error == T){
                
                generic_plot <- generic_plot +
                  geom_errorbar(aes(ymin = By - qnorm(0.975)*Byse, ymax = By + qnorm(0.975)*Byse), alpha = 0.3, colour = "blue") +
                  geom_errorbarh(aes(xmin = Bx - qnorm(0.975)*Bxse, xmax = Bx + qnorm(0.975)*Bxse), alpha = 0.3, colour = "blue")
                
              } else {
                generic_plot <- generic_plot
              }
              
              if(line == "egger"){
                generic_plot <- generic_plot + geom_abline(intercept = mr_egger(object)$Intercept, 
                                           slope = mr_egger(object)$Estimate, color = "blue")
                
              } else if (line == "ivw"){
                generic_plot <- generic_plot + geom_abline(intercept = 0, 
                                           slope = mr_ivw(object)$Estimate, color = "blue")
                
              } else {
                generic_plot <- generic_plot
              }
              
              if (interactive == FALSE) {
                if(labels == TRUE){
                  generic_plot + geom_text(data = NULL, aes(x = Bx, y = By, label = object@snps))
                } else {
                  generic_plot
                }
                
              } else {
                p <- plotly_build(generic_plot)
                
                p$x$data[[1]]$text <- object@snps
                
                if(error == TRUE) p$x$data[[3]]$text <- NULL
                if(line %in% c("ivw", "egger")) p$x$data[[2]]$text <- NULL
               
                p$x$layout$xaxis$tickfont$size <- 13
                p$x$layout$xaxis$titlefont$size <- 15
                p$x$layout$yaxis$tickfont$size <- 13
                p$x$layout$yaxis$titlefont$size <- 15
                
                p
              }

})

#' @docType methods
#' @rdname mr_plot

setMethod(f = "mr_plot",
          signature = "MRAll",
          function(object){
            
            df <- slot(object, "Values")
            n <- nrow(df)
            
            new.df <- data.frame(matrix(NA, nrow = n, ncol = 3))
            
            if(object@Method == "median"){
              new.df <- df[,1:2]
              new.df$Intercept <- rep(0, 3)
            } else if (object@Method == "ivw" ) { 
              new.df <- df[1:4,1:2]
              new.df$Intercept <- rep(0, 4)
            } else if (object@Method == "main"){
              new.df <- df[1:4,1:2]
              new.df$Intercept <- c(rep(0, 3), df[5,2])
            } else if (object@Method == "egger"){
              new.df <- df[c(1,3,5,7),1:2]
              new.df$Intercept <- df[c(2,4,6,8),2]
            } else {
              new.df <- df[c(1:8, 10, 12, 14), 1:2]
              new.df$Intercept <- c(rep(0,7), df[c( 9, 11, 13, 15), 2])
            }
            
            
            if ( object@Method!="all" ) {
              with(data=NULL, {        
                ggplot(data = NULL, aes(x = object@Data@betaX, y = object@Data@betaY)) +
                  geom_point() +
                  geom_hline(yintercept = 0, color = "red", alpha = 0.2) +
                  geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
                  geom_abline(data = new.df, 
                              aes(intercept = Intercept, slope = Estimate, color = Method), 
                              linetype = "solid",
                              show.legend = T, size = 1) +
                  
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
              with(data=NULL, {           
                ggplot(data = NULL, aes(x = object@Data@betaX, y = object@Data@betaY)) +
                  geom_point() +
                  geom_hline(yintercept = 0, color = "red", alpha = 0.2) +
                  geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
                  
                  geom_abline(data = new.df, aes(intercept = Intercept, slope = Estimate, 
                                                 color = Method, linetype = Method),
                              show.legend = TRUE, size = 1) +
                  
                  scale_colour_manual(name="Method", labels = new.df$Method, 
                                      breaks = c("Simple median", "Weighted median", "Penalized weighted median", 
                                                 "IVW", "Penalized IVW", "Robust IVW", "Penalized robust IVW",
                                                 "MR-Egger", "Penalized MR-Egger", "Robust MR-Egger", "Penalized robust MR-Egger"),
                                      values = c("#F8766D",   "#F8766D",    "#7CAE00", "#7CAE00",  
                                                 "#C77CFF", "#C77CFF", "#00BFC4",   "#00BFC4",  
                                                 "#00BFC4",  "#F8766D",    "#7CAE00")) +
                  # IVW      MR-Egger  PenIVW   PenEgger  PenRobIVW PenRobEgg PenWMed   RobIVW   RobEgg    SimpleMed WeightedMed
                  scale_linetype_manual(name="Method", labels = new.df$Method,
                                        breaks = c("Simple median", "Weighted median", "Penalized weighted median", 
                                                   "IVW", "Penalized IVW", "Robust IVW", "Penalized robust IVW",
                                                   "MR-Egger", "Penalized MR-Egger", "Robust MR-Egger", "Penalized robust MR-Egger"),
                                        values = c("solid", "dashed", "solid", "dashed", 
                                                   "solid",  "dashed", "dotted", "solid", 
                                                   "dashed", "dotted", "dotted")) +
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


