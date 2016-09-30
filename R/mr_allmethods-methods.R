# Calculating the estimate (and relevant statistics) using all the methods specified

#' @docType methods
#' @rdname mr_allmethods

setMethod("mr_allmethods",
          "MRInput",
          function(object,
                   method = "all", ...){

            headings <- c("Method", "Estimate", "Std Error", "95% CI ", " ", "P-value")

            if(method%in%c("all", "ivw", "median", "egger", "main")) {

            if(method == "all"){

              # all functions, robust and non-robust, penalised and non-penalised

              df <- data.frame(matrix(NA, 15, 6))
              colnames(df) <- headings
              Methods <- c("Simple median",
                                "Weighted median",
                                "Penalized weighted median",

                                "IVW",
                                "Penalized IVW",
                                "Robust IVW",
                                "Penalized robust IVW",

                                "MR-Egger",
                                "(intercept)",
                                "Penalized MR-Egger",
                                "(intercept)",
                                "Robust MR-Egger",
                                "(intercept)",
                                "Penalized robust MR-Egger",
                                "(intercept)")

              df[,1] <- Methods

              # Medians
              df[1,2:6] <- values(mr_median(object, weighting = "simple", ...))
              df[2,2:6] <- values(mr_median(object, weighting = "weighted", ...))
              df[3,2:6] <- values(mr_median(object, weighting = "penalized", ...))

              # IVW
              df[4,2:6] <- values(mr_ivw(object, model = "default"))
              df[5,2:6] <- values(mr_ivw(object, model = "default", penalized = T))
              df[6,2:6] <- values(mr_ivw(object, model = "default", robust = T))
              df[7,2:6] <- values(mr_ivw(object, model = "default", robust = T, penalized = T))

              # Egger
              df[8:9,2:6] <- values(mr_egger(object))
              df[10:11,2:6] <- values(mr_egger(object, penalized = T))
              df[12:13,2:6] <- values(mr_egger(object, robust = T))
              df[14:15,2:6] <- values(mr_egger(object, robust = T, penalized = T))

            } else if (method == "median") {

              # all functions, robust and non-robust, penalised and non-penalised

              df <- data.frame(matrix(NA, nrow = 3, ncol = 6))
              colnames(df) <- headings
              Methods <- c("Simple median",
                           "Weighted median",
                           "Penalized weighted median")

              df[,1] <- Methods

              # Medians
              df[1,2:6] <- values(mr_median(object, weighting = "simple", ...)) # simple median
              df[2,2:6] <- values(mr_median(object, weighting = "weighted", ...)) # weighted median
              df[3,2:6] <- values(mr_median(object, weighting = "penalized", ...)) # penalised weighted median

            } else if (method == "egger"){

              df <- data.frame(matrix(NA, nrow = 8, ncol = 6))
              colnames(df) <- headings

              Methods <-   c("MR-Egger",
                             "(intercept)",
                             "Penalized MR-Egger",
                             "(intercept)",
                             "Robust MR-Egger",
                             "(intercept)",
                             "Penalized robust MR-Egger",
                             "(intercept)")

              df[,1] <- Methods

              # Egger
              df[1:2,2:6] <- values(mr_egger(object)) # simple egger
              df[3:4,2:6] <- values(mr_egger(object, penalized = T)) # penalised egger
              df[5:6,2:6] <- values(mr_egger(object, robust = T)) # robust egger
              df[7:8,2:6] <- values(mr_egger(object, robust = T, penalized = T)) # penalised robust egger

            } else if (method == "ivw") {

              df <- data.frame(matrix(0, nrow = 4, ncol = 6))
              colnames(df) <- headings

              Methods <- c("IVW",
                           "Penalized IVW",
                           "Robust IVW",
                           "Penalized robust IVW")
              df[,1] <- Methods

              # IVW
              df[1,2:6] <- values(mr_ivw(object, model = "default")) # random simple ivw
              df[2,2:6] <- values(mr_ivw(object, model = "default", penalized = T)) # random penalised ivw
              df[3,2:6] <- values(mr_ivw(object, model = "default", robust = T)) # random robust ivw
              df[4,2:6] <- values(mr_ivw(object, model = "default", robust = T, penalized = T)) # random penalised robust ivw
            }
             else if (method == "main") {

              df <- data.frame(matrix(0, nrow = 5, ncol = 6))
              colnames(df) <- headings

              Methods <- c("Simple median",
                           "Weighted median",
                           "IVW",
                           "MR-Egger",
                           "(intercept)")
              df[,1] <- Methods

              df[1,2:6] <- values(mr_median(object, weighting = "simple", ...)) # simple median
              df[2,2:6] <- values(mr_median(object, weighting = "weighted", ...)) # weighted median
              df[3,2:6] <- values(mr_ivw(object, model = "default")) # random simple ivw
              df[4:5,2:6] <- values(mr_egger(object)) # random simple MR-Egger
            }

            return(new("MRAll", Data = object, Values = df, Method = method))
          }
 else {  cat("Method must be one of: all, ivw, median, egger, or main.\n"); return(NULL) }
 }
)
