
# This is an internal function
setGeneric(name = "values",
           def = function(object){standardGeneric("values")})

#--------------------------------------------------------------------------------------------

#' Median-based method
#'
#' @description The \code{mr_median} function implements the weighted median (default) or simple median method introduced by Bowden et al (2016) to calculate
#' the median of the ratio instrumental variable estimates evaluated using each genetic variant individually.
#'
#' @param object An \code{MRInput} object.
#' @param weighting The type of weighting applied. The default option is to calculate the weighted median (\code{"weighted"}); other options are \code{"simple"} and \code{"penalized"}.
#' @param distribution The type of distribution to use to calculate the 95\% confidence intervals, can be \code{"normal"} or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#' @param iterations The number of bootstrap samples to generate when calculating the estimated standard error. The default value is 10000.
#' @param seed The random seed to use when generating the bootstrap samples (for reproducibility). The default value is 314159265. If set to \code{NA}, the random seed will not be set (for example, if the function is used as part of a larger simulation).
#'
#' @details The median-based methods have greater robustness to individual genetic variants with strongly outlying causal estimates
#' compared with the inverse-variance weighted and MR-Egger methods. Formally, the simple median method gives a
#' consistent estimate of the causal effect when at least 50\% of the genetic variants are valid instrumental variables
#' (for the weighted median method, when 50\% of the weight comes from valid instrumental variables).
#'
#' When the weighting is \code{"simple"}, the estimate is obtained by calculating the ratio causal estimates
#' from each genetic variants theta = betaY/betaX, and finding the median estimate.
#'
#' When the weighting is \code{"weighted"}, the estimate is obtained by:
#'
#' \enumerate{
#'  \item Calculating the ratio causal estimates and ordering the genetic variants according to the magnitude of their estimates, i.e. \deqn{\theta_1 < \theta_2 < ... < \theta_J}{theta_1 < theta_2 < ... < theta_J}
#'  \item Calculate normalized inverse-variance weights for each genetic variant \eqn{w_1, w_2, ..., w_J}, as:
#'    \deqn{w_j = \frac{\beta_{Xj}^2}{se(\beta_{Yj})^2} / \sum_{i=1}^{J} \frac{\beta_{Xi}^2}{se(\beta_{Yi})^2}}{w_j = frac{betaX_j^2}{betaYse_j^2} / sum_i frac{betaX_i^2}{betaYse_l^2}}
#'  \item Find k such that
#'    \deqn{s_k = \sum_{i = 1}^k w_i < 0.5}{s_k = w_1 + w_2 + ... + w_k < 0.5}
#'    and
#'    \deqn{s_{k+1}  = \sum_{i = 1}^{k+1} w_i > 0.5}{s_(k+1) = w_1 + w_2 +... + w_k + w_(k+1) > 0.5}
#'  \item Calculate the weighted median estimate by extrapolation as:
#'    \deqn{\theta_{WM} = \theta_k + (\theta_{k+1} - \theta_k) \times \frac{0.5 - s_k}{s_{k+1} - s_k}}{theta_(WM) = theta_k + (theta_(k+1) - theta_k) * frac{0.5 - s_k}{s_(k+1) - s_k}} }
#'
#' The simple median estimate is the same as the weighted median estimate when all the weights are equal. Standard errors
#' for both the simple and weighted median methods are calculated through bootstrapping.
#'
#' When the weighting is \code{"penalized"}, the weighted method is used, but the contribution of genetic variants with outlying (heterogeneous) ratio estimates to the analysis is downweighted.
#'
#' @return The output from the function is a \code{WeightedMedian} object containing:
#'
#'  \item{Type}{The type of weights used: \code{"weighted"}, \code{"simple"}, or \code{"penalized"}.}
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Estimate}{The value of the causal estimate.}
#'  \item{StdError}{Standard error of the causal estimate calculated using bootstrapping.}
#'  \item{CILower}{The lower bound for the causal estimate based on the estimated bootstrapped standard error and the significance level provided.}
#'  \item{CIUpper}{The upper bound for the causal estimate based on the estimated bootstrapped standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-value associated with the estimate (calculated using \code{Estimate/StdError} as per a Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'
#' @examples mr_median(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   weighting = "weighted", iterations = 1000)
#'   # iterations is set to 1000 to reduce runtime for the mr_median method,
#'   # 10000 iterations are recommended in practice
#' mr_median(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   weighting = "simple", iterations = 1000)
#' mr_median(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   weighting = "penalized", iterations = 1000)
#'
#' @references Jack Bowden, George Davey Smith, Philip C Haycock, Stephen Burgess. Consistent estimation in Mendelian randomization with
#' some invalid instruments using a weighted median estimator. Genetic Epidemiology 2016; 40(4):304-314. doi: 10.1002/gepi.21965.
#'
#' @export

setGeneric(name = "mr_median",
           def = function(object, weighting = "weighted",
                          distribution = "normal", alpha = 0.05,
                          iterations = 10000, seed = 314159265
                          )
           {standardGeneric("mr_median")})

#--------------------------------------------------------------------------------------------

#' Inverse-variance weighted method
#'
#' @description The \code{mr_ivw} function implements the inverse-variance method, informally known as the "Toby Johnson" method. With a single
#' genetic variant, this is simply the ratio method.
#'
#' @param object An \code{MRInput} object.
#' @param model What type of model should be used: \code{"default"}, \code{"random"} or \code{"fixed"}. The random-effects model (\code{"random"}) is a multiplicative random-effects model, allowing overdispersion in the weighted linear regression (the residual standard error is not fixed to be 1, but is not allowed to take values below 1). The fixed-effect model (\code{"fixed"}) sets the residual standard error to be 1. The \code{"default"} setting is to use a fixed-effect model with 3 genetic variants or fewer, and otherwise to use a random-effects model.
#' @param robust Indicates whether robust regression using the \code{lmrob()} function from the package \code{robustbase} should be used in the method rather than standard linear regression (\code{lm}).
#' @param penalized Indicates whether a penalty should be applied to the weights to downweight the contribution of genetic variants with outlying ratio estimates to the analysis.
#' @param correl If the genetic variants are correlated, then this correlation can be accounted for. The matrix of correlations between must be provided in the \code{MRInput} object: the elements of this matrix are the correlations between the individual variants (diagonal elements are 1). If a correlation matrix is specified in the \code{MRInput} object, then \code{correl} is set to \code{TRUE}. If \code{correl} is set to \code{TRUE}, then the values of \code{robust} and \code{penalized} are taken as \code{FALSE}.
#' @param distribution The type of distribution used to calculate the confidence intervals. Options are \code{"normal"} (default) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#' @param ... Additional arguments to be passed to the regression method.
#'
#' @details With multiple uncorrelated genetic variants, this estimate can be thought of as: 1) the inverse-variance
#' weighted combination of the ratio estimates from a meta-analysis; 2) the ratio estimate from combining the
#' genetic variants into a weighted score and then using this score as an instrumental variable (the same estimate
#' is obtained from the two-stage least squares method using individual-level data); 3) the coefficient from
#' weighted linear regression of the associations with the outcome on the associations with the risk factor fixing
#' the intercept to zero and using the inverse-variance weights.
#'
#' Here, we implement the method using weighted linear regression. If the variants are correlated, the method is implemented using generalized weighted linear regression; this is hard coded using matrix algebra.
#'
#' The causal estimate is obtained by regression of the associations with the outcome on the associations with the risk factor, with the intercept set to zero and weights being the inverse-variances of the associations with the outcome.
#'
#' With a single genetic variant, the estimate is the ratio of coefficients betaY/betaX and the standard error is the first term of the delta method approximation betaYse/betaX.
#'
#' @return The output from the function is an \code{IVW} object containing:
#'
#'  \item{Model}{A character string giving the type of model used (\code{"fixed"}, \code{"random"}, or \code{"default"}).}
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Correlation}{The matrix of genetic correlations.}
#'  \item{Robust}{\code{TRUE} if robust regression has been used to calculate the estimate, \code{FALSE} otherwise.}
#'  \item{Penalized}{\code{TRUE} if weights have been penalized, \code{FALSE} otherwise.}
#'  \item{Estimate}{The value of the causal estimate.}
#'  \item{StdError}{Standard error of the causal estimate.}
#'  \item{CILower}{The lower bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper}{The upper bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-value associated with the estimate (calculated as Estimate/StdError as per Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{RSE}{The estimated residual standard error from the regression model.}
#'  \item{Heter.Stat}{Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.}
#'
#' @examples mr_ivw(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))
#' mr_ivw(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   robust = TRUE)
#' mr_ivw(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   penalized = TRUE)
#' mr_ivw(mr_input(calcium, calciumse, fastgluc, fastglucse, corr=calc.rho))
#'   ## correlated variants
#'
#' @references Original implementation: The International Consortium for Blood Pressure Genome-Wide Association Studies. Genetic variants in novel pathways influence blood pressure and cardiovascular disease risk. Nature 2011; 478:103-109. doi: 10.1038/nature10405.
#'
#' Detailed description of method: Stephen Burgess, Adam S Butterworth, Simon G Thompson. Mendelian randomization analysis with multiple genetic variants using summarized data. Genetic Epidemiology 2013; 37:658-665. doi: 10.1002/gepi.21758.
#'
#' Robust and penalized weights: Stephen Burgess, Jack Bowden, Frank Dudbridge, Simon G Thompson. Robust instrumental variable methods using multiple candidate instruments with application to Mendelian randomization. arXiv 2016; 1606.03729.
#'
#' Heterogeneity test: Fabiola del Greco, Cosetta Minelli, Nuala A Sheehan, John R Thompson. Detecting pleiotropy in Mendelian randomisation studies with summary data and a continuous outcome. Stat Med 2015; 34(21):2926-2940. doi: 10.1002/sim.6522.
#'
#' @export

setGeneric(name = "mr_ivw",
           def = function(object, model = "default",
                          robust = FALSE, penalized = FALSE, correl = FALSE,
                          distribution = "normal", alpha = 0.05, ...)
           {standardGeneric("mr_ivw")})

#--------------------------------------------------------------------------------------------

#' MR-Egger method
#'
#' @description The \code{mr_egger} function implements the MR-Egger method introduced by Bowden et al (2015).
#' 
#' This method provides: 1) a test of the for directional pleiotropy (the MR-Egger intercept test), 2) a test for a
#' causal effect, and 3) an estimate of the causal effect.
#' If the intercept term differs from zero, then the genetic variants are not all valid instrumental variables and
#' the standard (inverse-variance weighted) estimate is biased. If the InSIDE (Instrument Strength Independent of Direct Effect) assumption holds, then the MR-Egger slope parameter provides a test for a causal effect, and a consistent estimate of the causal effect even if the intercept differs from zero.
#'
#' @param object An \code{MRInput} object.
#' @param robust Indicates whether robust regression using the \code{lmrob()} function from the package \code{robustbase} should be used in the method.
#' @param penalized Indicates whether a penalty should be applied to the weights to downweight the contribution of genetic variants with outlying ratio estimates to the analysis.
#' @param correl If the genetic variants are correlated, then this correlation can be accounted for. The matrix of correlations between must be provided: the elements of this matrix are the correlations between the individual variants (diagonal elements are 1). If a correlation is specified, then the values of \code{"robust"} and \code{"penalized"} are taken as \code{FALSE}.
#' @param distribution The type of distribution used to calculate the confidence intervals, can be \code{"normal"} (the default option) or \code{"t-dist"}. If the distribution is \code{"t-dist"}, then a t-distribution is used in case of over-dispersion. In case of under-dispersion, the confidence interval is the wider of that using the estimated residual standard error and a t-distribution, or that using a residual standard error of 1 and a normal distribution. This ensures that under-dispersion is not "doubly penalized" by setting the residual standard error to 1 and using a t-distribution, and also that the random-effects analysis is no more precise than a fixed-effect analysis would be.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#' @param ... Additional arguments to be passed to the regression method.
#'
#' @details The causal estimate is obtained by regression of the associations with the outcome on the associations with the risk factor, with weights being the inverse-variances of the associations with the outcome. The intercept is estimated (in contrast with the inverse-variance weighted method, where the intercept is set to zero).
#'
#' As part of the analysis, the genetic variants are orientated so that all of the associations with the risk factor are positive (and signs of associations with the outcome are changed to keep the orientation consistent if required). Re-orientation of the genetic variants is performed automatically as part of the function.
#'
#' The MR-Egger model uses a random-effects model (\code{"random"}); a fixed-effect model does not make sense as pleiotropy leads to heterogeneity between the causal estimates targeted by the genetic variants. The (multiplicative) random-effects model allows over-dispersion in the regression model. Under-dispersion is not permitted (in case of under-dispersion, the residual standard error is set to 1).
#'
#' @return The output of the function is an \code{Egger} object containing:
#'
#'  \item{Model}{A character string giving the type of model used ("random").}
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Correlation}{The matrix of genetic correlations.}
#'  \item{Robust}{\code{TRUE} if robust estimate has been calculated, \code{FALSE} otherwise.}
#'  \item{Penalized}{\code{TRUE} if weights have been penalized, \code{FALSE} otherwise.}
#'  \item{Estimate}{The value of the causal estimate (slope coefficient).}
#'  \item{StdError.Est}{Standard error of the causal estimate.}
#'  \item{Pvalue.Est}{The p-value associated with the estimate (calculated as Estimate/StdError as per Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{CILower.Est}{The lower bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper.Est}{The upper bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{Intercept}{The value of the intercept estimate.}
#'  \item{StdError.Int}{Standard error of the intercept estimate.}
#'  \item{Pvalue.Int}{The p-value associated with the intercept.}
#'  \item{CILower.Int}{The lower bound of the intercept based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper.Int}{The upper bound of the intercept based on the estimated standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals (same as \code{alpha} above).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{Causal.pval}{The p-value for the MR-Egger causal estimate.}
#'  \item{Pleio.pval}{The p-value for the MR-Egger intercept test (a low p-value suggests either directional pleiotropy or failure of the InSIDE assumption, and indicates that the IVW estimate is biased).}
#'  \item{RSE}{The estimated residual standard error from the regression model.}
#'  \item{Heter.Stat}{Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that the regression model (including an intercept) fits the regression model with no additional variability. Rejection of the null hypothesis is expected if genetic variants are pleiotropic, and doesn't mean that the MR-Egger analysis or the InSIDE assumption is invalid.}
#' \item{I.sq}{A measure of heterogeneity between the genetic associations with the exposure (see Bowden IJE 2016). Low values of \code{I.sq} relate both to large differences in precision between MR-Egger and IVW estimates, and to more weak instrument bias (in a two-sample setting, this is attenuation of MR-Egger estimate towards the null).}
#'
#'
#' @examples mr_egger(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))
#' mr_egger(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   robust = TRUE)
#' mr_egger(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   penalized = TRUE)
#' mr_egger(mr_input(calcium, calciumse, fastgluc, fastglucse, corr=calc.rho))
#'   ## correlated variants
#'
#' @references Jack Bowden, George Davey Smith, Stephen Burgess. Mendelian randomization with invalid instruments: effect estimation and bias detection through Egger regression. International Journal of Epidemiology 2015; 44:512--525. doi: 10.1093/ije/dyv080.
#'
#' Confidence intervals, and robust and penalized weights: Stephen Burgess, Jack Bowden, Frank Dudbridge, Simon G Thompson. Robust instrumental variable methods using multiple candidate instruments with application to Mendelian randomization. arXiv 2016; 1606.03729.
#'
#' I-squared statistic: Jack Bowden and others. Assessing the suitability of summary data for Mendelian randomization analyses using MR-Egger regression: The role of the I2 statistic. Int J Epidemiol 2016 (to appear).
#'
#' @export

setGeneric(name = "mr_egger",
           def = function(object, 
                          robust = FALSE, penalized = FALSE, correl = FALSE,
                          distribution = "normal", alpha = 0.05, ...)
             {standardGeneric("mr_egger")})

#--------------------------------------------------------------------------------------------

#' Mendelian randomization estimation using all methods
#'
#' The function \code{mr_allmethods} implements Mendelian randomization analyses using summarized data to calculate estimates (as well as standard
#' errors and confidence interval limits) for all the methods included in the package (or alternatively for the group of methods chosen).
#'
#' @param object An \code{MRInput} object.
#' @param method Which estimation method should be included in the calculation. By default, all estimates are computed (\code{"all"}), but one can choose to show only the results of median-based, inverse-variance weighted, or MR-Egger methods separately through specifying \code{"median"}, \code{"ivw"}, \code{"egger"}, or \code{"main"} (gives main results only, that is simple and weighted median, IVW, and MR-Egger).
#' @param ... Additional arguments to be passed to other methods.
#'
#' @details See \code{mr_median}, \code{mr_egger}, and \code{mr_ivw} for details of how each of the methods is implemented.
#'
#' @return An object of type \code{MRAll} with the following slots :
#'
#'  \item{Data}{The MRInput object used to calculate the various values.}
#'  \item{Values}{A data.frame containing the various estimates.}
#'  \item{Method}{The choice of methods estimated (default is \code{"all"}).}
#'
#' @examples mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse, 
#'   by = chdlodds, byse = chdloddsse), method="all", iterations = 1000)
#'   # iterations is set to 1000 to reduce runtime for the mr_median method,
#'   # at least 10000 iterations are recommended in practice
#' mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse,
#'   by = chdlodds, byse = chdloddsse), method="ivw")
#' mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse,
#'   by = chdlodds, byse = chdloddsse), method="egger")
#' mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse,
#'   by = chdlodds, byse = chdloddsse), method="median", iterations = 1000)
#' mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse,
#'   by = chdlodds, byse = chdloddsse), method="main", iterations = 1000)
#'
#' @references See \code{mr_median}, \code{mr_egger}, and \code{mr_ivw}.
#'
#' @export

setGeneric(name = "mr_allmethods",
           def = function(object,
                          method = "all", ...){standardGeneric("mr_allmethods")})

#--------------------------------------------------------------------------------------------

#' Draw a scatter plot of the genetic associations and/or causal estimates
#'
#' The function \code{mr_plot} has two functionalities. It can generate a visual representation of both \code{MRInput} and \code{MRAll} objects.
#'
#' @param object An \code{MRInput} object or an \code{MRAll} object.
#' @param error When viewing an \code{MRInput} object, one can choose whether to include error bars (default is to include).
#' @param line When viewing an \code{MRInput} object, one can choose whether to include the IVW estimate (\code{line = "ivw"}) or the MR-Egger estimate (\code{line = "egger"}).
#' @param orientate When viewing an \code{MRInput} object, one can choose whether to orientate all genetic variants so that the associations with the risk factor are all positive. This is recommended particularly when plotting the MR-Egger estimate, although the default setting is \code{FALSE}.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @details The result is dependent on the type of object passed to \code{mr_plot}.
#' When the object is an \code{MRInput} object, the function uses \code{plotly} syntax to plot the association estimates against eachother. This plot is interactive and the user can hover over the various points to see the name of the associated genetic variant and its association estimates.
#' When the object is an \code{MRAll} object, the function generates a \code{ggplot} to compare the causal estimates proposed by different methods.
#'
#' @examples mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   line="ivw")
#' mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   error=FALSE, line="egger", orientate = TRUE)
#' mr_plot(mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse,
#'    by = chdlodds, byse = chdloddsse), method="all", iterations = 1000))
#'   # iterations is set to 1000 to reduce runtime for the mr_median method,
#'   # 10000 iterations are recommended in practice
#' mr_plot(mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse,
#'    by = chdlodds, byse = chdloddsse), method="ivw"))
#'
#' @export

setGeneric(name = "mr_plot",
           def = function(object, error = TRUE, line = "ivw", orientate=FALSE, ...){standardGeneric("mr_plot")})

#--------------------------------------------------------------------------------------------