
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
#'   weighting = "weighted", iterations = 100)
#'   # iterations is set to 100 to reduce runtime for the mr_median method,
#'   # 10000 iterations are recommended in practice
#' mr_median(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   weighting = "simple", iterations = 100)
#' mr_median(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   weighting = "penalized", iterations = 100)
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

#' Debiased inverse-variance weighted method
#'
#' @description The \code{mr_divw} function implements the debiased inverse-variance weighted method.
#'
#' @param object An \code{MRInput} object.
#' @param over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#' @param diagnostics Should the function returns the q-q plot for assumption diagnosis. Default is FALSE.
#'
#' @details The debiased inverse-variance weighted method (dIVW) removes the weak instrument bias of the IVW method and is more robust under many weak instruments.
#'
#' @return The output from the function is a \code{DIVW} object containing:
#'
#'  \item{Over.dispersion}{\code{TRUE} if the method has considered balanced horizontal pleiotropy, \code{FALSE} otherwise.}
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Estimate}{The value of the causal estimate.}
#'  \item{StdError}{Standard error of the causal estimate calculated using bootstrapping.}
#'  \item{CILower}{The lower bound for the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper}{The upper bound for the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-value associated with the estimate (calculated using \code{Estimate/StdError} as per a Wald test) using a normal distribution.}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{Condition}{A measure (average F-statistic -1)*sqrt(# snps) that needs to be large for reliable asymptotic approximation based on the dIVW estimator. It is recommended to be greater than 20.}
#'
#' @examples mr_divw(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))
#'
#' @references Ting Ye, Jun Shao, Hyunseung Kang (2021). Debiased Inverse-Variance Weighted Estimator in
#' Two-Sample Summary-Data Mendelian Randomization. The Annals of Statistics, 49(4), 2079-2100. Also available at \url{https://arxiv.org/abs/1911.09802}.
#'
#' @export

setGeneric(name = "mr_divw",
           def = function(object, over.dispersion = TRUE, alpha = 0.05, diagnostics=FALSE)
           {standardGeneric("mr_divw")})



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
#' @param weights Which weights to use in the weighted regression. If \code{"simple"} (the default option), then the IVW estimate is equivalent to meta-analysing the ratio estimates from each variant using inverse-variance weights based on the simplest expression of the variance for the ratio estimate (first-order term from the delta expansion - standard error of the association with the outcome divided by the association with the exposure). If \code{"delta"}, then the variance expression is the second-order term from the delta expansion. The second-order term incorporates uncertainty in the genetic association with the exposure -- this uncertainty is ignored using the simple weighting.
#' @param psi The correlation between the genetic associations with the exposure and the association with the outcome for each variant resulting from sample overlap. The default value is \code{0}, corresponding to a strict two-sample Mendelian randomization analysis (no overlap). If there is complete overlap between the samples, then the correlation should be set to the observational correlation between the exposure and the outcome. This correlation is only used in the calculation of standard errors if the option \code{weights} is set to \code{"delta"}.
#' @param correl If the genetic variants are correlated, then this correlation can be accounted for. The matrix of correlations between must be provided in the \code{MRInput} object: the elements of this matrix are the correlations between the individual variants (diagonal elements are 1). If a correlation matrix is specified in the \code{MRInput} object, then \code{correl} is set to \code{TRUE}. If \code{correl} is set to \code{TRUE}, then the values of \code{robust} and \code{penalized} are taken as \code{FALSE}, and \code{weights} is set to \code{"simple"}.
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
#' Simple versus delta weights (first-order versus second-order): Stephen Burgess, Jack Bowden. Integrating summarized data from multiple genetic variants in Mendelian randomization: bias and coverage properties of inverse-variance weighted methods. arXiv:1512.04486.
#'
#' @export

setGeneric(name = "mr_ivw",
           def = function(object, model = "default",
                          robust = FALSE, penalized = FALSE, weights = "simple", psi = 0, correl = FALSE,
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
#'  \item{Model}{A character string giving the type of model used (\code{"random"}).}
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
#' @description The function \code{mr_allmethods} implements Mendelian randomization analyses using summarized data to calculate estimates (as well as standard
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
#'   by = chdlodds, byse = chdloddsse), method="main", iterations = 100)
#'   # iterations is set to 100 to reduce runtime for the mr_median method,
#'   # at least 10000 iterations are recommended in practice
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
#' The function \code{mr_plot} has three functionalities. It can generate a visual representation of \code{MRInput}, \code{MRMVInput}  and \code{MRAll} objects.
#'
#' @param object An \code{MRInput} object or an \code{MRMVInput} object or an \code{MRAll} object.
#' @param error When viewing an \code{MRInput} or \code{MRMVInput} object, one can choose whether to include error bars (default is to include). For an \code{MRMVInput} object, the horizontal error bars only take into account uncertainty in the causal estimates.
#' @param line When viewing an \code{MRInput} object, one can choose whether to include the IVW estimate (\code{line = "ivw"}) or the MR-Egger estimate (\code{line = "egger"}). When viewing an \code{MRMVInput}, one can choose whether to include a line through the origin with gradient 1 (\code{line = TRUE}) or not.
#' @param orientate When viewing an \code{MRInput} or \code{MRMVInput} object, one can choose whether to orientate all genetic variants so that the associations with the risk factor are all positive. This is recommended particularly when plotting the MR-Egger estimate, although the default setting is \code{FALSE}.
#' @param interactive When viewing an \code{MRInput} or \code{MRMVInput} object, one can choose whether to produce an interactive graph using the \code{plotly} package, or a static graph using the regular \code{plot} command.
#' @param labels When viewing an \code{MRInput} or \code{MRMVInput} object with \code{interactive} set to \code{FALSE}, setting \code{labels} to \code{TRUE} means that the name of each genetic variants appears above the corresponding datapoint.
#'
#' @details The result is dependent on the type of object passed to \code{mr_plot}.
#' When the object is an \code{MRInput} object, the function uses either the \code{plot} command (if \code{interactive} is set to \code{FALSE}) or \code{plotly} syntax (if \code{interactive} is set to \code{TRUE}) to plot the association estimates against each other.
#' When the object is an \code{MRMVInput} object, functionality is similar except that we plot the estimated associations with the outcome on the y-axis, and fitted values of the associations with the outcome from the inverse-variance weighted method on the x-axis.
#' If \code{interactive} is set to \code{FALSE}, then a static graph is produced. By setting \code{labels} to \code{TRUE}, the names of the genetic variants appear above the points. This produces a less visually appealing graph, but one where it is easier to identify the individual genetic variants.
#' If \code{interactive} is set to \code{TRUE}, then the plot is interactive and the user can hover over the various points to see the name of the associated genetic variant and its association estimates.
#' When the object is an \code{MRAll} object, the function generates a \code{ggplot} to compare the causal estimates proposed by different methods.
#'
#' @examples mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   line="egger", orientate = TRUE)
#' mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   line="ivw", interactive=FALSE) # produces a static graph
#' mr_plot(mr_allmethods(mr_input(bx = ldlc, bxse = ldlcse,
#'    by = chdlodds, byse = chdloddsse), method="all", iterations = 50))
#'   # iterations is set to 50 to reduce runtime for the mr_median method,
#'   # 10000 iterations are recommended in practice
#'
#' @export

setGeneric(name = "mr_plot",
           def = function(object, error = TRUE, line = "ivw", orientate=FALSE, interactive = TRUE, labels = FALSE){standardGeneric("mr_plot")})

#--------------------------------------------------------------------------------------------

#' Maximum-likelihood method
#'
#' @description The \code{mr_maxlik} function implements the maximum-likelihood method introduced by Burgess et al (2013).
#'
#' @param object An \code{MRInput} object.
#' @param model What type of model should be used: \code{"default"}, \code{"random"} or \code{"fixed"}. The method naturally estimates a fixed-effect model, assuming that the same causal effect is estimated by each of the genetic variants. However, if there is heterogeneity in the causal estimates of the different variants, then confidence intervals under a fixed-effect model will be overly narrow. The random-effects model adds additional uncertainty by multiplying the standard error by the square-root of the likelihood ratio heterogeneity statistic divided by the number of genetic variants less one (unless this quantity is less than 1, in which case no modification to the standard error is made). This parallels the residual standard error in a regression model (the Cochran Q heterogeneity test statistic is equal to the square of the RSE multiplied by the number of genetic variants less one). The default setting (\code{"default"}) is to use a fixed-effect model with 3 genetic variants or fewer, and otherwise to use a random-effects model.
#' @param correl If the genetic variants are correlated, then this correlation can be accounted for. The matrix of correlations between must be provided in the \code{MRInput} object: the elements of this matrix are the correlations between the individual variants (diagonal elements are 1).
#' @param psi The correlation between the association with the exposure and the association with the outcome for each variant resulting from sample overlap.
#' @param distribution The type of distribution used to calculate the confidence intervals, can be \code{"normal"} (the default option) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#' @param ... Additional arguments to be passed to the optimization method.
#'
#' @details A likelihood function is defined by assuming that the summarized data for each genetic variant are normally distributed. A bivariate normal distribution is assumed for the associations of each genetic variant with the exposure and with the outcome. The mean of the association with the outcome is taken as the mean association with the exposure multiplied by the causal effect parameter.
#'
#' Thus, if there are \code{K} genetic variants, then \code{K+1} parameters are estimated by the method: one for each gene--exposure association, plus the causal parameter. If the number of genetic variants is large, then maximization of this function may be an issue. If the maximum likelihood estimate substantially differs from the inverse-variance weighted estimate, this may indicate that convergence has not occurred in the optimization algorithm.
#'
#' The variance-covariance matrices for the bivariate normal distributions are obtained from the standard error estimates provided. The correlation \code{psi} between genetic associations with the exposure and with the outcome due to sample overlap can be specified; its default value is zero.
#'
#' Two features why this method may be preferred over the inverse-variance weighted method are the incorporation in the model of uncertainty in the genetic associations with the exposure, and of correlation between the genetic association estimates with exposure and outcome for each variant. The method is implemented both for uncorrelated and correlated genetic variants. It can also be used for a single genetic variant.
#'
#' The original version of the maximum-likelihood method assumed that all genetic variants identify the same causal estimate; a fixed-effect model. The causal estimate may be overly precise if the fixed-effect model is incorrect and there is substantial heterogeneity in the causal estimates from the different variants. The random-effects analysis implemented here is an ad hoc solution to the problem of heterogeneity, but one that should result in reasonable confidence intervals that incorporate this heterogeneity.
#'
#' @return The output from the function is an \code{MaxLik} object containing:
#'
#'  \item{Model}{A character string giving the type of model used (\code{"fixed"}, \code{"random"}, or \code{"default"}).}
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Correlation}{The matrix of genetic correlations.}
#'  \item{Psi}{The correlation between genetic associations with the exposure and with the outcome.}
#'  \item{Estimate}{The value of the causal estimate.}
#'  \item{StdError}{Standard error of the causal estimate.}
#'  \item{CILower}{The lower bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper}{The upper bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-value associated with the estimate (calculated as Estimate/StdError as per Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{RSE}{The estimated residual standard error from the regression model (always equal to 1, as a fixed-effect model is required.}
#'  \item{Heter.Stat}{Heterogeneity statistic (likelihood ratio statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.}
#'
#' @examples mr_maxlik(mr_input(bx = ldlc[1:10], bxse = ldlcse[1:10],
#'   by = chdlodds[1:10], byse = chdloddsse[1:10]))
#' mr_maxlik(mr_input(bx = ldlc[1:10], bxse = ldlcse[1:10],
#'   by = chdlodds[1:10], byse = chdloddsse[1:10]), psi=0.2)
#' mr_maxlik(mr_input(calcium, calciumse, fastgluc, fastglucse, corr=calc.rho))
#'   ## correlated variants
#'
#' @references Stephen Burgess, Adam S Butterworth, Simon G Thompson. Mendelian randomization analysis with multiple genetic variants using summarized data. Genetic Epidemiology 2013; 37:658-665. doi: 10.1002/gepi.21758.
#'
#' @export

setGeneric(name = "mr_maxlik",
           def = function(object, model = "default",
                          correl = FALSE, psi = 0,
                          distribution = "normal", alpha = 0.05, ...)
             {standardGeneric("mr_maxlik")})

#--------------------------------------------------------------------------------------------

#' Multivariable inverse-variance weighted method
#'
#' @description The \code{mr_mvivw} function performs multivariable Mendelian randomization via the inverse-variance method. This is implemented by multivariable weighted linear regression.
#'
#' @param object An \code{MRMVInput} object.
#' @param model What type of model should be used: \code{"default"}, \code{"random"} or \code{"fixed"}. The random-effects model (\code{"random"}) is a multiplicative random-effects model, allowing overdispersion in the weighted linear regression (the residual standard error is not fixed to be 1, but is not allowed to take values below 1). The fixed-effect model (\code{"fixed"}) sets the residual standard error to be 1. The \code{"default"} setting is to use a fixed-effect model with 3 genetic variants or fewer, and otherwise to use a random-effects model.
#' @param robust Indicates whether robust regression using the \code{lmrob()} function from the package \code{robustbase} should be used in the method rather than standard linear regression (\code{lm}).
#' @param correl If the genetic variants are correlated, then this correlation can be accounted for. The matrix of correlations between must be provided in the \code{MRMVInput} object: the elements of this matrix are the correlations between the individual variants (diagonal elements are 1). If a correlation matrix is specified in the \code{MRMVInput} object, then \code{correl} is set to \code{TRUE}.
#' @param distribution The type of distribution used to calculate the confidence intervals. Options are \code{"normal"} (default) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#' @param ... Additional arguments to be passed to the regression method.
#'
#' @details Multivariable Mendelian randomization is an extension of Mendelian randomization to deal with genetic variants that are associated with multiple risk factors. Two scenarios are envisioned for its use: 1) risk factors that are biologically related, such as lipid fractions; and 2) risk factors where there is potentially a network of causal effects (mediation) from one risk factor to another. In both cases, under the extended assumptions of multivariable Mendelian randomization, coefficients represent the direct causal effects of each risk factor in turn with the other risk factors being fixed.
#'
#' We implement the method using multivariable weighted linear regression. If the variants are correlated, the method is implemented using generalized weighted linear regression; this is hard coded using matrix algebra.
#'
#' The causal estimate is obtained by regression of the associations with the outcome on the associations with the risk factors, with the intercept set to zero and weights being the inverse-variances of the associations with the outcome.
#'
#' @return The output from the function is an \code{MVIVW} object containing:
#'
#'  \item{Model}{A character string giving the type of model used (\code{"fixed"}, \code{"random"}, or \code{"default"}).}
#'  \item{Exposure}{A character vector with the names given to the exposure.}
#'  \item{Outcome}{A character string with the names given to the outcome.}
#'  \item{Robust}{\code{TRUE} if robust regression has been used to calculate the estimate, \code{FALSE} otherwise.}
#'  \item{Correlation}{The matrix of genetic correlations.}
#'  \item{Estimate}{A vector of causal estimates.}
#'  \item{StdError}{A vector of standard errors of the causal estimates.}
#'  \item{CILower}{The lower bounds of the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{CIUpper}{The upper bounds of the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-values associated with the estimates (calculated as Estimate/StdError as per Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{RSE}{The estimated residual standard error from the regression model.}
#'  \item{Heter.Stat}{Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.}
#'
#' @examples mr_mvivw(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse),
#'    by = chdlodds, byse = chdloddsse))
#'
#' @references Description of approach: Stephen Burgess, Simon G Thompson. Multivariable Mendelian Randomization: the use of pleiotropic genetic variants to estimate causal effects. American Journal of Epidemiology 2015; 181(4):251-260. doi: 10.1093/aje/kwu283.
#'
#' Description of inverse-variance weighted method: Stephen Burgess, Frank Dudbridge, Simon G Thompson. Re: "Multivariable Mendelian randomization: the use of pleiotropic genetic variants to estimate causal effects." American Journal of Epidemiology 2015; 181(4):290-291. doi: 10.1093/aje/kwv017.
#'
#' Use for mediation analysis: Stephen Burgess, Deborah J Thompson, Jessica MB Rees, Felix R Day, John R Perry, Ken K Ong. Dissecting causal pathways using Mendelian randomization with summarized genetic data: Application to age at menarche and risk of breast cancer. Genetics 2017; 207(2):481-487. doi: 10.1534/genetics.117.300191.
#'
#' @export

setGeneric(name = "mr_mvivw",
           def = function(object, model = "default", robust=FALSE,
                          correl = FALSE,
                          distribution = "normal", alpha = 0.05, ...)
           {standardGeneric("mr_mvivw")})

#--------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------

#' Multivariable MR-Egger method
#'
#' @description The \code{mr_mvegger} function performs multivariable Mendelian randomization via the MR-Egger method. This is implemented by multivariable weighted linear regression.
#'
#' @param object An \code{MRMVInput} object.
#' @param orientate The risk factor that genetic associations are orientated to. The univariable and multivariable versions of MR-Egger are both sensitive to the choice of parameterization of the genetic associations - which allele the associations are orientated with respect to (in other words, which allele is the effect allele). For univariable MR-Egger, this is resolved by setting the genetic associations with the exposure all to be positive. In multivariable MR-Egger, we have to choose which of the exposures to orientate the genetic associations to. The default option is \code{1}, meaning that genetic associations with the first exposure are set to be positive.
#' @param correl If the genetic variants are correlated, then this correlation can be accounted for. The matrix of correlations between must be provided in the \code{MRInput} object: the elements of this matrix are the correlations between the individual variants (diagonal elements are 1). If a correlation matrix is specified in the \code{MRInput} object, then \code{correl} is set to \code{TRUE}.
#' @param distribution The type of distribution used to calculate the confidence intervals. Options are \code{"normal"} (default) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#'
#' @details Multivariable MR-Egger is an extension of the MR-Egger method to deal with genetic variants that are associated with multiple risk factors.
#'
#' We implement the method using multivariable weighted linear regression. If the variants are correlated, the method is implemented using generalized weighted linear regression; this is hard coded using matrix algebra.
#'
#' The causal estimate is obtained by regression of the associations with the outcome on the associations with the risk factors, with the intercept estimated and weights being the inverse-variances of the associations with the outcome.
#'
#' @return The output from the function is an \code{MVEgger} object containing:
#'
#'  \item{Model}{A character string giving the type of model used (\code{"random"}).}
#'  \item{Orientate}{The number corresponding to the risk factor that the genetic associations are orientated to.}
#'  \item{Exposure}{A character vector with the names given to the exposure.}
#'  \item{Outcome}{A character string with the names given to the outcome.}
#'  \item{Correlation}{The matrix of genetic correlations.}
#'  \item{Estimate}{A vector of the causal estimates (slope coefficient).}
#'  \item{StdError.Est}{Standard errors of the causal estimates.}
#'  \item{Pvalue.Est}{The p-values associated with the estimates using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{CILower.Est}{The lower bound of the causal estimates based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper.Est}{The upper bound of the causal estimates based on the estimated standard error and the significance level provided.}
#'  \item{Intercept}{The value of the intercept estimate.}
#'  \item{StdError.Int}{Standard error of the intercept estimate.}
#'  \item{Pvalue.Int}{The p-value associated with the intercept.}
#'  \item{CILower.Int}{The lower bound of the intercept based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper.Int}{The upper bound of the intercept based on the estimated standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-values associated with the estimates (calculated as Estimate/StdError as per Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{RSE}{The estimated residual standard error from the regression model.}
#'  \item{Heter.Stat}{Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.}
#'
#' @examples mr_mvegger(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse),
#'    by = chdlodds, byse = chdloddsse), orientate = 1)
#'
#' @references Jessica Rees, Angela Wood, Stephen Burgess. Extending the MR?Egger method for multivariable Mendelian randomization to correct for both measured and unmeasured pleiotropy. Statistics in Medicine 2017; 36(29): 4705-4718. doi: 10.1002/sim.7492.
#'
#' @export

setGeneric(name = "mr_mvegger",
           def = function(object, orientate = 1,
                          correl = FALSE,
                          distribution = "normal", alpha = 0.05)
           {standardGeneric("mr_mvegger")})

#--------------------------------------------------------------------------------------------


#' Mode-based method of Hartwig
#'
#' @description The \code{mr_mbe} function implements the mode-based method introduced by Hartwig, Bowden and Davey Smith (2017).
#'
#' @param object An \code{MRInput} object.
#' @param weighting Whether the analysis should be \code{"weighted"} (the default option) or \code{"unweighted"}.
#' @param stderror Whether standard error estimates should be i) \code{"simple"} - calculated as the first-order term from the delta expansion - standard error of the association with the outcome divided by the association with the exposure), or ii) \code{"delta"} - calculated as the second-order term from the delta expansion (the default option). The second-order term incorporates uncertainty in the genetic association with the exposure -- this uncertainty is ignored using the simple weighting. The \code{"simple"} option is referred to by Hartwig et al as "assuming NOME", and the \code{"delta"} option as "not assuming NOME".
#' @param phi The choice of bandwidth in the kernel-smoothly density method. A value of 1 (the default value) represents the bandwidth value selected by the modified Silverman's bandwidth rule, as recommended by Hartwig et al. A value of 0.5 represents half that value, and so on.
#' @param seed The random seed to use when generating the bootstrap samples used to calculate the confidence intervals (for reproducibility). The default value is 314159265. If set to \code{NA}, the random seed will not be set (for example, if the function is used as part of a larger simulation).
#' @param iterations Number of iterations to use in the bootstrap procedure.
#' @param distribution The type of distribution used to calculate the confidence intervals, can be \code{"normal"} (the default option) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#'
#' @details The mode-based estimation (MBE) method takes the variant-specific ratio estimates from each genetic variant in turn, and calculates the modal estimate. This is implemented by constructing a kernel-smoothed density out of the ratio estimates, and taking the maximum value as the modal estimate. The standard error is calculated by a bootstrap procedure, and confidence intervals based on the estimate having a normal distribution.
#'
#' The method should give consistent estimates as the sample size increases if a plurality (or weighted plurality) of the genetic variants are valid instruments. This means that the largest group of variants with the same causal estimate in the asymptotic limit are the valid instruments.
#'
#' @return The output from the function is an \code{MRMBE} object containing:
#'
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Weighting}{A character string \code{"weighted"} or \code{"unweighted"}.}
#'  \item{StdErr}{A character string \code{"simple"} or \code{"delta"}.}
#'  \item{Phi}{The value of the bandwidth factor.}
#'  \item{Estimate}{The value of the causal estimate.}
#'  \item{StdError}{Standard error of the causal estimate.}
#'  \item{CILower}{The lower bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper}{The upper bound of the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-value associated with the estimate (calculated as Estimate/StdError as per Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'
#' @examples mr_mbe(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse), iterations=100)
#' mr_mbe(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'    phi=0.5, iterations=100)
#' mr_mbe(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'    weighting="weighted", stderror="delta", iterations=100)
#' # iterations set to 100 to reduce computational time,
#' #  more iterations are recommended in practice
#'
#' @references Fernando Pires Hartwig, George Davey Smith, Jack Bowden. Robust inference in summary data Mendelian randomization via the zero modal pleiotropy assumption. International Journal of Epidemiology 2017; 46(6): 1985-1998. doi: 10.1093/ije/dyx102.
#'
#' @export

setGeneric(name = "mr_mbe",
           def = function(object, weighting = "weighted", stderror = "simple", phi = 1, seed = 314159265, iterations = 10000,
                          distribution = "normal", alpha = 0.05)
             {standardGeneric("mr_mbe")})

#--------------------------------------------------------------------------------------------

#' Heterogeneity-penalized method
#'
#' @description Heterogeneity-penalized model-averaging method for efficient modal-based estimation.
#'
#' @param object An \code{MRInput} object.
#' @param prior The prior probability of a genetic variant being a valid instrument (default is 0.5).
#' @param CIMin The smallest value to use in the search to find the confidence interval (default is -1).
#' @param CIMax The largest value to use in the search to find the confidence interval (default is +1).
#' @param CIStep The step size to use in the search to find the confidence interval (default is 0.001). The confidence interval is determined by a grid search algorithm. Using the default settings, we calculate the likelihood at all values from -1 to +1 increasing in units of 0.001. If this range is too large or the step size is too small, then the grid search algorithm will take a long time to converge.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#'
#' @details This method was developed as a more efficient version of the mode-based estimation method of Hartwig et al. It proceeds by evaluating weights for all subsets of genetic variants (excluding the null set and singletons). Subsets receive greater weight if they include more variants, but are severely downweighted if the variants in the subset have heterogeneous causal estimates. As such, the method will identify the subset with the largest number (by weight) of variants having similar causal estimates.
#'
#' Confidence intervals are evaluated by calculating a log-likelihood function, and finding all points within a given vertical distance of the maximum of the log-likelihood function (which is the causal estimate). As such, if the log-likelihood function is multimodal, then the confidence interval may include multiple disjoint ranges. This may indicate the presence of multiple causal mechanisms by which the exposure may influence the outcome with different magnitudes of causal effect. As the confidence interval is determined by a grid search, care must be taken when chosing the minimum (\code{CIMin}) and maximum (\code{CIMax}) values in the search, as well as the step size (\code{CIStep}). The default values will not be suitable for all applications.
#'
#' The method should give consistent estimates as the sample size increases if a weighted plurality of the genetic variants are valid instruments. This means that the largest group of variants with the same causal estimate in the asymptotic limit are the valid instruments.
#'
#' The current implementation of the method evaluates a weight and an estimate for each of the subsets of genetic variants. This means that the method complexity doubles for each additional genetic variant included in the analysis. Currently, the method provides a warning message when used with 25+ variants, and fails to run with 30+.
#'
#' @return The output from the function is an \code{MRHetPen} object containing:
#'
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Prior}{The value of the bandwidth factor.}
#'  \item{Estimate}{The value of the causal estimate.}
#'  \item{CIRange}{The range of values in the confidence interval based on a grid search between the minimum and maximum values for the causal effect provided.}
#'  \item{CILower}{The lower limit of the confidence interval. If the confidence interval contains multiple ranges, then lower limits of all ranges will be reported.}
#'  \item{CIUpper}{The upper limit of the confidence interval. If the confidence interval contains multiple ranges, then upper limits of all ranges will be reported.}
#'  \item{CIMin}{The smallest value used in the search to find the confidence interval.}
#'  \item{CIMax}{The largest value used in the search to find the confidence interval.}
#'  \item{CIStep}{The step size used in the search to find the confidence interval.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'
#' @examples mr_hetpen(mr_input(bx = ldlc[1:10], bxse = ldlcse[1:10], by = chdlodds[1:10],
#'    byse = chdloddsse[1:10]), CIMin = -1, CIMax = 5, CIStep = 0.01)
#'
#' @references Stephen Burgess, Verena Zuber, Apostolos Gkatzionis, Christopher N Foley. Improving on a modal-based estimation method: model averaging for consistent and efficient estimation in Mendelian randomization when a plurality of candidate instruments are valid. bioRxiv 2017. doi: 10.1101/175372.
#'
#' @export

setGeneric(name = "mr_hetpen",
           def = function(object, prior=0.5, CIMin=-1, CIMax=1, CIStep=0.001, alpha = 0.05)
             {standardGeneric("mr_hetpen")})

#--------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------

#' Contamination mixture method
#'
#' @description Contamination mixture method for robust and efficient estimation under the 'plurality valid' assumption.
#'
#' @param object An \code{MRInput} object.
#' @param psi The value of the standard deviation of the distribution of invalid estimands (default value is 0, corresponding to 1.5 times the standard deviation of the ratio estimates).
#' @param CIMin The smallest value to use in the search to find the confidence interval. The default value is NA, which means that the method uses the smallest value of the lower bound of the 95\% confidence interval for the variant-specific ratio estimates as the smallest value.
#' @param CIMax The largest value to use in the search to find the confidence interval. The default value is NA, which means that the method uses the greatest value of the upper bound of the 95\% confidence interval for the variant-specific ratio estimates as the largest value.
#' @param CIStep The step size to use in the search to find the confidence interval (default is 0.01). The confidence interval is determined by a grid search algorithm. Using the default settings, we calculate the likelihood at all values from -1 to +1 increasing in units of 0.01. If this range is too large or the step size is too small, then the method will take a long time to run.
#' @param alpha The significance level used to calculate the confidence interval. The default value is 0.05.
#'
#' @details The contamination mixture method is implemented by constructing a likelihood function based on the variant-specific causal estimates. If a genetic variant is a valid instrument, then its causal estimate will be normally distributed about the true value of the causal effect. If a genetic variant is not a valid instrument, then its causal estimate will be normally distributed about some other value. We assume that the values estimated by invalid instruments are normally distributed about zero with a large standard deviation. This enables a likelihood function to be specified that is a product of two-component mixture distributions, with one mixture distribution for each variant. The computational time for maximizing this likelihood directly is exponential in the number of genetic variants. We use a profile likelihood approach to reduce the computational complexity to be linear in the number of variants.
#'
#' We consider different values of the causal effect in turn. For each value, we calculate the contribution to the likelihood for each genetic variant as a valid instrument and as an invalid instrument. If the contribution to the likelihood as a valid instrument is greater, then we take the variant's contribution as a valid instrument; if less, then its contribution is taken as an invalid instrument. This gives us the configuration of valid and invalid instruments that maximizes the likelihood for the given value of the causal effect. This is a profile likelihood, a one-dimensional function of the causal effect. The point estimate is then taken as the value of the causal effect that maximizes the profile likelihood.
#'
#' Confidence intervals are evaluated by calculating the log-likelihood function, and finding all points within a given vertical distance of the maximum of the log-likelihood function (which is the causal estimate). As such, if the log-likelihood function is multimodal, then the confidence interval may include multiple disjoint ranges. This may indicate the presence of multiple causal mechanisms by which the exposure may influence the outcome with different magnitudes of causal effect. As the confidence interval is determined by a grid search, care must be taken when chosing the minimum (\code{CIMin}) and maximum (\code{CIMax}) values in the search, as well as the step size (\code{CIStep}). The default values will not be suitable for all applications.
#'
#'
#' @return The output from the function is an \code{MRConMix} object containing:
#'
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Psi}{The value of the standard deviation parameter.}
#'  \item{Estimate}{The value of the causal estimate.}
#'  \item{CIRange}{The range of values in the confidence interval based on a grid search between the minimum and maximum values for the causal effect provided.}
#'  \item{CILower}{The lower limit of the confidence interval. If the confidence interval contains multiple ranges, then lower limits of all ranges will be reported.}
#'  \item{CIUpper}{The upper limit of the confidence interval. If the confidence interval contains multiple ranges, then upper limits of all ranges will be reported.}
#'  \item{CIMin}{The smallest value used in the search to find the confidence interval.}
#'  \item{CIMax}{The largest value used in the search to find the confidence interval.}
#'  \item{CIStep}{The step size used in the search to find the confidence interval.}
#'  \item{Pvalue}{The p-value associated with the estimate calculated using the likelihood function and a chi-squared distribution.}
#'  \item{Valid}{The numbers of genetic variants that were considered valid instruments at the causal estimate.}
#'  \item{ValidSNPs}{The names of genetic variants that were considered valid instruments at the causal estimate.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'
#' @examples mr_conmix(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds,
#'    byse = chdloddsse), psi = 3, CIMin = -1, CIMax = 5, CIStep = 0.01)
#'
#' @references Stephen Burgess, Christopher N Foley, Elias Allara, Joanna Howson. A robust and efficient method for Mendelian randomization with hundreds of genetic variants: unravelling mechanisms linking HDL-cholesterol and coronary heart disease. Nat Comms 2020. doi: 10.1038/s41467-019-14156-4.
#'
#' @export

setGeneric(name = "mr_conmix",
           def = function(object, psi=0, CIMin=NA, CIMax=NA, CIStep=0.01, alpha = 0.05)
             {standardGeneric("mr_conmix")})

#--------------------------------------------------------------------------------------------

#' Multivariable median-based method
#'
#' @description The \code{mr_mvmedian} function performs multivariable Mendelian randomization via the median method. This is implemented by multivariable weighted quantile regression, with the quantile set to 0.5 (the median).
#'
#' @param object An \code{MRMVInput} object.
#' @param distribution The type of distribution used to calculate the confidence intervals. Options are \code{"normal"} (default) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#' @param iterations The number of bootstrap samples to generate when calculating the estimated standard error. The default value is 10000.
#' @param seed The random seed to use when generating the bootstrap samples (for reproducibility). The default value is 314159265. If set to \code{NA}, the random seed will not be set (for example, if the function is used as part of a larger simulation).
#'
#' @details The multivariable median method is similar to the univariable weighted median method, except that it is implemented using quantile regression. The regression model is multivariable and weighted by the inverse of the variances of the variant-specific estimates. Confidence intervals are calculated by parametric bootstrap to estimate the standard error of the estimates, and then using quantiles of a normal or t-distribution (depending on the value of \code{distribution}).
#'
#' @return The output from the function is an \code{MVMedian} object containing:
#'
#'  \item{Exposure}{A character vector with the names given to the exposure.}
#'  \item{Outcome}{A character string with the names given to the outcome.}
#'  \item{Estimate}{A vector of causal estimates.}
#'  \item{StdError}{A vector of standard errors of the causal estimates.}
#'  \item{CILower}{The lower bounds of the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{CIUpper}{The upper bounds of the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-values associated with the estimates (calculated as Estimate/StdError as per Wald test) using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'
#'
#' @examples mr_mvmedian(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse),
#'    by = chdlodds, byse = chdloddsse), iterations = 100)
#'   # iterations is set to 100 to reduce runtime for the mr_mvmedian method,
#'   # 10000 iterations are recommended in practice
#'
#' @export

setGeneric(name = "mr_mvmedian",
            def = function(object,
                           distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265)
            {standardGeneric("mr_mvmedian")})

#--------------------------------------------------------------------------------------------

#' Draw a forest plot of causal estimates
#'
#' @description The \code{mr_forest} function draws a forest plot of causal estimates. The default option plots the variant-specific causal estimates (\code{by/bx}) and the estimate from the \code{mr_ivw} function using default settings (assuming variants are uncorrelated, random-effects for 4+ variants). Options allow users to plot estimates from a variety of different methods.
#'
#' @param object An \code{MRInput} object.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05, corresponding to 95\% confidence intervals.
#' @param snp_estimates Whether to plot the variant-specific estimates. Defaults to \code{TRUE}.
#' @param methods Takes a string of computation methods used to calculate estimates. Defaults to \code{"ivw"}. Options are: \code{"median"} (simple median estimate), \code{"wmedian"} (weighted median estimate), \code{"egger"} (MR-Egger estimate), \code{"mbe"} (mode-based estimate), \code{"conmix"} (contamination mixture estimate), and \code{"maxlik"} (maximum likelihood estimate).
#' @param ordered Determines by whether to arrange the variant-specific estimates in ascending order. Defaults to \code{FALSE}.
#'
#' @details As the function produces a \code{ggplot} object, graphical parameters can be changed by adding commands from the \code{ggplot2} package.
#'
#' @examples mr_forest(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   alpha = 0.01, ordered = TRUE)
#' mr_forest(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE)
#' forest = mr_forest(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))
#' # how to change x-axis limits
#' # library(ggplot2)
#' # forest2 = forest + coord_cartesian(xlim=c(-5,5))
#' # forest2
#'
#' @export



setGeneric(name = "mr_forest",
           def = function(object, alpha = 0.05, snp_estimates = TRUE, methods = "ivw", ordered = FALSE){standardGeneric("mr_forest")})

#--------------------------------------------------------------------------------------------

#' Leave-one-out estimates
#'
#' @description The \code{mr_loo} function draws a forest plot of causal estimates from the \code{mr_ivw} function using default settings (assuming variants are uncorrelated, random-effects for 4+ variants) omitting each variant in turn. So the estimate labelled \code{snp_1} includes all variants except the labelled variant, and so on. The \code{mr_ivw} estimate including all variants ("IVW estimate") is also provided for reference.
#'
#' @param object An \code{MRInput} object.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05, corresponding to 95\% confidence intervals.
#'
#' @details As the function produces a \code{ggplot} object, graphical parameters can be changed by adding commands from the \code{ggplot2} package.
#'
#' @examples mr_loo(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#'   alpha = 0.01)
#'
#' @export

setGeneric(name = "mr_loo",
           def = function(object, alpha = 0.05){standardGeneric("mr_loo")})

#--------------------------------------------------------------------------------------------

#' Draw a funnel plot of variant-specific estimates
#'
#' @description The \code{mr_funnel} function draws a funnel plot of variant-specific causal estimates. Estimates (\code{by/bx}) are plotted against the precision of the estimates (\code{abs(bx)/byse}). Precision is the reciprocal of standard error. A vertical dashed line is plotted at the estimate from the \code{mr_ivw} function.
#'
#' @param object An \code{MRInput} object.
#' @param CI A \code{logical} variable dicating as to whether to plot the confidence interval associated with each point. Default value is TRUE.
#'
#' @details As the function produces a \code{ggplot} object, graphical parameters can be changed by adding commands from the \code{ggplot2} package.
#'
#' @examples mr_funnel(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))
#'
#' @export

setGeneric(name = "mr_funnel",
           def = function(object, CI = TRUE){standardGeneric("mr_funnel")})

#--------------------------------------------------------------------------------------------

#' Multivariable MR-Lasso method
#'
#' @description The \code{mr_mvlasso} function performs the multivariable MR-Lasso method, which applies lasso-type penalization to the direct effects of genetic variants on the outcome.
#' The causal estimates are described as post-lasso estimates, and are obtained by performing the multivariable IVW method using only those genetic variants that are identified as valid by the lasso procedure.
#'
#' @param object An \code{MRMVInput} object.
#' @param orientate The risk factor that genetic associations are orientated to. The default option is \code{1}, meaning that genetic associations with the first risk factor are set to be positive.
#' @param distribution The type of distribution used to calculate the confidence intervals. Options are \code{"normal"} (default) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#' @param lambda The value of the tuning parameter used by the lasso procedure which controls the level of sparsity. If not specified, the tuning parameter will be calculated by the heterogeneity stopping rule.
#'
#' @details Multivariable MR-Lasso extends the multivariable IVW model to include an intercept term for each genetic variant. These intercept terms represent associations between the
#' genetic variants and the outcome which bypass the risk factors. The regularized regression model is estimated by multivariable weighted linear regression where the intercept terms are subject
#' to lasso-type penalization. The lasso penalization will tend to shrink the intercept terms corresponding to the valid instruments to zero.
#'
#' The lasso penalty relies on a tuning parameter which controls the level of sparsity. The default is to use a heterogeneity stopping rule, but a fixed value may be specified.
#'
#' As part of the analysis, the genetic variants are orientated so that all of the associations with one of the risk factors are positive (the first risk factor is used by default). Re-orientation
#' of the genetic variants is performed automatically as part of the function.
#'
#' The MR-Lasso method is performed in two steps. First, a regularized regression model is fitted, and some genetic variants are identified as valid instruments. Second, causal effects are estimated using standard multivariable IVW with only the valid genetic variants.
#' The post-lasso method will be performed as long as the number of genetic variants identified as valid instruments is greater than the number of risk factors.
#' The default heterogeneity stopping rule will always return more genetic variants as valid instruments than risk factors for identification.
#' The main estimates given by the method are the post-lasso estimates. However, parameter estimates from the regularized regression model used to identify invalid variants are also provided for completeness.
#'
#' If a substantial proportion of genetic variants are removed from the analysis, the multivariable MR-Lasso method may give a false impression of confidence in the causal estimate due to homogeneity of the variant-specific causal estimates amongst the remaining variants. However, it is not reasonable to claim that there is strong evidence for a causal effect after a large number of variants with heterogeneous estimates have been removed from the analysis.
#'
#' @return The output from the function is an \code{MVLasso} object containing:
#'
#'  \item{Exposure}{A character vector with the names given to the exposure.}
#'  \item{Outcome}{A character string with the names given to the outcome.}
#'  \item{Estimate}{A vector of causal estimates from the multivariable MR-Lasso method. These are the post-lasso estimates.}
#'  \item{StdError}{A vector of standard errors of the causal estimates from the multivariable MR-Lasso method.}
#'  \item{CILower}{The lower bounds of the confidence intervals for the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{CIUpper}{The upper bounds of the confidence intervals for the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-values associated with the (post-lasso) causal estimates using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{RegEstimate}{The estimates from the regularized regression model used in the multivariable MR-Lasso method.}
#'  \item{RegIntercept}{The intercept estimates estimates from the regularized regression model used in the multivariable MR-Lasso method.}
#'  \item{Valid}{The number of genetic variants that have been identified as valid instruments.}
#'  \item{ValidSNPs}{The names of genetic variants that have been identified as valid instruments.}
#'  \item{Lambda}{The value of the tuning parameter used to compute \code{RegEstimate}.}
#'
#'
#' @examples mr_mvlasso(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse),
#'    by = chdlodds, byse = chdloddsse))
#'
#' @references Andrew J Grant, Stephen Burgess. Pleiotropy robust methods for multivariable Mendelian randomization. arXiv 2020; 2008.11997
#'
#' @export

setGeneric(name = "mr_mvlasso",
           def = function(object,
                          orientate = 1, distribution = "normal", alpha = 0.05, lambda = numeric(0))
           {standardGeneric("mr_mvlasso")})

#--------------------------------------------------------------------------------------------

#' MR-Lasso method
#'
#' @description The \code{mr_lasso} function performs the MR-Lasso method, which applies lasso-type penalization to the direct effects of genetic variants on the outcome.
#' The causal estimate is described as a post-lasso estimate, and is obtained by performing the IVW method using only those genetic variants that are identified as valid by the lasso procedure.
#'
#' @param object An \code{MRInput} object.
#' @param distribution The type of distribution used to calculate the confidence intervals. Options are \code{"normal"} (default) or \code{"t-dist"}.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#' @param lambda The value of the tuning parameter used by the lasso procedure which controls the level of sparsity. If not specified, the tuning parameter will be calculated by the heterogeneity stopping rule.
#'
#' @details MR-Lasso extends the IVW model to include an intercept term for each genetic variant. These intercept terms represent associations between the
#' genetic variants and the outcome which bypass the risk factor. The causal effect estimates are estimated by weighted linear regression where the intercept terms are subject
#' to lasso-type penalization. The lasso penalization will tend to shrink the intercept terms corresponding to the valid instruments to zero.
#'
#' The lasso penalty relies on a tuning parameter which controls the level of sparsity. The default is to use a heterogeneity stopping rule, but a fixed value may be specified.
#'
#' As part of the analysis, the genetic variants are orientated so that all of the associations with the risk factor are positive (and signs of associations with the outcome are
#' changed to keep the orientation consistent if required). Re-orientation of the genetic variants is performed automatically as part of the function.
#'
#' The MR-Lasso method is performed in two steps. First, a regularized regression model is fitted, and some genetic variants are identified as valid instruments. Second, the causal effect is estimated using standard IVW with only the valid genetic variants.
#' The post-lasso method will be performed as long as at least two genetic variants are identified as valid instruments. The default heterogeneity stopping rule will always return at least two
#' genetic variants as valid instruments.
#' The main estimate given by the method is the post-lasso estimate. However, parameter estimates from the regularized regression model used to identify invalid variants are also provided for completeness.
#'
#' If a substantial proportion of genetic variants are removed from the analysis, the MR-Lasso method may give a false impression of confidence in the causal estimate due to homogeneity of the variant-specific causal estimates amongst the remaining variants. However, it is not reasonable to claim that there is strong evidence for a causal effect after a large number of variants with heterogeneous estimates have been removed from the analysis.
#'
#' @return The output from the function is an \code{MRLasso} object containing:
#'
#'  \item{Exposure}{A character vector with the names given to the exposure.}
#'  \item{Outcome}{A character string with the names given to the outcome.}
#'  \item{Estimate}{The causal estimate from the MR-Lasso method. This is the post-lasso estimate.}
#'  \item{StdError}{The standard error of the causal estimate from the MR-Lasso method.}
#'  \item{CILower}{The lower bound of the confidence interval for the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{CIUpper}{The upper bound of the confidence interval for the causal estimate based on the estimated standard error and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-value associated with the causal estimate using a normal or t-distribution (as specified in \code{distribution}).}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'  \item{RegEstimate}{The estimate from the regularized regression model used in the MR-Lasso method.}
#'  \item{RegIntercept}{The intercept estimates from the regularized regression model used in the MR-Lasso method.}
#'  \item{Valid}{The number of genetic variants that have been identified as valid instruments.}
#'  \item{ValidSNPs}{The names of genetic variants that have been identified as valid instruments.}
#'  \item{Lambda}{The value of the tuning parameter used to compute \code{RegEstimate}}
#'
#'
#' @examples mr_lasso(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))
#'
#' @references Jessica MB Rees, Angela M Wood, Frank Dudbridge, Stephen Burgess. Robust methods in Mendelian randomization via penalization of heterogeneous causal estimates. PLoS ONE 2019; 14(9):e0222362
#'
#' @export

setGeneric(name = "mr_lasso",
           def = function(object, distribution = "normal", alpha = 0.05, lambda = numeric(0))
           {standardGeneric("mr_lasso")})

#--------------------------------------------------------------------------------------------

#' Constrained maximum likelihood (cML) method
#'
#' @description Constrained maximum likelihood (cML) based Mendelian
#' Randomization method robust to both 
#' correlated and uncorrelated pleiotropy.
#'
#' @param object An \code{MRInput} object.
#' @param MA Whether model average is applied or not. Default is TRUE.
#' @param DP Whether data perturbation is applied or not. Default is TRUE. 
#' @param K_vec Set of candidate K's, the constraint parameter 
#' representing number of invalid IVs. Default is from 0 to (#IV - 2).
#' @param random_start Number of random starting points for cML, default is 0.
#' @param num_pert Number of perturbation when DP is TRUE, default is 200. 
#' @param random_start_pert Number of random start points for 
#' cML with data perturbation, default is 0.
#' @param maxit Maximum number of iterations for each optimization. 
#' Default is 100.
#' @param random_seed Random seed, default is 314. When \code{random_seed=NULL}, no random seed 
#' will be used and the results may not be reproducible.
#' @param n Sample size. When sample sizes of GWAS for exposure and outcome are different,
#' and/or when sample sizes of different SNPs are different, 
#' the smallest sample size is recommended to get conservative result and avoid type-I error.
#' See reference for more discussions.
#' @param Alpha Significance level for the confidence interval for estimate, default is 0.05.
#'
#' @details The MRcML method selects invalid IVs with correlated
#' and/or uncorrelated peliotropic effects using constrained maximum
#' likelihood. \code{cML-BIC} gives results of the selected model with 
#' original data, while \code{cML-MA-BIC} averages over all candidate models.
#' \code{cML-BIC-DP} and \code{cML-MA-BIC-DP} are the versions with 
#' data-perturbation to account for selection uncertainty when 
#' many invalid IVs have weak pleiotropic effects. 
#' 
#' When DP is performed,
#' two goodness-of-fit (GOF) tests are developed to check whether the 
#' model-based and DP- based variance estimates converge to the same estimate.
#' Small p-values of GOF tests indicate selection uncertainty is not ignorable, and 
#' results from DP is more reliable. See reference for more details.
#'
#' As the constrained maximum likelihood function is non-convex, multiple
#' random starting points could be used to find a global minimum. For some 
#' starting points the algorithm may not converge and a warning message will 
#' be prompted, typically this will not affect the results.
#' 
#'
#' @return The output from the function is an \code{MRcML} object containing:
#'
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Estimate}{Estimate of theta.}
#'  \item{StdError}{Standard error of estimate.}
#'  \item{Pvalue}{p-value of estimate.}
#'  \item{BIC_invalid}{Set of selected invalid IVs if cML-BIC is performed, i.e. without MA or DP.}
#'  \item{GOF1_p}{p-value of the first goodness-of-fit test.}
#'  \item{GOF2_p}{p-value of the second goodness-of-fit test.}
#'  \item{SNPs}{The number of SNPs that were used in the calculation.}
#'  \item{Alpha}{Significance level for the confidence interval for estimate, default is 0.05.}
#'  \item{CILower}{Lower bound of the confidence interval for estimate.}
#'  \item{CIUpper}{Upper bound of the confidence interval for estimate.}
#'  \item{MA}{Indicator of whether model average is applied.}
#'  \item{DP}{Indicator of whether data perturbation is applied.}
#'  
#'  
#' @examples # Perform cML-MA-BIC-DP:
#' mr_cML(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds,
#' byse = chdloddsse), num_pert=5, MA = TRUE, DP = TRUE, n = 17723)
#' # num_pert is set to 5 to reduce computational time
#' # the default value of 200 is recommended in practice
#' 
#' # Perform cML-BIC-DP:
#' mr_cML(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds,
#' byse = chdloddsse), MA = TRUE, DP = FALSE,, n = 17723)
#'    
#'
#' @references Xue, H., Shen, X., & Pan, W. (2021). 
#' Constrained maximum likelihood-based Mendelian randomization 
#' robust to both correlated and uncorrelated pleiotropic effects. 
#' The American Journal of Human Genetics, 108(7), 1251-1269.
#'
#' @export

setGeneric(name = "mr_cML",
           def = function(object, 
                          MA=TRUE, DP=TRUE,
                          K_vec = 0:(length(object@betaX) - 2),
                          random_start = 0,
                          num_pert = 200,
                          random_start_pert = 0,
                          maxit = 100,
                          random_seed = 314,
                          n,
                          Alpha = 0.05)
           {standardGeneric("mr_cML")})


#--------------------------------------------------------------------------------------------

#' Penalized inverse-variance weighted method
#'
#' @description The \code{mr_pivw} function implements the penalized inverse-variance weighted (pIVW) method.
#'
#' @param object An \code{MRInput} object.
#' @param lambda The penalty parameter in the pIVW estimator. It plays a role in the bias-variance trade-off of the estimator. It is recommended to choose \code{lambda=1} to achieve the smallest bias and valid inference. By default, \code{lambda=1}.
#' @param over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @param delta The z-score threshold for IV selection. \code{delta} should be greater than or equal to zero. By default, \code{delta=0} (i.e., no IV selection will be conducted).  See 'Details'.
#' @param sel.pval A numeric vector containing the P-values of the SNP effects on the exposure, which will be used for the IV selection. \code{sel.pval} should be provided when \code{delta} is not zero. See 'Details'.
#' @param Boot.Fieller If \code{Boot.Fieller=TRUE}, then the P-value and the confidence interval of the causal effect will be calculated based on the bootstrapping Fieller method. Otherwise, the P-value and the confidence interval of the causal effect will be calculated from the normal distribution. It is recommended to use the bootstrapping Fieller method when \code{Condition} is smaller than 10 (see 'Details'). By default, \code{Boot.Fieller=TRUE}.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#'
#' @details The penalized inverse-variance weighted (pIVW) estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously
#' in two-sample MR with summary statistics, i.e., an exposure sample (with IV-exposure effect \code{Bx} and standard error \code{Bxse}) and
#' an outcome sample (with IV-outcome effect \code{By} and standard error \code{Byse}).
#'
#' The pIVW estimator also allows for IV selection in three-sample MR, where weak IVs are screened out using
#' an extra sample (with IV-exposure effect \code{Bx*} and standard error \code{Bxse*}) independent of the exposure sample and outcome sample.
#' Generally, the P-value for \code{Bx*} can be computed by \code{sel.pval=2*pnorm(abs(Bx*/Bxse*), lower.tail = FALSE)},
#' Given \code{sel.pval} and a z-score threshold \code{delta}, the variants kept in the analysis will be those
#' with \code{sel.pval<2*pnorm(delta,lower.tail = FALSE)}.
#'
#' The \code{mr_pivw} function outputs a measure \code{Condition} that needs to be large for reliable asymptotic properties of the pIVW estimator.
#' We also refer to \code{Condition} as effective sample size, which is a function of a measure of IV strength and the number of IVs.
#' When \code{delta} is zero (i.e., no IV selection), \code{Condition = (average F-statistic -1)*sqrt(# snps)}. When \code{delta} is not zero
#' (i.e., IV selection is conducted), \code{Condition = [(average F-statistic -1)*sqrt(# snps)]/c},
#' where the numerator is computed using the selected variants, and the denominator \code{c} involves the selection probabilities
#' of all variants (see more details in the paper \url{https://onlinelibrary.wiley.com/doi/10.1111/biom.13732}). We suggest that \code{Condition} should be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.
#'
#'
#' @return The output from the function is a \code{PIVW} object containing:
#'
#'  \item{Over.dispersion}{\code{TRUE} if the method has considered balanced horizontal pleiotropy, \code{FALSE} otherwise.}
#'  \item{Boot.Fieller}{\code{TRUE} if the bootstrapping Fieller method is used to calculate the P-value and the confidence interval of the causal effect, \code{FALSE} otherwise.}
#'  \item{Lambda}{The penalty parameter in the pIVW estimator.}
#'  \item{Delta}{The z-score threshold for IV selection.}
#'  \item{Exposure}{A character string giving the name given to the exposure.}
#'  \item{Outcome}{A character string giving the name given to the outcome.}
#'  \item{Estimate}{The causal point estimate from the pIVW estimator.}
#'  \item{StdError}{The standard error associated with \code{Estimate}.}
#'  \item{CILower}{The lower bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then lower limits of all ranges will be reported.}
#'  \item{CIUpper}{The upper bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then upper limits of all ranges will be reported.}
#'  \item{Alpha}{The significance level used in constructing the confidence interval.}
#'  \item{Pvalue}{P-value associated with the causal estimate from the pIVW estimator, which is derived from the bootstrapping Fieller method or normal distribution.}
#'  \item{Tau2}{The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated by using all IVs in the data before conducting the IV selection.}
#'  \item{SNPs}{The number of SNPs after IV selection.}
#'  \item{Condition}{The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties. See 'Details'.}
#'
#' @examples mr_pivw(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))
#'
#' @references Xu S., Wang P., Fung W.K. and Liu Z. (2022). A Novel Penalized Inverse-Variance Weighted Estimator for Mendelian Randomization with Applications to COVID-19 Outcomes. Biometrics. doi: 10.1111/biom.13732.
#'
#' @export

setGeneric(name = "mr_pivw",
           def = function(object, lambda=1, over.dispersion=TRUE, delta=0, sel.pval=NULL, Boot.Fieller=TRUE, alpha=0.05)
           {standardGeneric("mr_pivw")})


#--------------------------------------------------------------------------------------------

#' Multivariable constrained maximum likelihood method
#'
#' @description The \code{mr_mvcML} function performs multivariable Mendelian randomization via the constrained maximum likelihood method, which is robust to both correlated and uncorrelated pleiotropy.
#'
#' @param object An \code{MRMVInput} object.
#' @param n Sample size. The smallest sample size among all (both exposures and outcome) GWAS used in the analysis is recommended.
#' @param DP Whether data perturbation is applied or not. Default is TRUE. 
#' @param rho_mat  The correlation matrix among the exposures and outcome GWAS estimates, which can be estimated by the intercept term from bivariate LDSC. See reference for more discussions. Default is the identify matrix, for example, in the absence of overlapping samples among GWAS datasets.
#' @param K_vec Set of candidate K's, the constraint parameter 
#' representing number of invalid IVs. It can range from 0 up to #IV - (#exposure + 1). Default is from 0 to (#IV/2).
#' @param random_start Number of random starting points for MVMRcML, default is 0.
#' @param num_pert Number of perturbation when DP is TRUE, default is 200. 
#' @param min_theta_range The lower bound of the uniform distribution for each initial value for theta generated from, default is -0.5. 
#' @param max_theta_range The uppder bound of the uniform distribution for each initial value for theta generated from, default is 0.5. 
#' @param maxit Maximum number of iterations for each optimization. Default is 100.
#' @param alpha Significance level for the confidence interval for estimate, default is 0.05.
#' @param seed The random seed to use when generating the perturbed samples (for reproducibility). The default value is 314159265. If set to \code{NA}, the random seed will not be set (for example, if the function is used as part of a larger simulation).
#'
#' @details Multivariable MRcML (MVMRcML) is an extension of MRcML to deal with multiple exposures of interest. It is robust to both correlated and uncorrelated pleiotropy as its univariable version.
#'
#' In practice, the data perturbation (DP) version is preferred in practice for a more robust inference as it can account for the uncertainty in model selection.
#' However, it may take a longer time especially when the number of IVs is large (so the range of \code{K_vec} can be large too).
#' One strategy is to try a small range of K (the number of invalid IVs) first (with a small \code{num_pert}), 
#' then adjust it if the number of selected invalid IVs fall close to the boundary.
#' You can also use other methods, e.g. \code{mr_mvlasso}, to get a rough sense of the number of invalid IVs.
#'
#' Similar to \code{mr_cML}, multiple random starting points could be used to find a global minimum.
#'
#' @return The output from the function is an \code{MVMRcML} object containing:
#'
#'  \item{Exposure}{A character vector with the names given to the exposure.}
#'  \item{Outcome}{A character string with the names given to the outcome.}
#'  \item{Estimate}{A vector of causal estimates.}
#'  \item{StdError}{A vector of standard errors of the causal estimates.}
#'  \item{CILower}{The lower bounds of the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{CIUpper}{The upper bounds of the causal estimates based on the estimated standard errors and the significance level provided.}
#'  \item{Alpha}{The significance level used when calculating the confidence intervals.}
#'  \item{Pvalue}{The p-values associated with the estimates (calculated as Estimate/StdError as per Wald test) using a normal distribution.}
#'  \item{BIC_invalid}{Set of selected invalid IVs by MVMRcML-BIC.}
#'  \item{K_hat}{The number of selected invalid IVs by MVMRcML-BIC, or a vector for each data perturbation in MVMRcML-DP.}
#'  \item{eff_DP_B}{The number of data perturbations with successful convergence in MVMRcML-DP.}
#'  \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
#'
#'
#' @examples 
#' # Perform MVMRcML-DP:
#' mr_mvcML(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse),
#'    by = chdlodds, byse = chdloddsse), n = 17723, num_pert = 5, random_start = 5)
#'   # num_pert is set to 5 to reduce runtime for the mr_mvcML method,
#'   # At least 100 perturbations should be used and more is preferred for a stable result.
#' 
#' rho_mat = matrix(c(1,-0.1,0.2,0,-0.1,1,-0.3,0,
#'  				  0.2,-0.3,1,0,0,0,0,1),ncol=4) ## Toy example of rho_mat
#' mr_mvcML(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse),
#'    by = chdlodds, byse = chdloddsse), n = 17723, num_pert = 100, rho_mat = rho_mat)
#'
#' # Perform MVMRcML-BIC:
#' mr_mvcML(mr_mvinput(bx = cbind(ldlc, hdlc, trig), bxse = cbind(ldlcse, hdlcse, trigse),
#'    by = chdlodds, byse = chdloddsse), n = 17723, DP = FALSE)
#'
#' @references Lin, Z., Xue, H., & Pan, W. (2023). Robust multivariable Mendelian randomization based on constrained maximum likelihood. The American Journal of Human Genetics, 110(4), 592-605.
#'
#' @export

setGeneric(name = "mr_mvcML",
           def = function(object,
                   n,
                   DP = TRUE,
                   rho_mat = diag(ncol(object@betaX)+1),
                   K_vec = 0:(ceiling(nrow(object@betaX)/2)),
                   random_start = 0,
                   num_pert = 100,
                   min_theta_range = -0.5,
                   max_theta_range = 0.5,
                   maxit = 100,
                   alpha = 0.05, 
                   seed = 314159265)
           {standardGeneric("mr_mvcML")})


