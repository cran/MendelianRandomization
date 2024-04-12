# MRInput class

#' MRInput Class
#'
#' @description An object containing the four vectors of summary statistics required to calculate Mendelian randomization estimates.
#'
#' @slot betaX A numeric vector of beta-coefficient values for genetic associations with the first variable (often referred to as the exposure, risk factor, or modifiable phenotype).
#' @slot betaY A numeric vector of beta-coefficient values for genetic associations with the second variable (often referred to as the outcome). For a disease outcome, the beta coefficients are log odds estimates from logistic regression analyses.
#' @slot betaXse The standard errors associated with the beta-coefficients in \code{betaX}.
#' @slot betaYse The standard errors associated with the beta-coefficients in \code{betaY}.
#' @slot correlation The matrix of correlations between genetic variants. If this variable is not provided, then we assume that genetic variants are uncorrelated.
#' @slot exposure The name of the exposure variable.
#' @slot outcome The name of the outcome variable.
#' @slot snps The names of the genetic variants (SNPs) included in the analysis. The slots \code{exposure}, \code{outcome}, and \code{snps} are not required, but may be useful for keeping track of various \code{MRInput} objects. They are also used by the \code{mr_plot} function.
#' @slot effect_allele The name of the effect allele for each SNP. The beta-coefficients are the associations with the exposure and outcome per additional copy of the effect allele.
#' @slot other_allele The name of the non-effect allele.
#' @slot eaf The expected allele frequencies (numeric). The slots \code{effect_allele}, \code{other_allele}, and \code{eaf} are neither required, nor currently used in the MendelianRandomization package. They are included for future compatibility with the MR-Base suite of functions.
#'
#' @details The beta-coefficients are assumed to be estimated for uncorrelated (independent) genetic variants, although a correlation matrix can be specified if the variants are correlated in their distributions. We also assume that the beta-coefficients for associations with the exposure and with the outcome are uncorrelated (corresponding to a two-sample Mendelian randomization analysis), although correlation between associations with the exposure and with the outcome generally have little impact on causal estimates or standard errors.
#' Estimates can either be specified by the user, or extracted from the PhenoScanner tool.
#'

setClass("MRInput",
         representation(betaX = "numeric",
                        betaY = "numeric",
                        betaXse = "numeric",
                        betaYse = "numeric",
                        exposure = "character",
                        outcome = "character",
                        snps = "character",
                        effect_allele = "character",
                        other_allele  = "character",
                        eaf           = "numeric",
correlation = "matrix"),

         prototype = prototype(betaX = ldlc,
                               betaY = chdlodds,
                               betaXse = ldlcse,
                               betaYse = chdloddsse,
                               exposure = "LDL-c",
                               outcome = "CHD",
                               snps = "snp",
                               effect_allele = lipid_effect,
                               other_allele  = lipid_other,
                               eaf = lipid_eaf,
                               correlation = calc.rho)
)

# Ensure the vectors are all of the same length
setValidity("MRInput",
            function(object)
            {if(!all(length(object@betaX) == length(object@betaY),
                     length(object@betaXse) == length(object@betaYse),
                     length(object@betaXse) == length(object@betaY)))
            { cat("Vectors do not all have the same length.")
            } else {}
            }
)

#--------------------------------------------------------------------------------------------

#' WeightedMedian Class
#'
#' @description An object containing the estimate produced using the median-based method as well as various statistics.
#'
#' @slot Type The type of median that has been calculated, \code{"simple"}, \code{"weighted"}, or \code{"penalized"}.
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate The causal point estimate from the median-based method.
#' @slot StdError The standard error associated with \code{Estimate} (obtained from bootstrapping).
#' @slot CILower The lower bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate from the Wald method.
#' @slot SNPs The number of SNPs that used in the calculation.

setClass("WeightedMedian",
         representation(Type = "character",
                        Exposure = "character",
                        Outcome = "character",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        SNPs = "numeric")
)


#--------------------------------------------------------------------------------------------

#' DIVW Class
#'
#' @description An object containing the estimate produced using the debiased inverse-variance weighted (dIVW) method as well as various statistics.
#'
#' @slot Over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate The causal point estimate from the median-based method.
#' @slot StdError The standard error associated with \code{Estimate} (obtained from bootstrapping).
#' @slot CILower The lower bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate from the Wald method.
#' @slot SNPs The number of SNPs that used in the calculation.
#' @slot Condition A measure defined as (average F-statistic -1)*sqrt(# snps) that needs to be large for reliable asymptotic approximation based on the dIVW estimator. It is recommended to be greater than 20.

setClass("DIVW",
         representation(Over.dispersion = "logical",
                        Exposure = "character",
                        Outcome = "character",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        SNPs = "numeric",
                        Condition = "numeric")
)


#--------------------------------------------------------------------------------------------

#' IVW Class
#'
#' @description An object containing the estimate produced using the inverse-variance weighted (IVW) method as well as various statistics.
#'
#' @slot Model The model used for estimation: random-effects (\code{"random"}) or fixed-effect (\code{"fixed"}). The default option (\code{"default"}) is to use a fixed-effect model when there are three or fewer genetic variants, and a random-effects model when there are four or more. The (multiplicative) random-effects model allows for heterogeneity between the causal estimates targeted by the genetic variants by allowing over-dispersion in the regression model. Under-dispersion is not permitted (in case of under-dispersion, the residual standard error is set to 1, as in a fixed-effect analysis).
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot Robust Whether robust regression was used in the regression model relating the genetic associations with the outcome and those with the exposure.
#' @slot Penalized Whether weights in the regression model were penalized for variants with heterogeneous causal estimates.
#' @slot Estimate The causal point estimate from the inverse-variance weighted method.
#' @slot StdError The standard error associated with \code{Estimate}.
#' @slot CILower The lower bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate.
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot RSE The estimated residual standard error from the regression model.
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.
#' @slot Fstat An approximation of the first-stage F statistic for all variants based on the summarized data.

setClass("IVW",
         representation(Model = "character",
                        Exposure = "character",
                        Outcome = "character",

                        Robust = "logical",
                        Penalized = "logical",
                        Correlation = "matrix",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        SNPs = "numeric",

                        RSE = "numeric",
                        Heter.Stat = "numeric",
                        Fstat = "numeric")
)

#--------------------------------------------------------------------------------------------

#' Egger Class
#'
#' @description An object containing the estimate produced using the MR-Egger method as well as various statistics.
#'
#' The MR-Egger model uses a random-effects model; a fixed-effect model does not make sense as pleiotropy leads to heterogeneity between the causal estimates targeted by the genetic variants. The (multiplicative) random-effects model allows over-dispersion in the regression model. Under-dispersion is not permitted (in case of under-dispersion, the residual standard error is set to 1).
#'
#' @slot Model Model always takes the value \code{random}, as only random-effects analyses are permitted.
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot Robust Whether robust regression was used in the regression model relating the genetic associations with the outcome and those with the exposure.
#' @slot Penalized Whether weights in the regression model were penalized for variants with heterogeneous causal estimates.
#' @slot Estimate The causal point estimate from the MR-Egger method.
#' @slot StdError.Est The standard error associated with \code{Estimate}.
#' @slot Pvalue.Est P-value associated with the causal estimate from the Wald method.
#' @slot CILower.Est The lower bound of the confidence interval for \code{Estimate} based on \code{StdError.Est}.
#' @slot CIUpper.Est The upper bound of the confidence interval for \code{Estimate} based on \code{StdError.Est}.
#' @slot Intercept The intercept estimate from the MR-Egger method. Under the InSIDE assumption, the intercept represents the average pleiotropic effect (average direct effect on the outcome) of a genetic variant. If the intercept differs from zero, this is evidence that the genetic variants are not all valid instruments; specifically, there is directional pleiotropy.
#' @slot StdError.Int The standard error associated with \code{Intercept}.
#' @slot Pvalue.Int P-value associated with the intercept from the Wald method.
#' @slot CILower.Int The lower bound of the confidence interval for \code{Intercept} based on \code{StdError.Int}.
#' @slot CIUpper.Int The upper bound of the confidence interval for \code{Estimate} based on \code{StdError.Int}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot Causal.pval P-value associated with the causal estimate.
#' @slot Pleio.pval P-value associated with the intercept (p-value for the MR-Egger intercept test of directional pleiotropy).
#' @slot RSE The estimated residual standard error from the regression model.
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that the MR-Egger regression model describes the associations with the outcome with no excess heterogeneity.
#' @slot I.sq A measure of heterogeneity between the genetic associations with the exposure (see Bowden IJE 2016: "Assessing the suitability of summary data for Mendelian randomization analyses using MR-Egger regression: The role of the I2 statistic."). Low values of \code{I.sq} relate both to large differences in precision between MR-Egger and IVW estimates, and to more weak instrument bias (in a two-sample setting, this is attenuation of MR-Egger estimate towards the null).

setClass("Egger",
         representation(Model = "character",
                        Exposure = "character",
                        Outcome = "character",

                        Robust = "logical",
                        Penalized = "logical",
                        Correlation = "matrix",

                        Estimate = "numeric",
                        StdError.Est = "numeric",
                        CILower.Est = "numeric",
                        CIUpper.Est = "numeric",
                        Pvalue.Est = "numeric",

                        Intercept = "numeric",
                        StdError.Int = "numeric",
                        CILower.Int = "numeric",
                        CIUpper.Int = "numeric",
                        Pvalue.Int = "numeric",

                        Pleio.pval = "numeric",
                        Causal.pval = "numeric",

                        Alpha = "numeric",

                        SNPs = "numeric",
                        RSE = "numeric",
                        Heter.Stat = "numeric",
                        I.sq = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MRAll Class
#'
#' @description An object containing the estimates produced using the \code{mr_allmethods} function.
#'
#' @slot Data The \code{mr_input} object that was used as an input to the \code{mr_allmethods} function. This includes the original data, so that a call to \code{mr_plot} can plot the original data and the various causal estimates.
#' @slot Values A data.frame object comprising estimates from the various methods called by the \code{mr_allmethods} function. The first column gives the names of the methods, then the causal estimates, standard errors, 95\% confidence intervals, and p-values.
#' @slot Method A string indicating whether all methods are implemented (\code{"all"}, the default option), or just main methods (\code{"main"}), or only a subset of methods (\code{"ivw"}, \code{"egger"}, or \code{"median"}).

setClass("MRAll",
         representation(Data = "MRInput",
                        Values = "data.frame", Method = "character")
)

#--------------------------------------------------------------------------------------------

#' MaxLik Class
#'
#' @description An object containing the estimate produced using the maximum-likelihood method as well as various statistics.
#'
#' @slot Model The model used for estimation: fixed-effect (\code{"fixed"}) or random-effects (\code{"random"}).
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot Psi The correlations between genetic associations with the exposure and with the outcome.
#' @slot Estimate The causal point estimate from the inverse-variance weighted method.
#' @slot StdError The standard error associated with \code{Estimate}.
#' @slot CILower The lower bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate.
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot RSE The estimated residual standard error from the regression model.
#' @slot Heter.Stat Heterogeneity statistic (likelihood ratio statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.

setClass("MaxLik",
         representation(Model = "character",
                        Exposure = "character",
                        Outcome = "character",

                        Correlation = "matrix",
                        Psi = "numeric",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        SNPs = "numeric",

                        RSE = "numeric",
                        Heter.Stat = "numeric")
)

#--------------------------------------------------------------------------------------------

# MRMVInput class

#' MRMVInput Class
#'
#' @description An object containing the summary statistics required to calculate multivariable Mendelian randomization estimates.
#'
#' @slot betaX A matrix of beta-coefficient values for genetic associations with the risk factor variables. These should be arranged so that column 1 are the beta-coefficients for risk factor 1, and row 1 are the beta-coefficients for genetic variant 1.
#' @slot betaY A numeric vector of beta-coefficient values for genetic associations with the second variable (often referred to as the outcome). For a disease outcome, the beta coefficients are log odds estimates from logistic regression analyses.
#' @slot betaXse The matrix of standard errors associated with the beta-coefficients in \code{betaX}.
#' @slot betaYse The vector of standard errors associated with the beta-coefficients in \code{betaY}.
#' @slot correlation The matrix of correlations between genetic variants. If this variable is not provided, then we assume that genetic variants are uncorrelated.
#' @slot exposure The names of the exposure variables.
#' @slot outcome The name of the outcome variable.
#' @slot snps The names of the genetic variants (SNPs) included in the analysis. The slots \code{exposure}, \code{outcome}, and \code{snps} are not required, but may be useful for keeping track of various \code{MRInput} objects. They are also used by the \code{mr_plot} function.
#' @slot effect_allele The name of the effect allele for each SNP. The beta-coefficients are the associations with the exposure and outcome per additional copy of the effect allele.
#' @slot other_allele The name of the non-effect allele.
#' @slot eaf The expected allele frequencies (numeric). The slots \code{effect_allele}, \code{other_allele}, and \code{eaf} are neither required, nor currently used in the MendelianRandomization package. They are included for future compatibility with the MR-Base suite of functions.
#'
#' @details The beta-coefficients are assumed to be estimated for uncorrelated (independent) genetic variants, although a correlation matrix can be specified if the variants are correlated in their distributions. We also assume that the beta-coefficients for associations with the exposure and with the outcome are uncorrelated (corresponding to a two-sample Mendelian randomization analysis), although correlation between associations with the exposure and with the outcome generally have little impact on causal estimates or standard errors.

setClass("MRMVInput",
         representation(betaX = "matrix",
                        betaY = "numeric",
                        betaXse = "matrix",
                        betaYse = "numeric",
                        exposure = "character",
                        outcome = "character",
                        snps = "character",
                        effect_allele = "character",
                        other_allele  = "character",
                        eaf           = "numeric",
correlation = "matrix"),

         prototype = prototype(betaX = cbind(ldlc, hdlc, trig),
                               betaY = chdlodds,
                               betaXse = cbind(ldlcse, hdlcse, trigse),
                               betaYse = chdloddsse,
                               exposure = c("LDL-c", "HDL-c", "triglycerides"),
                               outcome = "CHD",
                               snps = "snp",
                               effect_allele = lipid_effect,
                               other_allele  = lipid_other,
                               eaf = lipid_eaf,
                               correlation = calc.rho)
)

# Ensure the vectors are all of the same length
setValidity("MRMVInput",
            function(object)
            {if(!all(dim(object@betaX)[1] == length(object@betaY),
                     dim(object@betaXse)[1] == length(object@betaYse),
                     length(object@betaY) == length(object@betaYse),
                     dim(object@betaXse)[2] == dim(object@betaX)[2]))
            { cat("Inputs do not all have the same length.")
            } else {}
             if(dim(object@betaX)[1] < dim(object@betaX)[2])
            { cat("Multivariable MR requires more genetic variants than risk factors.")#
            } else {}

            }
)

#--------------------------------------------------------------------------------------------

#' MVIVW Class
#'
#' @description An object containing the estimates produced using the multivariable inverse-variance weighted (IVW) method as well as various statistics.
#'
#' @slot Model The model used for estimation: random-effects (\code{"random"}) or fixed-effect (\code{"fixed"}). The default option (\code{"default"}) is to use a fixed-effect model when there are three or fewer genetic variants, and a random-effects model when there are four or more. The (multiplicative) random-effects model allows for heterogeneity between the causal estimates targeted by the genetic variants by allowing over-dispersion in the regression model. Under-dispersion is not permitted (in case of under-dispersion, the residual standard error is set to 1, as in a fixed-effect analysis).
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Robust Whether robust regression was used in the regression model relating the genetic associations with the outcome and those with the exposure.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot Estimate The causal estimates from the inverse-variance weighted method.
#' @slot StdError The standard errors associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate.
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot RSE The estimated residual standard error from the regression model.
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.
#' @slot CondFstat Conditional F statistics: An approximation of the first-stage conditional F statistics for all variants based on the summarized data. This represents the instrument strength for each exposure conditional on other exposures in the model. This is only reported if the sample sizes for the genetic associations with the exposures are provided.

setClass("MVIVW",
         representation(Model = "character",
                        Exposure = "character",
                        Outcome = "character",

                        Robust = "logical",
                        Correlation = "matrix",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        SNPs = "numeric",

                        RSE = "numeric",
                        Heter.Stat = "numeric",
                        CondFstat = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MVIVWME Class
#'
#' @description An object containing the estimates produced using the multivariable inverse-variance weighted (IVW) method with measurement error, as well as various statistics.
#'
#' @slot Model The model used for estimation: random-effects (\code{"random"}) or fixed-effect (\code{"fixed"}). The default option (\code{"default"}) is to use a fixed-effect model when there are three or fewer genetic variants, and a random-effects model when there are four or more. The (multiplicative) random-effects model allows for heterogeneity between the causal estimates targeted by the genetic variants by allowing over-dispersion in the regression model. Under-dispersion is not permitted (in case of under-dispersion, the residual standard error is set to 1, as in a fixed-effect analysis).
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot Estimate The causal estimates from the inverse-variance weighted method.
#' @slot StdError The standard errors associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate.
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot RSE The estimated residual standard error from the regression model.
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.

setClass("MVIVWME",
         representation(Model = "character",
                        Exposure = "character",
                        Outcome = "character",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",
						
						Correlation = "matrix",

                        Pvalue = "numeric",
                        SNPs = "numeric",

                        RSE = "numeric",
                        Heter.Stat = "numeric"
						)
)


#--------------------------------------------------------------------------------------------

#' MVEgger Class
#'
#' @description An object containing the estimates produced using the multivariable MR-Egger method as well as various statistics.
#'
#' @slot Model Model always takes the value \code{random}, as only random-effects analyses are permitted.
#' @slot Orientate The number of the risk factor that genetic associations are orientated to. The default value is \code{1}, meaning that genetic associations with the first risk factor are set to be positive.
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot Estimate The causal estimates from the inverse-variance weighted method.
#' @slot StdError.Est The standard errors associated with \code{Estimate}.
#' @slot CILower.Est The lower bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper.Est The upper bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Pvalue.Est P-value associated with the causal estimate.
#' @slot Intercept The intercept estimate from the MR-Egger method. Under the InSIDE assumption, the intercept represents the average pleiotropic effect (average direct effect on the outcome) of a genetic variant. If the intercept differs from zero, this is evidence that the genetic variants are not all valid instruments; specifically, there is directional pleiotropy.
#' @slot StdError.Int The standard error associated with \code{Intercept}.
#' @slot CILower.Int The lower bound of the confidence interval for \code{Intercept} based on \code{StdError.Int}.
#' @slot CIUpper.Int The upper bound of the confidence interval for \code{Estimate} based on \code{StdError.Int}.
#' @slot Pvalue.Int P-value associated with the intercept.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot RSE The estimated residual standard error from the regression model.
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more variants may be pleiotropic.

setClass("MVEgger",
         representation(Model = "character",
                        Orientate = "numeric",
                        Exposure = "character",
                        Outcome = "character",

                        Correlation = "matrix",

                        Estimate = "numeric",
                        StdError.Est = "numeric",
                        CILower.Est = "numeric",
                        CIUpper.Est = "numeric",
                        Pvalue.Est = "numeric",

                        Intercept = "numeric",
                        StdError.Int = "numeric",
                        CILower.Int = "numeric",
                        CIUpper.Int = "numeric",
                        Pvalue.Int = "numeric",

                        Alpha = "numeric",
                        SNPs = "numeric",

                        RSE = "numeric",
                        Heter.Stat = "numeric")
)


#--------------------------------------------------------------------------------------------

#' MRMBE Class
#'
#' @description An object containing the estimate produced using the mode-based estimation method of Hartwig et al as well as various statistics.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Weighting Whether the analysis was \code{weighted} or \code{unweighted}.
#' @slot StdErr Whether the \code{simple} or \code{delta} version of the standard errors were used.
#' @slot Phi The value of the bandwidth factor.
#' @slot Estimate The causal estimate from the mode-based estimation method.
#' @slot StdError The standard errors associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate.
#' @slot SNPs The number of SNPs that were used in the calculation.

setClass("MRMBE",
         representation(Exposure = "character",
                        Outcome  = "character",
                        Weighting= "character",
                        StdErr   = "character",
                        Phi      = "numeric",
                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        SNPs = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MRHetPen Class
#'
#' @description An object containing the estimate produced using the heterogeneity-penalized model-averaging mode-based estimation method as well as various statistics.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Prior The value of the prior probability of a genetic variant being a valid instrument (default is 0.5).
#' @slot Estimate The causal estimate from the heterogeneity-penalized method.
#' @slot CIRange The confidence interval for \code{Estimate} based on a grid search.
#' @slot CILower The lower limit of the confidence interval. If the confidence interval contains multiple ranges, then lower limits of all ranges will be reported.
#' @slot CIUpper The upper limit of the confidence interval. If the confidence interval contains multiple ranges, then upper limits of all ranges will be reported.
#' @slot CIMin The smallest value used in the search to find the confidence interval.
#' @slot CIMax The largest value used in the search to find the confidence interval.
#' @slot CIStep The step size used in the search to find the confidence interval.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot SNPs The number of SNPs that were used in the calculation.

setClass("MRHetPen",
         representation(Exposure = "character",
                        Outcome  = "character",
                        Prior    = "numeric",
                        Estimate = "numeric",
                        CIRange  = "numeric",
                        CILower  = "numeric",
                        CIUpper  = "numeric",
                        CIMin    = "numeric",
                        CIMax    = "numeric",
                        CIStep   = "numeric",
                        Alpha    = "numeric",
                        SNPs = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MRConMix Class
#'
#' @description An object containing the estimate produced using the contamination mixture method as well as various statistics.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Psi The value of the standard deviation of the distribution of invalid estimands (default is 1.5 times the standard deviation of the ratio estimates).
#' @slot Estimate The causal estimate from the contamination mixture method.
#' @slot CIRange The confidence interval for \code{Estimate} based on a grid search.
#' @slot CILower The lower limit of the confidence interval. If the confidence interval contains multiple ranges, then lower limits of all ranges will be reported.
#' @slot CIUpper The upper limit of the confidence interval. If the confidence interval contains multiple ranges, then upper limits of all ranges will be reported.
#' @slot CIMin The smallest value used in the search to find the confidence interval.
#' @slot CIMax The largest value used in the search to find the confidence interval.
#' @slot CIStep The step size used in the search to find the confidence interval.
#' @slot Pvalue The p-value associated with the estimate calculated using the likelihood function and a chi-squared distribution.
#' @slot Valid The numbers of genetic variants that were considered valid instruments at the causal estimate.
#' @slot ValidSNPs The names of genetic variants that were considered valid instruments at the causal estimate.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot SNPs The number of SNPs that were used in the calculation.

setClass("MRConMix",
         representation(Exposure = "character",
                        Outcome  = "character",
                        Psi      = "numeric",
                        Estimate = "numeric",
                        CIRange  = "numeric",
                        CILower  = "numeric",
                        CIUpper  = "numeric",
                        CIMin    = "numeric",
                        CIMax    = "numeric",
                        CIStep   = "numeric",
                        Valid    = "numeric",
                        ValidSNPs= "character",
                        Pvalue   = "numeric",
                        Alpha    = "numeric",
                        SNPs = "numeric")
)

#' MRMVMedian class
#'
#' @description An object containing the estimates produced using the multivariable median method as well as various statistics.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate The causal estimates from the multivariable median method.
#' @slot StdError The standard errors associated with \code{Estimate} (obtained from bootstrapping).
#' @slot CILower The lower bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-values associated with the causal estimates from the Wald method.
#' @slot SNPs The number of SNPs that used in the calculation.

setClass("MVMedian",
         representation(Exposure = "character",
                        Outcome = "character",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        SNPs = "numeric")
)


#' MRMVLasso class
#'
#' @description An object containing the estimates produced using the multivariable MR-Lasso method as well as various statistics.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Orientate The number of the risk factor that genetic associations are orientated to. The default value is \code{1}, meaning that genetic associations with the first risk factor are set to be positive.
#' @slot Estimate The causal estimates from the multivariable MR-Lasso method.
#' @slot StdError The standard errors associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-values associated with the causal estimates from the multivariable MR-Lasso method.
#' @slot SNPs The number of SNPs used in the calculation.
#' @slot RegEstimate The estimates from the regularized regression model used in the multivariable MR-Lasso method.
#' @slot RegIntercept The intercept estimates from the regularized regression model used in the multivariable MR-Lasso method. An intercept estimate of zero identifies the corresponding genetic variant as a valid instrument. Genetic variants with non-zero intercept estimates will be excluded from the MR-Lasso method, which is obtained as a post-lasso estimator.
#' @slot Valid The number of genetic variants that have been identified as valid instruments.
#' @slot ValidSNPs The names of genetic variants that have been identified as valid instruments.
#' @slot Lambda The value of the tuning parameter used to compute \code{RegEstimate} (default is to calulate \code{Lambda} using the heterogeneity stopping rule).
#'

setClass("MVLasso",
         representation(Orientate = "numeric",
                        Exposure = "character",
                        Outcome = "character",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",

                        Alpha = "numeric",
                        Pvalue = "numeric",
                        SNPs = "numeric",

                        RegEstimate = "numeric",
                        RegIntercept = "numeric",
                        Valid = "numeric",
                        ValidSNPs = "character",
                        Lambda = "numeric")
)

#' MRLasso class
#'
#' @description An object containing the estimates produced using the MR-Lasso method as well as various statistics.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate The causal estimate from the MR-Lasso method.
#' @slot StdError The standard error associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate from the MR-Lasso method.
#' @slot SNPs The number of SNPs used in the calculation.
#' @slot RegEstimate The estimate from the regularized regression model used in the MR-Lasso method.
#' @slot RegIntercept The intercept estimates from the regularized regression model used in the MR-Lasso method. An intercept estimate of zero identifies the corresponding genetic variant as a valid instrument. Genetic variants with non-zero intercept estimates will be excluded from the post-lasso estimator.
#' @slot Valid The number of genetic variants that have been identified as valid instruments.
#' @slot ValidSNPs The names of genetic variants that have been identified as valid instruments.
#' @slot Lambda The value of the tuning parameter used to compute \code{RegEstimate} (default is to calulate \code{Lambda} using the heterogeneity stopping rule).
#'

setClass("MRLasso",
         representation(Exposure = "character",
                        Outcome = "character",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",

                        Alpha = "numeric",
                        Pvalue = "numeric",
                        SNPs = "numeric",

                        RegEstimate = "numeric",
                        RegIntercept = "numeric",
                        Valid = "numeric",
                        ValidSNPs = "character",
                        Lambda = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MRcML Class
#'
#' @description An object containing the results of MRcML.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate Estimate of theta.
#' @slot StdError Standard error of estimate.
#' @slot Pvalue p-value of estimate.
#' @slot BIC_invalid Set of selected invalid IVs if cML-BIC is performed, i.e. without MA or DP.
#' @slot GOF1_p p-value of the first goodness-of-fit test.
#' @slot GOF2_p p-value of the second goodness-of-fit test.
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot Alpha Significance level for the confidence interval for estimate, default is 0.05.
#' @slot CILower Lower bound of the confidence interval for estimate.
#' @slot CIUpper Upper bound of the confidence interval for estimate.
#' @slot MA Indicator of whether model average is applied.
#' @slot DP Indicator of whether data perturbation is applied.


setClass("MRcML",
         representation(Exposure = "character",
                        Outcome  = "character",
                        Estimate = "numeric",
                        StdError = "numeric",
                        Pvalue = "numeric",
                        BIC_invalid = "numeric",
                        GOF1_p = "numeric",
                        GOF2_p = "numeric",
                        SNPs = "numeric",
                        Alpha = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        MA = "logical",
                        DP = "logical")
)

#--------------------------------------------------------------------------------------------

#' PIVW Class
#'
#' @description An object containing the estimate produced using the penalized inverse-variance weighted (pIVW) method as well as various statistics.
#'
#' @slot Over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @slot Boot.Fieller If \code{Boot.Fieller=TRUE}, then the P-value and the confidence interval of the causal effect will be calculated based on the bootstrapping Fieller method. Otherwise, the P-value and the confidence interval of the causal effect will be calculated from the normal distribution. It is recommended to use the bootstrapping Fieller method when \code{Condition} (the estimated effective sample size) is smaller than 10. By default, \code{Boot.Fieller=TRUE}.
#' @slot Lambda The penalty parameter in the pIVW estimator. The penalty parameter plays a role in the bias-variance trade-off of the estimator. It is recommended to choose \code{lambda=1} to achieve the smallest bias and valid inference. By default, \code{lambda=1}.
#' @slot Delta The z-score threshold for IV selection. By default, \code{delta=0} (i.e., no IV selection will be conducted).
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate The causal point estimate from the pIVW estimator.
#' @slot StdError The standard error associated with \code{Estimate}.
#' @slot CILower The lower bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then lower limits of all ranges will be reported.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then upper limits of all ranges will be reported.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate from the pIVW estimator.
#' @slot Tau2 The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated by using all IVs in the data before conducting the IV selection.
#' @slot SNPs The number of SNPs after IV selection.
#' @slot Condition The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.

setClass("PIVW",
         representation(Over.dispersion = "logical",
                        Boot.Fieller = "logical",
                        Lambda = "numeric",
                        Delta = "numeric",
                        Exposure = "character",
                        Outcome = "character",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        Tau2 = "numeric",
                        SNPs = "numeric",
                        Condition = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MVMRcML Class
#'
#' @description An object containing the results of MVMRcML.
#'
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Estimate The causal estimates from the multivariable MRcML method.
#' @slot StdError The standard errors associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence intervals for \code{Estimate} based on \code{StdError}.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-values associated with the causal estimates from the multivariable MRcML method.
#' @slot BIC_invalid The index set of selected invalid IVs by MVMRcML-BIC.
#' @slot K_hat The number of selected invalid IVs by MVMRcML-BIC, or a vector for each data perturbation in MVMRcML-DP.
#' @slot eff_DP_B The number of data perturbations with successful convergence in MVMRcML-DP.
#' @slot SNPs The number of SNPs that were used in the calculation.
#' @slot DP Indicator of whether data perturbation is applied.


setClass("MVMRcML",
         representation(Exposure = "character",
                        Outcome  = "character",
                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",
                        Pvalue = "numeric",
                        BIC_invalid = "numeric",
                        K_hat = "numeric",
                        eff_DP_B = "numeric",
                        SNPs = "numeric",
                        DP = "logical")
)

#--------------------------------------------------------------------------------------------

#' PCGMM Class
#'
#' @description An object containing the estimates produced using the univariable principal components generalized method of methods (PC-GMM) method as well as various statistics.
#'
#' @slot robust Whether the robust model with overdispersion heterogeneity is estimated.
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot Estimate The causal estimate from the PC-GMM method.
#' @slot StdError The standard error associated with \code{Estimate}.
#' @slot CILower The lower bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Fstat The first-stage F statistic for all genetic principal components used as instruments.
#' @slot Overdispersion The estimate of the overdispersion parameter for the robust model.
#' @slot PCs The number of genetic principal components used to instrument the exposure.
#' @slot Pvalue P-value associated with the causal estimate.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value (for non-robust model): the null hypothesis is that all principal components estimate the same causal parameter; rejection of the null is an indication that one or more principal components may be pleiotropic.

setClass("PCGMM",
         representation(robust = "logical",
                        Exposure = "character",
                        Outcome = "character",
                        Correlation = "matrix",
                        
                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Fstat = "numeric",
                        Overdispersion = "numeric",
                        
                        PCs = "numeric",
                        Pvalue = "numeric",
                        
                        Alpha = "numeric",
                        Heter.Stat = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MVPCGMM Class
#'
#' @description An object containing the estimates produced using the multivariable principal components generalized method of methods (PC-GMM) method as well as various statistics.
#'
#' @slot robust Whether the robust model with overdispersion heterogeneity is estimated.
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot ExpCorrelation Whether an exposure correlation matrix was specified.
#' @slot CondFstat The conditional F-statistic for each exposure. 
#' @slot Estimate The causal estimates from the PC-GMM method.
#' @slot StdError The standard errors associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Overdispersion The estimate of the overdispersion parameter for the robust model.
#' @slot PCs The number of genetic principal components used to instrument the exposures.
#' @slot Pvalue P-value associated with the causal estimates.
#' @slot Alpha The significance level used in constructing confidence intervals (default is 0.05).
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value (for non-robust model): the null hypothesis is that all principal components estimate the same causal parameter; rejection of the null is an indication that one or more principal components may be pleiotropic.

setClass("MVPCGMM",
         representation(robust = "logical",
                        Exposure = "character",
                        Outcome = "character",
                        Correlation = "matrix",
                        ExpCorrelation = "logical",
                        
                        CondFstat = "numeric",
                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Overdispersion = "numeric",
                        
                        PCs = "numeric",
                        Pvalue = "numeric",
                        
                        Alpha = "numeric",
                        Heter.Stat = "numeric")
)

#--------------------------------------------------------------------------------------------

#' MVGMM Class
#'
#' @description An object containing the estimates produced using the multivariable generalized method of methods (GMM) method.
#'
#' @slot robust Whether the robust model with overdispersion heterogeneity is estimated.
#' @slot Exposure The names of the exposure variables.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants if specified. If not specified, an identity matrix will be returned.
#' @slot ExpCorrelation Whether an exposure correlation matrix was specified.
#' @slot CondFstat The conditional F-statistic for each exposure. 
#' @slot Estimate The causal estimates from the GMM method.
#' @slot StdError The standard errors associated with \code{Estimate}.
#' @slot CILower The lower bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot CIUpper The upper bounds of the confidence interval for \code{Estimate} based on \code{StdError}.
#' @slot Overdispersion The estimate of the overdispersion parameter for the robust model. If this is negative, then a value of zero is used in the method.
#' @slot Pvalue P-value associated with the causal estimates.
#' @slot Alpha The significance level used in constructing confidence intervals (default is 0.05).
#' @slot Heter.Stat Heterogeneity statistic (Cochran's Q statistic) and associated p-value (for non-robust model): the null hypothesis is that all genetic variants estimate the same causal parameter; rejection of the null is an indication that one or more genetic variants may be pleiotropic.

setClass("MVGMM",
         representation(robust = "logical",
                        Exposure = "character",
                        Outcome = "character",
                        Correlation = "matrix",
                        ExpCorrelation = "logical",
                        
                        CondFstat = "numeric",
                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Overdispersion = "numeric",
                        
                        Pvalue = "numeric",
                        
                        Alpha = "numeric",
                        Heter.Stat = "numeric")
)

#--------------------------------------------------------------------------------------------

#' CLR Class
#'
#' @description An object containing the confidence intervals produced using the conditional likelihood ratio method (CLR) and other identification-robust methods.
#'
#' @slot Exposure The name of the exposure variable.
#' @slot Outcome The name of the outcome variable.
#' @slot Correlation The matrix of correlations between genetic variants.
#' @slot ARlower The lower bounds of the causal estimate based on inverting Anderson and Rubin's test.
#' @slot ARupper The upper bounds of the causal estimate based on inverting Anderson and Rubin's test.
#' @slot Klower The lower bounds of the causal estimate based on inverting Kleibergen's test.
#' @slot Kupper The upper bounds of the causal estimate based on inverting Kleibergen's test.
#' @slot CLRlower The lower bounds of the causal estimate based on inverting Moreira's conditional likelihood ratio test.
#' @slot CLRupper The upper bounds of the causal estimate based on inverting Moreira's conditional likelihood ratio test.
#' @slot CIMin The smallest value used in the search to find the confidence interval.
#' @slot CIMax The largest value used in the search to find the confidence interval.
#' @slot CIStep The step size used in the search to find the confidence interval.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).

setClass("CLR",
         representation(Exposure = "character",
                        Outcome = "character",
                        Correlation = "matrix",
                        ARlower = "numeric",
                        ARupper = "numeric",
                        Klower = "numeric",
                        Kupper = "numeric",
                        CLRlower = "numeric",
                        CLRupper = "numeric",
                        CIMin = "numeric",
                        CIMax = "numeric",
                        CIStep = "numeric",
                        Alpha = "numeric")
)
