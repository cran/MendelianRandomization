% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/mr_lasso-methods.R
\docType{methods}
\name{mr_lasso}
\alias{mr_lasso}
\alias{mr_lasso,MRInput-method}
\title{MR-Lasso method}
\usage{
mr_lasso(object, distribution = "normal", alpha = 0.05, lambda = numeric(0))

\S4method{mr_lasso}{MRInput}(object, distribution = "normal", alpha = 0.05, lambda = numeric(0))
}
\arguments{
\item{object}{An \code{MRInput} object.}

\item{distribution}{The type of distribution used to calculate the confidence intervals. Options are \code{"normal"} (default) or \code{"t-dist"}.}

\item{alpha}{The significance level used to calculate the confidence intervals. The default value is 0.05.}

\item{lambda}{The value of the tuning parameter used by the lasso procedure which controls the level of sparsity. If not specified, the tuning parameter will be calculated by the heterogeneity stopping rule.}
}
\value{
The output from the function is an \code{MRLasso} object containing:

 \item{Exposure}{A character vector with the names given to the exposure.}
 \item{Outcome}{A character string with the names given to the outcome.}
 \item{Estimate}{The causal estimate from the MR-Lasso method. This is the post-lasso estimate.}
 \item{StdError}{The standard error of the causal estimate from the MR-Lasso method.}
 \item{CILower}{The lower bound of the confidence interval for the causal estimate based on the estimated standard error and the significance level provided.}
 \item{CIUpper}{The upper bound of the confidence interval for the causal estimate based on the estimated standard error and the significance level provided.}
 \item{Alpha}{The significance level used when calculating the confidence intervals.}
 \item{Pvalue}{The p-value associated with the causal estimate using a normal or t-distribution (as specified in \code{distribution}).}
 \item{SNPs}{The number of genetic variants (SNPs) included in the analysis.}
 \item{RegEstimate}{The estimate from the regularized regression model used in the MR-Lasso method.}
 \item{RegIntercept}{The intercept estimates from the regularized regression model used in the MR-Lasso method.}
 \item{Valid}{The number of genetic variants that have been identified as valid instruments.}
 \item{ValidSNPs}{The names of genetic variants that have been identified as valid instruments.}
 \item{Lambda}{The value of the tuning parameter used to compute \code{RegEstimate}}
}
\description{
The \code{mr_lasso} function performs the MR-Lasso method, which applies lasso-type penalization to the direct effects of genetic variants on the outcome.
The causal estimate is described as a post-lasso estimate, and is obtained by performing the IVW method using only those genetic variants that are identified as valid by the lasso procedure.
}
\details{
MR-Lasso extends the IVW model to include an intercept term for each genetic variant. These intercept terms represent associations between the
genetic variants and the outcome which bypass the risk factor. The causal effect estimates are estimated by weighted linear regression where the intercept terms are subject
to lasso-type penalization. The lasso penalization will tend to shrink the intercept terms corresponding to the valid instruments to zero.

The lasso penalty relies on a tuning parameter which controls the level of sparsity. The default is to use a heterogeneity stopping rule, but a fixed value may be specified.

As part of the analysis, the genetic variants are orientated so that all of the associations with the risk factor are positive (and signs of associations with the outcome are
changed to keep the orientation consistent if required). Re-orientation of the genetic variants is performed automatically as part of the function.

The MR-Lasso method is performed in two steps. First, a regularized regression model is fitted, and some genetic variants are identified as valid instruments. Second, the causal effect is estimated using standard IVW with only the valid genetic variants.
The post-lasso method will be performed as long as at least two genetic variants are identified as valid instruments. The default heterogeneity stopping rule will always return at least two
genetic variants as valid instruments.
The main estimate given by the method is the post-lasso estimate. However, parameter estimates from the regularized regression model used to identify invalid variants are also provided for completeness.

If a substantial proportion of genetic variants are removed from the analysis, the MR-Lasso method may give a false impression of confidence in the causal estimate due to homogeneity of the variant-specific causal estimates amongst the remaining variants. However, it is not reasonable to claim that there is strong evidence for a causal effect after a large number of variants with heterogeneous estimates have been removed from the analysis.
}
\examples{
mr_lasso(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse))

}
\references{
Jessica MB Rees, Angela M Wood, Frank Dudbridge, Stephen Burgess. Robust methods in Mendelian randomization via penalization of heterogeneous causal estimates. PLoS ONE 2019; 14(9):e0222362
}
