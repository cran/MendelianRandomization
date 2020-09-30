# Shortcut for creating an MRInput object

#' Inputting and formatting data for use in causal estimation
#'
#' @description The \code{mr_input} function is required for inputting and formatting data for use in any of the estimation functions provided in this package. The \code{MRInput} class outputted by the function can also be viewed graphically using the \code{mr_plot} function.
#'
#' @param bx A numeric vector of beta-coefficient values for genetic associations with the first variable (often referred to as the exposure, risk factor, or modifiable phenotype).
#' @param bxse The standard errors associated with the beta-coefficients \code{bx}.
#' @param by A numeric vector of beta-coefficient values for genetic associations with the second variable (often referred to as the outcome). For a disease outcome, the beta coefficients are log odds estimates from logistic regression analyses.
#' @param byse The standard errors associated with the beta-coefficients in \code{by}.
#' @param correlation The matrix of correlations between genetic variants. If this variable is not provided, then we assume that genetic variants are uncorrelated.
#' @param exposure The name of the exposure variable.
#' @param outcome The name of the outcome variable.
#' @param snps The names of the genetic variants (SNPs) included in the analysis. The inputs \code{exposure}, \code{outcome}, and \code{snps} are not required, but may be useful for keeping track of various \code{MRInput} objects. They are also used by the \code{mr_plot} function.
#' @param effect_allele The name of the effect allele for each SNP. The beta-coefficients are the associations with the exposure and outcome per additional copy of the effect allele.
#' @param other_allele The name of the non-effect allele.
#' @param eaf The expected allele frequencies (numeric). The slots \code{effect_allele}, \code{other_allele}, and \code{eaf} are neither required, nor currently used in the MendelianRandomization package. They are included for future compatibility with the MR-Base suite of functions.
#'
#' @details The beta-coefficients are assumed to be estimated for uncorrelated (independent) genetic variants, although a correlation matrix can be specified if the variants are correlated in their distributions. We also assume that the beta-coefficients for associations with the exposure and with the outcome are uncorrelated (corresponding to a two-sample Mendelian randomization analysis), although correlation between associations with the exposure and with the outcome generally have little impact on causal estimates or standard errors.
#' 
#' If the four variables are not all the same length, then an error message will be reported. The analyses will still try to run, but the output may be misleading. However, in some analyses (for example, the standard IVW and MR-Egger methods), the values of \code{bxse} are not used in the analysis, and can therefore safely be omitted (provided that the other variables are correctly labelled).
#'
#' @return An MRInput object containing:
#'
#' \item{betaX}{The genetic associations with the exposure.}
#' \item{betaXse}{The corresponding standard errors.}
#' \item{betaY}{The genetic associations with the outcome.}
#' \item{betaYse}{The corresponding standard errors.}
#' \item{correlation}{The matrix of genetic correlations.}
#' \item{exposure}{A character string giving the name given to the exposure.}
#' \item{outcome}{A character string giving the name given to the outcome.}
#' \item{snps}{A vector of character strings with the names of the genetic variants.}
#' \item{effect_allele}{A vector of character strings with the names of the effect alleles.}
#' \item{other_allele}{A vector of character strings with the names of the non-effect alleles.}
#' \item{eaf}{A numeric vector with the effect allele frequencies.}
#'
#' @seealso \code{extract.pheno.csv()} for a description of how an \code{MRInput} object can be extracted from PhenoScanner (\url{http://www.phenoscanner.medschl.cam.ac.uk/}).
#' 
#'


mr_input <- function(bx = 0, bxse = 0, by = 0, byse = 0,
                     correlation = matrix(),
                     exposure = "exposure",
                     outcome = "outcome",
                     snps = "snp",
                     effect_allele = NA,
                     other_allele  = NA,
                     eaf = NA){

  if("snp" %in% snps){
    snps <- paste("snp", 1:length(bx), sep = "_")
  } else {
    snps <- snps
  }

  new(Class = "MRInput",
      betaX = bx, betaY = by, betaXse = bxse, betaYse = byse,
      correlation = correlation,
      exposure = as.character(exposure), outcome = as.character(outcome), snps = as.character(snps),
      effect_allele = as.character(effect_allele), 
      other_allele  = as.character(other_allele), 
      eaf           = as.numeric(eaf))

  }


# Shortcut for creating an MRMVInput object

#' Inputting and formatting data for use in causal estimation
#'
#' @description The \code{mr_mvinput} function is required for inputting and formatting data for use in the multivariable Mendelian randomization functions provided in this package. 
#'
#' @param bx A matrix of beta-coefficient values for genetic associations with the risk factor variables. These should be arranged so that column 1 are the beta-coefficients for risk factor 1, and row 1 are the beta-coefficients for genetic variant 1.
#' @param bxse The matrix of standard errors associated with the beta-coefficients \code{bx}.
#' @param by A numeric vector of beta-coefficient values for genetic associations with the second variable (often referred to as the outcome). For a disease outcome, the beta coefficients are log odds estimates from logistic regression analyses.
#' @param byse The vector standard errors associated with the beta-coefficients in \code{by}.
#' @param correlation The matrix of correlations between genetic variants. If this variable is not provided, then we assume that genetic variants are uncorrelated.
#' @param exposure The names of the exposure variables.
#' @param outcome The name of the outcome variable.
#' @param snps The names of the genetic variants (SNPs) included in the analysis. The inputs \code{exposure}, \code{outcome}, and \code{snps} are not required, but may be useful for keeping track of various \code{MRInput} objects. They are also used by the \code{mr_plot} function.
#' @param effect_allele The name of the effect allele for each SNP. The beta-coefficients are the associations with the exposure and outcome per additional copy of the effect allele.
#' @param other_allele The name of the non-effect allele.
#' @param eaf The expected allele frequencies (numeric). The slots \code{effect_allele}, \code{other_allele}, and \code{eaf} are neither required, nor currently used in the MendelianRandomization package. They are included for future compatibility with the MR-Base suite of functions.
#'
#' @details The beta-coefficients are assumed to be estimated for uncorrelated (independent) genetic variants, although a correlation matrix can be specified if the variants are correlated in their distributions. We also assume that the beta-coefficients for associations with the exposure and with the outcome are uncorrelated (corresponding to a two-sample Mendelian randomization analysis), although correlation between associations with the exposure and with the outcome generally have little impact on causal estimates or standard errors.
#' 
#' If the variables are not all the same length, then an error message will be reported. The analyses will still try to run, but the output may be misleading. However, in some analyses (for example, the standard IVW and MR-Egger methods), the values of \code{bxse} are not used in the analysis, and can therefore safely be omitted (provided that the other variables are correctly labelled).
#'
#' @return An MRMVInput object containing:
#'
#' \item{betaX}{The genetic associations with the exposures.}
#' \item{betaXse}{The corresponding standard errors.}
#' \item{betaY}{The genetic associations with the outcome.}
#' \item{betaYse}{The corresponding standard errors.}
#' \item{correlation}{The matrix of genetic correlations.}
#' \item{exposure}{Character strings with the names given to the exposures.}
#' \item{outcome}{A character string giving the name given to the outcome.}
#' \item{snps}{A vector of character strings with the names of the genetic variants.}
#' \item{effect_allele}{A vector of character strings with the names of the effect alleles.}
#' \item{other_allele}{A vector of character strings with the names of the non-effect alleles.}
#' \item{eaf}{A numeric vector with the effect allele frequencies.}
#'
#'


mr_mvinput <- function(bx = matrix(), bxse = matrix(), by = 0, byse = 0,
                     correlation = matrix(),
                     exposure = "exposure",
                     outcome = "outcome",
                     snps = "snp",
                     effect_allele = NA,
                     other_allele  = NA,
                     eaf = NA){

  if("exposure" %in% exposure){
    exposure <- paste("exposure", 1:dim(bx)[2], sep = "_")
  } else {
    exposure <- exposure
  }

  if("snp" %in% snps){
    snps <- paste("snp", 1:dim(bx)[1], sep = "_")
  } else {
    snps <- snps
  }

  new(Class = "MRMVInput",
      betaX = bx, betaY = by, betaXse = bxse, betaYse = byse,
      correlation = correlation,
      exposure = as.character(exposure), outcome = as.character(outcome), snps = as.character(snps),
      effect_allele = as.character(effect_allele), 
      other_allele  = as.character(other_allele), 
      eaf           = as.numeric(eaf))

  }


