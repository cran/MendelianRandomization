#' Extract summarized data from a PhenoScanner .csv file (legacy)
#'
#' @description The function \code{extract.pheno.csv} extracts summarized data on associations with named exposure and outcome variables from a .csv file provided by PhenoScanner.
#'
#' @details Note that this function was written for a previous version of PhenoScanner. It has not been updated, as it has been overtaken by the \code{pheno_input} function that queries PhenoScanner directly from R. 
#'
#' The PhenoScanner bioinformatic tool (\url{http://phenoscanner.medschl.cam.ac.uk/}) is a curated database of publicly available results from large-scale genetic association studies. Queries can be made for individual genetic variants (SNPs and small indels), or for multiple variants in a single batch query. One of the output files is a .csv file containing all associations of variables with each of the SNPs. For commonly genotyped variants, associations with up to 200 variables may be reported. These association estimates and their standard errors can be used in Mendelian randomization analyses.
#'
#' The plan is to enable PhenoScanner to be queried directly from the MendelianRandomization package. However, this functionality is currently unavailable.
#'
#' The \code{extract.pheno.csv} function takes the output from the web version of PhenoScanner, and converts this into an \code{MRInput} object. PhenoScanner is still under development. This function is designed for output from PhenoScanner version 1.1 (Little Miss Sunshine).
#'
#' @param exposure The name of the exposure variable.
#' @param pmidE The PubMed ID (PMID) of the publication in which the genetic association estimates with the exposure were originally reported. Some variables are reported in multiple consortia (for example, associations with coronary artery disease by CARDIoGRAM in 2011 [PMID:21378990], by CARDIoGRAMplusC4D in 2013, and again by CARDIoGRAMplusC4D in 2015 [PMID:26343387]). Equally, some publications reported associations on multiple variables (for example, CARDIoGRAMplusC4D in 2015 [PMID:26343387] reported associations with coronary artery disease and with myocardial infarction). By providing the variable name and the PubMed ID, the set of associations is (almost) uniquely identified.
#' @param ancestryE The ancestry of individuals in which estimates were obtained. A small number of studies reported genetic association estimates for a single variable in a single publication for multiple ethnicities (for example, associations with log(eGFR creatinine) from CKD-Gen in 2016 [PMID:26831199] were reported for both Europeans and Africans). The combination of exposure name, PubMed ID, and ancestry uniquely defines the set of associations. Providing the ancestry also reminds analysts of the additional complication of conducting Mendelian randomization when associations with the exposure and with the outcome are in individuals of different ancestry. Most association estimates are obtained in \code{"European"} or \code{"Mixed"} populations, although some are obtained in \code{"African"}, \code{"Asian"}, or \code{"Hispanic"} populations.
#' @param outcome The name of the outcome variable.
#' @param pmidO The PubMed ID of the publication in which the genetic association estimates with the outcome were originally reported.
#' @param ancestryO The ancestry of individuals in which genetic association estimates with the outcome were obtained.
#' @param file The file path where the PhenoScanner .csv file can be found.
#' @param snps The names (rsIDs) of the genetic variants to be included in the analysis. The default option is \code{"all"}, indicating that all the genetic variants in the .csv file with beta-coefficients and standard errors for their associations with the risk factor and with the outcome should be used in the analysis. Otherwise, only variants whose names are included in the vector of character strings provided as \code{snps} will be included in the analysis.
#' @param rsq.proxy A proxy variant is a genetic variant in close correlation (high linkage disequilibrium) with the named variant. If PhenoScanner is run with proxies included, then proxies can be included in the analysis. In the second example below, with log(eGFR creatinine) as the exposure and Tanner stage as the outcome, the association of variant rs12785878 with the outcome is not reported. Instead, rs4944958 is used as a proxy for rs12785878. The association of rs4944958 with the outcome is used in the resulting \code{MRInput} object. The correlation between the two variants is reported as \code{R^2 = 1.000}. A message will always appear when a proxy variant is included in an analysis in place of the primary variant. The value of \code{rsq.proxy} is used as a threshold in the analysis; a variant is only included in the analysis if the value of \code{R^2} equals or excedes this threshold. The default option is \code{rsq.proxy = 1}, meaning that only perfect proxies are used in the analysis.
#'
#' @return The output of the \code{extract.pheno.csv} function is an \code{MRInput} object that can be used directly in any of the estimation functions (such as \code{mr_ivw}) or in the plotting function \code{mr_plot}. The output contains:
#'
#' \item{bx}{The genetic associations with the exposure.}
#' \item{bxse}{The corresponding standard errors.}
#' \item{by}{The genetic associations with the outcome.}
#' \item{byse}{The corresponding standard errors.}
#' \item{correlation}{The matrix of genetic correlations. Currently, this is set to the empty matrix (\code{matrix()}), meaning that only uncorrelated variants can be used in the \code{extract.pheno.csv} function.}
#' \item{exposure}{A character string giving the name of the exposure as provided in the PhenoScanner database.}
#' \item{outcome}{A character string giving the name of the outcome as provided in the PhenoScanner database.}
#' \item{snps}{A vector of character strings with the names of the genetic variants.}
#'
#' @examples path.noproxy <- system.file("extdata", "vitD_snps_PhenoScanner.csv",
#' package = "MendelianRandomization")
#' path.proxies <- system.file("extdata", "vitD_snps_PhenoScanner_proxies.csv",
#'  package = "MendelianRandomization")
#'  # these two files from PhenoScanner are provided
#'  # as part of the MendelianRandomization package
#'
#' extract.pheno.csv(
#'  exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
#'  outcome = "Tanner stage", pmidO = 24770850, ancestryO = "European", 
#'  file = path.noproxy)
#'
#' extract.pheno.csv(
#'  exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
#'  outcome = "Tanner stage", pmidO = 24770850, ancestryO = "European",
#'  rsq.proxy = 0.6, file = path.proxies)
#'
#' extract.pheno.csv(
#'  exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
#'  outcome = "Asthma", pmidO = 20860503, ancestryO = "European",
#'  rsq.proxy = 0.6, file = path.proxies)
#'
#' @references James R Staley, James Blackshow, Mihir A Kamat, Steve Ellis, Prvaeen Surendran, Benjamin B Sun, Dirk S Paul, Daniel Freitag, Stephen Burgess, John Danesh, Robin Young, and Adam S Butterworth. PhenoScanner: a database of human genotype--phenotype associations. Bioinformatics 2016. doi: 10.1093/bioinformatics/btw373.
#'
#' @export

extract.pheno.csv <- function(exposure, pmidE, ancestryE,
                              outcome, pmidO, ancestryO,
                              file,
                              rsq.proxy = 1, snps = "all"){

  dataTable <- read.csv(file, header = T, sep = ",", stringsAsFactors=FALSE)

if (length(dataTable$r2)>0) {
  snp.list.exposure = unique(dataTable[which(dataTable$Trait == exposure & dataTable$PMID == pmidE & dataTable$Ancestry == ancestryE & !is.na(dataTable$Beta) & !is.na(dataTable$SE) & dataTable$r2 >= rsq.proxy),1])
  snp.list.outcome  = unique(dataTable[which(dataTable$Trait == outcome & dataTable$PMID == pmidO & dataTable$Ancestry == ancestryO & !is.na(dataTable$Beta) & !is.na(dataTable$SE) & dataTable$r2 >= rsq.proxy),1])
  }
if (length(dataTable$r2)==0) {
  snp.list.exposure = unique(dataTable[which(dataTable$Trait == exposure & dataTable$PMID == pmidE & dataTable$Ancestry == ancestryE & !is.na(dataTable$Beta) & !is.na(dataTable$SE)),1])
  snp.list.outcome  = unique(dataTable[which(dataTable$Trait == outcome & dataTable$PMID == pmidO & dataTable$Ancestry == ancestryO & !is.na(dataTable$Beta) & !is.na(dataTable$SE)),1])
  }
  snp.list = intersect(snp.list.exposure, snp.list.outcome)
if (snps[1]!="all") { snp.list = intersect(snp.list, snps) }
if (length(snp.list) == 0) { cat("No variants found with beta-coefficients and standard errors for given risk factor and outcome combination. Please check spelling and PMIDs.\n"); return() }

row.exp = NULL; row.out = NULL
for (j in 1:length(snp.list)) {
 if (length(dataTable$r2)>0) {
  row.exp[j] = which(dataTable$Trait == exposure & dataTable$PMID == pmidE & dataTable$Ancestry == ancestryE & dataTable$r2 >= rsq.proxy & !is.na(dataTable$Beta) & !is.na(dataTable$SE) & dataTable$SNP == snp.list[j])[1]
  row.out[j] = which(dataTable$Trait == outcome  & dataTable$PMID == pmidO & dataTable$Ancestry == ancestryO & dataTable$r2 >= rsq.proxy & !is.na(dataTable$Beta) & !is.na(dataTable$SE) & dataTable$SNP == snp.list[j])[1]
  if (as.character(dataTable$Proxy.rsID[row.exp[j]]) != as.character(dataTable$SNP[row.exp[j]])) { cat("Variant ", as.character(dataTable$Proxy.rsID[row.exp[j]]), " used as a proxy for ", dataTable$SNP[row.exp[j]], " in association with the exposure. R^2 value: ", decimals(dataTable$r2[row.exp[j]],3), ".\n", sep="") }
  if (as.character(dataTable$Proxy.rsID[row.out[j]]) != as.character(dataTable$SNP[row.out[j]])) { cat("Variant ", as.character(dataTable$Proxy.rsID[row.out[j]]), " used as a proxy for ", dataTable$SNP[row.out[j]], " in association with the outcome. R^2 value: ", decimals(dataTable$r2[row.out[j]],3), ".\n", sep="") }
  }
 if (length(dataTable$r2)==0) {
  row.exp[j] = which(dataTable$Trait == exposure & dataTable$PMID == pmidE & dataTable$Ancestry == ancestryE & !is.na(dataTable$Beta) & !is.na(dataTable$SE) & dataTable$SNP == snp.list[j])[1]
  row.out[j] = which(dataTable$Trait == outcome  & dataTable$PMID == pmidO & dataTable$Ancestry == ancestryO & !is.na(dataTable$Beta) & !is.na(dataTable$SE) & dataTable$SNP == snp.list[j])[1]
  }
 }
  

  Bx. <- dataTable[row.exp, c("SNP", "Beta", "SE")]
  By. <- dataTable[row.out, c("SNP", "Beta", "SE")]

  dataSet <- merge(Bx., By., "SNP")

  return(mr_input(exposure = exposure, outcome = outcome, snps = as.character(dataSet[,1]),
                  bx=as.numeric(dataSet[,2]), bxse=as.numeric(dataSet[,3]), by=as.numeric(dataSet[,4]), byse=as.numeric(dataSet[,5]),
                  correlation = matrix()))
}

