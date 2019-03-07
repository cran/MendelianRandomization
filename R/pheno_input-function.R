#' PhenoScanner
#'
#' The \code{phenoscanner} function queries the PhenoScanner database of genotype-phenotype associations from inside R.
#' @param snpquery a vector of SNPs.
#' @param regionquery a vector of genomic regions.
#' @param genequery a vector of gene names.
#' @param catalogue the catalogue to be searched (options: None, GWAS, eQTL, pQTl, mQTL, methQTL).
#' @param pvalue the p-value threshold.
#' @param proxies the proxies database to be searched (options: None, AFR, AMR, EAS, EUR, SAS).
#' @param r2 the r2 threshold. 
#' @param build the genome build (options: 37, 38).
#' @import rjson
#' @return a list containing a data.frame of association results and a data.frame of SNP/Region/Gene information from PhenoScanner.
#' @examples
#' # SNP
#' res <- phenoscanner(snpquery="rs10840293")
#' head(res$results)
#' res$snps
#' 
#' # Gene
#' res <- phenoscanner(genequery="SWAP70")
#' head(res$results)
#' res$snps
#' 
#' # Region
#' res <- phenoscanner(regionquery="chr11:9685624-9774538")
#' head(res$results)
#' res$regions
#' @author PhenoScanner <phenoscanner@gmail.com>
#' @export

phenoscanner <- function(snpquery=NULL, genequery=NULL, regionquery=NULL, catalogue="GWAS", pvalue=1E-5, proxies="None", r2=0.8, build=37){
  if(is.null(snpquery) & is.null(regionquery) & is.null(genequery)) stop("no query has been requested")
  if((length(snpquery[1])+length(regionquery[1])+length(genequery[1]))>1) stop("only one query type allowed")
  if(!(catalogue=="None" | catalogue=="GWAS" | catalogue=="eQTL" | catalogue=="pQTL" | catalogue=="mQTL" | catalogue=="methQTL")) stop("catalogue has to be one of None, GWAS, eQTL, pQTL, mQTL or methQTL")
  if(!(proxies=="None" | proxies=="AFR" | proxies=="AMR" | proxies=="EAS" | proxies=="EUR" | proxies=="SAS")) stop("proxies has to be one of None, AFR, AMR, EAS, EUR or SAS")
  if(length(snpquery)>100) stop("a maximum of 100 SNP queries can be requested at one time")
  if(length(genequery)>10) stop("a maximum of 10 gene queries can be requested at one time")
  if(length(regionquery)>10) stop("a maximum of 10 region queries can be requested at one time")
  if(!(pvalue>0 & pvalue<=1)) stop("the p-value threshold has to be greater than 0 and less than or equal to 1")
  if(!(r2>=0.5 & r2<=1)) stop("the r2 threshold has to be greater than or equal to 0.5 and less than or equal to 1")
  if(!(build==37 | build==38)) stop("the build has to be equal to 37 or 38")
  if(!is.null(regionquery)){
    ub <- as.numeric(sub(".*-", "", sub(".*:", "",regionquery)))
    lb <- as.numeric(sub("-.*", "", sub(".*:", "",regionquery)))
    dist <- ub - lb
    if(any(dist>1000000)) stop("region query can be maximum of 1MB in size")
  }
  if(length(snpquery)>0){
    results <- data.frame()
    snps <- data.frame()
    n_queries <- length(snpquery) %/% 10
    if((length(snpquery) %% 10)>0){n_queries <- n_queries + 1}
    for(i in 1:n_queries){
      if(i < n_queries){qsnps <- paste0(snpquery[(1+10*(i-1)):(10*i)], collapse="+")}else{qsnps <- paste0(snpquery[(1+10*(i-1)):length(snpquery)], collapse="+")}
      json_file <- paste0("http://www.phenoscanner.medschl.cam.ac.uk/api/?snpquery=",qsnps,"&catalogue=",catalogue,"&p=",pvalue,"&proxies=",proxies,"&r2=",r2,"&build=",build)
      json_data <- fromJSON(file=json_file)
      if(length(json_data$results)==0 & length(json_data$snps)==0){
        print(paste0("Error: ",json_data$error))
        next
      }
      if(length(json_data$results)>0){
        fields <- json_data$results[[1]]; json_data$results[[1]] <- NULL
        if(length(json_data$results)>0){
          tables <- as.data.frame(matrix(unlist(json_data$results), ncol=length(fields), byrow=T), stringsAsFactors=F)
          names(tables) <- fields
          results <- rbind(results,tables)
          if(length(snpquery)==1){print(paste0(snpquery," -- queried"))}else{print(paste0(i," -- chunk of 10 SNPs queried"))}
        }else{if(length(snpquery)==1){print(paste0("Error: no results found for ",snpquery))}else{print(paste0("Error: no results found for chunk ",i))}}
      }
      if(length(json_data$snps)>0){
        fields_snps <- json_data$snps[[1]]; json_data$snps[[1]] <- NULL
        if(length(json_data$snps)>0){
          tables_snps <- as.data.frame(matrix(unlist(json_data$snps), ncol=length(fields_snps), byrow=T))
          names(tables_snps) <- fields_snps
          snps <- rbind(snps,tables_snps)
        }
      }
    }
    output <- list(snps=snps, results=results)
  }
  if(length(genequery)>0){
    results <- data.frame()
    genes <- data.frame()
    n_queries <- length(genequery)
    for(i in 1:n_queries){
      json_file <- paste0("http://www.phenoscanner.medschl.cam.ac.uk/api/?genequery=",genequery[i],"&catalogue=",catalogue,"&p=",pvalue,"&proxies=None&r2=1&build=",build)
      json_data <- fromJSON(file=json_file)
      if(length(json_data$results)==0 & length(json_data$genes)==0){
        print(paste0("Error: ", json_data$error))
        next
      }
      if(length(json_data$results)>0){
        fields <- json_data$results[[1]]; json_data$results[[1]] <- NULL
        if(length(json_data$results)>0){
          tables <- as.data.frame(matrix(unlist(json_data$results), ncol=length(fields), byrow=T), stringsAsFactors=F)
          names(tables) <- fields
          results <- rbind(results,tables)
          print(paste0(genequery," -- queried"))
        }else{print(paste0("Error: no results found for ",genequery))}
      }
      if(length(json_data$genes)>0){
        fields_genes <- json_data$genes[[1]]; json_data$genes[[1]] <- NULL
        if(length(json_data$genes)>0){
          tables_genes <- as.data.frame(matrix(unlist(json_data$genes), ncol=length(fields_genes), byrow=T))
          names(tables_genes) <- fields_genes
          genes <- rbind(genes,tables_genes)
        }
      }
    }
    output <- list(genes=genes, results=results)
  }
  if(length(regionquery)>0){
    results <- data.frame()
    regions <- data.frame()
    n_queries <- length(regionquery)
    for(i in 1:n_queries){
      json_file <- paste0("http://www.phenoscanner.medschl.cam.ac.uk/api/?regionquery=",regionquery[i],"&catalogue=",catalogue,"&p=",pvalue,"&proxies=None&r2=1&build=",build)
      json_data <- fromJSON(file=json_file)
      if(length(json_data$results)==0 & length(json_data$locations)==0){
        print(paste0("Error: ",json_data$error))
        next
      }
      if(length(json_data$results)>0){
        fields <- json_data$results[[1]]; json_data$results[[1]] <- NULL
        if(length(json_data$results)>0){
          tables <- as.data.frame(matrix(unlist(json_data$results), ncol=length(fields), byrow=T), stringsAsFactors=F)
          names(tables) <- fields
          results <- rbind(results,tables)
          print(paste0(regionquery," -- queried"))
        }else{print(paste0("Error: no results found for ",regionquery))}
      }
      if(length(json_data$locations)>0){
        fields_regions <- json_data$locations[[1]]; json_data$locations[[1]] <- NULL
        if(length(json_data$locations)>0){
          tables_regions <- as.data.frame(matrix(unlist(json_data$locations), ncol=length(fields_regions), byrow=T))
          names(tables_regions) <- fields_regions
          regions <- rbind(regions,tables_regions)
        }
      }
    }
    output <- list(regions=regions, results=results)
  }
  if(is.null(output)) stop("there is no output")
  return(output)
}

#' Extract summarized data from PhenoScanner
#'
#' @description The function \code{pheno_input} extracts summarized data on associations with named exposure and outcome variables from PhenoScanner.
#'
#' @details The PhenoScanner bioinformatic tool (\url{http://phenoscanner.medschl.cam.ac.uk}) is a curated database of publicly available results from large-scale genetic association studies. Queries can be made for individual genetic variants (SNPs and small indels), or for multiple variants in a single batch query. These association estimates and their standard errors can be used in Mendelian randomization analyses.
#'
#' The \code{phenoscanner} command is included in the \code{MendelianRandomization} package with permission of James Staley. The function is also available in a standalone package from github: \url{https://github.com/phenoscanner/phenoscanner}.
#'
#' @param snps The names (rsid) of the genetic variants to be included in the analysis.
#' @param exposure The name of the exposure variable.
#' @param pmidE The PubMed ID (PMID) of the publication in which the genetic association estimates with the exposure were originally reported. Some variables are reported in multiple consortia (for example, associations with coronary artery disease by CARDIoGRAM in 2011 [PMID:21378990], by CARDIoGRAMplusC4D in 2013, and again by CARDIoGRAMplusC4D in 2015 [PMID:26343387]). Equally, some publications reported associations on multiple variables (for example, CARDIoGRAMplusC4D in 2015 [PMID:26343387] reported associations with coronary artery disease and with myocardial infarction). By providing the variable name and the PubMed ID, the set of associations is (almost) uniquely identified.
#' @param ancestryE The ancestry of individuals in which estimates were obtained. A small number of studies reported genetic association estimates for a single variable in a single publication for multiple ethnicities (for example, associations with log(eGFR creatinine) from CKD-Gen in 2016 [PMID:26831199] were reported for both Europeans and Africans). The combination of exposure name, PubMed ID, and ancestry uniquely defines the set of associations. Providing the ancestry also reminds analysts of the additional complication of conducting Mendelian randomization when associations with the exposure and with the outcome are in individuals of different ancestry. Most association estimates are obtained in \code{"European"} or \code{"Mixed"} populations, although some are obtained in \code{"African"}, \code{"Asian"}, or \code{"Hispanic"} populations.
#' @param outcome The name of the outcome variable.
#' @param pmidO The PubMed ID of the publication in which the genetic association estimates with the outcome were originally reported.
#' @param ancestryO The ancestry of individuals in which genetic association estimates with the outcome were obtained.
#' @param correl The correlations between the genetic variants. If this is not specified, then the genetic variants are assumed to be uncorrelated. Note that for the correlations to reference the correct variants, the list of genetic variants needs to be in alphabetical order.
#'
#' @return The output of the \code{pheno_input} function is an \code{MRInput} object that can be used directly in any of the estimation functions (such as \code{mr_ivw}) or in the plotting function \code{mr_plot}. The output contains:
#'
#' \item{bx}{The genetic associations with the exposure.}
#' \item{bxse}{The corresponding standard errors.}
#' \item{by}{The genetic associations with the outcome.}
#' \item{byse}{The corresponding standard errors.}
#' \item{correlation}{The matrix of genetic correlations as specified by the user.}
#' \item{exposure}{A character string giving the name of the exposure as provided in the PhenoScanner database.}
#' \item{outcome}{A character string giving the name of the outcome as provided in the PhenoScanner database.}
#' \item{snps}{A vector of character strings with the names of the genetic variants.}
#'
#' @examples # pheno_input(snps=c("rs12916", "rs2479409", "rs217434", "rs1367117"),
#' # exposure = "Low density lipoprotein", pmidE = "24097068", ancestryE = "European",
#' # outcome = "Coronary artery disease", pmidO = "26343387", ancestryO = "Mixed")
#'
#' @references James R Staley, James Blackshow, Mihir A Kamat, Steve Ellis, Prvaeen Surendran, Benjamin B Sun, Dirk S Paul, Daniel Freitag, Stephen Burgess, John Danesh, Robin Young, and Adam S Butterworth. PhenoScanner: a database of human genotype--phenotype associations. Bioinformatics 2016. doi: 10.1093/bioinformatics/btw373.
#'
#' @export

pheno_input <- function(snps, exposure, pmidE, ancestryE, outcome, pmidO, ancestryO, correl=NULL) {

 dataTable <- phenoscanner(snpquery = snps, pvalue = 1)$results

 snp.list.exposure = unique(dataTable[which(dataTable$trait == exposure & dataTable$pmid == pmidE & dataTable$ancestry == ancestryE & !is.na(dataTable$beta) & !is.na(dataTable$se)),1])
 snp.list.outcome = unique(dataTable[which(dataTable$trait == outcome & dataTable$pmid == pmidO & dataTable$ancestry == ancestryO & !is.na(dataTable$beta) & !is.na(dataTable$se)),1])
 snp.list = intersect(snp.list.exposure, snp.list.outcome)
 if (length(snp.list) == 0) { cat("No variants found with beta-coefficients and standard errors for given risk factor and outcome combination. Please check spelling and PMIDs.\n"); return() }

 row.exp = NULL; row.out = NULL
 for (j in 1:length(snp.list)) {
  row.exp[j] = which(dataTable$trait == exposure & dataTable$pmid == pmidE & dataTable$ancestry == ancestryE & !is.na(dataTable$beta) & !is.na(dataTable$se) & dataTable$snp == snp.list[j])[1]
  row.out[j] = which(dataTable$trait == outcome & dataTable$pmid == pmidO & dataTable$ancestry == ancestryO & !is.na(dataTable$beta) & !is.na(dataTable$se) & dataTable$snp == snp.list[j])[1]
 }

 Bx. <- dataTable[row.exp, c("snp", "beta", "se")]
 By. <- dataTable[row.out, c("snp", "beta", "se")]

 dataSet <- merge(Bx., By., "snp")

if (is.matrix(correl)) {
  return(mr_input(exposure = exposure, outcome = outcome, snps = as.character(dataSet[,1]),
 bx=as.numeric(dataSet[,2]), bxse=as.numeric(dataSet[,3]), by=as.numeric(dataSet[,4]), byse=as.numeric(dataSet[,5]),
 correlation = correl)) }
 else {
 return(mr_input(exposure = exposure, outcome = outcome, snps = as.character(dataSet[,1]),
 bx=as.numeric(dataSet[,2]), bxse=as.numeric(dataSet[,3]), by=as.numeric(dataSet[,4]), byse=as.numeric(dataSet[,5]),
 correlation = matrix())) }
}
