## ------------------------------------------------------------------------
library(MendelianRandomization)

## ------------------------------------------------------------------------
MRInputObject <- mr_input(bx = ldlc, 
                          bxse = ldlcse, 
                          by = chdlodds, 
                          byse = chdloddsse)

MRInputObject  # example with uncorrelated variants

MRInputObject.cor <- mr_input(bx = calcium, 
                             bxse = calciumse, 
                             by = fastgluc, 
                             byse = fastglucse,
                             corr = calc.rho)

MRInputObject.cor  # example with correlated variants


## ---- eval=FALSE---------------------------------------------------------
#  MRInputObject <- mr_input(ldlc, ldlcse, chdlodds, chdloddsse)

## ------------------------------------------------------------------------
IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)

IVWObject <- mr_ivw(mr_input(bx = ldlc, bxse = ldlcse,
   by = chdlodds, byse = chdloddsse))

IVWObject

IVWObject.correl <- mr_ivw(MRInputObject.cor,
                    model = "default",
                    correl = TRUE,
                    distribution = "normal",
                    alpha = 0.05)

IVWObject.correl <- mr_ivw(mr_input(bx = calcium, bxse = calciumse,
   by = fastgluc, byse = fastglucse, corr = calc.rho))

IVWObject.correl

## ------------------------------------------------------------------------
WeightedMedianObject <- mr_median(MRInputObject, 
                                  weighting = "weighted", 
                                  distribution = "normal", 
                                  alpha = 0.05, 
                                  iterations = 10000, 
                                  seed = 314159265)

WeightedMedianObject <- mr_median(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

WeightedMedianObject 

SimpleMedianObject <- mr_median(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse), weighting = "simple")

SimpleMedianObject

## ------------------------------------------------------------------------
EggerObject <- mr_egger(MRInputObject, 
                        robust = FALSE,
                        penalized = FALSE,
                        correl = FALSE,
                        distribution = "normal",
                        alpha = 0.05)

EggerObject <- mr_egger(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

EggerObject

EggerObject.corr <- mr_egger(MRInputObject.cor, 
                        correl = TRUE,
                        distribution = "normal",
                        alpha = 0.05)

EggerObject.corr <- mr_egger(mr_input(bx = calcium, bxse = calciumse,
  by = fastgluc, byse = fastglucse, corr = calc.rho))

EggerObject.corr

## ------------------------------------------------------------------------
MaxLikObject <- mr_maxlik(MRInputObject, 
                          correl = FALSE,
                          psi = 0,
                          distribution = "normal",
                          alpha = 0.05)

MaxLikObject <- mr_maxlik(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

MaxLikObject

MaxLikObject.corr <- mr_maxlik(mr_input(bx = calcium, bxse = calciumse,
  by = fastgluc, byse = fastglucse, corr = calc.rho))

MaxLikObject.corr

## ------------------------------------------------------------------------
MRAllObject_all <- mr_allmethods(MRInputObject, method = "all")
MRAllObject_all

MRAllObject_egger <- mr_allmethods(MRInputObject, method = "egger")
MRAllObject_egger

MRAllObject_main <- mr_allmethods(MRInputObject, method = "main")
MRAllObject_main

## ---- eval = FALSE-------------------------------------------------------
#  mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#   error = TRUE, orientate = FALSE, line = "ivw")

## ------------------------------------------------------------------------
mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
 error = TRUE, orientate = FALSE, line = "ivw", interactive = FALSE)

## ------------------------------------------------------------------------
mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
 error = TRUE, orientate = FALSE, line = "ivw", interactive = FALSE, labels = TRUE)

## ------------------------------------------------------------------------
mr_plot(MRAllObject_all)

## ------------------------------------------------------------------------
mr_plot(MRAllObject_egger)

## ------------------------------------------------------------------------
mr_plot(mr_allmethods(mr_input(bx = hdlc, bxse = hdlcse,
  by = chdlodds, byse = chdloddsse)))

## ------------------------------------------------------------------------
path.noproxy <- system.file("extdata", "vitD_snps_PhenoScanner.csv",
  package = "MendelianRandomization")
path.proxies <- system.file("extdata", "vitD_snps_PhenoScanner_proxies.csv",
  package = "MendelianRandomization")


extract.pheno.csv(
  exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
  outcome = "Tanner stage", pmidO = 24770850, ancestryO = "European", 
  file = path.noproxy)

extract.pheno.csv(
  exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
  outcome = "Tanner stage", pmidO = 24770850, ancestryO = "European",
  rsq.proxy = 0.6, 
  file = path.proxies)

extract.pheno.csv(
  exposure = "log(eGFR creatinine)", pmidE = 26831199, ancestryE = "European",
  outcome = "Asthma", pmidO = 20860503, ancestryO = "European",
  rsq.proxy = 0.6, 
  file = path.proxies)

