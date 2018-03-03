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
                          model = "default",
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
MRMVInputObject <- mr_mvinput(bx = cbind(ldlc, hdlc, trig),
                              bxse = cbind(ldlcse, hdlcse, trigse),
                              by = chdlodds, 
                              byse = chdloddsse)

MRMVInputObject

## ------------------------------------------------------------------------
MRMVObject <- mr_mvivw(MRMVInputObject, 
                          model = "default",
						  correl = FALSE,
                          distribution = "normal",
                          alpha = 0.05)

MRMVObject <- mr_mvivw(MRMVInputObject)

MRMVObject

## ------------------------------------------------------------------------
MBEObject <- mr_mbe(MRInputObject, 
                    weighting = "weighted",
					stderror = "delta",
					phi = 1,
					seed = 314159265,
					iterations = 10000,
                    distribution = "normal",
                    alpha = 0.05)

MBEObject <- mr_mbe(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

MBEObject

## ---- eval = FALSE-------------------------------------------------------
#  HetPenObject <- mr_hetpen(MRInputObject,
#                            prior = 0.5,
#  					      CIMin = -1,
#  						  CIMax = 1,
#  						  CIStep = 0.001,
#                            alpha = 0.05)

## ------------------------------------------------------------------------
HetPenObject <- mr_hetpen(mr_input(bx = ldlc[1:10], bxse = ldlcse[1:10],
  by = chdlodds[1:10], byse = chdloddsse[1:10]), CIMin = -1, CIMax = 5, CIStep = 0.01)

HetPenObject

## ------------------------------------------------------------------------
bcrp    =c(0.160, 0.236, 0.149, 0.09, 0.079, 0.072, 0.047, 0.069)
bcrpse  =c(0.006, 0.009, 0.006, 0.005, 0.005, 0.005, 0.006, 0.011)
bchd    =c(0.0237903, -0.1121942, -0.0711906, -0.030848, 0.0479207, 0.0238895,
  0.005528, 0.0214852)
bchdse  =c(0.0149064, 0.0303084, 0.0150552, 0.0148339, 0.0143077, 0.0145478,
  0.0160765, 0.0255237)
  
HetPenObject.multimode <- mr_hetpen(mr_input(bx = bcrp, bxse = bcrpse,
  by = bchd, byse = bchdse))

HetPenObject.multimode

## ------------------------------------------------------------------------
MRInputObject <- mr_input(bx = ldlc, 
                          bxse = ldlcse, 
                          by = chdlodds, 
                          byse = chdloddsse)

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

