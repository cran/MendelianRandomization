## -----------------------------------------------------------------------------
library(MendelianRandomization) # loading the package

## -----------------------------------------------------------------------------
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


## ----eval=FALSE---------------------------------------------------------------
#  MRInputObject <- mr_input(ldlc, ldlcse, chdlodds, chdloddsse)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## ----eval = FALSE-------------------------------------------------------------
#  HetPenObject <- mr_hetpen(MRInputObject,
#                            prior = 0.5,
#                            CIMin = -1,
#                            CIMax = 1,
#                            CIStep = 0.001,
#                            alpha = 0.05)

## -----------------------------------------------------------------------------
HetPenObject <- mr_hetpen(mr_input(bx = ldlc[1:10], bxse = ldlcse[1:10],
  by = chdlodds[1:10], byse = chdloddsse[1:10]), CIMin = -1, CIMax = 5, CIStep = 0.01)

HetPenObject

## -----------------------------------------------------------------------------
bcrp    =c(0.160, 0.236, 0.149, 0.09, 0.079, 0.072, 0.047, 0.069)
bcrpse  =c(0.006, 0.009, 0.006, 0.005, 0.005, 0.005, 0.006, 0.011)
bchd    =c(0.0237903, -0.1121942, -0.0711906, -0.030848, 0.0479207, 0.0238895,
  0.005528, 0.0214852)
bchdse  =c(0.0149064, 0.0303084, 0.0150552, 0.0148339, 0.0143077, 0.0145478,
  0.0160765, 0.0255237)
  
HetPenObject.multimode <- mr_hetpen(mr_input(bx = bcrp, bxse = bcrpse,
  by = bchd, byse = bchdse))

HetPenObject.multimode

## -----------------------------------------------------------------------------
ConMixObject <- mr_conmix(MRInputObject, 
                          psi = 0,
                          CIMin = NA,
                          CIMax = NA,
                          CIStep = 0.01,
                          alpha = 0.05)

ConMixObject <- mr_conmix(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

ConMixObject

## -----------------------------------------------------------------------------
LassoObject <- mr_lasso(MRInputObject, 
                        distribution = "normal",
                        alpha = 0.05,
                        lambda = numeric(0))

LassoObject <- mr_lasso(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

LassoObject

## -----------------------------------------------------------------------------
DIVWObject <- mr_divw(MRInputObject,
                      over.dispersion = TRUE,
                      alpha = 0.05,
                      diagnostics = FALSE)

DIVWObject <- mr_divw(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

DIVWObject

## -----------------------------------------------------------------------------
PIVWObject <- mr_pivw(MRInputObject,
                      over.dispersion = TRUE,
                      delta = 0, 
                      sel.pval = NULL,
                      Boot.Fieller = TRUE,
                      alpha = 0.05)

PIVWObject <- mr_pivw(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse))

PIVWObject

## ----eval=FALSE---------------------------------------------------------------
#  cMLObject <- mr_cML(MRInputObject,
#                      MA = TRUE,
#                      DP = TRUE,
#                      K_vec = 0:(length(object@betaX)-2),
#                      random_start = 0,
#                      num_pert = 200,
#                      random_start_pert = 0,
#                      maxit = 100,
#                      random_seed = 314,
#                      n,
#                      Alpha = 0.05)

## -----------------------------------------------------------------------------
cMLObject <- mr_cML(mr_input(bx = ldlc, bxse = ldlcse,
  by = chdlodds, byse = chdloddsse), n = 17723)

cMLObject

## ----eval=FALSE---------------------------------------------------------------
#  pcGMMObject <- mr_pcgmm(MRInputObject.cor,
#                          nx,
#                          ny,
#                          r = NULL,
#                          thres = 0.999,
#                          robust = TRUE,
#                          alpha = 0.05)

## -----------------------------------------------------------------------------
pcGMMObject <- mr_pcgmm(mr_input(bx = calcium, bxse = calciumse,
   by = fastgluc, byse = fastglucse, corr = calc.rho),
   nx=6351, ny=133010)

pcGMMObject

## -----------------------------------------------------------------------------
MVMRInputObject <- mr_mvinput(bx = cbind(ldlc, hdlc, trig),
                              bxse = cbind(ldlcse, hdlcse, trigse),
                              by = chdlodds, 
                              byse = chdloddsse)

MVMRInputObject

MVMRInputObject.cor <- mr_mvinput(bx = cbind(ldlc, hdlc, trig),
                                  bxse = cbind(ldlcse, hdlcse, trigse),
                                  by = chdlodds, 
                                  byse = chdloddsse,
                                  correlation = diag(length(ldlc)))

## -----------------------------------------------------------------------------
MVIVWObject <- mr_mvivw(MVMRInputObject, 
                         model = "default",
                         correl = FALSE,
                         correl.x = NULL,
                         nx = NA,
                         distribution = "normal",
                         alpha = 0.05)

MVIVWObject <- mr_mvivw(MVMRInputObject)

MVIVWObject

## -----------------------------------------------------------------------------
MVIVWObject.condF <- mr_mvivw(MVMRInputObject, nx = 17723)

MVIVWObject.condF

## -----------------------------------------------------------------------------
MVEggerObject <- mr_mvegger(MVMRInputObject)

MVEggerObject

MVMedianObject <- mr_mvmedian(MVMRInputObject)

MVMedianObject

MVLassoObject <- mr_mvlasso(MVMRInputObject)

MVLassoObject

MVcMLObject <- mr_mvcML(MVMRInputObject, n = 17723)

MVcMLObject

MVGMMObject <- mr_mvgmm(MVMRInputObject, nx=rep(17723,3), ny=17723)

MVGMMObject

MVpcGMMObject <- mr_mvpcgmm(MVMRInputObject.cor, nx=rep(17723,3), ny=17723)

MVpcGMMObject

## -----------------------------------------------------------------------------
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

## ----eval = FALSE-------------------------------------------------------------
#  mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
#   error = TRUE, orientate = FALSE, line = "ivw")

## -----------------------------------------------------------------------------
mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
 error = TRUE, orientate = FALSE, line = "ivw", interactive = FALSE)

## -----------------------------------------------------------------------------
mr_plot(mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse),
 error = TRUE, orientate = FALSE, line = "ivw", interactive = FALSE, labels = TRUE)

## -----------------------------------------------------------------------------
mr_plot(MVMRInputObject, interactive = FALSE)

## -----------------------------------------------------------------------------
mr_plot(MRAllObject_all)

## -----------------------------------------------------------------------------
mr_plot(MRAllObject_egger)

## -----------------------------------------------------------------------------
mr_plot(mr_allmethods(mr_input(bx = hdlc, bxse = hdlcse,
  by = chdlodds, byse = chdloddsse)))

## -----------------------------------------------------------------------------
mr_forest(MRInputObject, ordered=TRUE)

## -----------------------------------------------------------------------------
mr_forest(MRInputObject,
          methods = c("ivw", "median", "wmedian", "egger", "maxlik", "mbe", "conmix"),
          snp_estimates = FALSE)

## -----------------------------------------------------------------------------
 mr_funnel(mr_input(bx = ldlc[1:8], bxse = ldlcse[1:8],
            by = chdlodds[1:8], byse = chdloddsse[1:8]))

## -----------------------------------------------------------------------------
mr_loo(MRInputObject)

## ----eval=FALSE---------------------------------------------------------------
#  library(ggplot2)
#  forest  = mr_forest(mr_input(ldlc, ldlcse, chdlodds, chdloddsse))
#  forest2 = forest + coord_cartesian(xlim=c(-5,5))
#  forest2

## ----eval = FALSE-------------------------------------------------------------
#  mr_ivw(pheno_input(snps=c("rs12916", "rs2479409",
#                            "rs217434", "rs1367117",
#                            "rs4299376", "rs629301",
#                            "rs4420638", "rs6511720"),
#                     exposure = "Low density lipoprotein",
#                     pmidE = "24097068",
#                     ancestryE = "European",
#                     outcome = "Coronary artery disease",
#                     pmidO = "26343387",
#                     ancestryO = "Mixed"))

