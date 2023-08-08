
#--------------------------------------------------------------------------------------------

#' Data on lipid effects on coronary artery disease (uncorrelated variants)
#'
#' @description Two sets of example data are included in the package: one illustrating uncorrelated variants, and the other correlated variants. These are the data on uncorrelated variants.
#'
#' The variables \code{ldlc}, \code{hdlc}, \code{trig}, and \code{chdlodds} are the genetic associations with (respectively) LDL-cholesterol, HDL-cholesterol, triglycerides, and coronary heart disease (CHD) risk for 28 genetic variants reported by Waterworth et al (2010). The respective standard errors of the associations are given as \code{ldlcse}, \code{hdlcse}, \code{trigse}, and \code{chdloddsse}.
#' 
#' These data can be used to test out the various functions in the package.
#'
#' @references Dawn Waterworth, Sally Ricketts, ..., Manj Sandhu: Genetic variants influencing circulating lipid levels and risk of coronary artery disease. Arterioscler Thromb Vasc Biol 2010; 30:2264-227. doi: 10.1161/atvbaha.109.201020.
#'
"ldlc"

#' @rdname ldlc
"hdlc"

#' @rdname ldlc
"hdlcse"

#' @rdname ldlc
"ldlcse"

#' @rdname ldlc
"trig"

#' @rdname ldlc
"trigse"

#' @rdname ldlc
"chdlodds"

#' @rdname ldlc
"chdloddsse"

#' @rdname ldlc
"lipid_effect"

#' @rdname ldlc
"lipid_other"

#' @rdname ldlc
"lipid_eaf"

#--------------------------------------------------------------------------------------------

#' Data on effect of calcium on fasting glucose (correlated variants)
#'
#' Two sets of example data are included in the package: one illustrating uncorrelated variants, and the other correlated variants. These are the data on correlated variants.
#'
#' The variables \code{calcium}, and \code{fastgluc} are the genetic associations with calcium and fasting glucose for 6 genetic variants reported by Burgess et al (2015). The respective standard errors of the associations are given as \code{calciumse} and \code{fastglucse}. The matrix of correlations between the genetic variants is given as \code{calc.rho}.
#' 
#' These data can be used to test out the various functions in the package.
#'
#' @references Stephen Burgess, Robert A Scott, Nic J Timpson, George Davey Smith, Simon G Thompson. Using published data in Mendelian randomization: a blueprint for efficient identification of causal risk factors. Eur J Epidemiol 2015; 30(7):543-552. doi: 10.1007/s10654-015-0011-z.
#'
"calcium"

#' @rdname calcium
"calciumse"

#' @rdname calcium
"fastgluc"

#' @rdname calcium
"fastglucse"

#' @rdname calcium
"calc.rho"

#---------------

#' Course data
#'
#' Data required for Mendelian randomization course.
#'
#' @keywords internal
#'
"coursedata"

