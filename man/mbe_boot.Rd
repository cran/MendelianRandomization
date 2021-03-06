% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_mbe-methods.R
\name{mbe_boot}
\alias{mbe_boot}
\title{Mode-based estimate (Hartwig) bootstrap function}
\usage{
mbe_boot(BetaIV.in, seBetaIV.in, weighting.in, iterations.in, phi.in)
}
\arguments{
\item{BetaIV.in}{Ratio causal estimates for each genetic variant considered individually.}

\item{seBetaIV.in}{Standard errors of ratio causal estimates.}

\item{iterations.in}{Number of bootstrap iterations.}

\item{phi.in}{Bandwidth multiplication factor.}
}
\value{
Bootstrapped mode-based estimates.
}
\description{
Internal function for calculating standard error of mode-based estimate.
}
\details{
None.
}
\examples{
mbe_boot(BetaIV.in = chdlodds/ldlc, seBetaIV.in = abs(chdloddsse/ldlc),
   weighting.in = "simple", iterations.in = 100)

}
\keyword{internal}
