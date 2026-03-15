\name{connectivity.surface}
\alias{connectivity.surface}
\title{Compute a Connectivity Surface from a Fitted Cost-Surface SCR Model}
\description{
  Computes a landscape connectivity surface using a fitted SCR model that
  includes an ecological distance (least-cost path) component. The connectivity
  surface represents the expected movement-weighted density across the state
  space, integrating over individual activity centers weighted by the
  cost-distance-based detection kernel.
}
\usage{
connectivity.surface(scr.fit, scrFrame, ssDF, costDF = NULL,
                     PROJ = NULL, directions = 8)
}
\arguments{
  \item{scr.fit}{A fitted model object of class \code{oSCR.fit}, fitted using
    \code{distmet = "ecol"} (ecological distance) so that cost-surface
    parameters are present in \code{scr.fit$outStats}.}

  \item{scrFrame}{The \code{scrFrame} object used to fit the model.}

  \item{ssDF}{The state-space object used to fit the model. A list of data
    frames, one per session, with columns \code{X} and \code{Y}.}

  \item{costDF}{A list of data frames, one per session, containing the cost
    surface covariate values. Must match the structure used when fitting
    the model.}

  \item{PROJ}{A projection string specifying the coordinate reference system,
    passed to \code{costDistance} in the \pkg{gdistance} package. Defaults
    to \code{"+proj=utm +zone=12 +datum=WGS84"} if \code{NULL}.}

  \item{directions}{Number of directions used in computing transition matrices
    for least-cost paths. Either 4, 8, or 16. Default is 8.}
}
\value{
  A list with two elements:
  \item{dwc}{A list of matrices, one per session. Each matrix contains the
    connectivity surface values — the movement-weighted expected density at
    each state-space pixel, integrating over activity center posteriors and
    the cost-distance-based detection kernel, separately for each sex.}
  \item{Ec}{Reserved for future use; currently returns \code{NULL}.}
}
\details{
  The connectivity surface is computed by constructing a cost-distance matrix
  between all state-space pixels using the fitted cost surface parameters,
  then weighting the predicted realized density surface (\code{predict.oSCR})
  by the movement kernel (half-normal based on ecological distance). For
  sex-structured models, separate kernels are computed for each sex using
  their respective sigma estimates, and the result is a sex-ratio-weighted
  combination.

  Requires the \pkg{raster} and \pkg{gdistance} packages.
}
\author{
  Andy Royle
}
\references{
  Royle, J.A., Chandler, R.B., Sun, C.C., and Fuller, A.K. (2013).
  Integrating resource selection information with spatial capture-recapture.
  \emph{Methods in Ecology and Evolution}, 4(6), 520--530.
}
\seealso{
  \code{\link{oSCR.fit}}, \code{\link{predict.oSCR}}, \code{\link{make.ssDF}}
}
\keyword{spatial capture-recapture}
\keyword{landscape connectivity}
