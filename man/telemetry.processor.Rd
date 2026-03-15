\name{telemetry.processor}
\alias{telemetry.processor}
\title{Process Telemetry Data for Integration with SCR Models}
\description{
  Prepares raw telemetry fix data for use in integrated SCR-telemetry models
  by assigning GPS fixes to state-space grid cells and computing per-individual
  fix frequency matrices. The output is formatted for direct use in the
  \code{telemetry} argument of \code{\link{make.scrFrame}}.
}
\usage{
telemetry.processor(ssDF, teldata)
}
\arguments{
  \item{ssDF}{A list of state-space data frames, one per session, each with
    columns \code{X} and \code{Y} giving pixel coordinates. Typically the
    output of \code{\link{make.ssDF}}.}

  \item{teldata}{A list of data frames, one per session. Each data frame must
    have at minimum three columns: individual ID (column 1), X coordinate
    (column 2), and Y coordinate (column 3) of each telemetry fix.}
}
\value{
  A list with three elements:

  \item{nfreq}{A list of matrices, one per session, with dimensions
    \eqn{n_{tel} \times n_G} (telemetered individuals by grid cells). Each
    entry gives the number of telemetry fixes for that individual in that
    grid cell. This is the \code{fixfreq} object required by
    \code{\link{make.scrFrame}}.}

  \item{dlst}{A list of matrices, one per session, with dimensions
    \eqn{n_G \times n_{tel}}. Each entry gives the Euclidean distance from
    each grid cell to the mean activity center of each telemetered individual.}

  \item{sbar}{A list of data frames, one per session, giving the mean
    X and Y coordinates (centroid) of each telemetered individual's fixes.}
}
\details{
  Each telemetry fix is assigned to its nearest state-space grid cell using
  k-nearest neighbour lookup via \code{get.knnx} from the \pkg{FNN} package.
  Fix counts are then aggregated per individual per grid cell to produce the
  \code{nfreq} matrix.

  The \code{nfreq} output corresponds directly to the \code{fixfreq} slot
  required in the \code{telemetry} argument of \code{\link{make.scrFrame}},
  and should be passed as \code{telemetry = list(fixfreq = nfreq, ...)}.

  Requires the \pkg{FNN} package.
}
\author{
  Andy Royle, Dan Linden
}
\references{
  Royle, J.A., Chandler, R.B., Sun, C.C., and Fuller, A.K. (2013).
  Integrating resource selection information with spatial capture-recapture.
  \emph{Methods in Ecology and Evolution}, 4(6), 520--530.
}
\seealso{
  \code{\link{make.scrFrame}}, \code{\link{oSCR.fit}}, \code{\link{telemetry}}
}
\examples{
## teldata should be a list of data frames with columns: id, X, Y
## one data frame per session
# tel.out <- telemetry.processor(ssDF = mySSDF, teldata = myTeldata)
# scrFrame <- make.scrFrame(caphist = y3d, traps = traplocs,
#               telemetry = list(fixfreq = tel.out$nfreq,
#                                indCovs = tel.indCovs))
}
\keyword{spatial capture-recapture}
\keyword{telemetry}
