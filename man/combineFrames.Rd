\name{combineFrames}
\alias{combineFrames}
\title{Combine Multiple scrFrame Objects into a Multi-Session scrFrame}
\description{
  This function is the inverse of \code{\link{subFrame}}: it takes a list of
  \code{scrFrame} objects and combines them into a single multi-session
  \code{scrFrame} suitable for use with \code{\link{oSCR.fit}}. Sessions are
  appended in the order they appear in \code{sfList}. Each input frame may
  itself be single- or multi-session.
}
\usage{
combineFrames(sfList)
}
\arguments{
  \item{sfList}{A list of objects of class \code{scrFrame}, one per session
    (or group of sessions). All elements must be of class \code{scrFrame}.}
}
\value{
  An object of class \code{scrFrame} containing all sessions from the input
  list, with the following components:

  \item{caphist}{List of \eqn{n \times J \times K} capture history arrays,
    one per session.}

  \item{traps}{List of trap coordinate data frames (columns \code{X} and
    \code{Y}), one per session.}

  \item{indCovs}{List of individual covariate data frames, one per session,
    or \code{NULL} if absent in all input frames.}

  \item{trapCovs}{List of lists of occasion-level trap covariate data frames,
    one per session, or \code{NULL} if absent in all input frames.}

  \item{trapOperation}{List of \eqn{J \times K} binary trap operation
    matrices, one per session, or \code{NULL} if absent in all input frames.}

  \item{telemetry}{Combined telemetry object (\code{fixfreq}, \code{indCovs},
    \code{cap.tel}), or \code{NULL} if absent in all input frames.}

  \item{occasions}{Integer vector of length equal to the number of sessions,
    giving the number of sampling occasions per session. Required internally
    by \code{do.trim} when \code{trimS} is specified in \code{oSCR.fit}.}

  \item{nocc}{Maximum number of occasions across all sessions.}

  \item{mmdm}{Mean maximum distance moved, pooled across all individuals and
    sessions. \code{NULL} if \code{mmdm} is absent from all input frames. Used
    by \code{oSCR.fit} as a starting value for \code{sigma}.}

  \item{type}{Inherited from the first element of \code{sfList}.}
}
\details{
  Components that are \code{NULL} in some input frames but present in others
  are padded with \code{NULL} for the missing sessions, consistent with the
  behaviour of \code{\link{make.scrFrame}}. A warning is issued when this
  occurs.

  The \code{occasions} slot is derived directly from the dimensions of the
  \code{caphist} arrays and is always set correctly, regardless of whether
  the input frames carry this slot. This ensures compatibility with
  \code{do.trim}, which is called internally by \code{oSCR.fit} whenever
  \code{trimS} is non-\code{NULL}.

  The \code{mmdm} slot is concatenated from the input frames if present. Note
  that \code{mmdm} as computed by \code{\link{make.scrFrame}} is a scalar
  pooled across all sessions; if a more accurate pooled \code{mmdm} is
  required for the combined frame it should be recomputed directly from the
  combined capture histories and trap locations.

  Telemetry objects are combined by appending the per-session elements of
  \code{fixfreq}, \code{indCovs}, and \code{cap.tel} across all input frames.
}
\author{
  Andy Royle
}
\seealso{
  \code{\link{subFrame}}, \code{\link{make.scrFrame}}, \code{\link{oSCR.fit}}
}
\examples{
## Round-trip: split a multi-session scrFrame then recombine
data(rbs)
sf1 <- subFrame(rbs$scrFrame, 1)
sf2 <- subFrame(rbs$scrFrame, 2)
sf.combined <- combineFrames(list(sf1, sf2))

## Fit a model to the combined frame
ss.combined <- list(rbs$ssDF[[1]], rbs$ssDF[[2]])
out <- oSCR.fit(model = list(D ~ 1, p0 ~ 1, sig ~ 1),
                scrFrame = sf.combined,
                ssDF = ss.combined,
                trimS = 2)
}
\keyword{spatial capture-recapture}
\keyword{scrFrame}
