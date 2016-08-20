\name{make.scrFrame}
\alias{make.scrFrame}
\title{
Make an scrFrame object.
}
\description{
A helper function for making the scrFrame object required for model fitting in \pkg{oSCR}
}
\usage{
make.scrFrame(caphist, traps, indCovs = NULL, trapCovs = NULL, 
              trapOperation = NULL, type = "scr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{caphist}{
A list of capture history arrays with dimensions \eqn{n x J x K} (individual-by-trap-by-occasion), one for each session. The arrays are either binary or counts.}
  \item{traps}{
A list of data frames, one for each session, containing the coordinates of all the traps/detectors. The coordinates must be named X and Y (note these are capitals). Traps that are not operational in specific occasions can be \emph{turned off} by supplying the appropriate \code{trapOperation} object.}
  \item{indCovs}{
A list of data frames, one for each session, containing individual-specific data. Currently \code{oSCR} deals only with  a two-class category such as sex, 
but can accommodate missing data (e.g., male=1, female=0, and unknown=NA).}
  \item{trapCovs}{
A two level list. A list for each session contains a further list of data frames, one for each occasion, i.e., trap covariates can be time-varying.
}
  \item{trapOperation}{
 A list of binary matrices with dimension \eqn{J x K} (trap-by-occasion), one for each session. The matrix denotes whether a trap was operational (1) or not (0), in each occasion. Traps with trapOperation equal to zero does not contribute to the likelihood.}
  \item{type}{
This should be set to SCR always. This option will, in the future, allow the user to provide data in the \pkg{secr} format.}
}
\details{
}
\value{
Returns an object of class scrFrame. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chris Sutherland}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(rbs)
caphist <- rbs.sf$caphist
traps <- rbs.sf$traps
trapCovs <- rbs.sf$trapCovs

sf <- make.scrFrame(caphist = caphist, traps = traps, trapCovs = trapCovs, type = "scr")
sf
}