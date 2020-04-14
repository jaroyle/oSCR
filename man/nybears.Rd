\name{nybears}
\alias{nybears}
\docType{data}
\title{
NY black bear data set 
}
\description{
NY black bear data set
}
\usage{data("nybears")}
\format{
  \describe{
    \item{\code{y2d}}{A matrix of capture data for 33 individuals at 103 traps.}
    \item{\code{traplocs}}{A matrix of UTM coordinates for 103 traps.}
    \item{\code{K}}{The number (5) of trapping occasions/surveys used to collect capture data.}
    \item{\code{ss}}{A data frame of UTM coordinates for the 4,329 state-space pixels, each 1 km x 1 km, buffering the trap array by ~10 km.}
    \item{\code{elevation}}{A vector of standardized elevation values for the 4,329 state-space pixels.}
    \item{\code{teldata}}{A data frame of GPS locational data for 3 individual bears, including the animal ID, date, time, and coordinates in both lat/lon and UTM.}
      }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Royle, J.A., Chandler, R.B., Sun, C.C. and Fuller, A.K., 2013. Integrating resource selection information with spatial capture-recapture. Methods in Ecology and Evolution, 4(6), pp.520-530.
}
\examples{
data(nybears)

?telemetry


}
\keyword{datasets}
