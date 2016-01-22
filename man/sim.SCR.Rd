\name{sim.SCR}
\alias{sim.SCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate SCR data
}
\description{
This function will simulate single session SCR encounter history data.
}
\usage{
sim.SCR(N = 100, K = 20, alpha0 = -2.5, sigma = 0.5, discard0 = TRUE, array3d = FALSE, ssRes = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Total number of individuals in the study area.}
  \item{K}{
Number of occasions.
}
  \item{alpha0}{
The baseline encounter probability.
}
  \item{sigma}{
The parameters describing the spatial scale of detection.
}
  \item{discard0}{
If \code{TRUE}, indivuduals with "all-zero" encounter histtories will be removed. If \code{FALSE}, the encounter histroy data will contain th e"all-zero" encounter histories and therefore have \code{N} rows.
}
  \item{array3d}{
The encounter history data can be returned as either an individual-bytrap-by-occasion 3-D array (\code{array3d=TRUE}), or a "flat" individual-by-trap encounter frequency matrix.

}
  \item{ssRes}{
The resolution of the state-space (pixel/grid size). This value is the length of the side of a square pixel.
}
}
\details{
This function is called by the more flexible \code{sim.SCR.ms} function that returns object formatted for direct model fitting.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chris Sutherland \email{csutherland@umass.edu}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
scr.data <- sim.SCR(N = 100, K = 20, alpha0 = -2.5, sigma = 0.5, discard0 = TRUE, array3d = FALSE, ssRes = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
