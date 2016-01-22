\name{sim.SCR.ms}
\alias{sim.SCR.ms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Multi-Session SCR data
}
\description{
This function will simulate multi-session SCR encounter history data and also allows for simulation of individual sex data.}
}
\usage{
sim.SCR.ms(sessions = 1, sex = F, sex.ratio = 0.5, N = 100, K = 5, alpha0 = -2.5, sigma = 0.5, discard0 = TRUE, array3d = FALSE, ssRes = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sessions}{
Number of session to simulate dats for
}
  \item{sex}{
Whether of not each individual is assigned a sex. The default is no sex data.
}
  \item{sex.ratio}{
If simulating sex data, you can specify the sex ration. This is either a single value, in which case all session have the same sex ratio, or it can be a vector containing session specific sex ratios.
}
  \item{N}{
Abundance. This is either a single value, in which case all sessions have the same number of individuals, or it can be a vector containing session specific abundances.
}
  \item{K}{
Number of sampling occasion. This is either a single value, in which case all session have the same number of visits, or it can be a vector containing session specific number of occasions.}
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
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
simdat <- sim.SCR.ms(sessions = 2, sex = F, sex.ratio = 0.5, N = 100, K = 5, alpha0 = -2.5, sigma = 0.5, discard0 = TRUE, array3d = FALSE, ssRes = 0.5)
sf <- simdat$sf
ssDF <- make.ssDF(sf,buffer=1.5,res=0.25)
out <- oSCR.fit(sf,list(~1,~1,~1,~1),ssDF)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
