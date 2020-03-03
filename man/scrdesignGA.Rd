\name{scrdesignGA}
\alias{scrdesignGA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for generating optimal SCR sampling designs.
}
\description{
Given any SCR-sensible objective function, \code{scrdesignGA()} uses implements a genetic algorithm to solve a \emph{k-of-n} problem. In the context of SCR, and for a given objective function and landscpae, the function finds the best subset of sampling locations (\emph{k})from all possible sampling locations (\emph{n}).    
}
\usage{
scrdesignGA(statespace = NULL, alltraps = NULL, ntraps = 9, beta0 = -0.6, sigma = 2, crit = 1, N = 100, verbose = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statespace}{
%%     ~~Describe \code{statespace} here~~
}
  \item{alltraps}{
%%     ~~Describe \code{alltraps} here~~
}
  \item{ntraps}{
%%     ~~Describe \code{ntraps} here~~
}
  \item{beta0}{
%%     ~~Describe \code{beta0} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{crit}{
%%     ~~Describe \code{crit} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
library(kofnGA)
data(pakistan)
plot(pakSS, asp=1, pch=16, cex=0.25, col="grey")
points(pakTT, pch=3, cex=1.2)

test_design <- scrdesignGA(statespace = pakSS, alltraps = pakTT,
                           ntraps=25, beta0 = 0.2, sigma=3, crit=2,
                           popsize = 10, ngen=50, keepbest=5)
par(mfrow=c(1,3))
plot(test_design, which=4)
}

\keyword{ ~methods }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~design }% __ONLY ONE__ keyword per line
