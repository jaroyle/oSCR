\name{make.mods}
\alias{make.mods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a table of oSCR models 
  %%  ~~function to do ... ~~
}
\description{
This function produces a table of oSCR models if you provide it the
basic models for each component. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.mods(density, detection, sigma, cost = c(~1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{density}{ density model set
%%     ~~Describe \code{density} here~~
}
  \item{detection}{ detection model set 
%%     ~~Describe \code{detection} here~~
}
  \item{sigma}{ sigma model set 
%%     ~~Describe \code{sigma} here~~
}
  \item{cost}{ cost function model set 
%%     ~~Describe \code{cost} here~~
}
}
\details{
Creates a table containing all possible models obtained by crossing
density, detection, sigma and cost models.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A table with 4 columns.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andy Royle %%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
oSCR.parfit 
}
\examples{


mods <- make.mods(density= c(~1,~session),
                  detection = c(~1,~session),
                  sigma = c(~1,~session),
                  cost = c(~1))
 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
