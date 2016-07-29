\name{SCRsmy}
\alias{SCRsmy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summarize SCR data %%  ~~function to do ... ~~
}
\description{
Products a simple summary of SCR data including
number of recaptures, MMDM, etc.. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCRsmy(y3d, traplocs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y3d}{  nind x ntraps x noccasions encounter array (binary)
}
  \item{traplocs}{
  ntraps x 2 matrix of trap coorindates 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
print out containing the following things:

1. Total capture events

2. Spatial recaptures. note: this is actually the total number of captures 
for individuals captured at more than 1 trap. Possibly not the most sensible
definition of spatial recaptures.

3. Ordinary capture events: number of captures if the data were reduced
to 2-d encounter histories

4. Captures lost in non-spatial model: difference between 1 and 3

5. Mean traps per individual: avg. number of unique traps per individual

6. MMDM -- mean maximum distance moved


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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
SCRsmy(beardata$bearArray, traps)

Total captures:  151
Spatial recaptures:  132
Ordinary capture events:  95
Captures lost in non-spatial model:  56
Mean traps/individual:  2.553191
MMDM:  4.647112
   (based on  24  individuals)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
