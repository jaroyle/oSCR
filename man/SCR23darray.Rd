\name{SCR23darray}
\alias{SCR23darray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert SCR encounter data file (EDF) into a 3-d array useful for
fitting models in BUGS/JAGS.

}
\description{
Convert SCR encounter data file (EDF) into a 3-d array useful for
fitting models in BUGS/JAGS. The output array is 'nind' x 'ntraps' x 'noccassions'.

Note: this was modified on 10/27/2013 , the previous version incorrectly
returned 'nind' x 'noccasions' x 'ntraps' 
}
\usage{
SCR23darray(edf, tdf) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edf}{
Encounter data file (EDF) with assumed order of columns: session,
individual, occasion, trap. All are integers. 
}
  \item{tdf}{
Trap deployment file (TDF) %%     ~~Describe \code{traps} here~~
}
}
\details{


  This was modified on 10/27/2013 , the previous version incorrectly
returned 'nind' x 'noccasions' x 'ntraps' 


  %%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a  3-d encounter array having dimensions 'nind' x 'ntraps' x
'noccasions' suitable for analysis in BUGS/JAGS. 
}
\references{
  SCR book, Section 5.9
  %% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  data(wolverine), wolvSCR0, array2secr
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library("scrbook")
data(wolverine)
traps<-wolverine$wtraps
y3d <-SCR23darray(wolverine$wcaps,wolverine$wtraps)
# produce 'nind' x 'ntraps' 2-d matrix
y2d <- apply(y3d,c(1,2),sum)
 
toad<-wolvSCR0(y3d,traps,nb=1000,ni=2000,buffer=1,M=100,coord.scale=10000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
