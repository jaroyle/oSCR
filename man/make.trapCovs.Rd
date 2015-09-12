\name{make.trapCovs}
\alias{make.trapCovs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format trap x occasion covariates for oSCR. 
  %%  ~~function to do ... ~~
}
\description{
This takes a simple  list-of-lists (lol) and restructures it into the
correct format for oSCR (as a list-of-list-of-data-frames).  
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.trapCovs(lol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lol}{
     A list of lists. e.g., lol$cov1 is a list of length nSessions where
     lol$cov[[s]] contains the nTraps x nOcc values of "cov1" for
     session s.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
