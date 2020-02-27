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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (statespace = NULL, alltraps = NULL, ntraps = 9, beta0 = -0.6, 
    sigma = 2, crit = 1, N = 100, verbose = 1, ...) 
{
    if (is.null(statespace)) 
        stop("Must supply a 'statespace' object (coords of the study area)")
    if (is.null(alltraps)) 
        stop("Must supply a 'alltraps' object (coords of all possible trap locations)")
    des <- kofnGA(n = nrow(alltraps), k = ntraps, OF = scrdesignOF, 
        verbose = verbose, ..., alltraps = alltraps, statespace = statespace, 
        beta0 = beta0, sigma = sigma, crit = crit)
    optimaltraps <- alltraps[des$bestsol, ]
    scrdesign <- list(des = des, statespace = statespace, alltraps = alltraps, 
        optimaltraps = optimaltraps, sigma = sigma, beta0 = beta0)
    class(scrdesign) <- "scrdesign"
    return(scrdesign)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
