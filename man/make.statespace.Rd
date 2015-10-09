\name{make.statespace}
\alias{make.statespace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make.statespace %%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.statespace(ll = NULL, minx = NA, maxx = NA, miny = NA, maxy = NA, nx = 40, ny = NULL, buffer = 0, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ll}{
%%     ~~Describe \code{ll} here~~
}
  \item{minx}{
%%     ~~Describe \code{minx} here~~
}
  \item{maxx}{
%%     ~~Describe \code{maxx} here~~
}
  \item{miny}{
%%     ~~Describe \code{miny} here~~
}
  \item{maxy}{
%%     ~~Describe \code{maxy} here~~
}
  \item{nx}{
%%     ~~Describe \code{nx} here~~
}
  \item{ny}{
%%     ~~Describe \code{ny} here~~
}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
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
function (ll = NULL, minx = NA, maxx = NA, miny = NA, maxy = NA, 
    nx = 40, ny = NULL, buffer = 0, delta) 
{
    if (is.null(ny)) 
        ny <- nx
    if (!is.null(ll)) {
        minx <- min(ll[, 1])
        maxx <- max(ll[, 1])
        miny <- min(ll[, 2])
        maxy <- max(ll[, 2])
        minx <- minx - buffer
        maxx <- maxx + buffer
        miny <- miny - buffer
        maxy <- maxy + buffer
    }
    xgr <- seq(minx + delta/2, maxx - delta/2, delta)
    ygr <- seq(miny + delta/2, maxy - delta/2, delta)
    gr <- as.matrix(expand.grid(xgr, ygr))
    dimnames(gr) <- list(NULL, c("X", "Y"))
    gr
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
