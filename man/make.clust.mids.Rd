\name{make.clust.mids}
\alias{make.clust.mids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute cluster mid points%%  ~~function to do ... ~~
}
\description{
blah blah blah  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.clust.mids(C, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
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
function (C, S) 
{
    mids <- data.frame(cl = unique(C[, 3]), Xcoord = c(tapply(C[, 
        1], C$cluster, mean)), Ycoord = c(tapply(C[, 2], C$cluster, 
        mean)))
    plot(S)
    points(C[, 1:2], pch = 21, bg = C[, 3], cex = 1.5)
    points(mids[, -1], pch = "+", col = mids[, 1], cex = 1.5)
    return(mids)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
