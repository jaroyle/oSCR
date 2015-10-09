\name{make.S}
\alias{make.S}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make.S%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.S(scrFrame, buffer, res, xy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
%%     ~~Describe \code{scrFrame} here~~
}
  \item{buffer}{
%%     ~~Describe \code{buffer} here~~
}
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{xy}{
%%     ~~Describe \code{xy} here~~
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
function (scrFrame, buffer, res, xy = NULL) 
{
    if (class(scrFrame) != "scrFrame") 
        stop("Need an object of class scrFrame to create the state space")
    if (is.null(xy)) 
        xy <- c(1, 2)
    if (res == NULL) 
        stop("You didnt provide a resolution value!")
    if (res == NULL) 
        stop("You didnt provide a buffer value!")
    S <- list()
    trpls <- scrFrame$traps
    for (i in 1:length(trpls)) {
        bl <- apply(trpls[[i]][, xy], 2, min)
        tr <- apply(trpls[[i]][, xy], 2, max)
        sxy <- expand.grid(seq(bl[1] - buffer, tr[1] + buffer, 
            res), seq(bl[2] - buffer, tr[2] + buffer, res))
        dd <- apply(e2dist(sxy, trpls[[i]][, xy]), 1, min)
        S[[i]] <- sxy[dd <= buffer, ]
    }
    return(S)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
