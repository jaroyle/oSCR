\name{make.trapCovs}
\alias{make.trapCovs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.trapCovs(lol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lol}{
%%     ~~Describe \code{lol} here~~
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
function (lol) 
{
    ncovs <- length(lol)
    nsessions <- length(lol[[1]])
    cnames <- names(lol)
    outlist <- list(rep(NULL, nsessions))
    for (s in 1:nsessions) {
        ntrps <- nrow(lol[[1]][[s]])
        tmp <- matrix(NA, nrow = ntrps, ncol = ncovs)
        for (c in 1:ncovs) {
            tmp[, c] <- lol[[c]][[s]][, 1]
        }
        tmp <- data.frame(tmp)
        colnames(tmp) <- cnames
        outlist[[s]] <- tmp
    }
    if (ncovs > 1) {
        for (s in 1:nsessions) {
            ntrps <- nrow(lol[[1]][[s]])
            nocc <- ncol(lol[[1]][[s]])
            for (k in 1:nocc) {
                tmp <- matrix(NA, nrow = ntrps, ncol = ncovs)
                for (c in 1:ncovs) {
                  tmp[, c] <- lol[[c]][[s]][, k]
                }
                tmp <- data.frame(tmp)
                colnames(tmp) <- cnames
                outlist[[s]][[k]] <- tmp
            }
        }
    }
    outlist
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
