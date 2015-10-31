\name{subFrame}
\alias{subFrame}
\title{
A function for efficiently subsetting a multi-session "scrFrame" object.}
\description{
This function returns an "scrFrame" object containing only those sessions identified by "subs".
}
\usage{
subFrame(scrFrame, subs)
}
\arguments{
  \item{scrFrame}{
An scrFrame which will be subsetted. Must be of class 'scrFrame'.
}
  \item{subs}{
A numeris vector indicating which of the sessions are to be retained.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An exact copy of the supplied scrFrame for only those session identified using the "subs" argument.
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
function (scrFrame, subs) 
{
    nms <- names(scrFrame)
    nms <- setdiff(nms, c("nocc", "type"))
    miss.nms <- setdiff(c("caphist", "traps", "trapCovs", "indCovs", 
        "trapOperation"), nms)
    scrFrame.sub <- list()
    for (i in nms) {
        scrFrame.sub[[i]] <- list()
        if (is.null(scrFrame[[i]])) {
            scrFrame.sub[[i]] <- NULL
        }
        else {
            for (j in 1:length(subs)) {
                scrFrame.sub[[i]][[j]] <- scrFrame[[i]][[subs[j]]]
            }
        }
    }
    for (i in miss.nms) {
        scrFrame.sub[[i]] <- NULL
    }
    new.nocc <- max(unlist(lapply(scrFrame.sum$caphist, function(x) dim(x)[3])))
    scrFrame.sub$nocc <- scrFrame$nocc
    scrFrame.sub$type <- scrFrame$type
    class(scrFrame.sub) <- "scrFrame"
    return(scrFrame.sub)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
