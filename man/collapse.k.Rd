\name{collapse.k}
\alias{collapse.k}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
collapse.k(scrFrame)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
%%     ~~Describe \code{scrFrame} here~~
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
data("rbs_ecography")
rbs.data <- data2oscr(edf = rbs.edf, #the EDF
                      sess.col = 1, #session column
                      id.col = 2, #individual column
                      occ.col = 3, #occasion column
                      trap.col = 4, #trap column
                      sex.col = 5, #sex column (optional)
                      tdf = list(rbs.tdf1, rbs.tdf2, #list of TDF
                                 rbs.tdf3, rbs.tdf4),
                      K = c(7,5,6,4), #occasion vector
                      ntraps = c(50,50,50,50), #no. traps vector
                      trapcov.names = c("jday","jday2"), #covariate names
                      tdf.sep = "/", #char used separator
                      sex.nacode = "U") #unknown sex code
scrFrame <- rbs.data$scrFrame
rbs2 <- collapse.k(scrFrame = scrFrame)
rbs.ss <- make.ssDF(rbs2,3,1)
f3 <- oSCR.fit(list(~1,~offset(log(effort)),~1), rbs2, rbs.ss, encmod = "P")

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
