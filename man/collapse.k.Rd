\name{collapse.k}
\alias{collapse.k}
\title{
Reduce capture histories to a single (k=1) count.
}
\description{
Reduce capture histories to a single (k=1) count.
}
\usage{
collapse.k(scrFrame)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
An scrFrame.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{scrFRame}{An edited scrFrame with counts and k=1.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chris Sutherland
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
reduced.rbs <- collapse.k(scrFrame = scrFrame)
rbs.ss <- make.ssDF(reduced.rbs,3,1)
fit <- oSCR.fit(list(~1,~offset(log(effort)),~1), reduced.rbs, rbs.ss, encmod = "P")
}