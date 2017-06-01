\name{data2oscr.new}
\alias{data2oscr.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
data2oscr.new(edf, sess.col = 1, id.col = 2, occ.col = 3, trap.col = 4, sex.col = NULL, tdf = NULL, K = NULL, ntraps = NULL, remove.zeros = FALSE, remove.extracaps = FALSE, sex.nacode = NULL, tdf.sep = "/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edf}{
%%     ~~Describe \code{edf} here~~
}
  \item{sess.col}{
%%     ~~Describe \code{sess.col} here~~
}
  \item{id.col}{
%%     ~~Describe \code{id.col} here~~
}
  \item{occ.col}{
%%     ~~Describe \code{occ.col} here~~
}
  \item{trap.col}{
%%     ~~Describe \code{trap.col} here~~
}
  \item{sex.col}{
%%     ~~Describe \code{sex.col} here~~
}
  \item{tdf}{
%%     ~~Describe \code{tdf} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{ntraps}{
%%     ~~Describe \code{ntraps} here~~
}
  \item{remove.zeros}{
%%     ~~Describe \code{remove.zeros} here~~
}
  \item{remove.extracaps}{
%%     ~~Describe \code{remove.extracaps} here~~
}
  \item{sex.nacode}{
%%     ~~Describe \code{sex.nacode} here~~
}
  \item{tdf.sep}{
%%     ~~Describe \code{tdf.sep} here~~
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

library(oSCR)

# could get EDF/TDF file from data(ocelot) but here we read in raw data files
###data(ocelot)

# ocelot data have 2 sessions
edf<- read.csv("ocelot_edf.csv")  ##[1:41,]
tdf1<- read.csv("ocelot_tdf1.csv")
tdf2<- read.csv("ocelot_tdf2.csv")

# convert basic EDF/TDF data into "3-d array" and list formats
data<-data2oscr(edf[1:41,],sex.col=5, K=c(44),   ntraps=c(23), 
             tdf=list(tdf1),
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 

# need to have occasions labeled same as in TDF. THIS IS DIFFERENT than old version
colnames(tdf1) <- c("Detector","X","Y", 1:44)
colnames(tdf2) <- c("Detector","X","Y", 1:96)

data.new<- data2oscr.new(edf[1:41,],sex.col=5, K=c(44),   ntraps=c(23), 
             tdf=list(tdf1),
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
