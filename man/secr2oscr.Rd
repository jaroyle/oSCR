\name{secr2oscr}
\alias{secr2oscr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert secr data file(s) to oSCR data objects
}
\description{
This will take input files for the R package 'secr' and convert them to list objects
used by oSCR. 
}
\usage{
secr2oscr(edf, sess.col = 1, id.col = 2, occ.col = 3, trap.col = 4, sex.col = NULL, tdf = NULL, K = NULL, ntraps = NULL, remove.zeros = FALSE, remove.extracaps = FALSE, sex.nacode = NULL, tdf.sep = "/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edf}{
NOT A LIST.
Encounter data file having basic information about session, individual, occasion,
trap and (optional) sex.   
}
  \item{sess.col}{
Column index (integer) of session info
}
  \item{id.col}{
Column index (integer) of individual ID
}
  \item{occ.col}{
Column index (integer) of occasion 
}
  \item{trap.col}{
Column index (integer) of trap 
}
  \item{sex.col}{
Optional column index (integer) of sex 
}
  \item{tdf}{
LIST: one for each session. 
Each element should be a trap deployment file. Necessary format:
Trap ID, Xcoord, Ycoord, trap deployment columns, character symbol, trap covariates

Optional columns: character symbol is any character symbol which indicates subsequent columns represent covariates. 

}
  \item{K}{
Optional vector number of occasions per session. If not provided this is deduced from the data. It is better to provide this.
}
  \item{ntraps}{
Optional vector number of traps per session. If not provided this is deduced from the data. It is better to provide it. 
}
  \item{remove.zeros}{
If TRUE then individuals that were never captured in a session are removed
from the session encounter history data. Why would you have these?
}
  \item{remove.extracaps}{
If TRUE then individuals captured > 1 time in a trap at an occasion have their
encounter frequency reduced to 1. 
}
  \item{sex.nacode}{
Which character codes represent unknown sex?
}
  \item{tdf.sep}{
Which character indicates that subsequent columns should be understood
as trap covariates?
}
}
\details{

This function takes one or more encounter data files and one or more
trap deployment files and creates list objects suitable for use by oSCR.fit.

It is best to make sure that you have a proper TDF for EACH EDF that is provided.

The EDF should be a single EDF.

The TDF should be provided in list form. One element per session.

}
\value{
Returns a list with elements
 \item{edf }{the encounter data file}
 \item{y3d }{a list of 3-d encounter history matrices, one for each session}
 \item{sex }{a list of sex information suitable for inputting into oSCR.fit}
 \item{traplocs }{a list of trap locations suitable for inputting into oSCR.fit}
 \item{trapopp }{a list of trap operation information suitable for inputting into oSCR.fit}
 \item{trapcovs }{a list of trap covariates suitable for inputting into oSCR.fit}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andy Royle
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
data(ocelot)

# This fails because one TDF is given and the data imply DIFFERENT # occasions
a<-secr2oscr(rbind(edf1,edf2),sex.col=5,   ntraps=c(23,23), 
             tdf=list(tdf1),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 

# This works
a<-secr2oscr(rbind(edf1,edf2),sex.col=5, K=c(44,96),   ntraps=c(23,23), 
             tdf=list(tdf1,tdf2),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 


 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
