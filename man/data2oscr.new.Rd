\name{data2oscr.new}
\alias{data2oscr.new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert standard SCR data file(s) to oSCR data objects
}
\description{
This will take input files similar to those used
for the R package 'secr' and convert them to list objects
used by oSCR. 
}
\usage{
data2oscr(edf, sess.col = 1, id.col = 2, occ.col = 3, trap.col = 4, sex.col = NULL, tdf = NULL, K = NULL, ntraps = NULL, remove.zeros = FALSE, remove.extracaps = FALSE, sex.nacode = NULL, tdf.sep = "/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edf}{
NOT A LIST.
Encounter data file having basic information about session, individual, occasion,
trap and (optional) sex.
Note: if the session column indicates multiple sessions then the result
will be to produce a list of appropriate length. e.g., if session has
values (1,2) then the output lists will all be of length 2. If you have
session indicating "sex" you probably don't want this to happen. 
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
a<-data2oscr(rbind(edf1,edf2),sex.col=5,   ntraps=c(23,23), 
             tdf=list(tdf1),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 

# This works. Always best to input K and ntraps!
a<-data2oscr(rbind(edf1,edf2),sex.col=5, K=c(44,96),   ntraps=c(23,23), 
             tdf=list(tdf1,tdf2),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 


# process the trap covariate info. This should be done internally but isn't yet.
dist1=rep(0,23)
dist2=as.vector(a$trapcovs[[2]])/1000
dist1<- matrix(dist1,nrow=23,ncol=44,byrow=FALSE)
dist2<- matrix(dist2, nrow=23,ncol=96, byrow=FALSE)
# oSCR wont work if there's only one covariate...so I pad this here with a dummy covariate
# "list of lists" format
tc<-  list(dist=list(dist1, dist2) , dist2=list(dist1,dist2) )
tc<- make.trapCovs(tc)


scrFrame  <- list(caphist=a$y3d, indCovs=NULL, 
                   type="scr", traps=a$traplocs,trapCovs=tc,
                   trapOperation=a$trapopp )
class(scrFrame) <- "scrFrame"

# make a state-space
ssDF <- make.ssDF(scrFrame, buffer=2000, res = 600)

##  slowish 
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame, ssDF=ssDF,plotit=FALSE , trimS=2000)

## No trim.... really slow!
out1b <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1), ssDF=ssDF,plotit=FALSE )

## Good starts, improves run-time (can make a huge difference in complex models)
out1c <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1), start.vals=c(-3.4, 6.4, -0.58),
                 ssDF=ssDF,plotit=FALSE, trimS=2000 )

# Fit a model with "distance to disturbance"
out1d <- oSCR.fit(model=list(D~1,p0~dist,sig~1), scrFrame,start.vals=c(-3.5, 6.4, 0, -1.5), ssDF=ssDF, plotit=FALSE, trimS=2000  )

# Session specific density
out1e <- oSCR.fit(scrFrame,model=list(D~session,p0~1,sig~1),start.vals=c(-3.5, 6.4, -1.5, -1.5), ssDF=ssDF,plotit=FALSE , trimS=2000)



## Make a density map
pred<- predict.oSCR(scrFrame,out1e,ssDF)

# Plot total density map
par(mfrow=c(1,2))
plot(pred$r[[1]])
plot(pred$r[[2]])

# Plot posterior of s for a specific individual
par(mfrow=c(1,1))
rs<- rasterFromXYZ(cbind(pred$ssDF[[1]][,1:2], pred$preds[[1]][1,]))
plot(rs)
#  library(scrbook)
# spatial.plot(pred$ssDF[[1]][,1:2],pred$preds[[1]][8,],cx=3)






 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
