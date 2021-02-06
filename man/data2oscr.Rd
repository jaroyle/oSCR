\name{data2oscr}
\alias{data2oscr}
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
data2oscr(edf, sess.col = 1, id.col = 2, occ.col = 3, trap.col = 4, sex.col = NULL, tdf = NULL, K = NULL, ntraps = NULL, remove.zeros = TRUE, trapcov.names=NULL,remove.extracaps = TRUE, sex.nacode = NULL, tdf.sep = "/")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edf}{
MATRIX -- NOT A LIST.
Encounter data file having basic information about session, individual, occasion,
trap and (optional) sex.
Note: if the session column indicates multiple sessions then the result
will be to produce a list of appropriate length. e.g., if session has
values (1,2) then the output lists will all be of length 2. If you have
session indicating "sex" you probably don't want this to happen so add a
new column for "sex" and keep that information independent of "session".
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

 IMPORTANT NOTE: The column names for the trap deployment columns (cols
 >= 4) must have the same NAMES as the traps referenced in the encounter
 data file (EDF).
}
  \item{K}{
Vector number of occasions per session. K = c(8) for a single season
 study with 8 occasions. 
}
  \item{ntraps}{
Vector number of traps per session. ntraps = c(38) for a single season
 study with 38 traps. 
}
  \item{remove.zeros}{
If TRUE then individuals that were never captured in a session are removed
from the session encounter history data. If FALSE then the resulting
 data structure preserves individual ID across sesions.  For example in
 a 2 session data set then output$y3d[[1]][7,,]  is the same individual
 as output$y3d[[2]][7,,].
For fitting models in oSCR you should have remove.zeros=TRUE but if  you
 want to fit open SCR models in JAGS you might use the remove.zeros =
 FALSE option so that the rows are alighned properly. 

}
\item{trapcov.names}{

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

The EDF should be a single EDF with sessions differentiated by a column
in the EDF.

The TDF should be provided in list form. One element per session.


}
\value{
Returns a list with elements
 \item{edf }{the encounter data file}
 \item{y3d }{a list of 3-d encounter history matrices, one for each session}
 \item{sex }{a list of sex information suitable for inputting into oSCR.fit}
 \item{traplocs }{a list of trap locations suitable for inputting into oSCR.fit}
 \item{trapopp }{a list of trap operation information suitable for inputting into oSCR.fit}
 \item{trapcovs }{a list of trap covariates suitable for inputting into
   oSCR.fit}
 \item{scrFrame}{ the scrFrame suitable for use in oSCR.fit. }

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

# Use data2oscr to process the black bear data (no covariates, no trap operation information)
library(oSCR)
data(beardata)
data<- data2oscr(beardata$edf,
                 sess.col = 1,
                 id.col = 2,
                 occ.col = 3,
                 trap.col = 4,
                 tdf = list(beardata$tdf[,1:3]),
                 K = c(8),
                 ntraps = c(38))
sf <- data$scrFrame
ssDF <- make.ssDF(sf, buffer=3, res = 0.5)
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame=sf, ssDF=ssDF,plotit=FALSE , trimS=8)
####pred<- predict.oSCR(out1)


# Approach 2
data <- data2oscr( edf = beardata$edf,
tdf = list(beardata$tdf),
sess.col = which(colnames(beardata$edf) %in% "Session"),
id.col = which(colnames(beardata$edf) %in% "ID"),
occ.col = which(colnames(beardata$edf) %in% "Occasion"),
trap.col = which(colnames(beardata$edf) %in% "trapID"),
sex.col = which(colnames(beardata$edf) %in% "Sex"),
sex.nacode = "NA",
K = 8,
ntraps = nrow(beardata$tdf))

# Ocelot data

data(ocelot)
 # need to have occasions labeled same as in TDF
colnames(tdf1) <- c("Detector","X","Y", 1:44)
colnames(tdf2) <- c("Detector","X","Y", 1:96)

# This fails because one TDF is given and the data imply DIFFERENT # occasions
data<-data2oscr(rbind(edf1,edf2),sess.col=1, id.col=2, occ.col=3, trap.col=4, sex.col=5,   
             ntraps=c(23,23),
             tdf=list(tdf1),
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") )

# This works. Always best to input K and ntraps!
data<-data2oscr(rbind(edf1,edf2),sess.col=1, id.col=2, occ.col=3, trap.col=4,sex.col=5, K=c(44,96),   ntraps=c(23,23),
             tdf=list(tdf1,tdf2),
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") )

# There is one covariate "distance to disturbance" which is not defined for session = 1
#    therefore I set it to 0 here. This is probably ok if we define the "exposure" covariate to be 1/dist
#    because then as "dist -> 0" (i.e., close to disturance) then exposure -> BIG whereas if dist -> infinity then
#    exposure -> 0.  Therefore 0 values are like "no exposure to disturbance" which is a good description of
#    pre-exposure.


# The covariate "distance" is in the input TDF file as a trap covariate. The function data2oscr processes that into a matrix
#  but it does NOT exist for session 1!
dist1 =  rep(0,23)
dist2 =  as.vector(data$trapcovs[[2]])/1000   # Convert to km
# Organize the covariates into a matrix that is J x K (traps x occasions)
dist1 <- matrix(dist1, nrow=23,ncol=44,byrow=FALSE)
dist2 <- matrix(dist2, nrow=23,ncol=96, byrow=FALSE)
# Transformation to exposure to disturbance. Scaling is pretty bad, some extreme values
exposure1<- dist1    # This is 0 for all traps
exposure2<- 1/dist2  # This goes to 0 as distance increase

# Package these into a suitable trapCovs format. Example of 2 covariates, defined for each session here:
tc<-  list(expose=list(exposure1, exposure2), expose2 = list(exposure1^2, exposure2^2) )

tc<-  make.trapCovs(tc)
# Example of 1 covariate
# "list of lists" format
tc <-  list(expose=list(exposure1, exposure2))
tc <-  make.trapCovs(tc)
# example of 2 covariates
#tc<-  list(dist=list(dist1, dist2) , dist2=list(dist1,dist2) )
#tc<- make.trapCovs(tc)

# Make the SCRframe. Sex is not included here. Either way is fine....
ocelot.sf  <-  make.scrFrame(caphist=data$y3d, indCovs=NULL,
                    traps=data$traplocs,trapCovs=tc ,
                    trapOperation=data$trapopp )


# Make a state-space
ocelot.ss <- make.ssDF(ocelot.sf, buffer=2000, res = 600)
plot(ocelot.ss[[1]], pch=20)
points(ocelot.sf$traps[[1]],pch=3,col="red",lwd=2)
spiderplot(y = data$y3d[[2]],traplocs = data$traplocs[[1]] )
spiderplot(y = data$y3d[[2]],traplocs = data$traplocs[[1]],add=TRUE)

##  fit a model
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), ocelot.sf, ssDF=ocelot.ss,plotit=FALSE , trimS=2000)

## Good starts, improves run-time (can make a huge difference in complex models)
out1c <- oSCR.fit(model=list(D~1,p0~1,sig~1), ocelot.sf, start.vals=c(-3.4, 6.4, -0.58),
                 ssDF=ocelot.ss,plotit=FALSE, trimS=2000 )

# Fit a model with "exposure to disturbance"
## default vals for some models are especially bad. Could be geographic scaling issue.
out5 <- oSCR.fit(model=list(D~1,p0~ expose ,sig~1), ocelot.sf, start.vals=c(-3.5, 6.4, 0.035, -1.5),
                 ssDF=ocelot.ss, plotit=FALSE, trimS=2000  )


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
