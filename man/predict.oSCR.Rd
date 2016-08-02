\name{predict.oSCR}
\alias{predict.oSCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Produces data objects necessary to create a map of realized density
  and posterior distributions of activity centers. 
}
\description{
Computes the estimated posterior distribution for each activity centers
and the total realized density.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.oSCR(scrFrame, scr.fit, ssDF, costDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
the scrFrame used in fitting the model for which predictions are required%%     ~~Describe \code{scrFrame} here~~
}
  \item{scr.fit}{
the oSCR fit object for which predictions are required %%     ~~Describe \code{scr.fit} here~~
}
  \item{ssDF}{
the state-space object used for obtaining the model fit for which
predictions are required %%     ~~Describe \code{ssDF} here~~
}
  \item{costDF}{
the cost surface dataframe for obtaining the model fit for which
predictions are required %%     ~~Describe \code{costDF} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{r}{ raster format of density predictions (sum of posterior
    probabilities
  for all activity centers including n+1 ) }
  \item{ssN}{ E(N) for each state-space pixel.  }
  \item{preds}{ Estimated posterior distribution
    of individual activity center locations }
  \item{ssDF}{ the ssDF used}
  \item{pbar}{probability of capture for each pixel }
      
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


###

# Load the beardata (also in oSCR now)
library("oSCR")
data(beardata)
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]
traplocs<- as.matrix(beardata$trapmat)
par(mfrow=c(1,1))
temp <- spiderplot(beardata$bearArray,beardata$trapmat)

 
# Fit a basic model
# Obtain posterior distributions of each individual
# Obtain a total density map
 

tdf<- cbind(1:38, beardata$trapmat, matrix(1, nrow=38, ncol=8))
edf<- beardata$flat
edf[,1]<- 1  # change all sessions to 1 , *important* !
data<- data2oscr(edf, sess.col = 1, id.col = 2, occ.col = 3, trap.col = 4, sex.col = NULL,
 tdf = list(tdf), K = c(8), ntraps = c(38), remove.zeros = TRUE, remove.extracaps = TRUE, sex.nacode = NULL, tdf.sep = "/")
# Make the SCRframe
 scrFrame  <- make.scrFrame(caphist=data$y3d,  traps=data$traplocs,trapCovs=NULL ,
                   trapOperation=data$trapopp )
sf <- data$scrFrame

# make a state-space. Units are strange, maybe "5 km"
ssDF <- make.ssDF(sf, buffer=3, res = 0.5)
dev.off()
plot(ssDF)
points(sf$traps[[1]],pch=3, lwd=2,col="red")

# Basic model SCR0
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame=sf, ssDF=ssDF,plotit=FALSE , trimS=8)

# Note: density intercept is log(per pixel) so total N = exp(intercept)*[# statespace points]
Nhat<- exp(out1$outStats[3,3])*nrow(ssDF[[1]])

# Obtain the predictions of density and activity center posteriors
pred<- predict.oSCR(sf, out1 ,ssDF)

#### predict.oSCR is just a wrapper. If you want you can get the raw materials here:
#### a<- oSCR.fit(scrFrame=scrFrame,ssDF=ssDF,start.vals=out1$rawOutput$estimate,predict=TRUE)
####

  

# Plot posterior of s for a specific individual
par(mfrow=c(1,1))
for(i in 1:47){
 rs<- rasterFromXYZ(cbind(pred$ssDF[[1]][,1:2], pred$preds[[1]][i,]))
 plot(rs)
 points(tdf[,2:3],pch=20)
 # Add the traps of captur
 caps<- data$y3d[[1]][i,,]
 caps<- apply(caps, 1, sum) # total captures per trap
 points(tdf[caps>0, 2:3], pch=as.character(caps[caps>0]), col="red",lwd=2)
 browser()
}
 
# posterior of s for an uncaptured individual
rs<- rasterFromXYZ(cbind(pred$ssDF[[1]][,1:2], pred$preds[[1]][48,]))
plot(rs)
points(tdf[,2:3],pch=20)
 
# Number of uncaptured individuals
n0<- sum(pred$preds[[1]][48,])
 
(Nhat<-  sum(pred$ssN[[1]],na.rm=TRUE)) 
# [1] 57.49587
  
# Plot total density map. Same as ssN except in "raster" format (with some missing values)
par(mfrow=c(1,1))
plot(pred$r[[1]])   # pred is a list, one for each session, $r is a raster 
points(tdf[,2:3],pch=20)
title("Realized density")
null <- spiderplot(beardata$bearArray,beardata$trapmat, add=TRUE)

# Conditional probability of capture (ESA is area-weighted sum of this)
par(mar=c(3,3,3,6))
 plot(pred$pbar[[1]])




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
