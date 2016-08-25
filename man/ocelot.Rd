\name{ocelot}
\alias{ocelot}
\docType{data}
\title{
Ocelot camera trapping data.
}
\description{
Joe Kolowski's ocelot data
}
\usage{data("ocelot")}
\format{
  The format is:
 chr "ocelot"
}
\details{
Data from a study of ocelots using a grid of 23 camera traps. The study
occurred twice at the same location and same trap locations, once with
44 nights of operation and then once again with 96 nights of operation
which occurred after some disturbance regime.  The two periods are 
naturally regarded as sessions. 

The data consist of 2 trap deployment files (tdf1 and tdf2) and 2 
encounter data files (edf1 and edf2).   These files should run
in 'secr' directly (examples to be provided....not done yet....)


}
\source{
Kolowski, J.M, and A. Alonso. 2010.
Density and activity patterns of ocelots (Leopardus pardalis) in northern Peru and the impact
of oil exploration activities. Biological Conservation 143:917-925. 
}
\references{
Kolowski, J.M, and A. Alonso. 2010.
Density and activity patterns of ocelots (Leopardus pardalis) in northern Peru and the impact
of oil exploration activities. Biological Conservation 143:917-925. 
}
\examples{

 


library(oSCR)
data(ocelot)

# Convert the data format
data<-data2oscr(rbind(edf1,edf2),sex.col=5, K=c(44,96),   ntraps=c(23,23), 
             tdf=list(tdf1,tdf2),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 


# Process the trap covariate info. This should be done internally but isn't yet.
# There is one covariate "distance to disturbance" which is not defined for session = 1
# therefore I set it to 0 here so that the linear predictor evaluates to 0 (i.e., no adjustment to p).

dist1 =  rep(0,23)
dist2 =  as.vector(data$trapcovs[[2]])/1000   # Convert to km
dist1 <- matrix(dist1, nrow=23,ncol=44,byrow=FALSE)
dist2 <- matrix(dist2, nrow=23,ncol=96, byrow=FALSE)
 
# Example of 1 covariate
tc <-  list(dist=list(dist1, dist2) )
tc <-  make.trapCovs(tc)

# Make the SCRframe
scrFrame  <- list(caphist=data$y3d, indCovs=NULL, 
                   type="scr", traps=data$traplocs,trapCovs=tc ,
                   trapOperation=data$trapopp )
class(scrFrame) <- "scrFrame"

# make a state-space
ssDF <- make.ssDF(scrFrame, buffer=2000, res = 600)

# Make a descriptive plot of everything
plot(ssDF[[1]], pch=20)
points(data$traplocs[[1]],pch=3,col="red",lwd=2)
spiderplot(data$y3d[[2]],data$traplocs[[1]],add=TRUE)

# Fit a basic model
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame, ssDF=ssDF,plotit=FALSE , trimS=2000)

## default vals for some models are especially bad. Could be geographic scaling issue.
out2 <- oSCR.fit(model=list(D~1,p0~dist,sig~1), scrFrame,start.vals=c(-3.5, 6.4, 0, -1.5), ssDF=ssDF, plotit=FALSE, trimS=2000  )

out3 <- oSCR.fit(scrFrame,model=list(D~session,p0~1,sig~1),start.vals=c(-3.5, 6.4, -1.5, -1.5), ssDF=ssDF,plotit=FALSE , trimS=2000)

## Summarize output
## Make a density map

pred<- predict.oSCR(scrFrame,out3,ssDF)

# Plot total density map
par(mfrow=c(1,2))
plot(pred$r[[1]])
plot(pred$r[[2]])

# Plot posterior of s for a specific individual
par(mfrow=c(1,1))
rs<- rasterFromXYZ(cbind(pred$ssDF[[1]][,1:2], pred$preds[[1]][1,]))
plot(rs)

# Compute distances from each point to seismic lines
# points should be traps or state-space grid points. 
library(rgeos)
ldist <- apply(gDistance(SpatialPoints(as.matrix(pts)),
                         seismic_lines,byid=T),2,min)





##
## 
## plot the cost-surface variable "disturb"
##
##

# Do this before fitting ASU models because all values should be 0 pre-disturb
  ocelot.cost[[1]][,"disturb"]<- 0

  plot(ocelot.cost[[2]][,c("X","Y")],pch=".")
  points(ocelot.cost[[2]][,c("X","Y")][ocelot.cost[[2]][,"disturb"]==1,],pch=20)
  ocelot.cost[[1]][,"disturb"]<- 0
  table(ocelot.cost[[1]][,"disturb"])
# plot(seismic_lines, add=TRUE)
  points(tdf1[,c("X","Y")],pch=3, lwd=3,col="red")
 



}
\keyword{datasets}
