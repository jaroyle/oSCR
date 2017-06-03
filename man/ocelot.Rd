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
 data(ocelot)
> ls()
[1] "edf1"          "edf2"          "seismic_lines" "tdf1"        "tdf2"         
}
\details{
Data from a study of ocelots using a grid of 23 camera traps. The study
occurred twice at the same location and same trap locations, once with
44 nights of operation and then once again with 96 nights of operation
which occurred after some disturbance regime.  The two periods are 
naturally regarded as sessions. 

The data consist of 2 trap deployment files (tdf1 and tdf2) and 2 
encounter data files (edf1 and edf2).   
The "seismic_lines" data file shows where geological testing lines were
established after the first session. A key question in this study is to
assess the effect of this disturbance on ocelot density or space usage. 
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
#library(rgeos)
#ldist <- apply(gDistance(SpatialPoints(as.matrix(pts)),
#                         seismic_lines,byid=T),2,min)





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
