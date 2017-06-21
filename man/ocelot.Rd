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
[1] "edf1"          "edf2"          "seismic_lines" "tdf1"        "tdf2"       "disturb"  
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

The object "disturb" has an index of disturbance computed for each 600 m x 600 m pixel in
the study area (state-space).
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
# Interested in using "disturb" as a covariate on density. There are values
# of disturbance computed for different distances from the center of the pixel...
library(raster)
plot(rasterFromXYZ(disturb[,c("X","Y","g1000")]))
 
data<- data2oscr(rbind(edf1,edf2),sex.col=5, K=c(44, 96),   ntraps=c(23,23), 
             tdf=list(tdf1,tdf2  ),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 
 
ssDF2<- disturb[,c("X","Y","g600", "g1000","g1200")]
ssDF1<- ssDF2
ssDF1[,3:5]<- 0
ocelot.ss<- list(ssDF1,ssDF2)

# Make the SCRframe. Sex is not included here. Either way is fine....
ocelot.sf  <-  make.scrFrame(caphist=data$y3d, indCovs=data$sex, 
                    traps=data$traplocs,trapCovs=NULL ,
                    trapOperation=data$trapopp )


# Will add some model fitting and summarization functions later

}
\keyword{datasets}
