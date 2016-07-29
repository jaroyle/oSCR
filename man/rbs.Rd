\name{rbs}
\alias{rbs}
\docType{data}
\title{
Red-backed salamander \emph{Plethodon cinereus} data. 
}
\description{
The rbs data is the data analyzed by Sutherland \emph{et al.}'s (2016) \emph{Herpetologia} paper.
}
\usage{data("rbs")}
\format{
rbs contains a \code{scrFrame object and an \code{ssDF object (see \code{?oSCR.fit})
}
\details{
This is a four-session data set where sessions are independent artificial cover object (ACO) plots in Ithaca, NY, USA.
}
\source{
Chris Sutherland, David J. Muñoz, David A.W. Miller, and Evan H. Campbell Grant (2016) \emph{Spatial Capture–Recapture: A Promising Method for Analyzing Data Collected Using Artificial Cover Objects}. Herpetologica: March 2016, Vol. 72, No. 1, pp. 6-12.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

################################################################################
### Analysis of the Ithaca (NY) salamander data
################################################################################
 
  library(oSCR)
  library(scrbook)
  library(gdata)
  library(raster)
  library(rgdal)
 
#######################
#Get the data
 
  data(rbs)

#######################
#Fit some models
 
  #NULL model
  rbs.0 <- oSCR.fit(rbs.sf, list(D~1, p0~1, a1~1), rbs.ss, trimS=2)
 
  #Variation in density
  rbs.Dsite <- oSCR.fit(rbs.sf, list(D~session, p0~1, a1~1), rbs.ss, trimS=2)
 
  #Variation in detectability by date and date^2
  rbs.Pday <- oSCR.fit(rbs.sf, list(D~session, p0~day, a1~1), rbs.ss, trimS=2)
  rbs.Pday <- oSCR.fit(rbs.sf, list(D~session, p0~day+day2, a1~1), rbs.ss, trimS=2)
 
 
##############################
#Summarize and compare models
 
  (fl <- fitList(list(rbs.0,rbs.Dsite,rbs.Pday,rbs.Pday),rename=T))
  (ms <- modSel.oSCR(fl))
  ms$aic
  ms$coef

}
\keyword{dataset}
