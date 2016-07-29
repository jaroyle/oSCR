\name{rbs}
\alias{rbs}
\docType{data}
\title{
Red-backed salamander \emph{Plethodon cinereus} data.
}
\description{
The rbs data is the data analyzed by Sutherland \emph{et al.} (2016) in their \emph{Herpetologia} paper.
}
\usage{data("rbs")}
\format{
  The format is:
 chr "rbs"
}
\details{
This is a four-session data set where sessions are independent artificial cover object (ACO) plots in Ithaca, NY, USA. The data object \code{rbs} contains a \code{scrFrame} object and an \code{ssDF} object (see \code{?oSCR.fit})
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
  rbs.0 <- oSCR.fit(model= list(D~1, p0~1, sig~1), 
                       scrFrame= rbs.sf, ssDF = rbs.ss, trimS=2)

  #Variation in density
  rbs.Dsite <- oSCR.fit(model= list(D~session, p0~1, sig~1), 
                       scrFrame= rbs.sf, ssDF = rbs.ss, trimS=2)

  #Variation in detectability by date and date^2
  rbs.Pday <- oSCR.fit(model= list(D~session, p0~day, sig~1), 
                       scrFrame= rbs.sf, ssDF = rbs.ss, trimS=2)
  rbs.Pday <- oSCR.fit(model= list(D~session, p0~day+day2, sig~1), 
                       scrFrame= rbs.sf, ssDF = rbs.ss, trimS=2)


##############################
#Summarize and compare models

  (fl <- fitList(list(rbs.0,rbs.Dsite,rbs.Pday,rbs.Pday),rename=T))
  (ms <- modSel.oSCR(fl))
  ms$aic
  ms$coef
}
\keyword{datasets}
