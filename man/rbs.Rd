\name{rbs_ecography}
\alias{rbs_ecography}
\docType{data}
\title{
Red-backed salamander \emph{Plethodon cinereus} data.
}
\description{
The rbs data is the data analyzed by Sutherland \emph{et al.} (2016) in their \emph{Herpetologia} paper and presented again in the official \code{oSCR} publication (Sutherland \emph{et al.} 2019). This is a four-session data set where sessions are independent artificial cover object (ACO) plots in Ithaca, NY, USA.
}
\usage{data("rbs_ecography")}
\format{
 The data contains an encounter data file (\code{rbs.edf}) and 4 trap deployment files ( \code{rbs.tdf1}, \code{rbs.tdf2}, \code{rbs.tdf3}, \code{rbs.tdf4}) that are used to build objects necessary for \code{oSCR.fit}.
}
\details{

}
\source{
Sutherland, C., Munoz, D. J., Miller, D. A., & Grant, E. H. C. (2016). Spatial capture-recapture: a promising method for analyzing data collected using artificial cover objects. \emph{Herpetologica}, 72(1), 6-12.

Sutherland, C., Royle, J. A., & Linden, D. W. (2019). oSCR: a spatial capture-recapture R package for inference about spatial ecological processes. \emph{Ecography}, 42(9), 1459-1469.

Ecography supplement: http://www.ecography.org/appendix/ecog-04551

}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library(oSCR)
data(rbs_ecography)
ls()

# trap deployment (session 4)
str(rbs.tdf4)

# encounter data
str(rbs.edf)

# prepare data for analysis
rbs.data <- data2oscr(
  edf = rbs.edf,                     #the EDF
  sess.col = 1,                      #session column
  id.col = 2,                        #individual column
  occ.col = 3,                       #occasion column
  trap.col = 4,                      #trap column
  sex.col = 5,                       #sex column (optional)
  tdf = list(rbs.tdf1, rbs.tdf2,     #list of TDFs
             rbs.tdf3, rbs.tdf4),
  K = c(7,5,6,4),                    #occasion vector
  ntraps = c(50,50,50,50),           #no. traps vector
  trapcov.names = c("jday","jday2"), #covariate names
  tdf.sep = "/",                     #char used separator
  sex.nacode = "U")                  #unknown sex code

names(rbs.data)
rbs.sf <- rbs.data$scrFrame
names(rbs.sf)
rbs.sf

# plot the encounters
par(mfrow=c(2,2), mar=c(2,2,2,2), oma=c(0,0,0,0))
plot(rbs.sf,jit = 2)

# make the state space
rbs.ss <- make.ssDF(scrFrame = rbs.sf, buffer = 4, res=0.5)
str(rbs.ss)

# plot the state space
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(ssDF = rbs.ss, scrFrame = rbs.sf)

# fit a model (ART ~100 min)
m1 <- oSCR.fit(
  model = list(D~1,                    #density
               p0~b + jday + jday2,    #detection
               sig~1),                 #space use
  scrFrame = rbs.sf, ssDF = rbs.ss)
m1

}
\keyword{datasets}
