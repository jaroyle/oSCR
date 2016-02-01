\name{rbs}
\alias{rbs}
\docType{data}
\title{
The slow worm data from Beni Schmidt
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("slowworm")}
\format{
has ssDF and scrFrame for the slowworm

}
\details{
Don't use the ssDF that comes with the data object. Rather, use the make.ssDF function to create a new one as below.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Sutherland and Royle book chapter. Meier thesis. Schmidt et al. manuscript. 
%%  ~~ possibly secondary sources and usages ~~
}
\examples{



########################################################
### The following analyiss produces what is used in the Sutherland and Royle chapter
###########################

load("slowworm.RData")
ssDF.new <- make.ssDF(scrFrame, res=10, buffer=75)
 
temp1 <-  oSCR.fit(scrFrame, model=list(D~1,p0~1,a1~1),ssDF=ssDF.new, plotit=FALSE, multicatch=TRUE)

p0.hat<- plogis(temp1$outStats[1,3])
sigma.hat<- sqrt(1/(2*exp(temp1$outStats[2,3])))

d0<- exp(temp1$outStats[3,3])
Nhat<- d0*nrow(ssDF.new[[1]])
Dhat<- (d0/100)*10000

y<- flatten(y1)
otis(y)

}
\keyword{datasets}
