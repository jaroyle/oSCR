\name{bcharea}
\alias{bcharea}
\docType{data}
\title{
buffer a convex hull and compute area
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(bcharea)}
\format{
  The format is:
 chr "bcharea"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library("rgeos")
library("scrbook")
library("maptools")
data(beardata)
trapmat<-beardata$trapmat

B.MMDM<- 2.373

aa<-chull(trapmat)
buffered.area<- area(as.matrix(trapmat[aa,]))

bcharea(B.MMDM,traplocs=trapmat)
cat("area of buffered convex hull: ",buffered.area,fill=TRUE)

}
\keyword{datasets}
