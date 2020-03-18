\name{pakistan}
\alias{pakistan}
\docType{data}
\title{
Landscape and detector coordinates for use in scrdesignGA().
}
\description{
Landscape and detector coordinates for use in scrdesignGA().
}
\usage{data("pakistan")}
\format{
pakSS is a two column dataframe of pixel coordinates for the entire landscape.
pakTT is a two column dataframe of coordinates of all possible sampling locations.
}
\details{
loads the pakistan lanscape and detector data.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

library(oSCR)
data("pakistan")
plot(pakSS, asp=1, pch=16, cex=0.25)
points(pakTT, pch=16, col="grey")

}
\keyword{datasets}
