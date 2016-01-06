\name{e2dist}
\alias{e2dist}
\title{
Compute the distance between trap locations and animal activity centers.
}
\description{
A fast alternative to using 'for' loops for computing a distance matrix. 
}
\usage{
e2dist(x, y)
}
\arguments{
  \item{x}{
    Coordinates of animal activity centers. A matrix with one row for
    each animal and 2 columns for easting and northing.
}
  \item{y}{
    Coordinates of trap locations. A matrix with one row for
    each trap and 2 columns for easting and northing.
}
}
\value{
  A distance matrix with one row per animal and one column per trap.
}
\author{
Andy Royle, aroyle@usgs.gov 
}
\examples{

nAnimals <- 10
nTraps <- 5
activity.centers <- cbind(runif(nAnimals), runif(nAnimals))
trap.locs <- cbind(seq(0, 1, length=nTraps), seq(0, 1, length=nTraps))
Dmat<- e2dist(activity.centers, trap.locs)


}
