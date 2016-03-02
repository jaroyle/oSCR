\name{extract.rast}
\alias{extract.rast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function for extracting raster data using an oSCR state space (ssDF) object.}
\description{
The function extracts data from a raster based uding X Y coordinated of the oSCR state space (ssDF) object using and appends the data to the existing data frame. The function will take multi session state space objects as input. The resolution of the raster summarires are defined by the discretization resolution of the state space.}
\usage{
extract.rast(ss, rast, mult = 1, cov.name = "val.1", func = median)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{An oSCR state space object, ssDF.}
  \item{rast}{A raster containing data from which density covariates are to be extracted.}
  \item{mult}{A multiplier used scale ssDF coordinated to the coordinate scale of the raster.}
  \item{cov.name}{The name that will be given to the generated covariate.}
  \item{func}{The function used to summarise the raster within the ssDF generated raster.}
}
\details{
}
\value{
}
\references{
}
\author{Chris Sutherland}
\note{}

\seealso{
}
\examples{

Mink - NLCD example to go here soon!
## The function is currently defined as
function (ss, rast, mult = 1, cov.name = "val.1", func = median) 
{
    for (i in 1:length(ss)) {
        tmpS <- ss[[i]][, c("X", "Y")] * mult
        id <- factor(1:nrow(tmpS))
        tmpR <- cbind(tmpS, id)
        aa <- rasterFromXYZ(tmpR, crs = projection(rast))
        bb <- rasterToPolygons(aa)
        r1 <- extract(rast, bb, fun = func)
        ss[[i]][, cov.name] <- r1
    }
    return(ss)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
