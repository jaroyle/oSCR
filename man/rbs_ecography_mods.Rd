\name{rbs_ecography_mods}
\alias{rbs_ecography_mods}
\docType{data}
\title{
Red-backed salamander \emph{Plethodon cinereus} fitted models.
}
\description{
These are the 16 fitted models used in the example in Sutherland \emph{et al.} (2019), the official \code{oSCR} publication.
}
\usage{data("rbs_ecography_mods")}
\format{
The 16 \code{oSCR.fit} model objects, from \code{m1} through \code{m16}.

}
\details{

}
\source{
Sutherland, C., Royle, J. A., & Linden, D. W. (2019). oSCR: a spatial capture-recapture R package for inference about spatial ecological processes. \emph{Ecography}, 42(9), 1459-1469.

Ecography supplement: http://www.ecography.org/appendix/ecog-04551
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library(oSCR)
data(rbs_ecography)
data(rbs_ecography_mods)
ls()[1:16]

# model 1 output
m1

# model selection
fl <- fitList.oSCR(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16),
drop.asu=T, rename=TRUE) #rename=T adds sensible model names
ms <- modSel.oSCR(fl)
ms
}
\keyword{datasets}
