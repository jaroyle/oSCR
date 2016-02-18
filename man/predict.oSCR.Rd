\name{predict.oSCR}
\alias{predict.oSCR}
\docType{data}
\title{
  spatial predictions from an oSCR.fit object
}
\description{
computes posterior distribution of each s and also the frequency of
uncaptured individuals. Adds them up to make a density map
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{}
\format{
this is not a data file, just an R script (see below)
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

library(scrbook)
library(oSCR)
# Prediction example

# Fit a model using oSCR.fit
out1<- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1),ssDF=ssDF,plotit=FALSE, 
multicatch=TRUE , trimS=3,start.vals=out1$outStats[,3])

# make a prediction by feeding in the MLEs and use predict=TRUE
out1p <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1),ssDF=ssDF,plotit=FALSE, 
multicatch=TRUE , trimS=3,start.vals=out1$outStats[,3],predict=TRUE)

nguys<- dim(out1p$preds[[1]])[1]
# Uses function spatial.plot from scrbook to make a graphical depicture of 
# the posterior of s for each individual including the all-0 encounter history
par(mar=c(3,3,3,8))

# This will loop over the number of individuals and make the posterior
#   plot, pausing after each until you hit <enter>
for(guy in 1:nguys){
post<- out1p$preds[[1]][guy,]
# can change the cut-off threshold here to smaller to get more posterior mass
post[post< 0.00001]<- NA   # set to NA to make the plot more informative
post<- log(post)
spatial.plot(out1p$ssDF[[1]],post,cx=3,col=terrain.colors(10))
traps<- out1p$traps[[1]]
points(out1p$traps[[1]],pch=20)

ncaps<- apply(out1p$data[[1]][guy,,],1,sum)
if(any(ncaps>0)){
points(traps[ncaps>0,],pch=20,col="white")
text(traps[ncaps>0,],as.character(ncaps[ncaps>0]))
}
browser()
}


## maybe str(predict.oSCR) ; plot(predict.oSCR) ...
}
\keyword{datasets}
