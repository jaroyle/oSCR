\name{predict.oSCR}
\alias{predict.oSCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Produces data objects necessary to create a map of realized density
  and posterior distributions of activity centers. 
}
\description{
Computes the estimated posterior distribution for each activity centers
and the total realized density.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.oSCR(scrFrame, scr.fit, ssDF, costDF = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
the scrFrame used in fitting the model for which predictions are required%%     ~~Describe \code{scrFrame} here~~
}
  \item{scr.fit}{
the oSCR fit object for which predictions are required %%     ~~Describe \code{scr.fit} here~~
}
  \item{ssDF}{
the state-space object used for obtaining the model fit for which
predictions are required %%     ~~Describe \code{ssDF} here~~
}
  \item{costDF}{
the cost surface dataframe for obtaining the model fit for which
predictions are required %%     ~~Describe \code{costDF} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andy Royle
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(oSCR)
data(ocelot)

# Convert the data format
data<-secr2oscr(rbind(edf1,edf2),sex.col=5, K=c(44,96),   ntraps=c(23,23), 
             tdf=list(tdf1,tdf2),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 


# Process the trap covariate info. This should be done internally but isn't yet.
# There is one covariate "distance to disturbance" which is not defined for session = 1
# therefore I set it to 0 here so that the linear predictor evaluates to 0 (i.e., no adjustment to p).

dist1 =  rep(0,23)
dist2 =  as.vector(data$trapcovs[[2]])/1000   # Convert to km
dist1 <- matrix(dist1, nrow=23,ncol=44,byrow=FALSE)
dist2 <- matrix(dist2, nrow=23,ncol=96, byrow=FALSE)
 
# Example of 1 covariate
tc <-  list(dist=list(dist1, dist2) )
tc <-  make.trapCovs(tc)

# Make the SCRframe
scrFrame  <- list(caphist=data$y3d, indCovs=NULL, 
                   type="scr", traps=data$traplocs,trapCovs=tc ,
                   trapOperation=data$trapopp )
class(scrFrame) <- "scrFrame"

# make a state-space
ssDF <- make.ssDF(scrFrame, buffer=2000, res = 600)

# Make a descriptive plot of everything
plot(ssDF[[1]], pch=20)
points(data$traplocs[[1]],pch=3,col="red",lwd=2)
spiderplot(data$y3d[[2]],data$traplocs[[1]],add=TRUE)

# Fit a basic model
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame, ssDF=ssDF,plotit=FALSE , trimS=2000)

## default vals for some models are especially bad. Could be geographic scaling issue.
out2 <- oSCR.fit(model=list(D~1,p0~dist,sig~1), scrFrame,start.vals=c(-3.5, 6.4, 0, -1.5), ssDF=ssDF, plotit=FALSE, trimS=2000  )

out3 <- oSCR.fit(scrFrame,model=list(D~session,p0~1,sig~1),start.vals=c(-3.5, 6.4, -1.5, -1.5), ssDF=ssDF,plotit=FALSE , trimS=2000)

## Summarize output
## Make a density map

pred<- predict.oSCR(scrFrame,out3,ssDF)

# Plot total density map
par(mfrow=c(1,2))
plot(pred$r[[1]])
plot(pred$r[[2]])

# Plot posterior of s for a specific individual
par(mfrow=c(1,1))
rs<- rasterFromXYZ(cbind(pred$ssDF[[1]][,1:2], pred$preds[[1]][1,]))
plot(rs)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
