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


predict.oSCR <- function(scrFrame, scr.fit, ssDF, costDF=NULL){
 
  mles <- scr.fit$rawOutput$estimate
  call <- scr.fit$call
  
  oSCR.fit2 <- oSCR.fit
  # modify the arguments to match the original call (skip the data)
  call.fix <- names(call)[!names(call) %in% c("","scrFrame","ssDF","costDF")]
  for(i in 1:length(call.fix)){formals(oSCR.fit2)[[call.fix[i]]]<- call[[call.fix[i]]]}
  prediction <- oSCR.fit2(scrFrame=scrFrame, ssDF=ssDF, costDF=costDF,
                          #critical arguments that differ for prediction
                          start.vals=mles, predict=TRUE)
  return(prediction)
}



library(scrbook)
library(oSCR)
# Prediction example

make.statespace <- 
function (ll = NULL, minx = NA, maxx = NA, miny = NA, maxy = NA, 
    nx = 40, ny = NULL, buffer = 0,delta) 
{
    if (is.null(ny)) 
        ny <- nx
    if (!is.null(ll)) {
        minx <- min(ll[, 1])
        maxx <- max(ll[, 1])
        miny <- min(ll[, 2])
        maxy <- max(ll[, 2])
        minx <- minx - buffer
        maxx <- maxx + buffer
        miny <- miny - buffer
        maxy <- maxy + buffer
    }
    xgr <- seq(minx + delta/2, maxx - delta/2, delta)
    ygr<-  seq(miny + delta/2, maxy - delta/2, delta)
    gr<- as.matrix(expand.grid(xgr,ygr))
    dimnames(gr)<-list(NULL,c("X","Y"))
gr
}

set.seed(2014)
 
 

 
parms<-list(N=100, alpha0= -.40, sigma=0.5,alpha2=0)
ssbuff<- 2

truth<- "multicatch"   # "binomial"
with.sex<- TRUE
 
# simulate some data

if(truth=="multicatch"){
 data<-simMnSCR(parms,K=7,ssbuff=2)
 traplocs<- data$X
 nind<-nrow(data$Ycat)
 Ycat<-data$Ycat
 # K = 7, ntraps = 25
 y3d<- array(0,dim=c(nrow(Ycat), 25, 7))
 
for(i in 1:nrow(Ycat)){
 for(k in 1:7){
   if(Ycat[i,k]<26)
   y3d[i,Ycat[i,k],k]<- 1
  }
 }
}
if(truth=="binomial"){
 # Careful here because I have a seed hard-wired into the code to produce the book data set
data<- simSCR0(N=parms$N, K=7, alpha0 = parms$alpha0, sigma = parms$sigma, array3d=TRUE, rnd=2014)
y3d<- data$Y
nind<- dim(y3d)[1]
traplocs<- data$traplocs
}
 

if(with.sex)
  sex<- rbinom(nind,1,0.5)

traplocs <- cbind(sort(rep(1:5, 5)), rep(1:5, 5))
 
## Convert data for oSCR

# trap list: dataframe with named X and Y columns
traplocs <- list(data.frame(X=traplocs[,1],Y=traplocs[,2]) )
  
if(with.sex){
  scrFrame <- list(caphist=list(y3d), indCovs=list(data.frame(sex=sex)),
                   type="scr", traps=traplocs,trapCovs=NULL, nocc= dim(y3d)[3] )
  class(scrFrame) <- "scrFrame"
}
if(!with.sex){
  scrFrame <- list(caphist=list(y3d), indCovs=NULL, 
                   type="scr", traps=traplocs,trapCovs=NULL, nocc= dim(y3d)[3] )
  class(scrFrame) <- "scrFrame"
}
ssDF <- list(     data.frame(      make.statespace(traplocs[[1]], delta =.5, buffer = 2) ) )


# Fit a model using oSCR.fit
out1<- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1),ssDF=ssDF,plotit=FALSE, 
multicatch=TRUE , trimS=3)

out1p<- predict.oSCR(scrFrame, scr.fit=out1, ssDF)

plot(out1p[[1]])
title("Estimated local density: E( N(s) | data )")

# Summing up the local density value = Nhat

 sum(values(out1p[[1]]))
#[1] 104.2319


# make a prediction by feeding in the MLEs and use predict=TRUE
#out1p <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1),ssDF=ssDF,plotit=FALSE, 
#multicatch=TRUE , trimS=3,start.vals=out1$outStats[,3],predict=TRUE)



###
## BELOW HERE IS OBSOLETE NEEDS REVISED to produce single-individual predictions
###

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
#browser()
}

Nhat<- sum(out1p$d.s)

# number of individuals not detected at each pixel
# this is wrong, should involve lik.condl ... ie Pr(y=0 | s)
n0<- out1p$d.s*out1p$lik.cond
# probability that a guy at pixel s is NOT captured
plot(rasterFromXYZ(cbind(out1p$ssDF[[1]],out1p$lik.cond)))
out1p$preds[[1]][nguys,]<- n0
# Does not include the "n0" weighting yet
library(raster)
total<- apply(out1p$preds[[1]],2,sum)
# check
sum(total)
r<- rasterFromXYZ(cbind(out1p$ssDF[[1]],total))
plot(r)
title("Density map")

## maybe str(predict.oSCR) ; plot(predict.oSCR) ...
}
\keyword{datasets}
