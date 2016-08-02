\name{spiderplot}
\alias{spiderplot}
\title{
Makes a spider plot of some SCR data.
}
\description{
Makes a spider plot and computes some summaries of individual captures.
}
\usage{
spiderplot(y3d,traplocs, add=FALSE)
}
\arguments{
  \item{y3d}{
 3-d encounter history array: nind x ntraps x nreps
}
  \item{traplocs}{
   ntraps x 2 matrix of trap coordinates.
}
\item{add=FALSE}{
 logical, if TRUE will add to current plot device.
 }
 
}
\value{
 Makes a spider plot. Returns a list with 

  xcent = distance from trap array centroid of each guy.

  avg.s = average capture location

  center = centroid of trap array
}
\author{
Andy Royle, aroyle@usgs.gov
}
\examples{

library("scrbook")

# make spiderplot for the beardata the hard way by converting the 
# encounter data file to a 3-d array:
data(beardata)
X<-as.matrix(cbind(id=1:nrow(beardata$trapmat),beardata$trapmat))
opps<-matrix(1,nrow=nrow(X),ncol=8)
dimnames(opps)<-list(NULL,1:8)
X<-cbind(X,opps)
a<-SCR23darray(beardata$flat,X)
toad<-spiderplot(a,beardata$trapmat)
title("Fort Drum hair snag array and mean capture locations")
# now grab the distance from centroid variable
xcent<-toad$xcent
# see Chapter 3 of the book

# the easy way
spiderplot(beardata$bearArray,beardata$trapmat)
title("Fort Drum hair snag array and mean capture locations")

# for the wolverine data:
data(wolverine)
y3d <- SCR23darray(wolverine$wcaps,wolverine$wtraps)
toad<- spiderplot(y3d,wolverine$wtraps[,2:3])
title("Wolverine camera trap array and mean capture locations")


##
## Here's a full script to fit the "individual covariate" model
## to the Fort Drum black bear data
## Provided by Bob Klaver (bklaver@iastate.edu)
##
library(scrbook)
library(R2WinBUGS)
data(beardata)
str(beardata)

trapmat<-beardata$trapmat
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]
M=200

nz<-M-nind

Yaug <- array(0, dim=c(M,ntraps,K))
Yaug[1:nind,,] <- beardata$bearArray

y<- apply(Yaug,c(1,3),sum) # summarize by ind x rep

y[y>1]<- 1                 # toss out multiple encounters b/c
                           #    traditional CR models ignore space

ytot <- apply(y, 1, sum)

bear.spider <- spiderplot(beardata$bearArray, beardata$trapmat)
xcent <- bear.spider$xcent

Xarg <- vector("numeric", length=nz)
xcent <- c(xcent, rep(NA,nz))

B <- 11.5

set.seed(2013)             # to obtain the same results each time

cat("model {
    p0 ~dunif(0, 1)
    alpha0 <- log(p0/(1-p0))
    alpha1 ~ dnorm(0, 0.01)
    psi ~dunif(0,1)
    
    for (i in 1:(nind+nz)) {
      xcent[i] ~dunif(0, B) 
      z[i] ~dbern(psi)       # zero inflated variables
      lp[i] <- alpha0 + alpha1*xcent[i]
      logit(p[i]) <- lp[i]
      mu[i] <- z[i]*p[i]
      y[i] ~dbin(mu[i], K)  # observation model
    }
    N<-sum(z[1:(nind+nz)])
} ",file="modelMcov.txt")


zst=c(rep(1,nind),rbinom(M-nind, 1, .5))

data2 <- list(y=ytot, nz=nz, nind = nind, K=K, xcent=xcent, B=B)
params2 <- c("p0", "psi", "N", "alpha0", "alpha1")
inits =  function() {  list(z=as.numeric(ytot >= 1), p0=runif(1),
                       alpha1=rnorm(1)) }

fit2 = bugs(data2, inits, params2, model.file="modelMcov.txt",n.chains=3,
            n.iter=110000, n.burnin=10000, n.thin=1,
            working.directory=getwd(), debug=TRUE)

print(fit2, digits=2)

plot(density(fit2$sims.list$N))

summary(fit2$sims.list$N/277.11)





}
