\name{telemetry}
\alias{telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integration of telemtry data with SCR data%%  ~~function to do ... ~~
}
\description{
This is an R script that demonstrates how to integrate telemetry data
with SCR models using oSCR...%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
telemetry()  # returns NULL
}
%- maybe also 'usage' for other objects documented here.
\details{
See below %%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Royle, J. Andrew, Richard B. Chandler, Catherine C. Sun, and Angela
K. Fuller. 2013. "Integrating resource selection information with spatial
capture-recapture." Methods in Ecology and Evolution 4(6):
520-530.

Linden, D. W., Sirén, A. P. K., and Pekins, P. J. 2018. "Integrating telemetry data into spatial capture–recapture modifies inferences on multi‐scale resource selection." Ecosphere 9(4):e02203. 10.1002/ecs2.2203 
}
\author{
Andy Royle, Dan Linden 
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
library(scrbook)

# test telemetry data
# based on code from supplement of Royle et al (2013) - MEE

## the following block of code makes up a covariate as a spatially correlated
## noise field, with an exponential spatial correlation function
set.seed(1234)
gr<-expand.grid(1:40,1:40)
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
# Below change to matrix multiplication
z<-as.vector( crossprod( t(chol(V)), rnorm(1600) ) )
spatial.plot(gr,z,cx=3)
library(raster)
r<- rasterFromXYZ(cbind(gr, z))
plot(r)


###
### Set some parameter values
###
alpha0 <- -2
sigma<- 2
alpha2<- 1   # This is the effect of the covariate on resource selection
Ntel<-4      # number of individuals with telemetry devices
nsim<-100
Nfixes<-20   # number of telemetry fixes per individual. Note: these
### should probably be regarded as indepedent locations of
### home range usage. 
N<- 100      # population size

#seed <- sample(1:999,1)
#set.seed(seed)

# Simulate activity centers of all N individuals in the population
Sid<- sample(1:1600,N,replace=TRUE)
# and coordinates
S<-gr[Sid,]

## Now lets simulate some SCR data:

# Make a trap array
X<- cbind( sort(rep( seq(5,35,5),7)), rep( seq(5,35,5),7))
ntraps<-nrow(X)
raster.point<-rep(NA,nrow(X))
# This just maps the trap locations to the raster cells
for(j in 1:nrow(X)){ 
  raster.point[j]<- (1:1600)[ (X[j,1]==gr[,1]) & (X[j,2] == gr[,2])]
}
points(X,pch=20,cex=2)

# Hazard model is used. This seems the most sensible. 
D<- e2dist(S,X) ## N x ntraps
Zmat<- matrix(z[raster.point],nrow=N,ncol=ntraps,byrow=TRUE) # note make dims the same
loglam<- alpha0 -(1/(2*sigma*sigma))*D*D + alpha2*Zmat
p<- 1-exp(-exp(loglam))
# Now simulate SCR data
K<- 10
y<-matrix(NA,nrow=N,ncol=ntraps)
for(i in 1:N){
  y[i,]<- rbinom(ntraps,K,p[i,])
}

# Subset data to captured individuals
cap<-apply(y,1,sum)>0
y<-y[cap,]
gr<-as.matrix(gr)
# Average capture location
sbar<- (n%*%gr)/as.vector(n%*%rep(1,nrow(gr)))

# Now draw centers of telemetered individuals
# We have to draw telemetry guys interior (i.e., the study area) or else make up more landscape --
# because we can't have truncated telemetry observations. 
poss.tel<- S[,1]>5 & S[,1]<35 & S[,2]>5 & S[,2]<35

# Need to account for capture so that dependent models can be explored
tel.guys.id <- sort(sample(which(cap & poss.tel),Ntel))
tel.guys <- Sid[tel.guys.id] #which s for telemetry guys (these 4 purposely selected)
cap.tel <- match(tel.guys.id,which(cap))   #which row in capture history for each telemetry guy
sid<-tel.guys
stel<-gr[sid,]

# Make a matrix to store RSF data (telemetry fixes)
n<-matrix(NA,nrow=Ntel,ncol=1600)

# for each telemetered guy simulate a number of fixes.
# note that nfix = 0 for most of the landscape pixels
par(mfrow=c(3,3))
lammat<-matrix(NA,nrow=Ntel,ncol=1600)
for(i in 1:Ntel){
  d<- Dmat[sid[i],]
  lam<- exp(1 - (1/(2*sigma*sigma))*d*d + alpha2* z)
  n[i,]<-rmultinom(1,Nfixes,lam/sum(lam))
  par(mar=c(3,3,3,6))
  lammat[i,]<-lam
  img<- matrix(lam,nrow=40,ncol=40,byrow=FALSE)
  image(1:40,1:40,(img),col=terrain.colors(10))
}



plot(S,pch=15,col="lightgray")
points(stel,pch=15,col=c("green","blue","purple","orange"))
points(X,pch=3,col="red")
for(j in 1:Ntel){
  points(X[which(y[cap.tel[j],]>0),,drop=FALSE],pch=22,col=c("green","blue","purple","orange")[j])
}

# Set up the data for oSCR
library(oSCR)
# The habitat variable "z" is used as a trapCov (as it would even if you
#   didn't have telemetry data)
trapCovs <- list(z=list(data.frame(z=matrix(z[raster.point],ntraps,K))))
trapCovs <- make.trapCovs(trapCovs)
# Set up the ssDF and the rsfDF
ssDF <- data.frame(X=gr[,1],Y=gr[,2],z=z)
rsfDF <- data.frame(X=gr[,1],Y=gr[,2],z=z)
# Simulate sex of individuals 
sex.cap <- rbinom(nrow(y),1,.5)

# 
# Distribute the binomial captures (n x ntraps) among the K surveys
# ...could also simply simulate binary captures (n x ntraps x K)
# Should  only be done when time effects are not important
y.arr <- array(0,dim=c(nrow(y),ntraps,K))
for(i in 1:nrow(y)){
  for(j in 1:ntraps){
    which.K <- sample(1:K,y[i,j])
    y.arr[i,j,which.K] <- 1
  }
}

# KEY STEP:
# Create the telemetry list (use telemetry.processor() for raw fix data,
# see below).  The inclusion of "cap.tel" indicates that some collared
# individuals were also captured -- this vector includes the row position
# in the capture history array for each collared individual.  If some
# collared individuals were not captured, they should be sorted last,
# possibly by assigning uniqe IDs with additional characters (e.g., 99 or XX).
telemetry <- list(fixfreq=list(n),indCovs=list(data.frame(sex=c(0,1,0,1))),
                  cap.tel=list(cap.tel))

# Create the scrFrame 
sftel <- make.scrFrame(caphist = list(y.arr),
                       traps = list(data.frame(X=X[,1],Y=X[,2])), 
                       #trapCovs = trapCovs,
                       telemetry = telemetry,
                       rsfDF = list(rsfDF),
                       indCovs = list(data.frame(sex=sex.cap))
)

###########################
#     Model fitting       #
###########################

# fit the SCR model with NO telemetry integration (z from traps only)
fit1 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 #rsfDF=list(rsfDF),RSF=TRUE,
                 trimS=7.5,
                 model=list(D~1,p0~z,sigma~1,path~1))

# Next we use the telemetry information to inform about both beta.z (RSF=TRUE) and
# also 'sigma', assuming independence between data (captures vs. collars)
fit2 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 rsfDF=list(rsfDF),RSF=TRUE,telemetry="ind",
                 trimS = 7.5,
                 model=list(D~1,p0~z,sigma~1,path~1))

# Now we use the telemetry information again (RSF = TRUE) but assume dependence
# between data since some collared guys were captured
fit3 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=TRUE,telemetry="dep",
                trimS = 7.5,
                model=list(D~1,p0~z,sigma~1,path~1))

# Here we fit the SCR model with RSF = FALSE, which only uses the
# telemetry data to inform about 'sigma' NOT the RSF parameters
fit4 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                  rsfDF=list(rsfDF),RSF=FALSE,telemetry="ind",
                  trimS = 7.5,
                  model=list(D~1,p0~z,sigma~1,path~1))


#
# Same models but no trapcovs. They are made internally from the rsfDF when provided
# to make.scrFrame.  Note: it is mandatory for rsfDF to be provided to
# oSCR.fit AND make.scrFrame (may be changed in the future).
#
sftel <- make.scrFrame(caphist = list(y.arr),
                       traps = list(data.frame(X=X[,1],Y=X[,2])), 
                       trapCovs = NULL, rsfDF=list(rsfDF),
                       telemetry = telemetry,
                       indCovs = list(data.frame(sex=sex.cap))
)

# fit the SCR model with telemetry integration
fit2b <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 rsfDF=list(rsfDF),RSF=TRUE,
                 model=list(D~1,p0~z,sigma~1,path~1))



#
#
# Now fit a model to the NY bear data
#
#
#


library(oSCR)
#
#---------------------------------------------------------------------------#
# 
data("nybears")
attach(nybears)
fixes <- teldata[,c("animalid","X_UTM","Y_UTM")]
colnames(fixes)<- c("ind","X","Y")

ntraps <- nrow(traplocs)
traplocs <- traplocs/10000
colnames(traplocs) <- c("X","Y")
trapCovs <- list(z=list(data.frame(z=matrix(elevation[trap2raster],ntraps,K))))
trapCovs <- make.trapCovs(trapCovs)
ssDF <- rsfDF <- data.frame(ssgrid,z=elevation)
sex.cap <- rbinom(nrow(y2d),1,.5)

# distribute the binomial captures (n x ntraps) among the K surveys
# could simply simulate a binary array of captures (n x ntraps x K)
y.arr <- array(0,dim=c(nrow(y2d),ntraps,K))
for(i in 1:nrow(y2d)){
  for(j in 1:ntraps){
    which.K <- sample(1:K,y2d[i,j])
    y.arr[i,j,which.K] <- 1
  }
}

# create telemetry list using precalculated fix frequencies (raw fix data requires telemetry.processor function)
telemetry <- list(fixfreq=list(ntel),indCovs=list(data.frame(sex=c(0,1,0))))
# create telemetry list by calculating fix frequences from raw fix data

# teldata2[,2:3]<- teldata2[,2:3]/10000
#telemetry <- list(fixfreq=telemetry.processor( list(rsfDF),list(teldata2))$nfreq,
#                 indCovs=list(data.frame(sex=c(0,1,0))))


# create the scrFrame
sftel <- make.scrFrame(caphist = list(y.arr),
                       traps = list(traplocs), 
                       trapCovs = trapCovs,
                       telemetry = telemetry,
                       indCovs = list(data.frame(sex=sex.cap))
)

# fit the SCR model with telemetry integration
fit <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=TRUE,telemetry="ind",
                model=list(D~1,p0~z,sigma~1,path~1))



#
# Same model but here converting from raw GPS coordinates into frequencies
# using the telemetry.processor function. NOTE: we thin the raw telemetry data
# 
kp<- seq(1,nrow(teldata2),length= 0.10*nrow(teldata2)) #keep x% of fixes
teldata2.thin <- teldata2
teldata2.thin[,2:3]<- teldata2[,2:3]/10000
teldata2.thin<- teldata2.thin[kp,]

# Note list structure needed for telemetry.processor, just like other oSCR functions
telemetry2 <- list(fixfreq=telemetry.processor(list(rsfDF),list(teldata2.thin))$nfreq,
                  indCovs=list(data.frame(sex=c(0,1,0))))

plot(ssDF[,1:2],pch=15,col="lightgray",main="telemetry")
points(ssDF[which(apply(telemetry$fixfreq[[1]],2,sum)>0),1:2],pch=15,col=rgb(1,0,0,1))
plot(ssDF[,1:2],pch=15,col="lightgray",main="telemetry2")
points(ssDF[which(apply(telemetry2$fixfreq[[1]],2,sum)>0),1:2],pch=15,col=rgb(1,0,0,1))

# create the scrFrame
sftel <- make.scrFrame(caphist = list(y.arr),
                       traps = list(traplocs), 
                       trapCovs = trapCovs, 
                       telemetry = telemetry2,
                       indCovs = list(data.frame(sex=sex.cap)) )

# fit the SCR model with telemetry integration
fit <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=TRUE,telemetry="ind",
                model=list(D~1,p0~z,sigma~1,path~1))








}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
