\name{telemetry}
\alias{telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integration of telemetry data with SCR data%%  ~~function to do ... ~~
}
\description{
This is an R script that demonstrates how to integrate telemetry data
with SCR models using oSCR...%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
?telemetry
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

Linden, D. W., Siren, A. P. K., and Pekins, P. J. 2018. "Integrating 
telemetry data into spatial capture-recapture modifies inferences on 
multi-scale resource selection." Ecosphere 9(4):e02203. 10.1002/ecs2.2203 
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
library(tidyverse)

# test telemetry data
# based on code from supplement of Royle et al (2013) - MEE

## the following block of code makes up a covariate as a spatially correlated
## noise field, with an exponential spatial correlation function
set.seed(1234)
gr <- expand.grid(X = 1:40, Y = 1:40)
Dmat <- as.matrix(dist(gr))
V <- exp(-Dmat / 5)
# below change to matrix multiplication
z <- as.vector(crossprod(t(chol(V)), rnorm(1600)))

# create the state-space dataframe and plot it
ssDF <- data.frame(gr,z=z)
ss.plot <- ggplot(ssDF,aes(x=X,y=Y)) + coord_equal() + theme_minimal()
ss.plot + geom_tile(aes(fill=z)) + scale_fill_viridis_c()

###
### set some parameter values
###
alpha0 <- -1.5  # log(encounter rate)
sigma<- 2       # movement scale
alpha2<- 1      # effect of covariate on resource selection
Ntel<-4         # number of individuals with telemetry devices
Nfixes<-40      # number of telemetry fixes per individual
N<- 50          # population size

#
# simulating multiple realizations of capture/telemetry data
#
# n.iter <- 5
# seeds <- sample(1:9999,n.iter)
#  fit.sim <- list()
#  fit.sim[[1]] <- fit.sim[[2]] <- data.frame(
#    iter=1:n.iter,seed=seeds,p0.Int=NA,sig.Int=NA,t.beta.z=NA,d0.Int=NA,psi=NA
#  )


## Now let's simulate some SCR data:

# Make a trap array
X <- expand.grid(X=seq(5,35,5),Y=seq(5,35,5))
ntraps <- nrow(X)
raster.point <- rep(NA, nrow(X))
# This just maps the trap locations to the raster cells
for (j in 1:nrow(X)) {
  raster.point[j] <- (1:1600)[(X[j, 1] == gr[, 1]) &
                                (X[j, 2] == gr[, 2])]
}

#for (iter in 1:n.iter){
#set.seed(seeds[iter])
set.seed(9933)
  
# Simulate activity centers of all N individuals in the population
Sid<- sample(1:1600,N,replace=TRUE)
# and coordinates
S<-gr[Sid,]
  
# Hazard model is used. This seems the most sensible. 
D<- oSCR::e2dist(S,X) ## N x ntraps
Zmat<- matrix(z[raster.point],nrow=N,ncol=ntraps,byrow=TRUE) # note make dims the same
loglam<- alpha0 -(1/(2*sigma*sigma))*D*D + alpha2*Zmat
p<- 1-exp(-exp(loglam))
# Now simulate SCR data
K<- 3
y<-array(NA,dim=c(N,ntraps)) #,K))
for(i in 1:N){
    y[i,]<- rbinom(ntraps,K,p[i,])
}

# Subset data to captured individuals
cap<-apply(y,1,sum)>0
y<-y[cap,]

# Now draw centers of telemetered individuals
# We have to draw telemetry guys interior (i.e., the study area) 
# or else make up more landscape because we can't have truncated 
# telemetry observations. 
poss.tel<- S[,1]>5 & S[,1]<35 & S[,2]>5 & S[,2]<35

# Need to account for capture so that dependent models can be explored
tel.guys.id <- sort(sample(which(cap & poss.tel),Ntel))
tel.guys <- Sid[tel.guys.id] #which s for telemetry guys (these 4 purposely selected)
cap.tel <- match(tel.guys.id,which(cap))   #which row in capture history for each telemetry guy
sid <- tel.guys
stel <- gr[sid,]

# Make a matrix to store RSF data (telemetry fixes)
nfix<-matrix(NA,nrow=Ntel,ncol=1600)

# for each telemetered guy simulate a number of fixes.
# note that nfix = 0 for most of the landscape pixels
lammat<-matrix(NA,nrow=Ntel,ncol=1600)
for(i in 1:Ntel){
  d<- Dmat[sid[i],]
  lam<- exp(1 - (1/(2*sigma*sigma))*d*d + alpha2* z)
  nfix[i,]<-rmultinom(1,Nfixes,lam/sum(lam))
  lammat[i,]<-lam
}
# Average fix location (matches activity centers)
(sbar <- (nfix \%*\% as.matrix(gr))/as.vector(nfix \%*\% rep(1,nrow(gr))))
stel

# expected encounter rates (lambda) for collared guys
lamDF <- pivot_longer(
  data.frame(gr,t(lammat)),-c(1:2),names_to="ind",values_to="lam")
# plot lambdas by collared individual
ss.plot + geom_tile(data=lamDF,aes(fill=lam)) +
  scale_fill_viridis_c(direction = -1) + facet_wrap(~ind)

# collection of capture locations and s for each n
caps.X <- data.frame(which(apply(y,c(1,2),sum)>0,arr.ind=T))
caps.X[,c("X","Y")] <- X[caps.X$col,]; caps.X <- caps.X[order(caps.X$row),]
caps.X <- data.frame(caps.X,s=S[cap,][caps.X$row,])
caps.X$tel <- match(caps.X$row,cap.tel)

# plot the capture data, emphasis on collared guys
ss.plot +
  # activity centers of all N
  geom_point(data=S,pch=15,col="gray",size=4) +
  # activity centers of collared Ntel
  geom_point(data=stel,aes(col=factor(1:4)),size=4,pch=15) +
  # trap of captures for n
  geom_point(data=caps.X,col="gray",size=4,pch=16) +
  geom_segment(data=caps.X,aes(x=X,y=Y,xend=s.X,yend=s.Y),col="gray",lwd=1.25) +
  # trap of captures for Ntel
  geom_point(data=caps.X,aes(color=factor(tel)),size=4,pch=16) +
  geom_segment(data=caps.X,aes(x=X,y=Y,xend=s.X,yend=s.Y,col=factor(tel)),lwd=1.25) +
  # trap locations
  geom_point(data=X,pch=3,size=3,alpha=0.55) +
  scale_color_brewer(palette="Set1") + labs(color="ID")


# Set up the data for oSCR
# Set up the rsfDF (matches ssDF)
rsfDF <- ssDF

# 
# Distribute the binomial captures (n x ntraps) among the K surveys
# ...could also simply simulate binary captures (n x ntraps x K)
# Should  only be done when survey effects are not important
y.arr <- array(0,dim=c(dim(y)[1],ntraps,K))
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
# possibly by assigning unique IDs with additional characters (e.g., 99 or XX).
telemetry <- list(fixfreq=list(nfix),
                  cap.tel=list(cap.tel))

# Create the scrFrame 
sftel <- make.scrFrame(caphist = list(y.arr),
                       traps = list(X), 
                       #trapCovs = trapCovs,
                       telemetry = telemetry,
                       rsfDF = list(rsfDF)
)
# spatial recaptures
table(apply(apply(y,1:2,sum),1,function(x){length(which(x>0))}))

###########################
#     Model fitting       #
###########################

# fit the SCR model with NO telemetry integration (z from traps only)
fit1 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 #rsfDF=list(rsfDF),RSF=TRUE,
                 trimS=sftel$mdm,
                 model=list(D~1,p0~z,sigma~1,path~1))

# Next we use the telemetry information to inform about both beta.z (RSF=TRUE) and
# also 'sigma', assuming independence between data (captures vs. collars)
fit2 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 rsfDF=list(rsfDF),RSF=TRUE,telemetry="ind",
                 trimS=sftel$mdm,
                 model=list(D~1,p0~z,sigma~1,path~1))

data.frame(fit1=fit1$outStats,fit2=fit2$outStats[,2:3],
           truth = c(alpha0,log(sigma),alpha2,log(N/1600)))

# storing the simulations comparing fits
#fit.sim[[1]][iter,] <- c(iter,seeds[iter],fit1$coef.mle[,2])
#fit.sim[[2]][iter,] <- c(iter,seeds[iter],fit2$coef.mle[,2])
#}

# Now we use the telemetry information again (RSF = TRUE) but assume dependence
# between data since some collared guys were captured
fit3 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=TRUE,telemetry="dep",
                trimS=sftel$mdm,
                model=list(D~1,p0~z,sigma~1,path~1))

# Here we fit the SCR model with RSF = FALSE, which only uses the
# telemetry data to inform about 'sigma' NOT the RSF parameters
fit4 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 rsfDF=list(rsfDF),RSF=FALSE,telemetry="ind",
                 trimS=sftel$mdm,
                 model=list(D~1,p0~z,sigma~1,path~1))

# compare predictions of s for collared individuals between the independent
# and dependent model fits
#

# telemetry = independent
pred2 <- predict.oSCR(fit2)
pred.s2 <- pivot_longer(
  data.frame(gr,t(pred2$preds[[1]][cap.tel,])),-c(1:2),names_to="ind",values_to="Pr(s)")

# plot Pr(s) by collared individual
ss.plot + geom_tile(data=pred.s2,aes(fill=`Pr(s)`)) +
  scale_fill_viridis_c(direction = -1) + facet_wrap(~ind)

# telemetry = dependent
pred3 <- predict.oSCR(fit3)
pred.s3 <- pivot_longer(
  data.frame(gr,t(pred3$preds[[1]][cap.tel,])),-c(1:2),names_to="ind",values_to="Pr(s)")

# plot Pr(s) by collared individual
ss.plot + geom_tile(data=pred.s3,aes(fill=`Pr(s)`)) +
  scale_fill_viridis_c(direction = -1) + facet_wrap(~ind)


#
#---------------------------------------------------------------------------#
# 
#
# Now fit a model to the NY bear data
#
#
library(oSCR)
library(tidyverse)

data("nybears")
ls(nybears)

# capture data
edf <- nybears$edf
K <- nybears$K
# rescale spatial coordinates so each unit = 10km
tdf <- nybears$tdf
tdf[,2:3] <- tdf[,2:3]/1e4
colnames(tdf)[2:3] <- c("X","Y")
ntraps <- nrow(tdf)

ss <- nybears$ss/1e4
colnames(ss) <- c("X","Y")

# telemetry data
fixes <- nybears$teldata[,c("animalid","X_UTM","Y_UTM")]
colnames(fixes)<- c("ind","X","Y")
fixes[,c("X","Y")] <- fixes[,c("X","Y")]/1e4

# create the state space and RSF surfaces with covariate
ssDF <- rsfDF <- data.frame(ss,z=nybears$elevation)


# thin the fixes, given that hourly fix rate contains too much dependence
kp<- seq(1,nrow(fixes),length=0.10*nrow(fixes)) #keep x% of fixes
fixes.thin <- fixes[kp,]

# create the individual by pixel frequencies of fixes
nfix <- telemetry.processor(list(rsfDF),list(fixes.thin))$nfreq

# create the telemetry object for oSCR (note, nfix already a list)
telemetry <- list(fixfreq=nfix, indCovs=list(data.frame(sex=c(0,1,0))))

# create the scrFrame from data2oscr
oSCR.dat <- data2oscr(edf = edf,
                      sess.col = 1,
                      id.col   = 2,
                      occ.col  = 3,
                      trap.col = 4,
                      sex.col  = 5,
                      tdf = list(tdf),
                      K = K,
                      ntraps = ntraps,
                      rsfDF = list(rsfDF),
                      telemetry = telemetry
                      )

sftel <- oSCR.dat$scrFrame

# fit the SCR model with telemetry integration; "SCR+RSF" in Royle et al. 2013)
fit <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=TRUE,telemetry="ind",
                model=list(D~1,p0~z,sigma~1,path~1))


}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
