\name{peromyscus}
\alias{peromyscus}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("peromyscus")}
\format{
  The format is:
List of 22
 $ Ycat       : num [1:1082, 1:10] 122 122 9 23 97 111 115 122 122 122 ...
 $ group.mem  : num [1:1082] 1 1 1 1 1 1 1 1 1 1 ...
 $ M          : int 1082
 $ n.sites    : num 8
 $ n.seasons  : num 3
 $ n.groups   : num 24
 $ n.traps    : num 121
 $ trap.locs  : int [1:121, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
 $ Yl         : num -1
 $ Yu         : num 13
 $ Xl         : num -1
 $ Xu         : num 13
 $ site       : int [1:24] 1 2 3 4 5 6 7 8 1 2 ...
 $ season     : num [1:24] 1 1 1 1 1 1 1 1 2 2 ...
 $ fire       : num [1:24] 0 0 0 0 0 0 0 0 0 0 ...
 $ thin       : num [1:24] 0 0 0 0 0 0 0 0 0 0 ...
 $ reencounter: num [1:1082, 1:10] 0 0 0 0 0 0 0 0 0 0 ...
 $ last.cap   : num [1:1082] 10 10 10 10 10 10 10 10 10 10 ...
 $ trap.space : num 25
 $ traps.avail: num [1:121, 1:24] 1 0 1 0 1 0 1 0 1 0 ...
 $ am.effect  : num [1:10] 1 0 1 0 1 0 1 0 1 0
 $ twotraps   : num [1:36] 1 3 5 7 9 11 23 25 27 29 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Converse, S. J., G. C. White, and W. M. Block. 2006. Small mammal responses to thinning and wildfire in ponderosa pine-dominated forsts of the southwestern United States. Journal of Wildlife Management 70:1711-1722. 
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Royle, J. A., and S. J. Converse. 2014. Hierarchical spatial capture-recapture models: modleing population density fr 
}
\examples{
data(peromyscus)

# This is a bit rough yet....

# grab the peromyscus data and process a few things out of it
# here we take group = 19 which has the most individuals (66) and
# a total of 185 total captures
data(peromyscus)
o<- peromyscus
y<- o$Ycat  # Data are categorical "trap of capture"
g<- o$group.mem
table(g)

# here we subset the data to trapping grid 19
y<- y[g==19,]
allna<- apply(is.na(y),1,sum)  # Gets rid of individuals not captured here
y<- y[!allna,]
table(o$traps.avail[,19])
traps<- o$trap.locs
K<- ncol(y)
J<- nrow(traps)

# here we create the 3-d array format which is conveneient for a lot
# of analysis in BUGS but also oSCR requires the  3-d array format
y3d<- array(0, dim=c(nrow(y), nrow(traps), ncol(y)) )
for(i in 1:nrow(y)){
  tmp<- y[i,]
  kid<- (1:K)[tmp<=J]
  jid<- tmp[tmp<=J]
  y3d[cbind(i,jid,kid)]<- 1
}

# load oSCR to do some analysis 
library(oSCR)

# We have an important covariate which is time of trap check, AM or PM
#   we put this into a matrix and treat it as a "trap*time covariate"
am<- matrix( o$am.effect, nrow=J,ncol=K,byrow=TRUE)

# To demo here is an example with 2 covariates:
weird<-  list(am=list(am),am2=list(am))
weird<- make.trapCovs(weird)

# "to" is the trap operation matrix
to<- matrix(1,nrow=121,ncol=K)
colnames(traps)<- c("X","Y")

# Create the scrFrame data structure which has basic information in it
scrFrame <- list(  
  caphist=list(y3d), 
  indCovs=NULL,
  traps=list(data.frame(traps)) ,
  trapCovs=weird,
  trapOperation=list(to),
  type="scr",
  nocc = K)
class(scrFrame)<- "scrFrame"

# Make the state-space data frame
ssDF.new <- make.ssDF(scrFrame, res=.5, buffer=2)

# plot this to see what it looks like
plot(ssDF.new[[1]],pch=20,cex=0.5)
points(traps,pch=3,col="red",cex=2,lwd=2)

# fit a basic SCR model
temp1 <-  oSCR.fit(scrFrame,  model=list(D~1,p0~1, sig~1),ssDF=ssDF.new, plotit=FALSE, multicatch=TRUE,trimS=3)

spacing<-  0.5*0.5*( (o$trap.space)^2 )   # square meters per pixel
p0.hat<- plogis(temp1$outStats[1,3])
sigma.hat<- exp(temp1$outStats[2,3])
#sigma.hat<- sqrt(1/(2*exp(temp1$outStats[2,3])))
d0<- exp(temp1$outStats[3,3])
Nhat<- d0*nrow(ssDF.new[[1]])
Dhat<- (d0/spacing)*10000   # now per ha

# fit SCRb
temp2 <-  oSCR.fit(scrFrame, model=list(D~1,p0~b,sig~1),ssDF=ssDF.new, plotit=FALSE, multicatch=TRUE, start =c(-3,2,-1,-2), trimS=3)

p0.hat<- plogis(temp2$outStats[1,3])
#sigma.hat<- sqrt(1/(2*exp(temp2$outStats[2,3])))
sigma.hat<- exp(temp2$outStats[2,3])
d0<- exp(temp2$outStats[3,3])
Nhat<- d0*nrow(ssDF.new[[1]])
Dhat<- (d0/spacing)*10000


# fit AM effect
 
temp3 <-  oSCR.fit(scrFrame, model=list(D~1,p0~ b + am,sig~1),ssDF=ssDF.new, plotit=FALSE, start=c(-3, 2 ,-1, 0, -2),multicatch=TRUE,trimS=3)


p0.hat<- plogis(temp3$outStats[1,3])
#sigma.hat<- sqrt(1/(2*exp(temp3$outStats[2,3])))
sigma.hat<- exp(temp3$outStats[2,3])
d0<- exp(temp3$outStats[3,3])
Nhat<- d0*nrow(ssDF.new[[1]])
Dhat<- (d0/spacing)*10000

fl<- fitList.oSCR(list(
"p(.)sig(.)" = temp1,
"p(b)sig(.)" = temp2,
"p(b+am)sig(.)" = temp3) )
scr<- modSel.oSCR(fl)
scr<- scr$coef.tab
scr[,2]<- sqrt(1/(2*exp(scr[,2])))
scr[,3]<- (exp(scr[,3])/(spacing))*10000  # per ha

# need to add "am" effect and also "b" to these functions
y<- flatten(y3d)

# Compute the standard capture-recapture models. 
(tab<- otis(y, cov=am[1,]))


## maybe str(peromyscus) ; plot(peromyscus) ...
}
\keyword{datasets}
