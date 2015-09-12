\name{oSCR}
\alias{oSCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oSCR(scrFrame, model = list(D ~ 1, p0 ~ 1, a1 ~ 1, path ~ 1), ssDF = NULL, costDF = NULL, distmet = c("euc", "user", "ecol")[1], sexmod = c("constant", "session")[1], DorN = c("D", "N")[1], directions = 8, Dmat = NULL, trimS = NULL, res = NULL, smallslow = FALSE, predict = FALSE, start.vals = NULL, PROJ = NULL, plotit = F, mycex = 0.5, tester = F, pl = 0, nlmgradtol = 1e-06, nlmstepmax = 10, multicatch = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
%%     ~~Describe \code{scrFrame} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{ssDF}{
%%     ~~Describe \code{ssDF} here~~
}
  \item{costDF}{
%%     ~~Describe \code{costDF} here~~
}
  \item{distmet}{
%%     ~~Describe \code{distmet} here~~
}
  \item{sexmod}{
%%     ~~Describe \code{sexmod} here~~
}
  \item{DorN}{
%%     ~~Describe \code{DorN} here~~
}
  \item{directions}{
%%     ~~Describe \code{directions} here~~
}
  \item{Dmat}{
%%     ~~Describe \code{Dmat} here~~
}
  \item{trimS}{
%%     ~~Describe \code{trimS} here~~
}
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{smallslow}{
%%     ~~Describe \code{smallslow} here~~
}
  \item{predict}{
%%     ~~Describe \code{predict} here~~
}
  \item{start.vals}{
%%     ~~Describe \code{start.vals} here~~
}
  \item{PROJ}{
%%     ~~Describe \code{PROJ} here~~
}
  \item{plotit}{
%%     ~~Describe \code{plotit} here~~
}
  \item{mycex}{
%%     ~~Describe \code{mycex} here~~
}
  \item{tester}{
%%     ~~Describe \code{tester} here~~
}
  \item{pl}{
%%     ~~Describe \code{pl} here~~
}
  \item{nlmgradtol}{
%%     ~~Describe \code{nlmgradtol} here~~
}
  \item{nlmstepmax}{
%%     ~~Describe \code{nlmstepmax} here~~
}
  \item{multicatch}{
%%     ~~Describe \code{multicatch} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


library(scrbook)
library(gdistance)
out<-make.EDcovariates()
covariate<-out$covariate.patchy
set.seed(2013)

N<-200
alpha0<- -2
sigma<- .5
K<- 5

alpha1<- 1/(2*sigma*sigma)
r<-raster(nrows=20,ncols=20)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84"
extent(r)<-c(.5,4.5,.5,4.5)
alpha2<-1
cost<- exp(alpha2*covariate)

tr1<-transition(cost,transitionFunction=function(x) 1/mean(x),directions=8)
tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)

# make up some trap locations
xg<-seq(1,4,1)
yg<-4:1
pts<-cbind( sort(rep(xg,4)),rep(yg,4))

traplocs<-pts
points(traplocs,pch=20,col="red")
ntraps<-nrow(traplocs)

S<-cbind(runif(N,.5,4.5),runif(N,.5,4.5))
D<-costDistance(tr1CorrC,S,traplocs)
probcap<-plogis(alpha0)*exp(-alpha1*D*D)
# now generate the encounters of every individual in every trap
Y<-matrix(NA,nrow=N,ncol=ntraps)
Y3d<- array(NA,dim=c(N,ntraps,K))
for(i in 1:nrow(Y)){
for(k in 1:K){
  ##Y[i,]<-rbinom(ntraps,K,probcap[i,])
 Y3d[i,1:ntraps,k]<- rbinom(ntraps,1,probcap[i,])
}
}
Y2d<- apply(Y3d,c(1,2),sum)
captured<- apply(Y2d,1,sum)>0
Y3d<-Y3d[captured,,]
Y<-Y2d[captured,]

##
## integrated likelihood for ecological distance from the SCR book
##
n0<- N-nrow(Y)
frog1<-nlm(intlik3ed,c(alpha0,alpha1,log(n0)),hessian=TRUE,y=Y,K=K,X=traplocs,
               distmet="euclid",covariate=covariate,alpha2=1)

frog2<-nlm(intlik3ed,c(alpha0,alpha1,log(n0),-.3),hessian=TRUE,y=Y,K=K,
               X=traplocs,distmet="ecol",covariate=covariate,alpha2=NA)

#
#
# Now run oSCR for comparison
#
#
traplocs<- data.frame(X=traplocs[,1],Y=traplocs[,2])

#this is from intlik we need to construct the state-space grid to be exactly the same
# or else estimates differ a little bit
#
  nc <- covariate@ncols
  nr <- covariate@nrows
  Xl <- covariate@extent@xmin
  Xu <- covariate@extent@xmax
  Yl <- covariate@extent@ymin
  Yu <- covariate@extent@ymax
  delta <- (Xu - Xl)/nc
  xg <- seq(Xl + delta/2, Xu - delta/2, delta)
  yg <- seq(Yl + delta/2, Yu - delta/2, delta)
  npix.x <- length(xg)
  npix.y <- length(yg)
  G <- cbind(rep(xg, npix.y), sort(rep(yg, npix.x)))
 
ssDFt<-list()
ssDFt[[1]] <- data.frame(x=G[,1], y=G[,2])
sex<- rbinom(nrow(Y),1,.5)

scrFrame <- list(caphist=list(Y3d), indCovs=list(sex=sex), type="scr", traps=list(traplocs), 
trapCovs=list(NULL), nocc=K)
  class(scrFrame) <- "scrFrame"

library(Formula)
test<-  oSCR(scrFrame, model=list(D~1, p0~1, a1~1, path~1), ssDF=ssDFt, costDF=NULL,
             distmet=c("euc","user","ecol")[1], directions = 16, Dmat=NULL,
             start.vals=NULL)

test$outStats
frog1$estimate

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
