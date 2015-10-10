\name{oSCR.fit}
\alias{oSCR.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
oSCR.fit%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oSCR.fit(scrFrame, model = list(D ~ 1, p0 ~ 1, a1 ~ 1, path ~ 1), ssDF = NULL, costDF = NULL, distmet = c("euc", "user", "ecol")[1], sexmod = c("constant", "session")[1], DorN = c("D", "N")[1], directions = 8, Dmat = NULL, trimS = NULL, start.vals = NULL, PROJ = NULL, pxArea = 1, plotit = F, mycex = 0.5, tester = F, pl = 0, nlmgradtol = 1e-06, nlmstepmax = 10, predict = FALSE, smallslow = FALSE, multicatch = FALSE, hessian = T)
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
  \item{start.vals}{
%%     ~~Describe \code{start.vals} here~~
}
  \item{PROJ}{
%%     ~~Describe \code{PROJ} here~~
}
  \item{pxArea}{
%%     ~~Describe \code{pxArea} here~~
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
  \item{predict}{
%%     ~~Describe \code{predict} here~~
}
  \item{smallslow}{
%%     ~~Describe \code{smallslow} here~~
}
  \item{multicatch}{
%%     ~~Describe \code{multicatch} here~~
}
  \item{hessian}{
%%     ~~Describe \code{hessian} here~~
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

# Here is a test script
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
library(scrbook)
# Code comes from simMnSCR in scrbook package
# hard-wired 5 x 5 grid with ssbuff = 2

outfile<- NULL
parms<-list(N=100, alpha0= -.40, sigma=0.5,alpha2=0)
ssbuff<- 2

truth<- "multicatch"
with.sex<- TRUE

for(dataset in 1:100){


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
 # Careful here because I haev a seed hard-wired into the code to produce the book data set
data<- simSCR0(N=parms$N, K=7, alpha0 = parms$alpha0, sigma = parms$sigma, array3d=TRUE, rnd=dataset*100 )
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

if(truth=="multicatch"){
out1 <- oSCR.fit(scrFrame,model=list(D~1,p0~1,a1~1),ssDF=ssDF,start.vals=c(-1, log(.5), 2) ,
plotit=FALSE, multicatch=TRUE )
exp(out1$outStats[3,3])*nrow(ssDF[[1]])
 ests1<- out1$outStats[,3]
}else{
 ests1<- rep(NA, ifelse(with.sex,4,3))
}

out1b <- oSCR.fit(scrFrame,model=list(D~1,p0~1,a1~1),ssDF=ssDF,start.vals=c(-1,log(.5), 2) ,
plotit=FALSE, multicatch=FALSE )
exp(out1b$outStats[3,3])*nrow(ssDF[[1]])
ests2<- out1b$outStats[,3]


outfile<- rbind(outfile, c(ests1,ests2))

}

if(with.sex)
colnames(outfile)<-c("logitp0","loga1","logD","lpsi","logitp0","loga1","logD","lpsi")
if(!with.sex)
colnames(outfile)<-c("logitp0","loga1","logD","logitp0","loga1","logD")


 binomial data, 100 simulated data sets, means
(tmp<- apply(outfile,2,mean) )
 
 
(mn<- round(
 c(
plogis(tmp[1]), sqrt(1/(2*exp(tmp[2]))), "N"=exp(tmp[3])*nrow(ssDF[[1]]), plogis(tmp[4]), 
plogis(tmp[5]), sqrt(1/(2*exp(tmp[6]))), "N"=exp(tmp[7])*nrow(ssDF[[1]]), plogis(tmp[8]) ), 3)   )




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
