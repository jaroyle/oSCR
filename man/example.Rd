\name{example}
\alias{example}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
example (old oSCR.fit man page)
}
\description{
 I'm just saving this for a period of time because the  R script at the end is useful. 
}
\usage{
oSCR.fit(scrFrame, model = list(D ~ 1, p0 ~ 1, sig ~ 1, path ~ 1), ssDF = NULL, costDF = NULL, distmet = c("euc", "user", "ecol")[1], sexmod = c("constant", "session")[1], DorN = c("D", "N")[1], directions = 8, Dmat = NULL, trimS = NULL, start.vals = NULL, PROJ = NULL, pxArea = 1, plotit = F, mycex = 0.5, tester = F, pl = 0, nlmgradtol = 1e-06, nlmstepmax = 10, predict = FALSE, smallslow = FALSE, multicatch = FALSE, hessian = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
  scrFrame <- list(caphist,
indCovs,
traps,
trapCovs,
trapOperation,
type,
nocc)

The object scrFrame is a list containing all of the named objects, which themseleves are lists containing a
data object for each of the R session in the analysis. For example, in the analysis in the main paper, we
considered four plots, i.e. R = 4 sessions.

 caphist: A list of 3-dimaensional n × J × K arrays, one for each or the R sessions, where n is the
number of individuals, J is the number of traps, and K is the number of sampling occasions. The array
is filled with 1’s or 0’s that denote whether individuals were captured (1) or not (0) in each trap in
each occasion.

 indCovs: A list of data frames, one for each of the R sessions, containing individual level data. Currently,
this data oject is only used for sex which can be male (1), female (0), or uknown (NA).

 traps: A list of data frames, one for each of the R sessions, containing the coordinates of the
traps/detectors. The data frame must have the headings X and Y. A data frame for each of R sessions
is

 trapCovs: A two level list object which is a list of lists of data frames. This list object is of length R,
and each of the R slots in the list is itself a list of data frames containing trap level covariates. The
objects requires a data frame for each sampling occasion which can be identical of differ by occasion,
which allows time varying trap level effects.

 trapOperation: A list of binary matrices, one for each of the R sessions, with dimension J × K. The
matrix denotes whether a trap was operational (1) or not (0), and if not, the assocaited cells in the the
encounter histories (i.e., in caphist) do not contiibute to the likelihood.

 type: This should be set to SCR always. This option will, in the future, allow the user to provide
data in the secr format.

 nocc: An integer denoting the largest number of occasions/visits for any of the R sessions.


IMPORTANT: While caphist, traps, type and nocc must be provided, indCovs, trapCovs and trapOperation can be
set to ‘NULL’.

}
  \item{model}{
 This is a list with 4 components describing the SCR model. The order of
  the components is: density model, detection model, "sigma" model and
  least-cost path model. For example, the null model is:   list(D ~ 1,
  p0 ~ 1, sig ~ 1,  path ~ 1)

  The density model can include "session" as an effect and then a
  session-specific density is fitted.  Otherwise variables are "within
  session" variables (i.e., covariates).

  p0 can include "sex" "b" for local (trap level) behavioral response,
    or b*sex or include covariates. 

  sig is the log(sigma) in the model p = exp(- (1/2*sigma*sigma)*dist(x,s)^2)
}
  \item{ssDF}{

A list of data frames, one for each of the R sessions, containing the
coordinates of the state-space pixels used to approximate the landscape of
interest. The data frame must have the headings X and Y. These data frames may
also have named columns that represent cell-specific (i.e., spacially varying)
covariate values.

 Although the fitting function in oSCR can generate a state-space
  object ‘on-the-fly’, we recommend that users specify the size and resolution
  of the stat-space in order avoid any errors, but which will also unsure a
  better understanding of the model structure. 

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
  This does likelihood analysis of sex-structured multi-session SCR
  models using either Euclidean or non-Euclidean (least-cost path)
  models.
  
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
  Borchers and Efford 2008;
  Royle et al. 2013;
  Royle et al. 2014; 
  Sutherland et al 2015;
  Fuller et al. 2015;
  
}
\author{
  Chris Sutherland
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

truth<- "multicatch"   # "binomial"
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
 # Careful here because I have a seed hard-wired into the code to produce the book data set
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
out1 <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1),ssDF=ssDF,
plotit=FALSE, multicatch=TRUE , trimS=3)
exp(out1$outStats[3,3])*nrow(ssDF[[1]])
 ests1<- out1$outStats[,3]
}else{
 ests1<- rep(NA, ifelse(with.sex,4,3))
}

out1b <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1),ssDF=ssDF,
plotit=FALSE, multicatch=FALSE, trimS=3)
exp(out1b$outStats[3,3])*nrow(ssDF[[1]])
ests2<- out1b$outStats[,3]


outfile<- rbind(outfile, c(ests1,ests2))

}

if(with.sex)
colnames(outfile)<-c("logitp0","logsig","logD","lpsi","logitp0","logsig","logD","lpsi")
if(!with.sex)
colnames(outfile)<-c("logitp0","logsig","logD","logitp0","logsig","logD")


# binomial data, 100 simulated data sets, means
(tmp<- apply(outfile,2,mean) )

# Summarize output 
if(with.sex){
(mn<- round( c(
plogis(tmp[1]), exp(tmp[2]),  "N"=exp(tmp[3])*nrow(ssDF[[1]]), plogis(tmp[4]), 
plogis(tmp[5]), exp(tmp[6]),  "N"=exp(tmp[7])*nrow(ssDF[[1]]),
 plogis(tmp[8]) ), 3) )
names(mn)<-c("p0","sigma","N","psi.sex","p0","sigma","N","psi.sex")
}
if(!with.sex){
(mn<- round( c(
plogis(tmp[1]), exp(tmp[2]),  "N"=exp(tmp[3])*nrow(ssDF[[1]]), 
plogis(tmp[4]), exp(tmp[5]),  "N"=exp(tmp[6])*nrow(ssDF[[1]])  ), 3) )
names(mn)<-c("p0","sigma","N","p0","sigma","N")

} 





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
