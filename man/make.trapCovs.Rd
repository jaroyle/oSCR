\name{make.trapCovs}
\alias{make.trapCovs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make.trapCovs%%  ~~function to do ... ~~
}
\description{
This function takes a sensible "list of lists" format for covariates and
converts it to the bizarre and unwiedly oSCR format.
This will be internalized at some point but for now we suggest to use
this function for convenience. 
}
\usage{
make.trapCovs(lol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lol}{
A list of lists.
}
}
\details{
  lol should be a list of lists where lol[[i]] is itself a LIST
  representing covariate i for each session. Thus lol[[i]][[s]] should
  be a J (traps) x K (occasions) matrix of covariate i for session s
  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  This is indescribable. But it can be input as the trapCovs argument to
  the scrFrame.  
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
Andy Royle %%  ~~who you are~~
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

# We have an important covariate which is time of trap check, AM or PM
#   we put this into a matrix and treat it as a "trap*time covariate"
am<- matrix( o$am.effect, nrow=J,ncol=K,byrow=TRUE)

# To demo here is an example with 2 covariates:
weird<-  list(am=list(am),am2=list(am))
weird<- make.trapCovs(weird)

# "to" is the trap operation matrix
to<- matrix(1,nrow=121,ncol=K)
colnames(traps)<- c("X","Y")

# load peromyscus data (not yet available)
# o$am.effect is a vector of length 10 , binary indicator

am<- matrix( o$am.effect, nrow=J,ncol=K,byrow=FALSE)
# Here is an example with two covariates although not meaningful ones
lol <-  list(am=list(am),am2=list(am*am))
lol <-  make.trapCovs(lol)

scrFrame <- make.scrFrame(  
  caphist=list(y3d), 
  indCovs=NULL,
  traps=list(data.frame(traps)) ,
  trapCovs=lol, 
  trapOperation=list(to))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
