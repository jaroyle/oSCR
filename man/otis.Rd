\name{otis}
\alias{otis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit basic capture-recapture models from Otis et al. (1978) %%  ~~function to do ... ~~
}
\description{
This functions fits model M0, model Mh, .... [more later]%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
otis(y2d,mhst=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y2d}{ nind x K (occasions) matrix of encounter histories. 
%%     ~~Describe \code{y2d} here~~
}
\item{mhst}{
(logit(p0), log(mu), log(n0)) where n0 = N - n.
If no starting values are provided, the function does a grid search to find
a decent starting value (this can take some time). 
}
}
\details{
Fits model M0 and two versions of model Mh: the logit-normal model and a
two-point finite mixture model.

Will add model Mb and model Mt and possibly others. 
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
Otis et al. (1978) 
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

library(scrbook)
data(beardata)
trapmat<-beardata$trapmat
nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]
bearArray<-beardata$bearArray

# Convert 3d array to 2d array
y<- flatten(bearArray)

# Fit some basic CR models 
otis(y)


# A simulation study:


N<- 100
K<- 8
nsims<- 20
for(sim in 1:nsims){
y<- matrix(y,nrow=N,ncol=K)
for(i in 1:N){
  p<- plogis(rnorm(1,0,-1.5))
  y[i,]<- rbinom(K, 1, p)
  }

otis(y)
}





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
