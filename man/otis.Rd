\name{otis}
\alias{otis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit some basic capture-recapture models based on Otis et al. (1978) %%  ~~function to do ... ~~
}
\description{
This functions fits models M0, Mt, Mh, and Mb and variations of those
models. 
}
\usage{
otis(y2d,mhst=NULL, cov=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y2d}{ nind x K (occasions) matrix of encounter histories.  A
    standard closed population capture-recapture data set. 
}
\item{mhst}{ Starting values for the logit-normal model Mh in order:
(logit(p0), log(mu), log(n0)) where n0 = N - n.
If no starting values are provided, the function does a grid search to find
a decent starting value (this can take some time). 
}
\item{cov}{  This is a K x 1 vector of some covariate that might explain
  variation in detection probability over time. e.g., a linear trend
  covariate or weather conditions.  If supplied then models which
  include this covariate will be fitted. 
  }
}
\details{
  Default setting (cov=NULL) will fit 4 models: M0, Mb, Mh
  (logit-normal) and Mh (2-point finite mixture).

  If cov= is specified then Mcov and Mb+cov are fitted

  Soon we will fit other models too. 
   
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an AIC table including other model summaries (Nhat) and a table of model coefficients.
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

# see also ?peromyscus 
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
nsims<- 1
for(sim in 1:nsims){
y<- matrix(0,nrow=N,ncol=K)
for(i in 1:N){
  p<- plogis(rnorm(1,-1.85, 0.5))
  y[i,]<- rbinom(K, 1, p)
  }

y<- y[apply(y,1,sum)>0,]
otis(y, cov=c(1,1,1,1,0,0,0,0))
}





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
