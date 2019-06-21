\name{SCRdesign}
\alias{SCRdesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes quasi-optimal (aka "good")  spatial capture-recapture designs %%  ~~function to do ... ~~
}
\description{
  This function finds the optimal design for a given state-space (grid
  of points) and potential trap locations (also a grid of points) by
  minimizing an objective function that is related to expected
  information provided by a design. For crit = 1, the
  designs maximize the expected encounter probability
  of an individual given the state-space and hence maximizes the
  expected sample size of individuals. For crit = 2 the p2bar criterion
  is maximized (probability of encounter in 2 or more traps). 
   %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SCRdesign(S = S, all.traps = C, clusters = NULL, fix = NULL, clust.mids = clust.mids, ntraps = 9, ndesigns = 10, nn = 19, beta0 = -0.6, sigma = 2, crit = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
    State-space object, a G x 2 matrix.
%%     ~~Describe \code{S} here~~
}
\item{all.traps}{
   Candidate trap locations, n x 2 matrix.
%%     ~~Describe \code{C} here~~
}
\item{clusters}{
   Logical. Preserve specified cluster structure of traps.
%%     ~~Describe \code{clusters} here~~
}
\item{fix}{
  A matrix of trap coordinates
  to include in the design. (these count
  toward the total number of traps)
%%     ~~Describe \code{fix} here~~p
}
 
\item{clust.mids}{
Definition of clusters. This is an ntraps x 3 matrix where each row
corresponds to the trap location and cluster ID of the trap.
 i.e., clust.mids[i,1:2] are the coordinates and clust.mids[i,3]
is the cluster ID.
If this is provided then the internal function
¡®make.clust.mids(C,S)¡¯ uses the mid point of the cluster to determine
exchange updates of the cluster. 
}
\item{ntraps}{
   The number of traps in the design.
}
\item{ndesigns}{
  The algorithm uses a swapping algorithm based on a random initial
  design. ndesigns is the number of random starts to use as initial
  designs. The algorithm does one optimization for each initial design.
%%     ~~Describe \code{ndesigns} here~~
}
  \item{nn}{
    The number of nearest-neighbors to consider in the exchange algorithm.
  }
  \item{beta0}{
     Intercept in the encounter model (log-encounter rate scale). 
%%     ~~Describe \code{beta0} here~~
}
\item{sigma}{
   Half normal detection model scale parameter.
%%     ~~Describe \code{sigma} here~~
}
\item{crit}{
   Integer 1-2 specifying which  criterion to optimize. criterion 1
   maximizes pbar (by minimizing 1-pbar). Criterion 2 maximizes p2bar. 
%%     ~~Describe \code{crit} here~~
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
  Royle et al. (2014), Chapter 10;
  Sutherland et al. (2018)
  }
\author{
Chris Sutherland and Andy Royle 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

 

# Define the state-space as 13 by 13 units
xlim <- c(0, 13)
ylim <- c(0, 13)
area <- diff(xlim)*diff(ylim)


rr <- 0.5
# Make the state space
myss <- expand.grid(X = seq(xlim[1]+rr/2, xlim[2]-rr/2,rr),
                        Y = seq(ylim[1]+rr/2, ylim[2]-rr/2,rr))
myss<- as.matrix(myss)

rr<- 0.33
xlim<-c(2,11)
ylim<-c(2,11)
mytraps<- as.matrix( expand.grid(X = seq(xlim[1]+rr/2, xlim[2]-rr/2,rr),
                        Y = seq(ylim[1]+rr/2, ylim[2]-rr/2,rr))   )

#CHANGE Qfn TO p2bar (temporary inconvenience here)

Qfn<- p2bar.fn

 tmp1<- SCRdesign(statespace = myss, 
                      all.traps = mytraps, 
                      ntraps = 36, nn = 23, ndesigns = 2,
                      sigma = 0.80, beta0 = log(0.6), crit = 1)

 tmp2<- SCRdesign(statespace = myss, 
                      all.traps = mytraps, 
                      ntraps = 36, nn = 23, ndesigns = 2,
                      sigma = 0.80, beta0 = log(0.6), crit = 2)

plot(myss)
points(mytraps,pch=20)
points(tmp1$Xlst[[1]],pch=20,col="red",cex=2)
points(tmp2$Xlst[[1]],pch=20,col="green",cex=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
