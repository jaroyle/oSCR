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
# Find optimal pbar and p2bar designs and simulate the performance of MLEs

# Set population size
N <- 50
D <- N/area

# Define parameters (should pass these as arguments)
 p0 <- 0.2                            #baseline encounter probability
 sigma <- 0.8                         #scale factor of detection function
 K <- 3                               # sampling occasions
 lam0<- K*p0
 # for the Poisson model this is approximately exp(beta0) = 0.2*3
 # Simulator function needs modified to allow for the Poisson model.
 # SCRdesign ONLY does the Poisson model.

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

#CHANGE TO p2bar
Qfn<- p2bar.fn

 tmp1<- SCRdesign(statespace = myss, 
                      all.traps = mytraps, 
                      ntraps = 36, nn = 23, ndesigns = 2,
                      sigma = 1.30, beta0 = log(lam0), crit = 1)

 tmp2<- SCRdesign(statespace = myss, 
                      all.traps = mytraps, 
                      ntraps = 36, nn = 23, ndesigns = 2,
                      sigma = 1.30, beta0 = log(lam0), crit = 2)

plot(myss)
points(mytraps,pch=20)
points(tmp1$Xlst[[1]],pch=20,col="red",cex=2)
points(tmp2$Xlst[[1]],pch=20,col="green",cex=2)




 
traps1<- tmp1$Xlst[[1]]
traps2<- tmp2$Xlst[[1]]




library(oSCR)
#Create a simulator function
simulator<- function(traps , nsim) {

  simout1<- matrix(NA,nrow=nsim,ncol=7) #create empty matrix for output
  colnames(simout1)<- c("p0","sig","d0","nind", "avg.caps","avg.spatial","mmdm")

  for(sim in 1:nsim){
    print(paste("Simulation Number", sim, sep = " ")) # keep track

    # Generate home range centers
    s <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1],ylim[2]) )
    rr <- 0.5

    # Make the state space
    myss <- expand.grid(X = seq(xlim[1]+rr/2, xlim[2]-rr/2,rr),
                        Y = seq(ylim[1]+rr/2, ylim[2]-rr/2,rr))
    cat("size of state-space: ", nrow(myss), " pixels", fill=TRUE)
    myss$Tr <- 1
    myss <- list(myss)
    class(myss) <- "ssDF"

    # individual-trap distance matrix
    D <- e2dist(s,traps)

    # Compute detection probabilities:
    pmat <- p0*exp(-D*D/(2*sigma*sigma)) #p for all inds-traps p_ij
    ntraps <- nrow(traps)
    y <- array(0, dim=c(N, ntraps, K)) # empty 3D array (inds by traps by occ)

    for(i in 1:N){# loop through each individual/activity center
      for(j in 1:ntraps){# loop through each trap
        y[i,j,1:K]<- rbinom(K, 1, pmat[i,j]) # y ~ binomial(p_ijk)

      }
    }

    ncap <- apply(y,c(1), sum)       # sum of captures for each individual
    y<- y[ncap>0,,]                  # reduce the y array to include only captured individuals

    # Some summary information, that is actually printed for you later with "print(scrFrame)"
    caps.per.ind.trap<- apply(y,c(1,2),sum) #shows # capts for each indv across all traps

    # Make the SCRframe
    colnames(traps)<- c("X","Y")
    sf <- make.scrFrame(caphist=list(y), traps=list(traps))

    plot(myss) # plot the state space
    spiderplot(y = sf$caphist[[1]], traps, add=TRUE)

    # Fit a basic model SCR0
    out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame = sf, ssDF=myss, trimS = 4*sigma)

    stats <- print(sf)[[1]]  # pulls avg caps, avg spatial caps, and mmdm
    est <- out1$outStats$mle       # pulls p0, sigma, and d0 estimates from the model
    simout1[sim,] <- c(plogis(est[1]), exp(est[2]), exp(est[3]), dim(y)[1], stats)

  }
  return(simout1)
}

# simulate pbar design
simout1<- simulator(traps1 ,nsim=50)   # runs simulation on the first design
# simulate p2bar design
simout2<- simulator(traps2 ,nsim=50)   # runs simulation on the first design

colMeans(simout1) #shows average values of the estimates
colMeans(simout2)
 
apply(simout1,2,sd) #calculates standard deviation of all the estimates
apply(simout2,2,sd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
