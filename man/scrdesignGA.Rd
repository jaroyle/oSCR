\name{scrdesignGA}
\alias{scrdesignGA}
%- Also NEED an '\alias' for EACH other topic documented here. % CHIRS: I'm not sure what to do here (gld)
\title{
Generating near-optimal sampling designs for spatial capture-recapture.
}
\description{
Given any SCR-sensible objective function, \code{scrdesignGA()} implements a genetic algorithm to solve this \emph{k-of-n} design problem. In the context of SCR, and for a given objective function and landscape, the function finds the best subset of sampling locations \emph{(k)} from all possible sampling locations \emph{(n)} . The arguments of this function that are related to the biology of the species should be parameterized using the values that are expected, preferrably informed by previous research.
}
\usage{
scrdesignGA(statespace = NULL, alltraps = NULL, ntraps = 9, 
            beta0 = -0.2*5, sigma = 2, N = 100, crit = 1, 
            verbose = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{statespace}{
A discretized representation of the study area (state-space), defined as a convex hull around the trapping extent with an additional buffer of approximately \code{3*sigma} to ensure the incorporation of the activity centers of all individuals that could be captured.
}
  \item{alltraps}{
A two-column dataframe, labelled X and Y accordingly, comprised of all possible trapping locations. These  locations should be user-defined after accounting for logistic constraints, primarily accessibility. Design   generation then relies on subsetting these locations.
}
  \item{ntraps}{
The number of traps available for constructing the design. More traps will allow for a larger sample size, providing more data to an SCR model, resulting in improvement in terms of bias and       precision of estimates.
}
  \item{beta0}{
The probability of encounter of an individual across occassions at its activtiy center, specified here as the logaritm of baseline detection probability (\code{p0}) multiplied by the number of sampling occasions, i.e., \code{log(p0*K)}.
}
  \item{sigma}{
The spatial scaling parameter in the half-normal encounter model, which relates the baseline encounter probability, \code{p0}, to the distance between the individual's activity center and the location of the trap where it was captured.
}
  \item{N}{
The number of individuals in the landscape, which is used to determine population density in the statespace.
}
  \item{crit}{
Integer. The design criteria used in the objective function for optimization. Available options are model-based criteria, which optimize:
    \itemize{
      \item{\code{crit = 1 }} The probability of capture.
      \item{\code{crit = 2 }} The probability of capture at more than one trap.
      \item{\code{crit = 3 }} An intermediate between criteria 1 and 2.
      \item{\code{crit = 4 }} The probability of more than one capture.
    }
Additionally, user-defined criteria can also be implement, see \code{oSCR::scrdesignOF} for more details.
}
  \item{\dots}{
Below are the descriptions used in the help file for \code{\link[kofnGA:kofnGA]{kofnGA::kofnGA()}} \emph{(credit to the authors of \code{\link[kofnGA]{kofnGA}})}
}
  \item{popsize}{
The size of the population; equivalently, the number of offspring produced each generation.
}
  \item{keepbest}{
The \code{keepbest} least fit offspring each generation are replaced by the \code{keepbest} most fit members of the previous generation. Used to implement elitism.
}
  \item{ngen}{
The number of generations to run.
}
  \item{cluster}{
If non-null, the objective function evaluations for each generation are done in
parallel. \code{cluster} can be either a cluster as produced by \code{\link[parallel]{makeCluster}},
or an integer number of parallel workers to use. If an integer, \code{makeCluster(cluster)} will
be called to create a cluster, which will be stopped on function exit.
}
  \item{verbose}{
An integer controlling the display of progress during search.  If \code{verbose}
takes positive value \code{v}, then the iteration number and best objective function value are
displayed at the console every \code{v} generations. Otherwise nothing is displayed. Default is
zero (no display).
}

}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of class "scrdesign" with the following elements:

\item{des}{
A list containing the resulting components of the design-finding process. For more information, see \code{\link[kofnGA:kofnGA]{kofnGA::kofnGA()}}.
}
\item{statespace}{
A 3-column dataframe containing the coordinates (X and Y) and the expected per-pixel density of the discretized state-space.
}
\item{alltraps}{
A 2-column dataframe containing the coordinates (X and Y) of all possible locations available for traps.
}
\item{optimaltraps}{
Most importantly, the optimal design produced by the algorithm. A 2-column dataframe containing the selected coordinates (X and Y) of each trap.
}
\item{sigma}{
The value of sigma, the spatial scaling parameter in the half-normal encounter model, originally specified by the user.
}
\item{beta0}{
The value of the encounter probability across sampling occasions, as originally specified by the user.
}

}
\references{
TBD 😄
}
\author{
%% For Chris to add
}
\note{
For best designs, \code{ngen} should be ~\code{1500}+. Further, the \code{kofnGA} package recommends \code{popsize = 200}, and \code{keepbest = 20} for best performance.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For more details on the particulars of the \code{kofnGA} genetic algorithm, see the help file for: \code{\link[kofnGA:kofnGA]{kofnGA::kofnGA()}}. The help file for \code{kofnGA} also provides information on specifying the parallel evaluation of that function, which can also be used here for \code{scrdesignGA}.

}
\examples{

#----Load the package----

library(oSCR)


#----Load data from Pakistan----

data(pakistan) 
# This loads the statespace: pakSS 
# and the possible trap locations: pakTT


#----Plot this data----

plot(pakSS, asp=1, pch=16, col="grey")
points(pakTT, pch=20, cex = 0.5)


#----Run the design-finding algo----

test_design <- scrdesignGA(
  statespace = pakSS, alltraps = pakTT,  # Study area components
  ntraps = 25,                           # Number of available traps
  beta0 = 0.2*5, sigma = 3, crit = 1,    # Expected data
  popsize = 50, keepbest= 5, ngen = 50   # GA settings 
  ) # FOR BEST DESIGNS: ngen should be ~1500


#----Plot the results----

par(mfrow=c(1,3))          # Setup plotting area
plot(test_design, which=4) # which = 4 plots all three together
par(mfrow=c(1,1))          # Reset plotting area

}

\keyword{ ~methods }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~design }% __ONLY ONE__ keyword per line
