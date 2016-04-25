\name{oSCR.fit}
\alias{oSCR.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
oSCR.fit%%  ~~function to do ... ~~
}
\description{
  This does likelihood analysis of sex-structured multi-session SCR
  models using either Euclidean or non-Euclidean (least-cost path)
  models.
}
\usage{
oSCR.fit(model = list(D ~ 1, p0 ~ 1, sig ~ 1, path ~ 1), scrFrame,  ssDF = NULL, costDF = NULL, 
         distmet = c("euc", "user", "ecol")[1], sexmod = c("constant", "session")[1], 
         DorN = c("D", "N")[1], directions = 8, Dmat = NULL, trimS = NULL, start.vals = NULL, 
         PROJ = NULL, pxArea = 1, plotit = F, mycex = 0.5, tester = F, pl = 0, nlmgradtol = 1e-06, 
         nlmstepmax = 10, predict = FALSE, smallslow = FALSE, multicatch = FALSE, hessian = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 This is a list with 4 components describing the SCR model. The order of
  the components is: density model, detection model, sigma model and
  asymmetric space use model. For example, the null model is wiould be written as \code{list(D ~ 1,
  p0 ~ 1, sig ~ 1,  path ~ 1)}.

  The density model can include "session" as an effect and then a
  session-specific density is fitted.  Otherwise variables are "within
  session" variables (i.e., covariates).

  p0 can include "sex" "b" for local (trap level) behavioral response,
    or b*sex or include covariates. 

  sig is the log(sigma) in the model p = exp(- (1/2*sigma*sigma)*dist(x,s)^2)
}
  \item{scrFrame}{
  An object of class scrFrame.
  }
  \item{ssDF}{

A list of data frames, one for each of the R sessions, containing the
coordinates of the state-space pixels used to approximate the landscape of
interest. The data frame must have the headings X and Y. These data frames may
also have named columns that represent cell-specific (i.e., spacially varying)
covariate values.

 Although the fitting function in oSCR can generate a state-space
  object ?on-the-fly?, we recommend that users specify the size and resolution
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
    Area of state-space pixels. Default value of pxArea=1 returns "Area"
    = number of pixels in the state-space.
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

data(rbs)
caphist <- rbs$scrFrame$caphist
traps <- rbs$scrFrame$traps
trapCovs <- rbs$scrFrame$trapCovs
ss <- rbs$ssDF
starts <- c(-2.9, -0.051, -0.920)
sf <- make.scrFrame(caphist = caphist, traps = traps, trapCovs = trapCovs, type = "scr")

m0 <- oSCR.fit(model = list(~1,~1,~1,~1), scrFrame = sf, ssDF = ss, trimS = 4, start.vals = starts)
m0

}

