\name{oSCR.parfit}
\alias{oSCR.parfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallel (multi-core) model fitting for oSCR %%  ~~function to do ... ~~
}
\description{
This is a wrapper that executes oSCR.fit for multiple cores using
doParallel. %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oSCR.parfit(mods, data, ncores = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mods}{
A table of models to fit, see make.mods() %%     ~~Describe \code{mods} here~~
}
  \item{data}{
A list with the required oSCR.fit objects: data$sf = scrFrame, data$ss =
ssDF, data$cs = costDF

For example, the function sim.SCR.ms produces the correct objects (see below)
%%     ~~Describe \code{data} here~~
}
\item{ncores}{
  How many cores should be used to fit the models?
%%     ~~Describe \code{ncores} here~~
}
}
\details{
Uses the doParallel package..... %%  ~~ If necessary, more details than the description above ~~
}
\value{
  an oSCR fitList object 
  
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
make.mods %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


mods <- make.mods(density= c(~1,~session),
                  detection = c(~1,~session),
                  sigma = c(~1,~session),
                  cost = c(~1))
 
data<- sim.SCR.ms(sessions=2, sex=TRUE, sex.ratio=0.5, N = 100, K = 3,
                       alpha0 = -1.5, sigma = 0.5, discard0 = TRUE,
                       array3d = FALSE, ssRes = 0.5) 
 
# Note: In practice, allowing for model-specific starting values will probably
# be necessary to avoid crashes. Not done yet.

fmL <- oSCR.parfit(mods=mods, scrFrame=data$sf, ssDF=data$ssDF, ncores=4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
