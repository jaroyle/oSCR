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
otis(y2d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y2d}{ nind x K (occasions) matrix of encounter histories. 
%%     ~~Describe \code{y2d} here~~
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
