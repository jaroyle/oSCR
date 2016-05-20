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


# load peromyscus data (not yet available)
# o$am.effect is a vector of length 10 , binary indicator

am<- matrix( o$am.effect, nrow=J,ncol=K,byrow=FALSE)
# here is an example with two covariates
weird<-  list(am=list(am),am2=list(am))
weird<-  make.trapCovs(weird)

scrFrame <- list(  
  caphist=list(y3d), 
  indCovs=NULL,
  traps=list(data.frame(traps)) ,
  trapCovs=weird,
  trapOperation=list(to),
  type="scr",
  nocc = K)
class(scrFrame)<- "scrFrame"

# See the help file ?peromyscus for some models that are run using this scrFrame






}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
