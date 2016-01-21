\name{fitList.oSCR}
\alias{fitList.oSCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a fitList for Convenient Summarization 
}
\description{
Organizes multiple fitted model objects into a list that can be used for summarizing multiple models, model selection or for model-averaged prediction. The latter is not yet implemented
}
\usage{
fitList.oSCR(x, rename = F, drop = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A list of fitted models. List elements can be named.
}
  \item{rename}{
An option that, if \code{TRUE}, will generate model-specific names using details from the fitted model object. Default is \code{FALSE} which use the existing names of the list. If the list un-named, models are numbered in the order they appear in the list \code{x}. If \code{rename = TRUE}, the model names are extracted from the model object and are displayed in the following order: Density, \code{D()}, baseline detection, \code{p(), sigma}, \code{sig()}, and assymetric space use model, \code{asu()}. The dot notation "." is the null model.    
}
  \item{drop}{
An option to drop the name of a particular model component. For example, \code{drop=1} will remove the \code{D()} part of the name, \code{drop=2} will remove the \code{p()} part of the name, etc. \strong{Not yet operational}.   
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
%% ~put references to the literature/web site here ~
}
\author{
Chris Sutherland \email{csutherland@umass.edu}
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
rbs.sf <- subFrame(rbs$scrFrame,c(1,2))
rbs.ss <- list(rbs$ssDF[[1]],rbs$ssDF[[2]])

rbs.0 <- oSCR.fit(rbs.sf, list(D~1, p0~1, a1~1), rbs.ss, trimS=2)
rbs.Dsite <- oSCR.fit(rbs.sf, list(D~session, p0~1, a1~1), rbs.ss, trimS=2)
rbs.Pday <- oSCR.fit(rbs.sf, list(D~session, p0~day, a1~1), rbs.ss, trimS=2)
rbs.Pday2 <- oSCR.fit(rbs.sf, list(D~session, p0~day+day2, a1~1), rbs.ss, trimS=2)

fl <- fitList.oSCR(list(rbs.0,rbs.Dsite,rbs.Pday,rbs.Pday),rename=T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
