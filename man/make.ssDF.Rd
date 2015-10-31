\name{make.ssDF}
\alias{make.ssDF}
\title{
A function for creating a discretized state-space data object of class 'ssDF'.  
}
\description{
A function for creating a discretized state-space data object of class "ssDF". The object is a list of data frames, one dataframe for each session and each has columns "X", and "Y" which are pixel coordinates. The function generates session specific state-space objects which are created based on the session specific trap locations contained in the supplied 'scrFrame'. The buffer region and the resolution of the discrete surface are user specified. The function also allows random spatially referenced covariates to be added to the ssDF; either a continous, a categorical,  "cont" and  "fact", respectively, or both).  
}
\usage{
make.ssDF(scrFrame, buffer, res, cont.cov = F, fact.cov = F)
}
\arguments{
  \item{scrFrame}{
An scrFrame which must be of class 'scrFrame'.
}
  \item{buffer}{
The required buffer distance. This should be in the units of the trap loaction coordinates.
}
  \item{res}{
The required resolution of the state space pixels. This should also be in the units of the trap loaction coordinates.
}
  \item{cont.cov}{
If TRUE, a random continuos covariate is returned with name "cont".
}
  \item{fact.cov}{
If TRUE, a random categorical covariate is returned with name "cont".
}
}
\details{
TBA
}
\value{
A list of dataframes. One dataframe per session. Each dataframe has at minimum an "X" and "Y", column which are pixel coordinates. If selected, continous ("cont") and/or categorical ("fact") covariates are added. The list is of class "ssDF""}
\references{
None
}
\author{
Chris Sutherland
}
\note{
NA
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# simulate an scrFrame and create an associated ssDF object with no covariates:

ssDF <- make.ssDF(scrFrame=sf, buffer=2, res=0.25, cont.cov = F, fact.cov = F)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
