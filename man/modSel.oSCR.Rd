\name{modSel.oSCR}
\alias{modSel.oSCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize Multiple Fitted Models
}
\description{
A method for producing an AIC ordered model table and a table of model specific coefficients values.
}
\usage{
modSel.oSCR(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "oSCR.fitList" created by the function fitList.oSCR}
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
%%  ~~who you are~~
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
ms <- modSel.oSCR(fl)
ms
ms$aic
ms$coef
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
