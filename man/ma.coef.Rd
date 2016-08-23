\name{ma.coef}
\alias{ma.coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model average coefficient estimates %%  ~~function to do ... ~~
}
\description{
This function will computed model-averaged parameters from a bunch of
SCR models organized into an oSCR.fitList 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ma.coef(ms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{
%%     ~~Describe \code{ms} here~~
}
}
\details{

General info about model averaging:
 
If a coefficient doenst appear in the model there are two ways to model average:
 
1)       Reweight using weights from only models that include the effect

2)       Weight normally but set the coefficient value to 0 if the effect is not present.

 
If you do all subsets, then both are fine and they should be the same. If you don’t fit all subsets, 
then they are not the same.
 
From two of the model averaging R package documentation:
 
From the package MuMIn:

The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. 
An alternative, the ‘full’ average assumes that a variable is included in every model, but in 
some models the corresponding coefficient (and its respective variance) is set to zero. Unlike 
the ‘subset average’, it does not have a tendency of biasing the value away from zero. The ‘full’ 
average is a type of shrinkage estimator and for variables with a weak relationship to the response 
they are smaller than ‘subset’ estimators."
 
AICcmodavg:

"an alternative version of model-averaging parameter estimates that consists in shrinking estimates 
toward 0 to reduce model selection bias as in Burnham and Anderson (2002, p. 152), 
Anderson (2008, pp. 130-132) and Lukacs et al. (2010). Specifically, models without the parameter of 
interest have an estimate and variance of 0."
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ms) 
{
    ests <- ms$coef.tab[, -1]
    ests[is.na(ests)] <- 0
    wts <- ms$aic.tab$weight
    vi <- t(ests != 0) \%*\% wts
    ma.beta.shrink <- colSums(ests * wts)
    tmp <- matrix(ma.beta.shrink, dim(ests)[1], dim(ests)[2], 
        byrow = T)
    if (ms$se) {
        vars <- ms$se.tab[, -1]
        vars[is.na(vars)] <- 0
        ma.se.shrink <- colSums(wts * sqrt(vars^2 + (ests - tmp)^2))
    }
    else {
        ma.se.shrink <- NA
    }
    ma.beta <- numeric(ncol(ests))
    ma.se <- numeric(ncol(ests))
    for (i in 1:ncol(ests)) {
        rmv <- which(ests[, i] != 0)
        rescale.wt <- wts[rmv]/sum(wts[rmv])
        ma.beta[i] <- sum(ests[rmv, i] * rescale.wt)
        if (ms$se) {
            ma.se[i] <- sum(rescale.wt * sqrt(vars[rmv, i]^2 + 
                (ests[rmv, i] - ma.beta[i])^2))
        }
        else {
            ma.se[i] <- NA
        }
    }
    ma.coef <- data.frame(colnames(ests), cbind(ma.beta, ma.se, 
        ma.beta.shrink, ma.se.shrink, vi))
    rownames(ma.coef) <- NULL
    colnames(ma.coef) <- c("Parameter", "Estimate", paste("Std. Error", 
        sep = ""), "Estimate*", paste("Std. Error*", sep = ""), 
        "RVI")
    class(ma.coef) <- c("ma.coef")
    return(ma.coef)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
