\name{test.oSCR}
\alias{test.oSCR}
\title{Run Validation Tests for oSCR Model Fitting}
\description{
  Fits one or more benchmark models to built-in datasets and compares the
  output to known reference values, allowing the user to verify that
  \code{\link{oSCR.fit}} is producing correct results in their environment.
}
\usage{
test.oSCR(test = 1)
}
\arguments{
  \item{test}{Integer specifying which test(s) to run. Options are:
    \describe{
      \item{\code{1}}{Fits a null SCR model to the red-backed salamander
        data (\code{rbs}) for a single session and compares output to
        known reference estimates.}
      \item{\code{2}}{Fits a sex-structured null SCR model to the ocelot
        camera trapping data (\code{ocelot}) and compares output to known
        reference estimates.}
      \item{\code{3}}{Runs both test 1 and test 2.}
    }
  }
}
\value{
  Prints two output tables to the console for each test: the known reference
  parameter estimates and the estimates from the current model fit. Returns
  \code{NULL} invisibly. The user should visually compare the two tables to
  confirm they match.
}
\details{
  Test 1 uses \code{subFrame} to extract a single session from the
  red-backed salamander dataset, builds a state space with
  \code{make.ssDF}, and fits a null model (\code{D ~ 1, p0 ~ 1, sig ~ 1})
  using provided starting values. The resulting parameter table is printed
  alongside the stored reference output \code{rbs.known}.

  Test 2 fits a null sex-structured model to the ocelot dataset using
  stored \code{ocelot.sf} and \code{ocelot.ss} objects and compares to
  \code{ocelot.known}.

  These tests are useful after installing or updating oSCR, or when
  switching R versions, to confirm numerical reproducibility.
}
\author{
  Chris Sutherland, Andy Royle, Dan Linden
}
\seealso{
  \code{\link{oSCR.fit}}, \code{\link{rbs_ecography}}, \code{\link{ocelot}}
}
\examples{
## Run the salamander test
test.oSCR(test = 1)

## Run the ocelot test
test.oSCR(test = 2)

## Run both tests
test.oSCR(test = 3)
}
\keyword{spatial capture-recapture}
\keyword{validation}
