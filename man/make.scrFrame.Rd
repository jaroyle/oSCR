\name{make.scrFrame}

\alias{make.scrFrame}

\title{Create a formatted oSCR data object}

\description{
A helper function for creating the primary \pkg{oSCR} data object, an scrFrame, used for fitting spatial capture-recapture models with the fitting function \code{\link{oSCR.fit}}.}

\usage{
make.scrFrame(caphist, traps, indCovs = NULL, trapCovs = NULL, trapOperation = NULL, telemetry = NULL, rsfDF = NULL, type = "scr")}

\arguments{

  \item{caphist}{A list of capture history arrays with dimensions \eqn{n x J x K} (individual-by-trap-by-occasion), one array for each session. The arrays are either binary detection/non-detection data or encounter frequencies (counts).}

  \item{traps}{A list of trap coordinates, one for each session. The coordinate columns must be named `X' and `Y' (note these are capitals). All traps used in the study must be provided here and should coincide with the columns of the \code{caphist} array. Traps that are not operational in specific occasions can be `turned off' by supplying the appropriate \code{trapOperation} object.}

  \item{indCovs}{A list, one for each session, containing individual-specific data. Currently two types of individual covariates are used:
  
  \emph{sex}, a binary (two-class) covariate codes as 1 for males, 0 for females, and NA for individuals with unknown sex.
  
  \emph{removed}, a integer value detailing on which occasion an individual was removed from the study (e.g., a mortality event). If there were no removals then this information is not required.\\
  
  If supplied, the columns of the \code{indCovs} must be named \emph{sex} and/or \emph{removed}.}

  \item{trapCovs}{
  A two level list (a list of lists, lol). A list for each session contains a list of data frames, one for each occasion, i.e., trap covariates can be time-varying. The occasion level data has a row for each trap matching the \code{traps} object, and a named column for each covariate.}

  \item{trapOperation}{A list of binary matrices with dimension \eqn{J x K} (trap-by-occasion), one for each session. The matrix denotes whether a trap was operational (1) or not (0) in each occasion. Traps with trapOperation equal to zero does not contribute to the likelihood.}
  
  \item{telemetry}{A list with 2 objects containing data from collared individuals:
  
  1) \code{fixfreq}, a list of matrices with dimension \eqn{n.tel x nG} (collared individual-by-grid cell), one for each session.  The values are counts of the number of telemetry/GPS fixes within each grid cell of the rsfDF, as returned by the \code{telemetry.processor()} function.
  
  2) \code{indCovs}, a list containing individual-specific data, one for each session.  This is structured the same as the \code{indCovs} object for captured individuals, and currently only accommodates sex.
  }
  \item{rsfDF}{
  Similar to the \code{ssDF}, a list of data frames for each session containing the coordinates of the state-space pixels (labeled "X" and "Y") and covariate values used for fitting the resource selection function.  This can be exactly the same as the \code{ssDF}, or represent a grid with an even finer resolution (coarser resolution NOT recommended).

  By specifying the \code{rsfDF} when making the \code{scrFrame}, covariates within the \code{rsfDF} are automatically added to \code{trapCovs} if missing. 
  }
  \item{type}{This should be set to "SCR".}
  }

\details{The scrFrame is the formatted data object required for fitting SCR models using \code{\link{oSCR.fit}}. A variety of error checks are conducted to ensure that the individual data objects are consistent.}

\value{An object of class scrFrame}

\author{Chris Sutherland}

\examples{
data(rbs)
caphist <- rbs.sf$caphist
traps <- rbs.sf$traps
trapCovs <- rbs.sf$trapCovs
sf <- make.scrFrame(caphist = caphist, traps = traps, trapCovs = trapCovs, type = "scr")
}
