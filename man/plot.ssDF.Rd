\name{plot.ssDF}

\alias{plot.ssDF}

\title{
Function for plotting oSCR state space objects (ssDF)
}

\description{
Plot the ssDF ojects. The function also allows superimposing traps and individual detections.
}

\usage{
plot.ssDF(ssDF, scrFrame = NULL, collapse = FALSE, spider = FALSE)
}

\arguments{
  \item{ssDF}{
  The oSCR formatted state space object (an ssDF) 
  }
  
  \item{scrFrame}{
  The oSCR formatted encounter and effort object (an scrFrame). If supplied, traps are autmatically added   to the statespace plot. 
  }
  \item{collapse}{
  Logical. Option to collapse a multi-session ssDF onto a single plot.
  }
  \item{spider}{
  Logical. Option to add spatial average detection locations and lines to traps in which those indidivuals   were observed.
  }
}

\details{
Returns either a plot of the discrete state space, the state space and traps (if an scrFrame is provded), and the state space, traps, and detections (if an scrFrame is provded and spider = TRUE).
}

\value{}

\references{}

\author{
Chris Sutherland
}

\note{}

\seealso{}
\examples{
data("rbs_ecography_mods")
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(1,1,1,1))
plot.ssDF(m1$ssDF)                            #just the statespace
plot.ssDF(m1$ssDF,m1$scrFrame)                #statespace and traps
plot.ssDF(m1$ssDF,m1$scrFrame, spider = TRUE) #statespace, traps and detections
}
