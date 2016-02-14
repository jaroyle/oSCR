\name{beardata}
\alias{beardata}
\docType{data}
\title{
bear data %%   ~~ data name/kind ... ~~
}
\description{
comes from gardner et al. 2010... 
}
\usage{data("beardata")}
\format{
  The format is:

  List of 4

  $ trapmat  :'data.frame':	38 obs. of  2 variables:

  ..$ V1: num [1:38] 448 439 439 442 442 ...

  ..$ V2: num [1:38] 4886 4881 4879 4884 4881 ...

  $ bearArray: num [1:47, 1:38, 1:8] 0 0 0 0 0 0 0 0 0 0 ...

  $ flat     : num [1:151, 1:4] 1 1 1 1 1 1 1 1 1 1 ...

  ..- attr(*, "dimnames")=List of 2

  .. ..$ : NULL

  .. ..$ : chr [1:4] "Session" "ID" "Occasion" "trapID"

  $ sex      : num [1:47] 1 1 2 1 1 1 1 2 1 2 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

# Here is a comparison of oSCR and secr
data(beardata)
 
library("secr")
library("scrbook")
data("beardata")
 
#oSCR data
y3d <- beardata$bearArray
traps.df <-data.frame(X=beardata$trapmat[,1],Y=beardata$trapmat[,2])
sex.df <- data.frame(sex=beardata$sex-1)
bear.sf <- list(caphist=list(y3d),
                traps = list(traps.df),
                indCovs = list(sex.df),
                trapCovs = NULL,
                type="scr")
class(bear.sf) <- "scrFrame"
 
bear.ss <- make.ssDF(bear.sf, buffer=20, res=2)
 
 
#secr data
secr.traps <- as.matrix(cbind(c(1:dim(beardata$trapmat)[1]), beardata$trapmat * 1000))
colnames(secr.traps) <- c("trapID", "x", "y")
traps1 <- as.data.frame(secr.traps[, 1:3])
trapfile1 <- read.traps(data = traps1, detector = "proximity")
bear.cap <- make.capthist(as.data.frame(beardata$flat), trapfile1, fmt = "trapID", noccasions = 8)
secr.mask <- read.mask(data = bear.ss[[1]]*1000, spacing=2000)
 
# fit a selection of models in secr:
s.bear.0 = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ 1, sigma ~ 1), mask = secr.mask)
#s.bear.t = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ t, sigma ~ 1), mask = secr.mask)
s.bear.b = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ bk, sigma ~ 1), mask = secr.mask)
s.bear.sex = secr.fit(bear.cap, model = list(D ~ session, g0 ~ session, sigma ~ session), mask = secr.mask)
 
#s.bear.B = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ b, sigma ~ 1), buffer = 20000)
#s.bear.Bt = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ b + t, sigma ~ 1), buffer = 20000)
#s.bear.h2 = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ h2, sigma ~ h2), buffer = 20000)
 
# fit a selection of models in oSCR NEW:
 
o.bear.0 =   oSCR.fit(bear.sf, model = list(D ~ 1, p0 ~ 1,       a1 ~ 1), ssDF=bear.ss)
#o.bear.t =   oSCR.fit(bear.sf, model = list(D ~ 1, p0 ~ t,       a1 ~ 1), ssDF=bear.ss)
o.bear.b =   oSCR.fit(bear.sf, model = list(D ~ 1, p0 ~ b,       a1 ~ 1), ssDF=bear.ss)
o.bear.sex = oSCR.fit(bear.sf, model = list(D ~ 1, p0 ~ sex,     a1 ~ sex), ssDF=bear.ss)


## maybe str(beardata) ; plot(beardata) ...
}
\keyword{datasets}
