\name{beardata}
\alias{beardata}
\docType{data}
\title{
Fort Drum bear data %%   ~~ data name/kind ... ~~
}
\description{
  Data from the Fort Drum black bear hair snare study, see Gardner et
  al. (2010). 
}
\usage{data("beardata")}
\format{
  The format is:

  beardata is a List of 4

  $ trapmat  :'data.frame':	38 obs. of  2 variables:

  ..$ V1: num [1:38] 448 439 439 442 442 ...

  ..$ V2: num [1:38] 4886 4881 4879 4884 4881 ...

  $ y3d: num [1:47, 1:38, 1:8] 0 0 0 0 0 0 0 0 0 0 ...

  $ edf     : num [1:151, 1:4] 1 1 1 1 1 1 1 1 1 1 ...

  ..- attr(*, "dimnames")=List of 2

  .. ..$ : NULL

  .. ..$ : chr [1:4] "Session" "ID" "Occasion" "trapID"

  $ sex      : num [1:47] 1 1 2 1 1 1 1 2 1 2 ...

  bearfood is a raster with [X,Y,bearfood]

}
\details{
loads the Fort Drum beardata.
A list with "trapmat", "y3d" (3 d array), "edf" (2 d matrix, secr
format, with sessions coding for sex) and "sex" (sex of the individuals,
1=female, 2=male)

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

#
#
# Here is an example of using secr with the "hcov" capability to handle
#    sex. This is equivalent to how oSCR deals with sex
#
library("secr")
library(oSCR)
data("beardata")
 
#oSCR data
y3d <- beardata$y3d
traps.df <-data.frame(X=beardata$trapmat[,1],Y=beardata$trapmat[,2])
sex.df <- data.frame(sex=beardata$sex-1)
bear.sf <- make.scrFrame(caphist=list(y3d),
                traps = list(traps.df),
                indCovs = NULL, 
                trapCovs = NULL)
 
 
bear.ss <- make.ssDF(bear.sf, buffer=20, res=2)

#secr data
secr.traps <- as.matrix(cbind(c(1:dim(beardata$trapmat)[1]), beardata$trapmat * 1000))
colnames(secr.traps) <- c("trapID", "x", "y")
traps1 <- as.data.frame(secr.traps[, 1:3])
trapfile1 <- read.traps(data = traps1, detector = "proximity")
bear.cap <- make.capthist(as.data.frame(beardata$edf), trapfile1, fmt = "trapID", noccasions = 8)
secr.mask <- read.mask(data = bear.ss[[1]]*1000, spacing=2000)
 
# fit a selection of models in secr:
s.bear.0 = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ 1, sigma ~ 1), mask = secr.mask)
#s.bear.t = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ t, sigma ~ 1), mask = secr.mask)
s.bear.b = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ bk, sigma ~ 1), mask = secr.mask)
 
#s.bear.B = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ b, sigma ~ 1), buffer = 20000)
#s.bear.Bt = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ b + t, sigma ~ 1), buffer = 20000)
#s.bear.h2 = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ h2, sigma ~ h2), buffer = 20000)
 
# fit a selection of models in oSCR NEW:
o.bear.0 =   oSCR.fit(model = list(D ~ 1, p0 ~ 1, sig ~ 1),bear.sf, ssDF=bear.ss)
# o.bear.t =   oSCR.fit( model = list(D ~ 1, p0 ~ t, sig ~ 1), scrFrame=bear.sf, ssDF=bear.ss)
o.bear.b =   oSCR.fit(model = list(D ~ 1, p0 ~ b, sig ~ 1), scrFrame=bear.sf, ssDF=bear.ss)






# Models with sex

y3d <- beardata$y3d
traps.df <-data.frame(X=beardata$trapmat[,1],Y=beardata$trapmat[,2])
sex.df <- data.frame(sex=beardata$sex-1)
bear.sf <- make.scrFrame(caphist=list(y3d),
                traps = list(traps.df),
                indCovs = list(sex.df),
                trapCovs = NULL)
 
bear.ss <- make.ssDF(bear.sf, buffer=20, res=2)

# Equivalent models treating session == sex for secr
#s.bear.sex = secr.fit(bear.cap, model = list(D ~ session, g0 ~ session, sigma ~ session), mask = secr.mask)
#o.bear.sex = oSCR.fit(model = list(D ~ 1, p0 ~ sex, sig ~ sex), scrFrame=bear.sf, ssDF=bear.ss)


# Now reformulate things in secr to use the "hcov" approach which is the
#     same as oSCR 
secr.traps <- as.matrix(cbind(c(1:dim(beardata$trapmat)[1]), beardata$trapmat * 1000))
colnames(secr.traps) <- c("trapID", "x", "y")
traps1 <- as.data.frame(secr.traps[, 1:3])
trapfile1 <- read.traps(data = traps1, detector = "proximity")
bear.df <- as.data.frame(beardata$edf)
bear.df$Session <- rep(1,nrow(bear.df))
bear.cap <- make.capthist(bear.df, trapfile1, fmt = "trapID", noccasions = 8)
indCovs <- data.frame(sex = beardata$sex)
covariates(bear.cap) <- indCovs
secr.mask <- read.mask(data = bear.ss[[1]]*1000, spacing=2000)
 
#### try sex effect on p0:
s.bear.0 = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ 1, sigma ~ 1), hcov= "sex", mask = secr.mask)
o.bear.0 = oSCR.fit( model = list(D ~ 1, p0 ~ 1, sig ~ 1), bear.sf, ssDF=bear.ss)
 
s.bear.Sex = secr.fit(bear.cap, model = list(D ~ 1, g0 ~ h2, sigma ~ 1), hcov= "sex", mask = secr.mask)
o.bear.Sex = oSCR.fit(model = list(D ~ 1, p0 ~ sex, sig ~ 1), bear.sf, ssDF=bear.ss)



##
##
## beardata has a habitat raster called "bearfood"
##
##

#library(raster)
#bearfood <-rasterFromXYZ(cbind(Sgrid, bearfood))
plot(bearfood)
points(traps.df,pch=20)
title("Index of bear food")
  


## maybe str(beardata) ; plot(beardata) ...
}
\keyword{datasets}
