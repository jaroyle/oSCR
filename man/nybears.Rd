\name{nybears}
\alias{nybears}
\docType{data}
\title{
NY black bear data set 
}
\description{
NY black bear data set
}
\usage{data("nybears")}
\format{
  The format is:
List of 10
 $ y2d        : num [1:33, 1:103] 0 0 0 0 0 0 0 0 0 0 ...
 $ K          : num 5
 $ traplocs   :List of 2
  ..$ :'data.frame':	25 obs. of  2 variables:
  .. ..$ X: int [1:25] 1 1 1 1 1 2 2 2 2 2 ...
  .. ..$ Y: int [1:25] 1 2 3 4 5 1 2 3 4 5 ...
  ..$ :'data.frame':	25 obs. of  2 variables:
  .. ..$ X: int [1:25] 1 1 1 1 1 2 2 2 2 2 ...
  .. ..$ Y: int [1:25] 1 2 3 4 5 1 2 3 4 5 ...
 $ elevation  : num [1:4293] 1.246 1.356 1.152 1.091 0.681 ...
 $ trap2raster: int [1:103] 3753 3559 3363 3690 2608 2887 2746 2674 1914 1844 ...
 $ ssgrid     : num [1:4293, 1:2] 28.3 28.4 28.5 28.7 28.8 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "X_UTM" "Y_UTM"
 $ ntel       : num [1:3, 1:4293] 0 0 0 0 0 0 0 0 0 0 ...
 $ ss         : num [1:4293, 1:2] 283448 284466 285483 286501 287519 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "X_UTM" "Y_UTM"
 $ teldata    :'data.frame':	1948 obs. of  8 variables:
  ..$ animalid: Factor w/ 3 levels "CU803","CU818",..: 1 1 1 1 1 1 1 1 1 1 ...
  ..$ fixnum  : int [1:1948] 746 748 749 750 751 752 753 754 755 756 ...
  ..$ Date    : Factor w/ 34 levels "01.08.2011","02.08.2011",..: 18 18 18 18 18 18 18 18 18 18 ...
  ..$ Time    : Factor w/ 1264 levels "0:00:14","0:00:16",..: 55 688 938 973 1036 1075 1154 1197 1239 136 ...
  ..$ Lat     : num [1:1948] 42.3 42.3 42.3 42.3 42.3 ...
  ..$ Long    : num [1:1948] -77.3 -77.3 -77.3 -77.3 -77.3 ...
  ..$ X_UTM   : num [1:1948] 312256 312326 312390 312328 312390 ...
  ..$ Y_UTM   : num [1:1948] 4685583 4685295 4685257 4685283 4685284 ...
 $ teldata2   :'data.frame':	1948 obs. of  3 variables:
  ..$ animalid: Factor w/ 3 levels "CU803","CU818",..: 1 1 1 1 1 1 1 1 1 1 ...
  ..$ X_UTM   : num [1:1948] 312256 312326 312390 312328 312390 ...
  ..$ Y_UTM   : num [1:1948] 4685583 4685295 4685257 4685283 4685284 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Royle, J.A., Chandler, R.B., Sun, C.C. and Fuller, A.K., 2013. Integrating resource selection information with spatial capture-recapture. Methods in Ecology and Evolution, 4(6), pp.520-530.
}
\examples{
data(nybears)



plot (nybears$ss, col="red",pch=".")
points(nybears$traps,pch=5,col="blue")
points(nybears$teldata2[,2:3],pch=".")

colnames(teldata2)<- c("ind","X","Y")
colnames(ss)<- c("X","Y")

teldata<- list(teldata2)
ssDF<- list(ss)

telemetry.processor<- function(ssDF, teldata){
library(FNN)
nsess<- length(teldata)
sbar<- list()
dlst<- list()
nfreq<- list()
for(s in 1:nsess){
s.grid <- as.vector(get.knnx(ssDF[[s]][,c("X","Y")],teldata[[s]][,c("X","Y")],1)$nn.index)
sbar[[s]]<- aggregate(teldata[[s]][,2:3],list(teldata[[s]][,1]),mean)
#  Group.1    X_UTM   Y_UTM
#1   CU803 310627.1 4680376
#2   CU818 279953.5 4699489
#3   CU905 277557.5 4690543
  tmp<-table(s.grid, teldata[[s]][,1])
 tmp2<- matrix(0,nrow=nrow(ssDF[[s]]), ncol=nrow(sbar[[s]]) )
 tmp2[as.numeric(rownames(tmp)),]<- tmp
nfreq[[s]]<- tmp2
dlst[[s]]<-t( e2dist(sbar[[s]][,2:3], ssDF[[s]][,c("X","Y")])  )
#lambda<- exp(-(1/(2*sigma*sigma))*d*d)

}

list("dlst"=dlst, "nfreq" = nfreq,  "sbar"=sbar)

}


data<- telemetry.processor(ssDF,teldata)

## maybe str(nybears) ; plot(nybears) ...
}
\keyword{datasets}
