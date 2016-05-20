\name{ocelot}
\alias{ocelot}
\docType{data}
\title{
Ocelot camera trapping data.
}
\description{
Joe Kolowski's ocelot data
}
\usage{data("ocelot")}
\format{
  The format is:
 chr "ocelot"
}
\details{
Data from a study of ocelots using a grid of 23 camera traps. The study
occurred twice at the same location and same trap locations, once with
44 nights of operation and then once again with 96 nights of operation
which occurred after some disturbance regime.  The two periods are 
naturally regarded as sessions. 

The data consist of 2 trap deployment files (tdf1 and tdf2) and 2 
encounter data files (edf1 and edf2).   These files should run
in 'secr' directly (examples to be provided....not done yet....)


}
\source{
Kolowski, J.M, and A. Alonso. 2010.
Density and activity patterns of ocelots (Leopardus pardalis) in northern Peru and the impact
of oil exploration activities. Biological Conservation 143:917-925. 
}
\references{
Kolowski, J.M, and A. Alonso. 2010.
Density and activity patterns of ocelots (Leopardus pardalis) in northern Peru and the impact
of oil exploration activities. Biological Conservation 143:917-925. 
}
\examples{


library(oSCR)
data(ocelot)

# Create the oSCR data file
a<-secr2oscr(rbind(edf1,edf2),sex.col=5, K=c(44,96),   ntraps=c(23,23), 
             tdf=list(tdf1,tdf2),   
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") ) 


# process the trap covariate info. This should be done internally but isn't yet.
dist1=rep(0,23)
dist2=as.vector(a$trapcovs[[2]])/1000
dist1<- matrix(dist1,nrow=23,ncol=44,byrow=FALSE)
dist2<- matrix(dist2, nrow=23,ncol=96, byrow=FALSE)
# oSCR wont work if there's only one covariate...so I pad this here with a dummy covariate
# "list of lists" format
tc<-  list(dist=list(dist1, dist2) , dist2=list(dist1,dist2) )
tc<- make.trapCovs(tc)


#Examples
#data(rbs)
#caphist <- rbs$scrFrame$caphist
#traps <- rbs$scrFrame$traps
#trapCovs <- rbs$scrFrame$trapCovs
#sf <- make.scrFrame(caphist = caphist, traps = traps, trapCovs = trapCovs, type = "scr")





scrFrame  <- list(caphist=a$y3d, indCovs=NULL, 
                   type="scr", traps=a$traplocs,trapCovs=tc,
                   trapOperation=a$trapopp )
class(scrFrame) <- "scrFrame"

# make a state-space

 ssDF <- make.ssDF(scrFrame, buffer=2000, res = 600)

##  slowish 
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame, ssDF=ssDF,plotit=FALSE , trimS=2000)
# doesn't work  

## No trim.... really slow!
out1b <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1), ssDF=ssDF,plotit=FALSE )
# works

## Good starts, improves run-time (can make a huge difference in complex models)
out1c <- oSCR.fit(scrFrame,model=list(D~1,p0~1,sig~1), start.vals=c(-3.4, 6.4, -0.58),
                 ssDF=ssDF,plotit=FALSE, trimS=2000 )

## default vals for some models are especially bad. Could be geographic scaling issue.
out1d <- oSCR.fit(model=list(D~1,p0~dist,sig~1), scrFrame,start.vals=c(-3.5, 6.4, 0, -1.5), ssDF=ssDF, plotit=FALSE, trimS=2000  )

out1e <- oSCR.fit(scrFrame,model=list(D~session,p0~1,sig~1),start.vals=c(-3.5, 6.4, -1.5, -1.5), ssDF=ssDF,plotit=FALSE , trimS=2000)



## Summarize output
## Make a density map

pred<- predict.oSCR(scrFrame,out1e,ssDF)

# Plot total density map
par(mfrow=c(1,2))
plot(pred$r[[1]])
plot(pred$r[[2]])

# Plot posterior of s for a specific individual
par(mfrow=c(1,1))
rs<- rasterFromXYZ(cbind(pred$ssDF[[1]][,1:2], pred$preds[[1]][1,]))
plot(rs)
#  library(scrbook)
# spatial.plot(pred$ssDF[[1]][,1:2],pred$preds[[1]][8,],cx=3)



## maybe str(ocelot) ; plot(ocelot) ...
}
\keyword{datasets}
