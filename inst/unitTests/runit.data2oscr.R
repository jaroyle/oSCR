test.data2oscr <- function()
{
library(oSCR)
data(beardata)
data<- data2oscr(beardata$edf,
                 sess.col = 1,
                 id.col = 2,
                 occ.col = 3,
                 trap.col = 4,
                 tdf = list(beardata$tdf[,1:3]),
                 K = c(8),
                 ntraps = c(38))
sf <- data$scrFrame
ssDF <- make.ssDF(sf, buffer=3, res = 0.5)

out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), scrFrame=sf, ssDF=ssDF,plotit=FALSE , trimS=8)
####pred<- predict.oSCR(out1)

checkEquals(out1$coef[,"mle"],c(-2.1635078,  0.6782806 ,-3.0474904), tol=1e-6)
# checks that a fit error is generated because habitat doesn't exist
checkException( oSCR.fit(model=list(D~habitat,p0~1,sig~1), scrFrame=sf, ssDF=ssDF,plotit=FALSE , trimS=8) )




data(ocelot)
 # need to have occasions labeled same as in TDF
colnames(tdf1) <- c("Detector","X","Y", 1:44)
colnames(tdf2) <- c("Detector","X","Y", 1:96)
 
# This works. Always best to input K and ntraps!
data<-data2oscr(rbind(edf1,edf2),sess.col=1, id.col=2, occ.col=3, trap.col=4,sex.col=5, K=c(44,96),   ntraps=c(23,23),
             tdf=list(tdf1,tdf2),
             remove.zeros=TRUE, remove.extracaps=TRUE,sex.nacode=c("U") )

checkEquals(apply(data$y3d[[1]],c(1),sum), c(1, 1, 1, 2, 4, 1 ,2 ,2 ,1 ,2, 6, 1, 2, 1, 1 ,5, 1, 2 ,1 ,1, 1, 1))
checkEquals(apply(data$y3d[[1]],c(2),sum), c(1 ,0 ,1, 0, 0 ,2 ,3, 2, 6 ,3 ,2 ,2 ,1 ,3 ,3 ,1, 1, 1 ,2, 1 ,2 ,2, 1))
checkEquals( apply(data$y3d[[2]],c(1),sum), c(1,  1  ,1 , 2,  6,  3,  3 , 2,  2,  6 ,10,  4,  1,  1 , 2  ,2 , 1,
  3, 16,  2  ,2,  3 , 1 , 3 , 1, 1))
checkEquals( apply(data$y3d[[2]],c(2),sum),c(4 , 0 , 1 , 1 , 2 , 5 , 4  ,1 , 4 , 5,  4 , 5  ,5 , 9, 10 , 1  ,4 ,
 4 , 4 , 2,  2  ,3 , 0))




# The covariate "distance" is in the input TDF file as a trap covariate. The function data2oscr processes that into a matrix
#  but it does NOT exist for session 1!
dist1 =  rep(0,23)
dist2 =  as.vector(data$trapcovs[[2]])/1000   # Convert to km
# Organize the covariates into a matrix that is J x K (traps x occasions)
dist1 <- matrix(dist1, nrow=23,ncol=44,byrow=FALSE)
dist2 <- matrix(dist2, nrow=23,ncol=96, byrow=FALSE)
# Transformation to exposure to disturbance. Scaling is pretty bad, some extreme values
exposure1<- dist1    # This is 0 for all traps
exposure2<- 1/dist2  # This goes to 0 as distance increase

# Package these into a suitable trapCovs format. Example of 2 covariates, defined for each session here:
tc<-  list(expose=list(exposure1, exposure2), expose2 = list(exposure1^2, exposure2^2) )
tc<-  make.trapCovs(tc)
# Example of 1 covariate
# "list of lists" format
tc <-  list(expose=list(exposure1, exposure2))
tc <-  make.trapCovs(tc)
# example of 2 covariates
#tc<-  list(dist=list(dist1, dist2) , dist2=list(dist1,dist2) )
#tc<- make.trapCovs(tc)
# Make the SCRframe. Sex is not included here. Either way is fine....
ocelot.sf  <-  make.scrFrame(caphist=data$y3d, indCovs=NULL,
                    traps=data$traplocs,trapCovs=tc ,
                    trapOperation=data$trapopp )
checkEquals(ocelot.sf$mmdm, c(1761.706) , tolerance=0.001)



# Make a state-space
ocelot.ss <- make.ssDF(ocelot.sf, buffer=2000, res = 600)
 

##  fit a model
out1 <- oSCR.fit(model=list(D~1,p0~1,sig~1), ocelot.sf, ssDF=ocelot.ss,plotit=FALSE , trimS=2000)
checkEquals( out1$coef[,"mle"], c(-3.342607  ,6.413376, -1.580025), tolerance=0.00001)
 

# Fit a model with "exposure to disturbance"
## default vals for some models are especially bad. Could be geographic scaling issue.
out3 <- oSCR.fit(model=list(D~1,p0~ expose ,sig~1), ocelot.sf, start.vals=c(-3.5, 6.4, 0.035, -1.5),
                 ssDF=ocelot.ss, plotit=FALSE, trimS=2000  )
checkEquals( out3$coef[,"mle"], c( -3.38571156  ,6.38633420  ,0.06387413 ,-1.56412735),tolerance=0.00001)
   

}

