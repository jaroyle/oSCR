#library(oSCR)
library(kofnGA)
source("R/scrdesignOF.R")
source("R/scrdesignGA.R")
source("R/plot.scrdesign.R")

SS <- expand.grid(seq(0,30,0.5),seq(0,30,0.5))
TT <- expand.grid(seq(4,26,1),seq(4,26,1))
plot(SS, pch=16,cex=0.2,asp=1)
points(TT, pch=3)

par(mfrow=c(3,3)) 
des1 <- scrdesignGA(statespace = SS, alltraps = TT, ntraps = 20, beta0 = -0.6,
                    sigma = 1.2, crit = 1, popsize=50, keepbest=25, ngen=100)
plot(des1,which=4)

des2 <- scrdesignGA(statespace = SS, alltraps = TT, ntraps = 20, beta0 = -0.6,
                    sigma = 1.2, crit = 2, popsize=50, keepbest=25, ngen=100)
plot(des2,which=4)

des3 <- scrdesignGA(statespace = SS, alltraps = TT, ntraps = 20, beta0 = -0.6,
                    sigma = 1.2, crit = 3, popsize=50, keepbest=25, ngen=100)
plot(des3,which=4)

