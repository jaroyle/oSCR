
################################################################################
## A function that creates an state-space object of class 'ssDF'

  make.ssDF <- function(scrFrame,buffer,res,cont.cov=F, fact.cov=F){
                    ssDF <- list()
                   for(i in 1:length(scrFrame$traps)){
                     trpls <- scrFrame$traps[[i]]
                     bl <- apply(trpls[,c("X","Y")],2,min)
                     tr <- apply(trpls[,c("X","Y")],2,max)
                     sxy <- expand.grid(seq(bl[1]-buffer,tr[1]+buffer,res),
                                        seq(bl[2]-buffer,tr[2]+buffer,res))
                     dd <- apply(e2dist(sxy,trpls[,c("X","Y")]),1,min)
                     ssDF[[i]] <- sxy[dd<=buffer,]
                     colnames(ssDF[[i]]) <- c("X","Y")
                    if(cont.cov==T)
                      ssDF[[i]]$cont <- runif(nrow(ssDF[[i]]), -1, 1)
                    if(fact.cov==T){
                      cat <- c("good","bad")
                      prob <- c(0.6,0.4)
                      tmp <- sample(cat,nrow(ssDF[[i]]),replace=T, prob=prob)
                      ssDF[[i]]$fact <- factor(tmp)
                    }
                   }
                    class(ssDF) <- "ssDF"
                    return(ssDF)
                  }

################################################################################
