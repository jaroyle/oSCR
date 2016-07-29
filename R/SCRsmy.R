SCRsmy <-
function (y3d, traplocs=NULL) 
{
    nind <- dim(y3d)[1]
ntraps<-    totcaps <- nperiods <- sprecaps <- rep(NA, nind)
    mmdm<- rep(NA,nind)
    for (i in 1:nind) {
        x <- y3d[i, , ]
        if(!is.null(traplocs) ){
            
         intraps<- apply(x,1,sum)
         if(sum(intraps>0) > 1){
           tmp<- traplocs[intraps>0,]
           dmat<- e2dist(tmp,tmp)
           mmdm[i]<- max(dmat)
         }
       }
        ntraps[i] <- sum(apply(x, 2, sum) > 0)
        ncaps <- sum(x)
        nperiods[i] <- sum(apply(x, 1, sum) > 0)
        sprecaps[i] <- ifelse(ntraps[i] > 1, 1, 0) * ncaps
        totcaps[i] <- sum(x)
    }
    cat("Total captures: ", sum(totcaps), fill = TRUE)
    cat("Spatial recaptures: ", sum(sprecaps), fill = TRUE)
    cat("Ordinary capture events: ", sum(nperiods), fill = TRUE)
    cat("Captures lost in non-spatial model: ", sum(totcaps) - 
        sum(nperiods), fill = TRUE)
    cat("Mean traps/individual: ", mean(ntraps), fill=TRUE)
    if(!is.null(traplocs)){
    cat("MMDM: ", mean(mmdm,na.rm=TRUE),fill=TRUE)
    cat("   (based on ", sum(!is.na(mmdm)), " individuals)", fill=TRUE)
        }
        }
