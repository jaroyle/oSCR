# ststespace: state space points
# all.traps:  set of trap locations
# fix:        a set of points that are fixed (i.e., always in the design)
# clust.mids: centroid of a cluster (for selction purposes?)
# ntraps:     number of location we WANT
# ndesigns:   number of different 'optimal' designs to return
# nn:         number of neighbours for the swapping algorithm
# beta0:      estimated bld
# sigma:      estimated sigma
# crit: which criterion to use

SCRdesign <- function(statespace=NULL, all.traps=NULL,
                      clusters=NULL, fix=NULL, clust.mids=clust.mids,
                      ntraps=9, ndesigns=10, nn=19, beta0=-0.6, sigma=2,
                      crit=2, N=100){

  if(is.null(statespace)) stop("Must supply a 'statespace' object (coordinates of the study area)")
  if(is.null(all.traps)) stop("Must supply a 'statespace' object (coordinates of the study area)")

  ngrid <- nrow(statespace)
  Dlist <- list()
  Qhistory <- NULL
  Qdesign <- NULL

  #single locations
  if(is.null(clusters)){
    Cd <- round(e2dist(all.traps,all.traps),8)
    NN2 <- NN <- matrix(0,nrow=nrow(Cd),ncol=ncol(Cd))

    for(i in 1:nrow(Cd)){
      xx<-Cd[i,]
      NN[i,] <- (xx>0 & xx<= sort(xx)[nn]) # select nearest nn neighbors
      NN2[i,] <- (xx>0 & xx<= sort(xx)[3])  # select nearest 3 neighbors
    }
  }
  #clusters
  if(!is.null(clusters)){
    Cd <- round(e2dist(clust.mids[,-1],clust.mids[,-1]),8)
    NN2 <- NN <- matrix(0,nrow=nrow(Cd),ncol=ncol(Cd))

    for(i in 1:nrow(Cd)){
      xx <- Cd[i,]
      NN[i,] <- (xx>0 & xx<= sort(xx)[nn])
      NN2[i,] <- (xx>0 & xx<= sort(xx)[3])
    }
  }

  if(is.null(fix)){ X <- all.traps }
  if(!is.null(fix)){ X <- rbind(all.traps,fix) }

  for(m in 1:ndesigns){

    Qbest <- 10^10
    if(is.null(clusters)){
      X.current <- sample(1:nrow(all.traps),ntraps)
      if(is.null(fix)){
        X <- all.traps[X.current,]
      }else{
        X <- rbind(all.traps[X.current,],fix)
      }
    }else{
      X.current <- sample(1:nrow(clusters),ntraps)         # pick 'ntraps' clusters
      which.sites <- apply(clusters[X.current,],2,sum)>0   # which sites are is the design
      if(is.null(fix)){
        X <- all.traps[which.sites,]
      }else{
        X <- rbind(all.traps[which.sites,],fix)
      }
    }

    Q <- Qfn(X=X, S=statespace, N=N, sigma=sigma, beta0=beta0)[crit]     # evaluate Qfn|X
    Qhistory <- c(Qhistory,Q)         # store Q
    Qdesign <- c(Qdesign,m)           # store the design ID
    cat("Initial Q: ", Q, fill=TRUE)  # print starting Q

    if(is.nan(Q)){
      Dlist[[m]]<- list(Q=NA, X=X, X.current=X.current)
      next # go back to the top and start the 'm' loop again
    }
    repeat{ # repeat until 'break' is satisfied i.e. [Qbest == Q]
      for(i in 1:ntraps){                      # here we will replace each site with [nn] alternatives
        chk <- NN[X.current[i],]               # chk is a vector of [1 = a near neighbour]
        chk[X.current] <- 0                    # remove any sites already in the design
        x.consider <- (1:ncol(NN))[chk==1]     # these are the alternatives to consider
        qtest <- rep(10^10,length(x.consider))
        if(length(x.consider)>0){
          for(j in 1:length(x.consider)){
            Xtest.current <- X.current         # this is just a list of clusters thats all.
            Xtest.current[i] <- x.consider[j]  # switch focal with alternative site
            if(!is.null(clusters)){
              which.test <- apply(clusters[Xtest.current,],2,sum)>0 # new design
              if(is.null(fix)){
                Xtest <- all.traps[which.test,]
              }else{
                Xtest <- rbind(all.traps[which.test,],fix)
              }
            }else{
              if(is.null(fix)){
                Xtest <- all.traps[Xtest.current,]
              }else{
                Xtest <- rbind(all.traps[Xtest.current,],fix)
              }
            }
            getQs <- Qfn(X=Xtest,S=statespace,N=N,sigma=sigma,beta0=beta0)
            qtest[j] <- getQs[crit]
          }
        }else{
          qtest <- NaN
        }
        if(any(is.nan(qtest))){
          Dlist[[m]] <- list(Q=NA, X=X, X.current=X.current)
          next
        }
        if(min(qtest) < Q){
          Q <- min(qtest)                    # best switch
          kp <- qtest==min(qtest)            # which switch resulted in the best swith
          X.current[i] <- x.consider[kp][1]  # make the switch permanent
          if(is.null(clusters)){
            X <- all.traps[X.current,]
          }else{
            which.sites <- apply(clusters[X.current,],2,sum)>0
            X <- all.traps[which.sites,]
          }
        }
      }
      if(Qbest == Q){                     #
        break                             #
      }                                   #
      if(Q<Qbest) Qbest <- Q              #
      if(Q>Qbest) cat("ERROR",fill=TRUE)  ## checks whether and better switches can be made
      Qhistory <- c(Qhistory,Q)
      Qdesign <- c(Qdesign,m)
    }

    Dlist[[m]]<- list(Q = Qbest,            # final Q value
                      X = X,                # starting desgin
                      X.current = X.current # final/best design
                      )
    m <- m+1
  }

  Qvec <- rep(NA,length(Dlist))
  Xid <- matrix(NA,nrow=ntraps,ncol=length(Dlist))
  Xlst <- list()

  for(i in 1:length(Dlist)){
    Qvec[i] <- Dlist[[i]]$Q
    Xid[,i] <- Dlist[[i]]$X.current
    Xlst[[i]] <- Dlist[[i]]$X
  }
  design.rank <- order(Qvec)
  tmp.xl <- list()
  tmp.Qh <- NULL
  tmp.Qi <- NULL
  for(i in 1:ndesigns){
    tmp.xl[[i]] <- Xlst[[design.rank[i]]]
    tmp.Qh <- c(tmp.Qh, Qhistory[Qdesign == which(design.rank==i)])
    tmp.Qi <- c(tmp.Qi, Qdesign[Qdesign == which(design.rank==i)])
  }
  Qvec <- Qvec[design.rank]
  Xid <- Xid[,design.rank]
  Xlst <- tmp.xl
  Qdata <- data.frame(score = tmp.Qh, design = tmp.Qi)
  output <- list(Qvec = Qvec, Xid = Xid, Xlst = Xlst, Qdata = Qdata,
                 all.traps = all.traps, statespace = statespace)
  return(output)
}
