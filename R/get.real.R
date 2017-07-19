get.real <- function(model, type = c("dens", "det", "sig", "all")[1], newdata = NULL, 
                     d.factor=1){

  
  pp <- model$outStats$mle
  vcv <- solve(model$rawOutput$hessian)
  names(pp) <- paste(model$outStats$parameters)
  names(pp)[grep("p0.(Intercept)",names(pp),fixed=T)] <- "p0"
  names(pp)[grep("d0.(Intercept)",names(pp),fixed=T)] <- "d0"
  names(pp)[grep("sig.(Intercept)",names(pp),fixed=T)] <- "sig0"
  
  mods <- model$model
  
  if(type == "dens"){
    d.mod <- update.formula(mods[[1]], NULL ~ .) 
    tmp.dm <- model.matrix(d.mod,model$ssDF[[1]])[,,drop=FALSE]
    nms <- paste0("d.beta.",colnames(tmp.dm))
    nms[1] <- "d0"
    id <- match(nms,names(pp))
    nms <- names(pp)[id]
    
    if(is.null(newdata)){
      pred.list <- list()
      for(i in 1:length(model$ssDF)){
        if(!any(c("psi.constant","psi.1") %in% names(pp))){
          session <- rep(paste0("session.",i),nrow(model$ssDF[[i]]))
          tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=1))
          pred <- do.call(rbind, tmp.ls)
          colnames(pred) <- c("estimate", "se", "lwr","upr")
          X <- model$ssDF[[i]]$X
          Y <- model$ssDF[[i]]$Y
          pred.list[[i]] <- data.frame(session,pred,X,Y)
        }else{
          if("psi.constant" %*% names(pp)){
            session <- rep(paste0("session.",i),nrow(model$ssDF[[i]]))
            sex <- rep(c("f","m"),each=nrow(model$ssDF[[i]]))
            tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=3))
            pred1 <- do.call(rbind, tmp.ls)          
            tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=2))
            pred2 <- do.call(rbind, tmp.ls)
            pred <- rbind(pred1,pred2)
            colnames(pred) <- c("estimate", "se", "lwr","upr")
            X <- rep(model$ssDF[[i]]$X,2)
            Y <- rep(model$ssDF[[i]]$Y,2)
            pred.list[[i]] <- data.frame(session,sex,pred,X,Y)
          }
          if("psi.1" %*% names(pp)){
            session <- rep(paste0("session.",i),nrow(model$ssDF[[i]]))
            sex <- rep(c("f","m"),each=nrow(model$ssDF[[i]]))
            tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=5))
            pred1 <- do.call(rbind, tmp.ls)          
            tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=4))
            pred2 <- do.call(rbind, tmp.ls)
            pred <- rbind(pred1,pred2)
            colnames(pred) <- c("estimate", "se", "lwr","upr")
            X <- rep(model$ssDF[[i]]$X,2)
            Y <- rep(model$ssDF[[i]]$Y,2)
            pred.list[[i]] <- data.frame(session,sex,pred,X,Y)
          }
        }
     }
     return(pred.list)
    }else{
      tmp.dm <- model.matrix(d.mod,newdata)[,,drop=FALSE]
      nms <- paste0("d.beta.",colnames(tmp.dm))
      nms[1] <- "d0"
      id <- match(nms,names(pp))
      if(!any(c("psi.constant","psi.1") %in% names(pp))){
        tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=1))
        pred <- do.call(rbind, tmp.ls)
        colnames(pred) <- c("estimate", "se", "lwr","upr")
      }else{
        if("psi.constant" %*% names(pp)){
          sex <- rep(c("f","m"),each=nrow(model$ssDF[[i]]))
          tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=3))
          pred1 <- do.call(rbind, tmp.ls)          
          tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=2))
          pred2 <- do.call(rbind, tmp.ls)
          pred <- rbind(pred1,pred2)
          colnames(pred) <- c("estimate", "se", "lwr","upr")
        }
        if("psi.1" %*% names(pp)){
          sex <- rep(c("f","m"),each=nrow(model$ssDF[[i]]))
          tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=5))
          pred1 <- do.call(rbind, tmp.ls)          
          tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=4))
          pred2 <- do.call(rbind, tmp.ls)
          pred <- rbind(pred1,pred2)
          colnames(pred) <- c("estimate", "se", "lwr","upr")
        }
      }
      return(pred)
    }
  }
  if(type == "det"){
    p.mod <- update.formula(mods[[2]], NULL ~ .)

    if(is.null(newdata)){
      pred.list <- list()
      nt <- nrow(model$scrFrame$traps[[1]])
      tmp.df <- data.frame(session=factor(rep(1:nt,2)),
                           sex=factor(rep(c("0","1"),each=nt)))
      if(!is.null(model$scrFrame$trapCovs)){
        tmp.df <- data.frame(tmp.df,model$scrFrame$trapCovs[[1]][[1]])
      }
      tmp.dm <- model.matrix(p.mod,tmp.df)[,,drop=FALSE]
      nms <- paste0("t.beta.",colnames(tmp.dm))
      nms[1] <- "p0"
      if(any(grepl("t.beta.sex1",nms))){
        nms[grep("t.beta.sex1",nms)] <- "p0.male"
      }
      if(any(grepl("t.beta.session",nms))){
        id <- grep("t.beta.session",nms)
        for(i in 1:length(grep("t.beta.session",nms))){
          nms[id[i]] <- paste0("p0.session",i+1)
        }
      }
      id <- match(nms,names(pp))

      for(i in 1:length(model$scrFrame$caphist)){
        pred.list[[i]] <- list()
        nt <- nrow(model$scrFrame$traps[[i]])
        for(k in 1:dim(model$scrFrame$caphist[[i]])[3]){
          tmp.df <- data.frame(session=factor(rep(1:nt,2)),
                               sex=factor(rep(c("0","1"),each=nt)),
                               occasion=factor(rep(k,nt*2)))
          if(!is.null(model$scrFrame$trapCovs)){
            tmp.df <- data.frame(tmp.df,model$scrFrame$trapCovs[[i]][[k]])
          }
          tmp.dm <- model.matrix(p.mod,tmp.df)[,,drop=FALSE]
          tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=6))
          pred <- do.call(rbind, tmp.ls)
          colnames(pred) <- c("estimate", "se", "lwr","upr")
          pred.list[[i]][[k]] <- data.frame(session=tmp.df$session,
                                       occasion=tmp.df$occasion,
                                       sex=tmp.df$sex,
                                       pred)
        }
      }
    }else{
      tmp.dm <- model.matrix(p.mod,newdata)[,,drop=FALSE]
      nms <- paste0("t.beta.",colnames(tmp.dm))
      nms[1] <- "p0"
      if(any(grepl("t.beta.sex1",nms))){
        nms[grep("t.beta.sex1",nms)] <- "p0.male"
      }
      if(any(grepl("t.beta.session",nms))){
        id <- grep("t.beta.session",nms)
        for(i in 1:length(grep("t.beta.session",nms))){
          nms[id[i]] <- paste0("p0.session",i+1)
        }
      }
      id <- match(nms,names(pp))
      tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=6))
      pred <- do.call(rbind, tmp.ls)
      colnames(pred) <- c("estimate", "se", "lwr","upr")
      return(pred)
    }
  }
  if(type == "sig"){
    s.mod <- update.formula(mods[[3]], NULL ~ .)
    nses <- length(model$scrFrame$caphist)
    if(is.null(newdata)){
      tmp.df <- model$scrFrame$sigCovs
      tmp.df <- rbind(tmp.df,tmp.df)
      tmp.df$sex <- factor(rep(c("0","1"),each=nrow(tmp.df)/2))
      tmp.dm <- model.matrix(s.mod,tmp.df)
      nms <- paste0("sig.",colnames(tmp.dm))
      if(any(grepl("sig.sex1",nms))){
        nms[grep("sig.sex1",nms)] <- "sig.male"
      }
      nms[1] <- "sig0"
      tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=7))
      pred <- do.call(rbind, tmp.ls)
      tmp.df2 <- data.frame(Session=factor(rep(1:nses,each=2)),
                            Sex = factor(rep(c("f","m"))))
      colnames(pred) <- c("estimate", "se", "lwr","upr")
      pred <- data.frame(tmp.df2,pred)
      return(pred)
    }else{
      tmp.df <- newdata
      tmp.dm <- model.matrix(s.mod,tmp.df)
      nms <- paste0("sig.",colnames(tmp.dm))
      if(any(grepl("sig.sex1",nms))){
        nms[grep("sig.sex1",nms)] <- "sig.male"
      }
      nms[1] <- "sig0"
      tmp.ls <- apply(tmp.dm,1,function(x) get.err(x=x,p=pp,vcv=vcv,nms=nms,err=7))
      pred <- do.call(rbind, tmp.ls)
      colnames(pred) <- c("estimate", "se", "lwr","upr")
      return(pred)
    }  
  }
}