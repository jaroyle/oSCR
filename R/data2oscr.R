data2oscr <- 
  function(edf, sess.col = NULL, id.col = NULL, occ.col = NULL, trap.col = NULL,
           sex.col = NULL, tdf = NULL, K = NULL, ntraps = NULL, remove.zeros = TRUE,
           trapcov.names = NULL, remove.extracaps = TRUE, sex.nacode = NULL, 
           tdf.sep = "/"){
    
    #tdf must be [NAME X Y] [optional: trap operation] [sep] [named trap covariates]  
    
    ## Some safety checks
    if(is.null(sess.col) | is.null(id.col) | is.null(occ.col) | is.null(trap.col)){
      cat("required information missing: sess, id, occ or trap",fill=TRUE)
      return(NULL)
    }
    if (is.null(K) | is.null(ntraps)) {
      cat("Must specify K and ntraps", fill = TRUE)
      return(NULL)
    }
    if (!is.null(sex.col) & is.null(sex.nacode)) {
      ux <- length(unique(edf[, sex.col]))
      if (ux > 2) {
        cat("error: more than 2 sex codes, no sex.nacode specified", fill = TRUE)
        return(NULL)
      }
    }
    if (!is.null(sex.col) & !is.null(sex.nacode)) {
      edf[edf[, sex.col] %in% sex.nacode, sex.col] <- NA
      ux <- length(unique(edf[, sex.col][!is.na(edf[, sex.col])]))
      if (ux > 2) {
        cat("error: more than 2 sex codes, no sex.nacode specified", fill = TRUE)
        return(NULL)
      }
    }
    ##
    ## End of safety checks
    ##
    
    
    Xid <- factor(edf[, id.col])
    if (!is.numeric(Xid)) {
      Xid <- as.numeric(as.factor(as.character(Xid)))
    }
    nind<- max(Xid)
    
    # convert to integer
    edf[,id.col]<- Xid
    
    Xsess <- edf[, sess.col]
    # Note this assumes all sessions contained in EDF....possibly not?
    if (!is.numeric(Xsess)) {
      Xsess <- as.numeric(as.factor(as.character(Xsess)))
    }
    # account for traps with no captures and thus not appearing in EDF
    # this has to be done "by session"
    new.edf<- list()
    
    nsess<- max(Xsess)
    if(nsess>1 & length(K) ==1){
      cat("Must have length(K) = nsessions",fill=TRUE)
      return(NULL)
    }
    if(!is.list(tdf)){
      cat("tdf must be a list",fill=TRUE)
      return(NULL)
    }
    if(length(tdf)!=nsess){
      cat("length(tdf) != apparent number of sessions", fill=TRUE)
      return(NULL)
    }
    
    
    ### process the TDF information here
    all.tcnames <- NULL
    traplocs <- list(NULL)
    trapopp <- list(NULL)
    trapcovs <- list(NULL)
    trapnames<- list(NULL)
    occnames<-list(NULL)
    
    for (s in 1:nsess){
      allnames <- tdf[[s]][, 1]
      trapnames[[s]]<- allnames
      if (any(is.na(match(edf[Xsess==s,trap.col],  allnames)))) {
        cat("some trap names in EDF not in TDF", fill = TRUE)
        return(NULL)
      }
      traplocs[[s]] <- as.matrix(tdf[[s]][, 2:3])
      colnames(traplocs[[s]]) <- c("X", "Y")
      if(ncol(tdf[[s]])>3){
        xx <- tdf[[s]][, 4:ncol(tdf[[s]])]
        is.trapcovs <- any(xx[1, ] == tdf.sep, na.rm = TRUE)
        if(is.trapcovs){
          xx.check <- which(xx[1, ] %in% tdf.sep)
          tc.nams <- dimnames(xx)[[2]][(xx.check + 1):ncol(xx)]
          trapcovs[[s]] <- as.matrix(xx[, (xx.check + 1):ncol(xx)])
          colnames(trapcovs[[s]]) <- tc.nams
          if(xx[1, 1] == tdf.sep){
            trapopp[[s]] <- matrix(1,nrow(xx),K[s])
          }else{
            trapopp[[s]] <- as.matrix(xx[, 1:(xx.check - 1)])
          }
          all.tcnames <- c(all.tcnames, tc.nams)  # Not used?
        }
        occnames[[s]]<- 1:K[s]#1:ncol(trapopp[[s]])#dimnames(trapopp[[s]])[[2]] <- not general
        trapopp[[s]] <- matrix(1,nrow(xx),K[s])
      }else{
        trapcovs[[s]] <- NULL
        trapopp[[s]] <- matrix(1,nrow(tdf[[s]]),K[s])
        occnames[[s]]<- 1:K[s]#1:ncol(trapopp[[s]])
      }
    }
    if(all(sapply(trapcovs,is.null))){
      trapcovs <- NULL
    }
    
    
    
    #####
    ## Process the EDF information here
    ##
    if (!is.null(sex.col)) {
      Xsex <- edf[,sex.col]
      if (!is.numeric(Xsex)) {
        Xsex <- as.numeric(as.factor(as.character(Xsex))) -     1
      }
      xx<- cbind(Xid, Xsex)
      usex.all <- xx[!duplicated(xx),]
    }
    
    
    caphist<- list()
    nn<-list()
    usex<- list()
    for(s in 1:nsess){
      
      
      if (any(is.na(match(edf[Xsess==s,occ.col],  occnames[[s]])))) {
        cat("some occassion names in EDF not in TDF", fill = TRUE)
        return(NULL)
      }
      
      new.edf[[s]]<- data.frame(edf[Xsess==s,])
      
      trapid<-  match(as.character(new.edf[[s]][,trap.col]), as.character(trapnames[[s]]) )
      ### levels(new.edf[[s]][,trap.col])<- unique(trapnames[[s]])
      #### levels(new.edf[[s]][,occ.col])<- unique(occnames[[s]])
      occid<- match(as.character(new.edf[[s]][,occ.col]), as.character(occnames[[s]]))
      ntraps[s]<- length(unique(tdf[[s]][,1]))
      
      y3d<- array(0, c(nind, ntraps[s], K[s]) )
      
      xx<- cbind( "individual" = new.edf[[s]][,id.col],
                  "occasion" = occid,
                  "trap" = trapid )
      for (obs in 1:nrow(xx)) {
        y3d[xx[obs, "individual"], 
            xx[obs, "trap"], 
            xx[obs, "occasion"]] <- y3d[xx[obs, "individual"], xx[obs,
                                                                  "trap"], xx[obs, "occasion"]] + 1
      }
      caphist[[s]] <- y3d
      nn[[s]] <- apply(y3d, c(1), sum)
      
      if(remove.zeros==TRUE){
        if (!is.null(sex.col)) {
          Xsex <- new.edf[[s]][, sex.col]
          if (!is.numeric(Xsex)) {
            Xsex <- as.numeric(as.factor(as.character(Xsex))) - 1
          }
          xx<- cbind(xx, "sex" = Xsex)
          usex[[s]] <- xx[!duplicated(xx[, c("individual",  "sex")]), c("individual", "sex")]
        }
      }else{
        if(!is.null(sex.col)) usex[[s]]<- usex.all
      }
      
      
    }
    
    
    
    
    
    for (s in 1:nsess) {
      if (remove.zeros)
        caphist[[s]] <- caphist[[s]][nn[[s]] > 0, , ,drop=F]
      if (remove.extracaps)
        caphist[[s]][caphist[[s]] > 1] <- 1
      nn[[s]] <- apply(caphist[[s]], c(1), sum)
      if (any(nn[[s]] == 0))
        cat("Some individuals in session", s, " have 0 captures",
            fill = TRUE)
      if (max(caphist[[s]] > 1))
        cat("Some individuals in session", s, " captured > 1 time in a trap/occasion",
            fill = TRUE)
    }
    
    
    
    if (!is.null(sex.col)) {
      sex.oscr <- list()
      for (s in 1:nsess) {
        sex.tmp <- usex[[s]][, 2]
        if (is.factor(sex.tmp))
          sex.tmp <- as.numeric(sex.tmp) - 1
        if (!is.factor(sex.tmp))
          sex.tmp <- as.numeric(as.factor(sex.tmp)) - 1
        sex.oscr[[s]] <- data.frame(sex = sex.tmp)
      }
    }
    else sex.oscr = NULL
    
    
    # make the trapCovs
    if(!is.null(trapcov.names)){
      trapCovs <- list()
      trapCov.list <- list()
      levels.list <- list()
      
      for(cov in trapcov.names){#covariate
        #create named list and track whether covariates are levels
        trapCov.list[[cov]] <- list()
        lev.tracker <- NULL
        
        for(i in 1:length(tdf)){#session
          tmp.tdf <- tdf[[i]]
          tmp.ind <- which(colnames(tmp.tdf) %in% paste0(cov,".",1:K[i]))
          tmp.covs <- tmp.tdf[,tmp.ind,drop=FALSE]
          
          #keep track of factor levels
          if(is.character(unlist(tmp.covs)))
            lev.tracker <- sort(unique(c(lev.tracker,unique(unlist(tmp.covs)))))
          if(is.factor(unlist(tmp.covs)))
            lev.tracker <- sort(unique(c(lev.tracker,levels(unlist(tmp.covs)))))
          if(is.numeric(unlist(tmp.covs)))
            lev.tracker <- NA
          levels.list[[cov]] <- lev.tracker
          
          #some checks:
          # make sure columns are named correctly:
          if(any(duplicated(colnames(tmp.covs))))
            stop("Duplicate column names not allowed \n add .k to index occasions")
          
          #make sure there are enough columns (single columns treated as time invariant)
          if(!(ncol(tmp.covs) %in% c(1,K))){
            stop("Each covariate must have either 1 or K columns")
          }else{
            if(ncol(tmp.covs)==1){
              #expand single covariates
              tmp.covs <- matrix(rep(tmp.covs[,paste0(cov,".1")],K[i]),nrow=nrow(tmp.covs),ncol=K[i],byrow=F)
              colnames(tmp.covs) <- paste0(cov,".",1:K[i])
            }else{
              tmp.covs <- as.matrix(tmp.covs)
            }
          }        
          trapCov.list[[cov]][[i]] <- tmp.covs 
        }
      }
      
      nms <- names(trapCov.list)
      
      for(s in 1:length(tdf)){
        trapCovs[[s]] <- list()
        for(k in 1:K[s]){
          tmp.df <- data.frame(Session=rep(factor(s,levels=c(1:length(tdf))),nrow(tdf[[s]])))
          for(cov in names(trapCov.list)){
            if(any(is.na(levels.list[[cov]]))){
              tmp.df[,cov] <- trapCov.list[[cov]][[s]][,k]
            }else{
              tmp.df[,cov] <- factor(trapCov.list[[cov]][[s]][,k],levels=levels.list[[cov]])
            }
          }
          #need to think about how to add a check for auto vs. generated session values
          #          if(any("session" %in% names(trapCov.list))){
          #            colnames(tmp.df) <- c(nms) 
          #          }else{
          #            colnames(tmp.df) <- c("Session",nms)
          #          }
          colnames(tmp.df) <- c("Session",nms)
          trapCovs[[s]][[k]] <- tmp.df
        }
      }
    }else{
      trapCovs <- NULL
    }
    indcovs <- sex.oscr
    caphist <- caphist
    scrFrame <- make.scrFrame(caphist = caphist, 
                              indCovs = sex.oscr,
                              traps = traplocs, 
                              trapCovs = trapCovs, 
                              trapOperation = trapopp)
    
    list(edf = new.edf, y3d = caphist, sex = sex.oscr, traplocs = traplocs,
         trapopp = trapopp, trapcovs = trapcovs, scrFrame = scrFrame)
  }
