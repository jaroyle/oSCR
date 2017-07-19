data2oscr <-
  function (edf, sess.col = NULL, id.col = NULL, occ.col = NULL, trap.col = NULL,
            sex.col = NULL, tdf = NULL, K = NULL, ntraps = NULL,
            remove.zeros = TRUE, trapcov.names = NULL,
            remove.extracaps = TRUE, sex.nacode = NULL, tdf.sep = "/")
  {
    
    ## Some safety checks
    ##
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
        cat("error: more than 2 sex codes, no sex.nacode specified",
            fill = TRUE)
        return(NULL)
      }
    }
    if (!is.null(sex.col) & !is.null(sex.nacode)) {
      edf[edf[, sex.col] %in% sex.nacode, sex.col] <- NA
      ux <- length(unique(edf[, sex.col][!is.na(edf[, sex.col])]))
      if (ux > 2) {
        cat("error: more than 2 sex codes, no sex.nacode specified",
            fill = TRUE)
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
      xx <- tdf[[s]][, 4:ncol(tdf[[s]])]
      is.trapcovs <- any(xx[1, ] == tdf.sep)
      if (is.trapcovs){
        xx.check <- (1:ncol(xx))[xx[1, ] == tdf.sep]
        tc.nams <- dimnames(xx)[[2]][(xx.check + 1):ncol(xx)]
        trapcovs[[s]] <- as.matrix(xx[, (xx.check + 1):ncol(xx)])
        colnames(trapcovs[[s]]) <- tc.nams
        trapopp[[s]] <- as.matrix(xx[, 1:(xx.check - 1)])
        all.tcnames <- c(all.tcnames, tc.nams)  # Not used?
      }
      else {
        trapopp[[s]] <- as.matrix(xx[, 1:ncol(xx)])
      }
      occnames[[s]]<- 1:ncol(trapopp[[s]])#dimnames(trapopp[[s]])[[2]] <- not general
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
        y3d[xx[obs, "individual"], xx[obs, "trap"], xx[obs,
                                                       "occasion"]] <- y3d[xx[obs, "individual"], xx[obs,
                                                                                                     "trap"], xx[obs, "occasion"]] + 1
      }
      caphist[[s]] <- y3d
      nn[[s]] <- apply(y3d, c(1), sum)
      
      if(remove.zeros==TRUE){
        if (!is.null(sex.col)) {
          Xsex <- new.edf[[s]][, sex.col]
          if (!is.numeric(Xsex)) {
            Xsex <- as.numeric(as.factor(as.character(Xsex))) -     1
          }
          xx<- cbind(xx, "sex" = Xsex)
          usex[[s]] <- xx[!duplicated(xx[, c("individual",  "sex")]), c("individual", "sex")]
        }
      }else{
        usex[[s]]<- usex.all
      }
      
      
    }
    
    
    
    
    
    for (s in 1:nsess) {
      if (remove.zeros)
        caphist[[s]] <- caphist[[s]][nn[[s]] > 0, , ]
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
    
    if (!is.null(sex.col))
      scrFrame <- make.scrFrame(caphist = caphist, indCovs = sex.oscr,
                                traps = traplocs, trapCovs = NULL, trapOperation = trapopp)
    if (is.null(sex.col))
      scrFrame <- make.scrFrame(caphist = caphist, indCovs = NULL,
                                traps = traplocs, trapCovs = NULL, trapOperation = trapopp)
    
    list(edf = new.edf, y3d = caphist, sex = sex.oscr, traplocs = traplocs,
         trapopp = trapopp, trapcovs = trapcovs, scrFrame = scrFrame)
  }
