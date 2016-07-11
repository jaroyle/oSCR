oSCR.fit.new <-
function (model = list(D ~ 1, p0 ~ 1, sig ~ 1, asu ~1), scrFrame,
          ssDF = NULL, costDF = NULL, distmet = c("euc", "user", "ecol")[1],
          sexmod = c("constant", "session")[1], encmod = c("B", "P", "CLOG")[1],
          DorN = c("D", "N")[1], directions = 8, Dmat = NULL,
          trimS = NULL, start.vals = NULL, PROJ = NULL, pxArea = 1,
          plotit = F, mycex = 1, tester = F, pl = 0, nlmgradtol = 1e-06,
          nlmstepmax = 10, predict = FALSE, smallslow = FALSE, multicatch = FALSE,
          se = TRUE, print.level = 0, getStarts = FALSE){
  my.model.matrix <- function(form, data){
    cont.args <- lapply(data.frame(data[,sapply(data.frame(data), is.factor)]),
                        contrasts, contrasts = FALSE)
    mdm <- suppressWarnings(model.matrix(form, data, contrasts.arg = cont.args))
    return(mdm)
  }
  hessian <- ifelse(se, TRUE, FALSE)
  max.dist <- NULL
  for(i in 1:length(scrFrame$caphist)){
   for(j in 1:nrow(scrFrame$caphist[[i]])) {
     where <- apply(scrFrame$caphist[[i]][j, , ], 1, sum)>0
     if(sum(where) > 1)
     max.dist <- c(max.dist,
                   max(0,dist(scrFrame$traps[[i]][where,c("X", "Y")]),
                   na.rm = T))
   }
  }
  mmdm <- mean(max.dist[max.dist > 0], na.rm = T)
  ptm <- proc.time()
  starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
  cl <- match.call(expand.dots = TRUE)
  if(!require(abind))
    stop("need to install package 'abind'")

  if(!require(Formula))
    stop("need to load package 'Formula'")

  if(distmet == "ecol"){
   if(!require(raster))
     stop("need to install package 'raster'")
   if(!require(gdistance))
     stop("need to install package 'gdistance'")
  }

  if(!inherits(scrFrame, "scrFrame"))
    stop("Data must be of class 'scrFrame'")

  if(encmod %in% c("B","CLOG") & max(unlist(lapply(scrFrame$caphist, max))) > 1)
    stop("Data in caphist must be Binary")

  if(distmet == "ecol" & is.null(PROJ))
     message("Projection not provided, using default: '+proj=utm +zone=12 +datum=WGS84'")

  if(!is.null(ssDF) & length(ssDF) != length(scrFrame$caphist))
    stop("A 'state space' object must be provided for EACH session.")

  if(multicatch){
    for(s in 1:length(scrFrame$caphist)){
      captures <- apply(scrFrame$caphist[[s]], c(1, 3),sum)
      if(any(captures > 1))
        stop("error: multicatch system cannot have > 1 capture.")
    }
  }

  if(predict & is.null(start.vals))
    stop("Starting values required to predict: use MLEs)")

  maxY <- unlist(lapply(scrFrame$caphist, max))
  if(any(maxY > 1) & encmod %in% c("B","CLOG"))
    stop("caphist must be binary when using the Binomial/Cloglog encounter model")

  if(all(maxY == 1) & encmod == "P")
    message("caphist looks binary but Poisson encounter model is selected")
  pars.p0 <- NULL
  names.p0 <- NULL
  pars.sig <- NULL
  names.sig <- NULL
  pars.beta.trap <- NULL
  names.beta.trap <- NULL
  pars.beta.den <- NULL
  names.beta.den <- NULL
  pars.n0 <- NULL
  names.n0 <- NULL
  pars.beta.den <- NULL
  names.beta.den <- NULL
  pars.dist <- NULL
  names.dist <- NULL
  pars.dist <- NULL
  names.dist <- NULL
  pars.sex <- NULL
  names.sex <- NULL
  singleS <- NULL
  singleT <- NULL
  singleG <- NULL
  D <- list()
  YY <- list()
  dm.den <- list()
  tmp.dm <- list()
  dm.trap <- list()
  dm.cost <- list()
  posterior <- list()
  dHPP <- FALSE
  dIPP <- FALSE
  n0Session <- FALSE
  trap.covs <- FALSE
  pDot <- FALSE
  pTime <- FALSE
  pJustsex <- FALSE
  pJustsesh <- FALSE
  pBothsexnsesh <- FALSE
  pBehave <- FALSE
  anySex <- FALSE
  aDot <- FALSE
  aJustsex <- FALSE
  aJustsesh <- FALSE
  aBothsexnsesh <- FALSE
  bDot <- FALSE
  bJustsex <- FALSE
  bJustsesh <- FALSE
  bBothsexnsesh <- FALSE
  #warnings <- list()

  if(length(model) == 3)
    model[[4]] <- formula(~1)

  if((length(all.vars(model[[4]]))>0) & distmet=="euc"){
    stop("non-euclidean 'distmet' must be selected when asu model is specified.")

  for(i in 1:4){
    model[[i]] <- update.formula(model[[i]], NULL ~ .)
  }
  if(is.null(ssDF)) {
    message("Generating a state space based on traps")
    dHPP <- TRUE
    ssDF <- make.ssDF(scrFrame, buffer, res)
  }

  ns <- length(scrFrame$caphist)
  nt <- length(scrFrame$traps)
  nK <- unlist(lapply(scrFrame$caphist, function(x) dim(x)[3]))
  hiK <- max(nK)
  nG <- unlist(lapply(ssDF, nrow))
  nnn <- all(unlist(lapply(ssDF, function(x){"session" %in% names(x) })))
  areaS <- NULL
  if("session" %in% all.vars(model[[1]]) & (!nnn)) {
   for(s in 1:ns){
     ssDF[[s]]$session <- factor(rep(s, nrow(ssDF[[s]])), levels = 1:ns)
   }
  }
  tmp.fn <- function(x){"session" %in% names(x)}
  nnnn <- all(unlist(lapply(scrFrame$trapCovs[[1]][[1]], tmp.fn)))
  if("session" %in% all.vars(model[[2]]) & (!nnnn)){
   for(s in 1:ns){
    for(m in 1:length(scrFrame$trapCovs[[s]])){
      tmp.f <- factor(rep(s, nrow(scrFrame$trapCovs[[s]][[m]])), levels = 1:ns)
      scrFrame$trapCovs[[s]][[m]]$session <- tmp.f
    }
   }
  }

  allvars.D <- all.vars(model[[1]])
  dens.fx <- allvars.D[!allvars.D %in% c("D","session")]
  allvars.T <- all.vars(model[[2]])
  trap.fx <- allvars.T[!allvars.T %in% c("p0","session","sex","t", "T", "b")]
  allvars.sig <- all.vars(model[[3]])
  allvars.dist <- all.vars(model[[4]])
  var.p0.1 <- "sex" %in% allvars.T
  var.p0.2 <- "session" %in% allvars.T
  var.p0.3 <- "t" %in% allvars.T
  var.p0.4 <- any(c("sex:session", "session:sex") %in%
                  attr(terms(model[[2]]),"term.labels"))
  var.sig.1 <- "sex" %in% allvars.sig
  var.sig.2 <- "session" %in% allvars.sig
  var.sig.3 <- any(c("sex:session", "session:sex") %in%
                   attr(terms(model[[3]]), "term.labels"))
  var.b.1 <- "b" %in% attr(terms(model[[2]]), "term.labels")
  var.b.2 <- any(c("b:sex", "sex:b") %in%
                 attr(terms(model[[2]]),"term.labels"))
  var.b.3 <- any(c("b:session", "session:b") %in%
                 attr(terms(model[[2]]),"term.labels"))
  var.b.4 <- any(c("b:session:sex","b:sex:session","sex:session:b",
                   "sex:b:session","session:b:sex","session:sex:b") %in%
                 attr(terms(model[[2]]), "term.labels"))
  pBehave <- any(c(var.b.1, var.b.2, var.b.3, var.b.4))

  for(s in 1:ns){
   if(!is.null(trimS)){
     pixels.prior <- rep(T, nG[s])
     pixels.post <- apply(e2dist(scrFrame$traps[[s]][, c("X", "Y")],
                          ssDF[[s]][, c("X", "Y")]), 2, min) <= trimS
     pixels <- (pixels.prior & pixels.post)
     pixels <- ifelse(pixels, 1, 0)
   }else{
     pixels <- rep(1, nG[s])
   }
   areaS <- c(areaS, sum(pixels) * pxArea)
  }

  if(length(trap.fx) > 0){
    trap.covs <- TRUE
    tcovnms <- colnames(scrFrame$trapCovs[[1]][[1]])
    tCovMissing <- trap.fx[which(!trap.fx %in% tcovnms)]
    if(length(tCovMissing) > 0){
      stop("I cant find these covariates in 'scrFrame$trapCovs':",
      for (i in tCovMissing) print(i))
    }
    mod2 <- update(model[[2]], ~. - sex - session - t - b -b:sex - sex:b -
                   b:session - session:b - b:session:sex - b:sex:session -
                   sex:session:b - sex:b:session - session:b:sex - session:sex:b)
    if(any(c("session") %in% allvars.T))
      tSession <- TRUE
    for(s in 1:ns){
      tmp.dm <- list()
      for(k in 1:nK[s]){
        tmp.dm[[k]] <- model.matrix(mod2, scrFrame$trapCovs[[s]][[k]])[,-1,drop=FALSE]
        if(s == 1 && k == 1)
          t.nms <- colnames(tmp.dm[[k]])
        if(nrow(tmp.dm[[k]]) != nrow(scrFrame$trapCovs[[s]][[k]])){
          mis <- setdiff(rownames(scrFrame$trapCovs[[s]][[k]]),
                         rownames(my.model.matrix(mod2, scrFrame$trapCovs[[s]][[k]])))
          tmp.insert <- matrix(NA, length(mis), ncol(tmp.dm[[k]]))
          row.names(tmp.insert) <- mis
          tmp.dm[[k]] <- rbind(tmp.dm[[k]], tmp.insert)
          tmp.dm[[k]] <- tmp.dm[[k]][order(as.numeric(row.names(tmp.dm[[k]]))),]
        }
      }
      dm.trap[[s]] <- tmp.dm
    }
    if(any(paste0("session:", t.nms) %in% attr(terms(model[[2]]),"term.labels"))){
      t.nms.sess <- t.nms[which(paste0("session:", t.nms) %in%
                                attr(terms(model[[2]]), "term.labels"))]
      tmpTsess <- rep(1:ns, each = length(t.nms.sess))
      tmpTcovs <- rep(t.nms.sess, ns)
      names.beta.trap1 <- paste("t.beta.",tmpTcovs,".session", tmpTsess, sep = "")
      if(length(t.nms) > length(t.nms.sess)){
        t.nms.nosess <- t.nms[!t.nms %in% t.nms.sess]
        names.beta.trap <- c(names.beta.trap1,paste("t.beta.",t.nms.nosess,sep=""))
      }else{
        names.beta.trap <- names.beta.trap1
      }
      pars.beta.trap <- rep(0,length(names.beta.trap))
    }else{
      names.beta.trap <- paste("t.beta.", t.nms, sep = "")
      pars.beta.trap <- rep(0,length(names.beta.trap))
    }
  }

  if(length(dens.fx) > 0){
    dcovnms <- colnames(ssDF[[1]])
    dCovMissing <- dens.fx[which(!dens.fx %in% dcovnms)]
    if(length(dCovMissing) > 0){
      stop("I can't find these covariates in 'ssDF'", for (i in dCovMissing) print(i))
    }
  }
  if(DorN == "N"){
   if(length(dens.fx) > 0){
     dIPP <- TRUE
     mod1 <- update(model[[1]], ~. - sex - session)
     for(s in 1:ns){
       dm.den[[s]] <- model.matrix(mod1, ssDF[[s]])[,-1,drop=FALSE]
       if (s == 1)
        d.nms <- colnames(dm.den[[s]])
     }
     if("Session" %in% all.vars(model[[1]])){
       tmpDsess <- rep(1:ns, each = length(d.nms))
       tmpDcovs <- rep(d.nms, ns)
       names.beta.den <- paste("d.beta.",tmpDcovs,".session",tmpDsess, sep = "")
       pars.beta.den <- rep(0,length(names.beta.den))
       names.n0 <- paste0("n0.session", 1:ns)
       pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
     }
     if("session" %in% all.vars(model[[1]])){
       names.beta.den <- paste0("d.beta.", d.nms, sep = "")
       pars.beta.den <- rep(0,length(names.beta.den))
       names.n0 <- paste("n0.session", 1:ns)
       pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
     }
     if(!("session" %in% all.vars(model[[1]])) &
        !("Session" %in% all.vars(model[[1]]))){
        names.beta.den <- paste0("d.beta.", d.nms, sep = "")
        pars.beta.den <- rep(0,length(names.beta.den))
        names.n0 <- paste0("n0.")
        pars.n0 <- log(mean(unlist(lapply(scrFrame$caphist,nrow))))
     }
   }else{
     dHPP <- TRUE
     names.beta.den <- NULL
     pars.beta.den <- NULL
     for(s in 1:ns){
       dm.den[[s]] <- my.model.matrix(~1, ssDF[[s]])
     }
     if("Session" %in% all.vars(model[[1]])){
       n0Session <- TRUE
       names.n0 <- paste0("n0.session", 1:ns)
       pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
     }
     if("session" %in% all.vars(model[[1]])){
       n0Session <- TRUE
       names.n0 <- paste0("n0.session", 1:ns)
       pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
     }
     if(!("session" %in% all.vars(model[[1]])) &
        !("Session" %in% all.vars(model[[1]]))){
        names.n0 <- paste0("n0.")
        pars.n0 <- log(mean(unlist(lapply(scrFrame$caphist,nrow))))
     }
   }
  }
  if(DorN == "D"){
    mod1 <- update(model[[1]], ~. - sex)
    for(s in 1:ns){
      dm.den[[s]] <- model.matrix(mod1, as.data.frame(ssDF[[s]]))
    }
    d.nms <- colnames(dm.den[[1]])
    names.beta.den <- paste("d.beta", d.nms, sep = ".")
    chx <- grep(fixed=TRUE,"Intercept", names.beta.den)
    if(length(chx) > 0)
      names.beta.den[chx] <- "d0.(Intercept)"
    pars.d0 <- log(mean((unlist(lapply(scrFrame$caphist,nrow)))/unlist(lapply(ssDF,nrow))))
    pars.beta.den <- c(pars.d0, rep(0.1, length(names.beta.den) - 1))
    pars.n0 <- NULL
    names.n0 <- NULL
  }
  if(distmet == "ecol" && length(allvars.dist) == 0){
    message("distmet='ecol' without a cost surface.\n
             Euclidean distance model will be used.")
    distmet <- "euc"
  }
  if(length(allvars.dist) > 0){
    ccovnms <- colnames(costDF[[1]])
    cCovMissing <- allvars.dist[which(!allvars.dist %in% ccovnms)]
    if(length(cCovMissing) > 0){
      stop("I cant find theses covariates in 'costDF'",
            for (i in cCovMissing) print(i))
    }
  }

  mod4 <- update(model[[4]], ~. - sex - session - 1)
  if(distmet == "ecol"){
   for(s in 1:ns){
     dm.cost[[s]] <- my.model.matrix(mod4, costDF[[s]])
     names.dist <- paste("c.beta.", allvars.dist, sep = "")
     pars.dist <- rep(0, length(names.dist))
   }
  }
  if(!smallslow){
   if(distmet == "euc"){
    for(s in 1:ns){
      D[[s]] <- e2dist(scrFrame$traps[[s]][, c("X","Y")], ssDF[[s]][, c("X", "Y")])
    }
   }
  }
  if("indCovs" %in% names(scrFrame)){
   if("sex" %in% names(scrFrame$indCovs[[1]])){
     anySex <- TRUE
   }
  }

  tmp.p0.names <- "p0.(Intercept)"
  if(sum(var.p0.1, var.p0.2, var.p0.3, var.p0.4) == 0){
    tmp.p0.names <- "p0.(Intercept)"
    pDot <- TRUE
  }
  if(var.p0.1 & !var.p0.2){
    tmp.p0.names <- c("p0.(Intercept)", "p0.male")
    pJustsex <- TRUE
  }
  if(!var.p0.1 & var.p0.2){
   if(ns > 1){
     tmp.p0.names <- c("p0.(Intercept)", paste0("p0.session", 2:ns))
     pJustsesh <- TRUE
   }else{
     tmp.p0.names <- "p0.(Intercept)"
     pDot <- TRUE
   }
  }
  if(var.p0.1 & var.p0.2){
    tmp.p0.names <- c("p0.(Intercept)", "p0.male")
    if(ns > 1){
      tmp.p0.names <- c(tmp.p0.names, paste0("p0.session", 2:ns))
      pJustsesh <- TRUE
    }else{
      tmp.p0.names <- tmp.p0.names
      pJustsex <- TRUE
    }
    pJustsex <- TRUE
  }
  if(var.p0.3){
    tmp.p0.names <- c(tmp.p0.names, paste0("p0.t", 2:hiK))
    pTime <- TRUE
  }
  if(var.p0.4){
   if(ns > 1){
     tmp.p0.names <- c("p0.(Intercept)",
                       paste0("p0.f.session", 2:ns),
                       paste0("p0.m.session", 1:ns))
     pBothsexnsesh <- TRUE
   }else{
     tmp.p0.names <- c("p0.(Intercept)", "p0.male")
     pJustsex <- TRUE
   }
  }
  names.p0 <- tmp.p0.names
  pars.p0 <- rep(0, length(names.p0))
  pars.p0[1] <- -1.5

  if(any(var.p0.1, var.p0.1, var.sig.1) && !anySex)
    stop("Sex defined in a model but no sex data provided.")

  if(var.b.1 & !var.b.2 & !var.b.3 & !var.b.4){
    pars.p0 <- c(pars.p0, 0)
    names.p0 <- c(names.p0, "p.behav")
    bDot <- TRUE
  }
  if(var.b.2 & !var.b.4){
    pars.p0 <- c(pars.p0, 0, 0)
    names.p0 <- c(names.p0, "p.behav.f", "p.behav.m")
    bJustsex <- TRUE
  }
  if(var.b.3 & !var.b.4){
    pars.p0 <- c(pars.p0, rep(0, ns))
    names.p0 <- c(names.p0, paste0("p.behav.session", 1:ns))
    bJustsesh <- TRUE
  }
  if(var.b.4){
    pars.p0 <- c(pars.p0, rep(0, 2 * ns))
    names.p0 <- c(names.p0,
                  paste0("p.behav.f.session", 1:ns),
                  paste0("p.behav.m.session", 1:ns))
    bBothsexnsesh <- TRUE
  }
  tmp.sig.names <- "sig.(Intercept)"
  if(sum(var.sig.1, var.sig.2, var.sig.3) == 0){
    aDot <- TRUE
  }
  if(var.sig.2 & !var.sig.3){
   if (ns > 1){
     tmp.sig.names <- c(tmp.sig.names, paste0("sig.session.",2:ns))
     aJustsesh <- TRUE
   }else{
     aDot <- TRUE
   }
  }
  if(var.sig.1 & !var.sig.3){
    tmp.sig.names <- c(tmp.sig.names, "sig.male")
    aJustsex <- TRUE
  }
  if(var.sig.3){
   if(ns > 1){
     tmp.sig.names <- c(tmp.sig.names,
                        paste0("sig.f.session",:ns),
                        paste0("sig.m.session", 1:ns))
     aBothsexnsesh <- TRUE
   }else{
     aJustsex <- TRUE
   }
  }
  names.sig <- tmp.sig.names
  pars.sig <- rep(0.1, length(names.sig))
  pars.sig[1] <- log(0.5 * mmdm)
  YY <- scrFrame$caphist
  if(anySex){
   if(sexmod == "constant"){
     pars.sex <- 0
     names.sex <- "psi.constant"
   }
   if(sexmod == "session"){
     pars.sex <- rep(0, ns)
     names.sex <- paste("psi", 1:ns, sep = "")
   }
  }else{
    pars.sex <- NULL
    names.sex <- NULL
  }

  pv <- round(c(pars.p0, pars.sig, pars.beta.trap, pars.beta.den, pars.dist,
                pars.n0, pars.sex), 4)
  pn <- c(names.p0, names.sig, names.beta.trap, names.beta.den, names.dist,
          names.n0, names.sex)

  if(!is.null(start.vals)){
   if(length(pv) == length(start.vals)){
     pv <- start.vals
   }else{
     message("Number of starting values provided doesnt match the model.")
     message("Using internally generated values (see getStarts = T output).")
   }
  }
  if(pBehave){
    prevcap <- list()
    for(s in 1:length(YY)){
      Ys <- YY[[s]]
      prevcap[[s]] <- array(0, dim = c(dim(Ys)[1], dim(Ys)[2],dim(Ys)[3]))
      first <- matrix(0, dim(Ys)[1], dim(Ys)[2])
      for(i in 1:dim(Ys)[1]){
       for(j in 1:dim(Ys)[2]){
        if(sum(Ys[i, j, ]) > 0){
          first[i, j] <- min((1:(dim(Ys)[3]))[Ys[i, j, ] > 0])
          prevcap[[s]][i, j, 1:first[i, j]] <- 0
          if(first[i, j] < dim(Ys)[3])
            prevcap[[s]][i, j, (first[i, j] + 1):(dim(Ys)[3])] <- 1
        }
       }
      }
      zeros <- array(0, c(1, dim(prevcap[[s]])[2], dim(prevcap[[s]])[3]))
            prevcap[[s]] <- abind(prevcap[[s]], zeros, along = 1)
    }
  }

  if(getStarts == TRUE){
    oSCR.start <- list(parameters = pn, values = pv)
    return(oSCR.start)
  }

  nR<- nC<- list()
  trimR <- trimC <- list()
  for(s in 1:length(YY)){
    Ys <- YY[[s]]
    if(!multicatch){
      zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
      Ys <- abind(Ys, zeros, along = 1)
    }
    if(multicatch){
      zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
      Ys <- abind(Ys, zeros, along = 1)
    }
    trimR[[s]] <- list()
    trimC[[s]] <- list()
    nR[[s]]<- nC[[s]]<- list()
    for(i in 1:nrow(Ys)){
      trimR[[s]][[i]]<- list()
      nR[[s]][[i]]<- list()
      if(is.null(trimS)){
        pp <- rep(T, ncol(Ys))
        trimC[[s]][[i]] <- rep(T, nG[s])
        for(k in 1:nK[s]){
          trimR[[s]][[i]][[k]] <- pp
        }
      }else{
      if(i < nrow(Ys)){
        pp <- apply(Ys[i, , ], 1, sum) > 0
        for(k in 1:nK[s]){
         if (!is.null(scrFrame$trapOperation)){
           trimR[[s]][[i]][[k]] <- (apply(rbind(rep(trimS*3+2, nrow(scrFrame$traps[[s]])),
             e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X", "Y")]), sum(pp), 2),
             scrFrame$traps[[s]][,c("X", "Y")])), 2, min) <= (2 * trimS)) &
             (scrFrame$trapOperation[[s]][,k]==1)
         }else{
           trimR[[s]][[i]][[k]] <- apply(rbind(rep(trimS*3+2, nrow(scrFrame$traps[[s]])),
             e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X", "Y")]), sum(pp), 2),
             scrFrame$traps[[s]][, c("X", "Y")])), 2, min) <= (2 * trimS)
         }
         nR[[s]][[i]][[k]]<- sum(trimR[[s]][[i]][[k]])
        }
        trimC[[s]][[i]] <- apply(rbind(rep(trimS+2,nG[s]),
          e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X", "Y")]), sum(pp), 2),
          ssDF[[s]][, c("X","Y")])), 2, min, na.rm = T) <= trimS
      }else{
        pp <- rep(T, ncol(Ys))
        timC[[s]][[i]] <- apply(rbind(rep(trimS + 2, nG[s]),
          e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X", "Y")]), sum(pp), 2),
          ssDF[[s]][, c("X","Y")])), 2, min, na.rm = T) <= trimS
        for(k in 1:nK[s]){
         if(!is.null(scrFrame$trapOperation)){
           trimR[[s]][[i]][[k]]<- pp & (scrFrame$trapOperation[[s]][,k]==1)
         }else{
           trimR[[s]][[i]][[k]]<- pp
         }
         nR[[s]][[i]][[k]]<- sum(trimR[[s]][[i]][[k]])
        }
      }}
      if(multicatch){
       for(k in 1:nK[s]){
         trimR[[s]][[i]][[k]] <- rep(T, length(trimR[[s]][[i]][[k]]))
         nR[[s]][[i]][[k]]<- sum(trimR[[s]][[i]][[k]])
       }
      }
      nC[[s]][[i]]<- sum(trimC[[s]][[i]])
    }
  }


###############################
## NO SEX! Likelihood function

  msLL.nosex <- function(pv = pv, pn = pn, YY = YY, D = D, hiK = hiK, nG = nG,
                         nK = nK, dm.den = dm.den, dm.trap = dm.trap){

  #p0 parameters:
  alpha0 <- array(0, dim = c(ns, hiK, 2))
  tmpP <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.(Intercept)", names.p0)]]
  if(pDot & !pTime){
    alpha0[, , ] <- tmpP
  }
  if(pDot & pTime){
    tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
    for(s in 1:ns){
      alpha0[s, , 1] <- tmpP + tmpT
  }}
  if(pJustsesh & !pTime){
    tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",names.p0)]])
    for(s in 1:ns){
      alpha0[s, , 1] <- tmpP + tmpSS[s]
  }}
  if(pJustsesh & pTime){
    tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
    tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session", names.p0)]])
    for(s in 1:ns){
     alpha0[s, , 1] <- tmpP + tmpSS[s] + tmpT
  }}
  BRmat <- array(0, c(ns, hiK, 1))
  if(bDot){
    BRmat[, , 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav",names.p0)]]
  }
  if(bJustsesh){
   for(k in 1:hiK){
     BRmat[, k, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.session",names.p0)]]
  }}
  alpha0[, , 2] <- alpha0[, , 1] + BRmat[, , 1]

  alphsig <- numeric(ns)
  if(aDot){
    tmpA <- pv[pn %in% names.sig[grep(fixed=TRUE,"sig.(Intercept)", names.sig)]]
    for(s in 1:ns){
     alphsig[s] <- tmpA
  }}
  if(aJustsesh){
    tmpA <- pv[pn %in% names.sig[grep(fixed=TRUE,"sig.(Intercept)", names.sig)]]
    tmpSS <- c(0, pv[pn %in% names.sig[grep(fixed=TRUE,"sig.session",names.sig)]])
    for(s in 1:ns){
      alphsig[s] <- tmpA + tmpSS[s]
    }}
  alphsig <- 1/(2 * exp(alphsig)^2)

  if(trap.covs){
    t.beta <- matrix(NA, ns, length(t.nms))
    if(any(paste0("session:", t.nms)%in%attr(terms(model[[2]]),"term.labels"))){
     for(s in 1:ns){
      if(length(t.nms) > length(t.nms.sess)){
        t.beta[s, ] <- pv[pn %in%
          c(names.beta.trap[grep(fixed=TRUE,paste("session",s, sep = ""),
          names.beta.trap)], names.beta.trap[as.vector(unlist(sapply(t.nms.nosess,
          function(x){grep(fixed=TRUE,x, names.beta.trap)})))])]
      }else{
        t.beta[s, ] <- pv[pn %in% c(names.beta.trap[grep(fixed=TRUE,
          paste("session", s, sep = ""), names.beta.trap)])]
     }}
    }else{
      for(s in 1:ns){
        t.beta[s, ] <- pv[pn %in% names.beta.trap]
  }}}
  if(DorN == "N"){
   if(dIPP){
     d.beta <- matrix(NA, ns, length(d.nms))
     if("session" %in% all.vars(model[[1]])){
       for(s in 1:ns){
         d.beta[s, ] <- pv[pn %in% names.beta.den[grep(fixed=TRUE,
           paste("session",s, sep = ""), names.beta.den)]]
       }
     }else{
       for(s in 1:ns){
         d.beta[s, ] <- pv[pn %in% names.beta.den]
     }}}
  }else{
    d.beta <- pv[pn %in% names.beta.den]
  }
  if(distmet == "ecol"){
    dist.beta <- pv[pn %in% names.dist]
  }
  if(n0Session)
    n0 <- exp(pv[pn %in% names.n0])
  if(!n0Session)
    n0 <- rep(exp(pv[pn %in% names.n0]), ns)

  outLik <- 0
  if(predict){
    preds <- list()
    lik.bits <- list()
    ss.bits <- list()
  }
  for(s in 1:length(YY)){
    Ys <- YY[[s]]
    if(predict)
      preds[[s]] <- matrix(NA, nrow = nrow(Ys) + 1, ncol = nrow(ssDF[[s]]))
    if(!multicatch){
      zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
      Ys <- abind(Ys, zeros, along = 1)
    }
    if(multicatch){
      zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
      Ys <- abind(Ys, zeros, along = 1)
    }
    if(distmet == "ecol"){
      cost <- exp(dm.cost[[s]] %*% exp(dist.beta))
      costR <- rasterFromXYZ(cbind(costDF[[s]][, c(1,2)], cost))
      if(is.null(PROJ)){
        projection(costR) <- "+proj=utm +zone=12 +datum=WGS84"
      }else{
        projection(costR) <- PROJ
      }
      trFn <- function(x) (1/(mean(x)))
      tr <- transition(costR, transitionFunction = trFn,direction = directions)
      trLayer <- geoCorrection(tr, scl = F)
      D[[s]] <- costDistance(trLayer,
                             as.matrix(scrFrame$traps[[s]][,c("X", "Y")]),
                             as.matrix(ssDF[[s]][, c(c("X","Y"))]))
    }
    if(smallslow){
     if(distmet == "euc"){
       D <- list()
       D[[s]]<-e2dist(scrFrame$traps[[s]][,c("X","Y")], ssDF[[s]][, c("X", "Y")])
     }
    }
    lik.marg <- rep(NA, nrow(Ys))
    if(!is.null(trimS)){
      pixels.prior <- rep(T, nG[s])
      pixels.post <- apply(D[[s]], 2, min) <= trimS
      pixels <- (pixels.prior & pixels.post)
      pixels <- ifelse(pixels, 1, 0)
    }else{
      pixels <- rep(1, nG[s])
    }
    if(DorN == "N"){
     if(dIPP){
       d.s <- exp(dm.den[[s]] %*% d.beta[s, ])
       pi.s <- (d.s * pixels)/sum(d.s * pixels)
     }
     if(dHPP){
      pi.s <- pixels/(sum(pixels))
     }
    }else{
      d.s <- exp(dm.den[[s]] %*% d.beta)
      pi.s <- (d.s * pixels)/sum(d.s * pixels)
    }
    Kern <- exp(-alphsig[s] * D[[s]]^2)
    for(i in 1:nrow(Ys)){
      if(plotit){
        pp <- sum(trimR[[s]][[i]][[k]])
        plot(ssDF[[s]][,c("X", "Y")], pch=16, col="grey", cex=0.5, asp=1,
             main=paste("Session:",s," Individual: ",i," traps: ",sum(pp),sep=" "))
        points(ssDF[[s]][trimC[[s]][[i]],c("X", "Y")], pch=16, col=2, cex=mycex)
        points(scrFrame$traps[[s]][,c("X", "Y")], pch=15, col="grey", cex=1)
        points(scrFrame$traps[[s]][trimR[[s]][[i]][[k]], c("X")],
               scrFrame$traps[[s]][trimR[[s]][[i]][[k]],c("Y")], pch=15, col=1, cex=1)
        pp <- apply(Ys[i,,],1,sum)>0
        points(scrFrame$traps[[s]][pp,c("X", "Y")], pch=15, col=3, cex=mycex)
      }
      lik.cond <- numeric(nG[s])
      for(k in 1:nK[s]){
       if(pBehave){
         a0 <- alpha0[s,k,1] * (1 - c(prevcap[[s]][i,,k])) +
               alpha0[s,k,2] * c(prevcap[[s]][i,,k])
       }else{
         a0 <- rep(alpha0[s, k, 1], nrow(D[[s]]))
       }
       if(trap.covs){
         a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
       }
       if(encmod == "B")
         probcap <- c(plogis(a0[trimR[[s]][[i]][[k]]])) *
                    Kern[trimR[[s]][[i]][[k]], trimC[[s]][[i]] ]
       if(encmod %in% c("P","CLOG"))
         probcap <- c(exp(a0[trimR[[s]][[i]][[k]]])) *
                      Kern[trimR[[s]][[i]][[k]], trimC[[s]][[i]]]
       if(!multicatch){
        if(encmod == "B"){
          probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]][[k]], k],
            sum(trimC[[s]][[i]])), 1, probcap[1:length(probcap)], log = TRUE))
        }
        if(encmod == "P"){
          probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]][[k]], k],
            sum(trimC[[s]][[i]])), probcap[1:length(probcap)], log = TRUE))
        }
        if(encmod == "CLOG"){
          probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i, trimR[[s]][[i]][[k]], k],
            sum(trimC[[s]][[i]])),1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
        }
       }else{
         probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
         probcap <- t(t(probcap)/colSums(probcap))
         vvv <- rep(c(Ys[i, trimR[[s]][[i]][[k]], k],
                      1 - any(Ys[i, trimR[[s]][[i]][[k]], k] > 0)),
                    sum(trimC[[s]][[i]]))
         vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv == 1])
         probcap[1:length(probcap)] <- vvv
       }
       if(is.matrix(probcap))
         lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + colSums(probcap)
       if(!is.matrix(probcap))
         lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + probcap
      }

      lik.cond[trimC[[s]][[i]]] <- exp(lik.cond[trimC[[s]][[i]]])
      lik.marg[i] <- sum(lik.cond * pi.s)
      if(predict)
        preds[[s]][i, ] <- (lik.cond * pi.s)/lik.marg[i]
    }

    if(!predict){
     if(DorN == "N"){
       nv <- c(rep(1, length(lik.marg) - 1), n0[s])
       part1 <- lgamma((nrow(Ys) - 1) + n0[s] + 1) -  gamma(n0[s] + 1)
       part2 <- sum(nv * log(lik.marg))
     }
     if(DorN == "D"){
       nv <- c(rep(1, length(lik.marg) - 1), 1)
       atheta <- 1 - lik.marg[nrow(Ys)]
       nind <- nrow(Ys) - 1
       part1 <- nind * log(sum(d.s * pixels)) - sum(d.s * pixels) * atheta
       part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
     }
     ll <- -1 * (part1 + part2)
     outLik <- outLik + ll
    }
    if(predict){
      lik.bits[[s]] <- cbind(lik.mar = lik.marg)
      ss.bits[[s]] <- cbind(pi.s, d.s, lik.cond)
      colnames(ss.bits[[s]]) <- c("pi.s", "d.s", "lik.cond")
    }
  }
  if(!predict){
    out <- outLik
    return(out)
  }
  if(predict){
   return(list(preds=preds, lik.bits=lik.bits, ss.bits=ss.bits, ssDF=ssDF,
               data=YY, traps=scrFrame$traps))
  }
  }#END OF LIKELIHOOD FUNCTION



###############################
## SEX! Likelihood function

  msLL.sex <- function(pv, pn, YY, D, Y, nG, nK, hiK, dm.den, dm.trap){

  alpha0 <- array(0, c(ns, hiK, 2, 2))
  tmpP <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.(Intercept)", names.p0)]]
  if(pDot & !pTime){
    alpha0[,,,1] <- tmpP
  }
  if(pTime){
    tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
    for(s in 1:ns){
      alpha0[s, , 1, 1] <- tmpP + tmpT
      alpha0[s, , 2, 1] <- tmpP + tmpT
  }}
  if(pJustsex & !pTime & !pJustsesh & !pBothsexnsesh){
   tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
   for(s in 1:ns){
     alpha0[s, , 1, 1] <- tmpP
     alpha0[s, , 2, 1] <- tmpP + tmpS
  }}
  if(pJustsex & pTime & !pJustsesh & !pBothsexnsesh){
    tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
    tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
    for(s in 1:ns){
      alpha0[s, , 1, 1] <- tmpP + tmpT
      alpha0[s, , 2, 1] <- tmpP + tmpT + tmpS
  }}
  if(pJustsesh & !pTime & !pJustsex & !pBothsexnsesh){
    tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",names.p0)]])
    for(s in 1:ns){
      alpha0[s, , 1, 1] <- tmpP + tmpSS[s]
      alpha0[s, , 2, 1] <- tmpP + tmpSS[s]
    }}
    if(pJustsesh & pTime & !pJustsex & !pBothsexnsesh){
      tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
      tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",names.p0)]])
      for(s in 1:ns){
        alpha0[s, , 1, 1] <- tmpP + tmpT + tmpSS[s]
        alpha0[s, , 2, 1] <- tmpP + tmpT + tmpSS[s]
    }}
    if(pJustsesh & pJustsex & !pTime & !pBothsexnsesh){
      tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
      tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",names.p0)]])
      for(s in 1:ns){
        alpha0[s, , 1, 1] <- tmpP + tmpSS[s]
        alpha0[s, , 2, 1] <- tmpP + tmpSS[s] + tmpS
    }}
    if(pJustsesh & pJustsex & pTime & !pBothsexnsesh){
      tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
      tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
      tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",names.p0)]])
      for(s in 1:ns){
        alpha0[s, , 1, 1] <- tmpP + tmpT + tmpSS[s]
        alpha0[s, , 2, 1] <- tmpP + tmpT + tmpSS[s] + tmpS
    }}
    if(pBothsexnsesh & !pTime){
      tmpSSF <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.f.session",names.p0)]])
      tmpSSM <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.m.session", names.p0)]]
      for(k in 1:hiK){
        alpha0[, k, 1, 1] <- tmpP + tmpSSF
        alpha0[, k, 2, 1] <- tmpP + tmpSSM
    }}
    if(pBothsexnsesh & pTime){
      stop("model with time varying parameters AND a sex-session interaction not implemented")
      tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",names.p0)]])
      tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
      tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
      for(s in 1:ns){
        alpha0[s, , 1, 1] <- tmpP + tmpT + tmpSS[s]
        alpha0[s, , 2, 1] <- tmpP + tmpT + tmpSS[s] + tmpS
    }}
    BRmat <- array(0, c(ns, hiK, 2, 1))
    if(bDot)
      BRmat[,,,1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav",names.p0)]]

    if(bJustsex){
      BRmat[,,1,1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.f",names.p0)]]
      BRmat[,,2,1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.m",names.p0)]]
    }
    if(bJustsesh){
     for(k in 1:hiK){
       BRmat[,k,1,1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.session",names.p0)]]
       BRmat[, k, 2, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.session",names.p0)]]
    }}
    if(bBothsexnsesh){
     for(k in 1:hiK){
       BRmat[, k, 1, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.f.session",names.p0)]]
       BRmat[, k, 2, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.m.session",names.p0)]]
    }}
    alpha0[, , 1, 2] <- alpha0[, , 1, 1] + BRmat[, , 1, 1]
    alpha0[, , 2, 2] <- alpha0[, , 2, 1] + BRmat[, , 2, 1]

    alphsig <- matrix(0, ns, 2)
    tmpA <- pv[pn %in% names.sig[grep(fixed=TRUE,"sig.(Intercept)", names.sig)]]
    if(aDot)
      alphsig[, ] <- tmpA

    if(aJustsex & !aJustsesh){
      tmpSex <- pv[pn %in% names.sig[grep(fixed=TRUE,"sig.male", names.sig)]]
      alphsig[, 1] <- tmpA
      alphsig[, 2] <- tmpA + tmpSex
    }
    if(aJustsesh & !aJustsex){
      tmpSS <- c(0, pv[pn %in% names.sig[grep(fixed=TRUE,"sig.session",names.sig)]])
      for(s in 1:ns){
        alphsig[s, 1] <- tmpA + tmpSS[s]
        alphsig[s, 2] <- tmpA + tmpSS[s]
    }}

    if(aJustsesh & aJustsex){
      tmpSS <- c(0, pv[pn %in% names.sig[grep(fixed=TRUE,"sig.session",names.sig)]])
      tmpSex <- pv[pn %in% names.sig[grep(fixed=TRUE,"sig.male", names.sig)]]
      for(s in 1:ns){
        alphsig[s, 1] <- tmpA + tmpSS[s]
        alphsig[s, 2] <- tmpA + tmpSS[s] + tmpSex
    }}
    if(aBothsexnsesh){
      tmpSF <- c(0, pv[pn %in% names.sig[grep(fixed=TRUE,"sig.f.session",names.sig)]])
      tmpSM <- pv[pn %in% names.sig[grep(fixed=TRUE,"sig.m.session",names.sig)]]
      alphsig[, 1] <- tmpA + tmpSF
      alphsig[, 2] <- tmpA + tmpSM
    }
    alphsig <- 1/(2 * exp(alphsig)^2)

    if(trap.covs){
      t.beta <- matrix(NA, ns, length(t.nms))
      if(any(paste0("session:", t.nms) %in% attr(terms(model[[2]]),"term.labels"))){
       for(s in 1:ns){
        if(length(t.nms) > length(t.nms.sess)){
          t.beta[s, ] <- pv[pn %in% c(names.beta.trap[grep(fixed=TRUE, paste("session",
            s, sep = ""), names.beta.trap)], names.beta.trap[as.vector(unlist(
            sapply(t.nms.nosess, function(x){grep(fixed=TRUE,x, names.beta.trap)})))])]
        }else{
          t.beta[s, ] <- pv[pn %in% c(names.beta.trap[grep(fixed=TRUE,
                            paste("session",s, sep = ""), names.beta.trap)])]
       }}
      }else{
        for(s in 1:ns){
          t.beta[s, ] <- pv[pn %in% names.beta.trap]
    }}}
    if(DorN == "N"){
     if(dIPP){
       d.beta <- matrix(NA, ns, length(d.nms))
       if("session" %in% all.vars(model[[1]])){
        for(s in 1:ns){
          d.beta[s, ] <- pv[pn %in% names.beta.den[grep(fixed=TRUE,
                            paste("session",s, sep = ""), names.beta.den)]]
       }}else{
         for(s in 1:ns){
           d.beta[s, ] <- pv[pn %in% names.beta.den]
     }}}
    }else{
      d.beta <- pv[pn %in% names.beta.den]
    }
    if(distmet == "ecol")
      dist.beta <- pv[pn %in% names.dist]
    if(n0Session)
      n0 <- exp(pv[pn %in% names.n0])
    if(!n0Session)
      n0 <- rep(exp(pv[pn %in% names.n0]), ns)
    if(sexmod == "constant")
      psi.sex <- rep(plogis(pv[grep(fixed=TRUE,"psi", pn)]), ns)
    if(sexmod == "session")
      psi.sex <- plogis(pv[grep(fixed=TRUE,"psi", pn)])

    outLik <- 0
    if(predict){
      preds <- list()
      lik.bits <- list()
      ss.bits <- list()
    }
    for(s in 1:length(YY)){
      Ys <- YY[[s]]
      if(predict)
        preds[[s]] <- matrix(NA, nrow = nrow(Ys) + 1, ncol = nrow(ssDF[[s]]))
      if(!multicatch){
        zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
        Ys <- abind(Ys, zeros, along = 1)
      }
      if(multicatch){
        zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
        Ys <- abind(Ys, zeros, along = 1)
      }
      sx <- c(scrFrame$indCovs[[s]]$sex + 1, NA)
      if(distmet == "ecol"){
        cost <- exp(dm.cost[[s]] %*% exp(dist.beta))
        costR <- rasterFromXYZ(cbind(costDF[[s]][, c(1,2)], cost))
        if(is.null(PROJ)){
          projection(costR) <- "+proj=utm +zone=12 +datum=WGS84"
        }else{
          projection(costR) <- PROJ
        }
        trFn <- function(x) (1/(mean(x)))
        tr <- transition(costR, transitionFunction = trFn, direction = directions)
        trLayer <- geoCorrection(tr, scl = F)
        D[[s]] <- costDistance(trLayer,
                               as.matrix(scrFrame$traps[[s]][,c("X", "Y")]),
                               as.matrix(ssDF[[s]][, c("X","Y")]))
      }
      if(smallslow){
       if(distmet == "euc"){
         D <- list()
         D[[s]] <- e2dist(scrFrame$traps[[s]][, c("X","Y")], ssDF[[s]][, c("X", "Y")])
      }}
      lik.marg <- lik.marg1 <- lik.marg2 <- rep(NA, nrow(Ys))
      if(!is.null(trimS)){
        pixels.prior <- rep(T, nG[s])
        pixels.post <- apply(D[[s]], 2, min) <= trimS
        pixels <- (pixels.prior & pixels.post)
        pixels <- ifelse(pixels, 1, 0)
      }else{
        pixels <- rep(1, nG[s])
      }
      if(DorN == "N"){
       if(dIPP){
         d.s <- exp(dm.den[[s]] %*% d.beta[s, ])
         pi.s <- (d.s * pixels)/sum(d.s * pixels)
       }
       if(dHPP){
         pi.s <- pixels/(sum(pixels))
       }
      }else{
        d.s <- exp(dm.den[[s]] %*% d.beta)
        pi.s <- (d.s * pixels)/sum(d.s * pixels)
      }
      for(i in 1:nrow(Ys)){
      if(plotit){
        pp <- sum(trimR[[s]][[i]][[k]])
        plot(ssDF[[s]][,c("X", "Y")], pch=16, col="grey", cex=0.5, asp=1,
             main=paste("Session:",s," Individual: ",i," traps: ",sum(pp),sep=" "))
        points(ssDF[[s]][trimC[[s]][[i]],c("X", "Y")], pch=16, col=2, cex=mycex)
        points(scrFrame$traps[[s]][,c("X", "Y")], pch=15, col="grey", cex=1)
        points(scrFrame$traps[[s]][trimR[[s]][[i]][[k]], c("X")],
               scrFrame$traps[[s]][trimR[[s]][[i]][[k]],c("Y")], pch=15, col=1, cex=1)
        pp <- apply(Ys[i,,],1,sum)>0
        points(scrFrame$traps[[s]][pp,c("X", "Y")], pch=15, col=3, cex=mycex)
      }
      if(!is.na(sx[i])){
        lik.cond <- numeric(nG[s])
        for(k in 1:nK[s]){
          if(pBehave){
            a0 <- alpha0[s, k, sx[i], 1] * (1 - c(prevcap[[s]][i,, k])) +
                  alpha0[s, k, sx[i], 2] * c(prevcap[[s]][i,, k])
          }else{
            a0 <- rep(alpha0[s, k, sx[i], 1], nrow(D[[s]]))
          }
          if(trap.covs){
            a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
          }
          if(encmod == "B")
            probcap <- c(plogis(a0[trimR[[s]][[i]][[k]]])) * exp(-alphsig[s, sx[i]] *
                         D[[s]][trimR[[s]][[i]][[k]], trimC[[s]][[i]]]^2)
          if(encmod %in% c("P","CLOG"))
            probcap <- c(exp(a0[trimR[[s]][[i]][[k]]])) * exp(-alphsig[s, sx[i]] *
                         D[[s]][trimR[[s]][[i]][[k]], trimC[[s]][[i]]]^2)
          if(!multicatch){
           if(encmod == "B"){
             probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]][[k]], k],
               sum(trimC[[s]][[i]])),1, probcap[1:length(probcap)], log = TRUE))
           }
           if(encmod == "P"){
             probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]][[k]], k],
               sum(trimC[[s]][[i]])), probcap[1:length(probcap)], log = TRUE))
           }
           if(encmod == "CLOG"){
             probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]][[k]], k],
               sum(trimC[[s]][[i]])),1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
           }
          }else{
            probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
            probcap <- t(t(probcap)/colSums(probcap))
            vvv <- rep(c(Ys[i, trimR[[s]][[i]][[k]], k],
                         1 - any(Ys[i, trimR[[s]][[i]][[k]], k] > 0)),
                        sum(trimC[[s]][[i]]))
            vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv ==1])
            probcap[1:length(probcap)] <- vvv
          }
          if(is.matrix(probcap))
            lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + colSums(probcap)
          if(!is.matrix(probcap))
            lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + probcap
        }
        lik.cond[trimC[[s]][[i]]] <- exp(lik.cond[trimC[[s]][[i]]])
        tmpPsi <- (sx[i] == 1) * (1 - psi.sex[s]) + (sx[i] == 2) * psi.sex[s]
        lik.marg[i] <- sum(lik.cond * pi.s) * tmpPsi
        if(predict){
         preds[[s]][i, ] <- lik.cond * pi.s * tmpPsi/lik.marg[i]
        }
      }else{
        lik.cond1<- lik.cond2 <- numeric(nG[s])
        for(k in 1:nK[s]){
         if(pBehave){
           a0.1 <- alpha0[s, k, 1, 1] * (1 - c(prevcap[[s]][i,,k])) +
                   alpha0[s, k, 1, 2] * c(prevcap[[s]][i,,k])
           a0.2 <- alpha0[s, k, 2, 1] * (1 - c(prevcap[[s]][i,,k])) +
                   alpha0[s, k, 2, 2] * c(prevcap[[s]][i,,k])
         }else{
           a0.1 <- rep(alpha0[s, k, 1, 1], nrow(D[[s]]))
           a0.2 <- rep(alpha0[s, k, 2, 1], nrow(D[[s]]))
         }
         if(trap.covs){
           a0.1 <- a0.1 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
           a0.2 <- a0.2 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
         }
         if(encmod == "B")
           probcap <- c(plogis(a0.1[trimR[[s]][[i]][[k]]])) * exp(-alphsig[s, 1] *
             D[[s]][trimR[[s]][[i]][[k]], trimC[[s]][[i]]]^2)
         if(encmod %in% c("P","CLOG"))
           probcap <- c(exp(a0.1[trimR[[s]][[i]][[k]]])) * exp(-alphsig[s, 1] *
             D[[s]][trimR[[s]][[i]][[k]], trimC[[s]][[i]]]^2)
         if(!multicatch){
          if(encmod == "B"){
            probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]][[k]], k],
              sum(trimC[[s]][[i]])),1, probcap[1:length(probcap)], log = TRUE))
          }
          if(encmod == "P"){
            probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]][[k]], k],
              sum(trimC[[s]][[i]])), probcap[1:length(probcap)], log = TRUE))
          }
          if(encmod == "CLOG"){
            probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]][[k]], k],
             sum(trimC[[s]][[i]])),1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
          }
         }else{
           probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
           probcap <- t(t(probcap)/colSums(probcap))
           vvv <- rep(c(Ys[i,trimR[[s]][[i]][[k]], k],
                        1-any(Ys[i, trimR[[s]][[i]][[k]],k] > 0)),
                        sum(trimC[[s]][[i]]))
           vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv == 1])
           probcap[1:length(probcap)] <- vvv
         }
         if(is.matrix(probcap))
           lik.cond1[trimC[[s]][[i]]]<- lik.cond1[trimC[[s]][[i]]] + colSums(probcap)
         if(!is.matrix(probcap))
           lik.cond1[trimC[[s]][[i]]]<- lik.cond1[trimC[[s]][[i]]] + probcap

         if(encmod == "B")
           probcap <- c(plogis(a0.2[trimR[[s]][[i]][[k]]])) * exp(-alphsig[s, 2] *
             D[[s]][trimR[[s]][[i]][[k]], trimC[[s]][[i]]]^2)
         if(encmod %in% c("P","CLOG"))
           probcap <- c(exp(a0.2[trimR[[s]][[i]][[k]]])) * exp(-alphsig[s, 2] *
             D[[s]][trimR[[s]][[i]][[k]], trimC[[s]][[i]]]^2)
         if(!multicatch){
          if(encmod == "B"){
            probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]][[k]], k],
              sum(trimC[[s]][[i]])),1, probcap[1:length(probcap)], log = TRUE))
          }
          if(encmod == "P"){
            probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]][[k]], k],
              sum(trimC[[s]][[i]])), probcap[1:length(probcap)], log = TRUE))
          }
          if(encmod == "CLOG"){
            probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]][[k]], k],
             sum(trimC[[s]][[i]])),1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
          }
         }else{
           probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
           probcap <- t(t(probcap)/colSums(probcap))
           vvv <- rep(c(Ys[i,trimR[[s]][[i]][[k]], k],
                        1-any(Ys[i, trimR[[s]][[i]][[k]],k] > 0)),
                        sum(trimC[[s]][[i]]))
           vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv == 1])
           probcap[1:length(probcap)] <- vvv
         }
         if(is.matrix(probcap))
           lik.cond2[trimC[[s]][[i]]]<- lik.cond1[trimC[[s]][[i]]] + colSums(probcap)
         if(!is.matrix(probcap))
           lik.cond2[trimC[[s]][[i]]]<- lik.cond1[trimC[[s]][[i]]] + probcap
        }

        lik.cond1[trimC[[s]][[i]]] <- exp(lik.cond1[trimC[[s]][[i]]])
        lik.cond2[trimC[[s]][[i]]] <- exp(lik.cond2[trimC[[s]][[i]]])
        lik.marg1[i] <- sum(lik.cond1 * pi.s)
        lik.marg2[i] <- sum(lik.cond2 * pi.s)
        lik.marg[i] <- lik.marg1[i]*(1 - psi.sex[s]) + lik.marg2[i]*psi.sex[s]
        if(predict){
          lik.cond <- (lik.cond1 * (1 - psi.sex[s]) + lik.cond2 * psi.sex[s])
          preds[[s]][i, ] <- (lik.cond * pi.s)/lik.marg[i]
      }}}
      if(!predict){
       if(DorN == "N"){
         nv <- c(rep(1, length(lik.marg) - 1), n0[s])
         part1 <- lgamma((nrow(Ys) - 1) + n0[s] + 1) -  lgamma(n0[s] + 1)
         part2 <- sum(nv * log(lik.marg))
       }
       if(DorN == "D"){
         nv <- c(rep(1, length(lik.marg) - 1), 1)
         atheta <- 1 - lik.marg[nrow(Ys)]
         nind <- nrow(Ys) - 1
         part1 <- nind * log(sum(d.s * pixels)) - sum(d.s * pixels) * atheta
         part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
       }
       ll <- -1 * (part1 + part2)
       outLik <- outLik + ll
      }
      if(predict){
        lik.bits[[s]] <- cbind(lik.mar = lik.marg)
        ss.bits[[s]] <- cbind(pi.s, d.s, lik.cond)
        colnames(ss.bits[[s]]) <- c("pi.s", "d.s", "lik.cond")
      }
    }
    if(!predict){
      out <- outLik
      return(out)
    }
    if(predict){
      return(list(preds = preds, ss.bits = ss.bits, lik.bits = lik.bits,
                  ssDF = ssDF, data = YY, traps = scrFrame$traps))
    }
  }#END OF LIKELIHOOD FUNCTION

  if(getStarts == FALSE){
   if(!predict){
    message("Fitting model: D",paste(model)[1],", p0",paste(model)[2],", sigma",
             paste(model)[3],", asu",paste(model)[4], sep = " ")
    if(!anySex){
     message("Using 'no sex' likelihood.")
     message(Sys.time())
     message(paste(pn, " ", sep = " | "))
     message(" ")
     myfit <- suppressWarnings(nlm(msLL.nosex, p=pv,pn=pn, YY=YY, D=D, nG=nG,
               nK=nK, hiK=hiK, dm.den=dm.den, dm.trap=dm.trap, hessian=hessian,
               print.level=print.level, iterlim=200))
    }else{
     message("Using 'sex' likelihood.")
     message(Sys.time())
     message(paste(pn, " ", sep = " | "))
     message(" ")
     myfit <- suppressWarnings(nlm(msLL.sex, p=pv, pn=pn, YY=YY, D=D, nG=nG,
                nK=nK, hiK=hiK, dm.den=dm.den, dm.trap=dm.trap, hessian=hessian,
                print.level=print.level, iterlim=200))
    }
    links <- rep(NA, length(pn))
    pars <- myfit$estimate
    if(encmod == "B"){
      links[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- "(logit)"
    }else{
      links[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- "(log)"
    }
    links[grep(fixed=TRUE,"sig.(Intercept)", pn)] <- "(log)"
    links[grep(fixed=TRUE,"n0.", pn)] <- "(log)"
    links[grep(fixed=TRUE,"d0.(Intercept)", pn)] <- "(log)"
    links[grep(fixed=TRUE,"psi", pn)] <- "(logit)"
    links[grep(fixed=TRUE,"beta", pn)] <- "(Identity)"
    trans.mle <- rep(0, length(pv))
    if(encmod == "B"){
      trans.mle[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- plogis(pars[grep(fixed=TRUE,"p0.(Intercept)", pn)])
    }else{
      trans.mle[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- exp(pars[grep(fixed=TRUE,"p0.(Intercept)", pn)])
    }
    trans.mle[grep(fixed=TRUE,"sig.(Intercept)", pn)] <- exp(pars[grep(fixed=TRUE,"sig.(Intercept)", pn)])
    trans.mle[grep(fixed=TRUE,"n0.", pn)] <- exp(pars[grep(fixed=TRUE,"n0.", pn)])
    trans.mle[grep(fixed=TRUE,"d0.(Intercept)", pn)] <- exp(pars[grep(fixed=TRUE,"d0.(Intercept)", pn)])
    trans.mle[grep(fixed=TRUE,"psi", pn)] <- plogis(pars[grep(fixed=TRUE,"psi", pn)])
    trans.mle[grep(fixed=TRUE,"beta", pn)] <- pars[grep(fixed=TRUE,"beta", pn)]
    if(pBehave){
      links[grep(fixed=TRUE,"pBehav", pn)] <- "(Identity)"
      trans.mle[grep(fixed=TRUE,"pBehav", pn)] <- pars[grep(fixed=TRUE,"pBehav", pn)]
    }
    std.err <- rep(rep(NA, length(pv)))
    trans.se <- rep(NA, length(pv))
    if("hessian" %in% names(myfit)){
     if(sum(myfit$hessian) != 0){
       #Need a check for this error and return mles and a warning
       #Error in solve.default(myfit$hessian) :
       #Lapack routine dgesv: system is exactly singular: U[1,1] = 0
       std.err <- sqrt(diag(solve(myfit$hessian)))
     }else{
       warning("Something went wrong! Try better starting values.")
    }}

#    outStats <- data.frame(parameters = pn, link = links,
#                mle = myfit$estimate, std.er = std.err, mle.tr = trans.mle,
#                se.tr = trans.se)
    outStats <- matrix(cbind(pn,myfit$estimate,std.err,myfit$estimate/std.err,
                             2*pnorm(myfit$estimate/std.err)))
    rownames(outStats) <- pn
    colnames(outStats) <-c("Estimate","SE","z","P(>|z|)")
    VcV <- NULL
    if(DorN == "N"){
      ED <- (exp(pars[grep(fixed=TRUE,"n0.", pn)]) +
             unlist(lapply(scrFrame$caphist,nrow)))/areaS
    }else{
      ED <- NULL
    }
    endtime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    output <- list(call=cl, rawOutput=myfit, outStats=outStats,
                   coef.mle=data.frame(param=pn, mle=myfit$estimate),Area=areaS,
                   ED=ED, nObs=unlist(lapply(scrFrame$caphist,nrow)),mmdm=mmdm,
                   nll=myfit$minimum, AIC=2*myfit$minimum+2*length(myfit$estimate),
                   started=starttime, ended=endtime, proctime=(proc.time()-ptm)[3]/60)
    class(output) <- "oSCR.fit"
    return(output)
   }
   if(predict){
    message("Predicting form model: D", paste(model)[1], ", p0", paste(model)[2],
            ", sigma", paste(model)[3],", cost", paste(model)[4], sep = " ")
    if(!anySex){
     message("Using 'no sex' likelihood.")
     message(Sys.time())
     message(paste(pn, " ", sep = " | "))
     message(" ")
     myfit <- msLL.nosex(p = start.vals, pn = pn, YY = YY, D = D, hiK = hiK,
                         nG = nG, nK = nK, dm.den = dm.den, dm.trap = dm.trap)
    }else{
     message("Using 'sex' likelihood.")
     message(Sys.time())
     message(paste(pn, " ", sep = " | "))
     message(" ")
     myfit <- msLL.sex(p = start.vals, pn = pn, YY = YY, D = D, nK = nK, nG = nG,
                       hiK = hiK, dm.den = dm.den, dm.trap = dm.trap)
    }
    return(myfit)
   }
}}
