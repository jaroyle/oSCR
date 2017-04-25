oSCR.fit <-
function (model = list(D ~ 1, p0 ~ 1, sig ~ 1, asu ~1), scrFrame, ssDF = NULL, 
          costDF = NULL, rsfDF = NULL, distmet = c("euc", "user", "ecol")[1],
          sexmod = c("constant", "session")[1], encmod = c("B", "P", "CLOG")[1],
          DorN = c("D", "N")[1], directions = 8, Dmat = NULL, trimS = NULL, 
          start.vals = NULL, PROJ = NULL, pxArea = 1, plotit = F, mycex = 1, 
          tester = F, pl = 0, nlmgradtol = 1e-06, nlmstepmax = 10, predict = FALSE, 
          smallslow = FALSE, multicatch = FALSE, se = TRUE, print.level = 0, 
          getStarts = FALSE, theta = 2, RSF = FALSE, telemetry = c("none","ind","dep")[1]){
  
  ptm <- proc.time()
  starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
  my.model.matrix <- function(form, data) {
    cont.arg <- lapply(data.frame(data[,sapply(data.frame(data), is.factor)]), 
                       contrasts, contrasts = FALSE) 
    mdm <- suppressWarnings(model.matrix(form, data, contrasts.arg = cont.arg))
    return(mdm)
  }
  hessian <- ifelse(se, TRUE, FALSE)
  mmdm <- scrFrame$mmdm
  mmdm[is.na(mmdm)] <- mean(mmdm,na.rm=T)
  
  #ADD A CHECK FOR WHETHER TRIMS IS TOO SMALL
  #if((!is.null(trimS)) & (trimS < (0.6*mdm)))
  #  warning("The trimS value is smaller than half the max observed 
  #           distance moved and is probably too small.")
    
  cl <- match.call(expand.dots = TRUE)
  model.call <- as.list(paste(model))
    
  if (!require(abind)) stop("need to install package 'abind'")
  if (!require(Formula)) stop("need to load package 'Formula'")
  if (distmet == "ecol") {
    if (!require(raster)) stop("need to install package 'raster'")
    if (!require(gdistance)) stop("need to install package 'gdistance'")
  }
  if (!inherits(scrFrame, "scrFrame")) stop("Data must be of class 'scrFrame'")
  if ((encmod %in% c("B","CLOG")) & (max(unlist(lapply(scrFrame$caphist, max))) > 1)) {
    stop("Error: Data in caphist must be Binary")
  }
  if (distmet == "ecol" & is.null(PROJ))
    message("Projection not provided, using default: '+proj=utm +zone=12 +datum=WGS84'")
  if (!is.null(ssDF) & length(ssDF) != length(scrFrame$caphist))
    stop("Error: A 'state space' object must be provided for EACH session.")
  if (multicatch) {
    for (s in 1:length(scrFrame$caphist)) {
      captures <- apply(scrFrame$caphist[[s]], c(1, 3), sum)
      if (any(captures > 1)) 
        stop("Error: multicatch system cannot have > 1 capture.")
    }
  }
  if (predict & is.null(start.vals))
    stop("Starting values required to predict (hint: use estimated MLEs)")
  maxY <- unlist(lapply(scrFrame$caphist, max))
  if (any(maxY > 1) & encmod %in% c("B","CLOG"))
    stop("caphist must be binary when using the Binomial/Cloglog encounter model")
  if (all(maxY == 1) & encmod == "P")
    warning("caphist looks binary but Poisson encounter model is selected")
  if (theta >2 | theta <1)
    warning("theta should be between 1 (exponential) and 2 (half-normal) for 
            power model distance function")
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
  pars.sex <- NULL
  names.sex <- NULL
  D <- list()
  Drsf <- list()
  dm.den <- list()
  tmp.dm <- list()
  dm.trap <- list()
  dm.cost <- list()
  dm.rsf <- list()
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
  bDot <- FALSE
  bJustsex <- FALSE
  bJustsesh <- FALSE
  bBothsexnsesh <- FALSE
  telem <- FALSE
  warnings <- list()
  if (length(model) == 3) {
    model[[4]] <- formula(~1)
  }
  if((length(labels(terms(model[[4]])))>0) & distmet=="euc"){
    stop("Error: asu model specified but no 'dismet'. Use distmet=ecol.)")
  }
  for (i in 1:4) {
    model[[i]] <- update.formula(model[[i]], NULL ~ .)
  }
  if (is.null(ssDF)) {
    message("Generating a state space based on traps")
    dHPP <- TRUE
    ssDF <- make.ssDF(scrFrame, buffer, res)
  }
  if (RSF) {
    if (is.null(rsfDF))
      stop("Error: Cannot fit RSF without rsfDF!")
  }
  if (telemetry %in% c("ind","dep")) {
    if (is.null(scrFrame$telemetry)){
      stop("Error: No telemetry data in scrFrame!")
    }
    telem <- TRUE
    YYtel <- scrFrame$telemetry$fixfreq
    
    if (is.null(rsfDF)){
      rsfDF <- ssDF
    }
    if (ncol(YYtel[[1]]) != nrow(rsfDF[[1]]))
      stop("Error: Grid cells for telemetry fixes do not match rsfDF")
  }
  ns <- length(scrFrame$caphist)
  nt <- length(scrFrame$traps)
  nK <- unlist(lapply(scrFrame$caphist, function(x) dim(x)[3]))
  hiK <- max(nK)
  nG <- unlist(lapply(ssDF, nrow))
  nnn <- all(unlist(lapply(ssDF, function(x) {"session" %in% names(x)})))
  areaS <- NULL
  
  if ("session" %in% all.vars(model[[1]]) & (!nnn)) {
    for (s in 1:ns) {
      ssDF[[s]]$session <- factor(rep(s, nrow(ssDF[[s]])), levels = 1:ns)
    }
  }
  nnnn <- all(unlist(lapply(scrFrame$trapCovs[[1]][[1]], 
                            function(x) {"session" %in% names(x) })))
  if ("session" %in% all.vars(model[[2]]) & (!nnnn)) {
    for(s in 1:ns){
      for(m in 1:length(scrFrame$trapCovs[[s]])){
        sesf <- factor(rep(s, nrow(scrFrame$trapCovs[[s]][[m]])),levels = 1:ns)
        scrFrame$trapCovs[[s]][[m]]$session <- sesf 
      }
    }
  }
  
  allvars.D <- all.vars(model[[1]])
  dens.fx <- allvars.D[!allvars.D %in% c("D", "session")]
  allvars.T <- all.vars(model[[2]])
  trap.fx <- allvars.T[!allvars.T %in% c("p0", "session", "sex", "t", "T", "b")]
  allvars.sig <- all.vars(model[[3]])
  allvars.dist <- all.vars(model[[4]])
  var.p0.1 <- "sex" %in% allvars.T
  var.p0.2 <- "session" %in% allvars.T
  var.p0.3 <- "t" %in% allvars.T
  var.p0.4 <- any(c("sex:session", "session:sex") %in% 
                    attr(terms(model[[2]]), "term.labels"))
  var.sig.1 <- "sex" %in% allvars.sig
  var.sig.2 <- "session" %in% allvars.sig
  var.sig.3 <- any(c("sex:session", "session:sex") %in% 
                     attr(terms(model[[3]]), "term.labels"))
  var.b.1 <- "b" %in% attr(terms(model[[2]]), "term.labels")
  var.b.2 <- any(c("b:sex", "sex:b") %in% 
                   attr(terms(model[[2]]), "term.labels"))
  var.b.3 <- any(c("b:session", "session:b") %in% 
                   attr(terms(model[[2]]), "term.labels"))
  var.b.4 <- any(c("b:session:sex", "b:sex:session", "sex:session:b", 
                   "sex:b:session", "session:b:sex", "session:sex:b") %in% 
                   attr(terms(model[[2]]), "term.labels"))
  pBehave <- any(c(var.b.1, var.b.2, var.b.3, var.b.4))
  
  for (s in 1:ns) {
    if (!is.null(trimS)){
      pixels.prior <- rep(T, nG[s])
        pixels.post <- apply(e2dist(scrFrame$traps[[s]][ , c("X", "Y")], 
                                    ssDF[[s]][, c("X", "Y")]), 2, min) <= trimS
        pixels <- (pixels.prior & pixels.post)
        pixels <- ifelse(pixels, 1, 0)
    }
    else {
      pixels <- rep(1, nG[s])
    }
    areaS <- c(areaS, sum(pixels) * pxArea) #should be ssDF$res
  }
  #turn off RSF if no spatial trap covariates (i.e., p0~1)
  if (RSF & length(trap.fx) == 0){
    RSF <- FALSE
  }
  
  #trap covariates
  #can be altered to have session, sex, and b in the DM
  if (length(trap.fx) > 0) {
    trap.covs <- TRUE
    tcovnms <- colnames(scrFrame$trapCovs[[1]][[1]])
    tCovMissing <- trap.fx[which(!trap.fx %in% tcovnms)]
    if (length(tCovMissing) > 0) {
      stop("I cant find these covariates in 'scrFrame$trapCovs'",
      for (i in tCovMissing) print(i))
    }
    mod2 <- update(model[[2]], ~. - sex - session - t - b - b:sex - sex:b - 
                   b:session - session:b - sex:session - session:sex - 
                   b:session:sex - b:sex:session - sex:session:b - sex:b:session - 
                   session:b:sex - session:sex:b)
    if (any(c("session") %in% allvars.T)) tSession <- TRUE
    for (s in 1:ns) {
      tmp.dm <- list()
      for (k in 1:nK[s]) {
        #why the -1 here?
        tmp.dm[[k]] <- model.matrix(mod2, scrFrame$trapCovs[[s]][[k]])[,-1,drop=FALSE]
        if (s == 1 && k == 1) t.nms <- colnames(tmp.dm[[k]])
        if (nrow(tmp.dm[[k]]) != nrow(scrFrame$trapCovs[[s]][[k]])) {
          mis <- setdiff(rownames(scrFrame$trapCovs[[s]][[k]]),
                         rownames(my.model.matrix(mod2, scrFrame$trapCovs[[s]][[k]])))
          tmp.insert <- matrix(NA, length(mis), ncol(tmp.dm[[k]]))
          row.names(tmp.insert) <- mis
          tmp.dm[[k]] <- rbind(tmp.dm[[k]], tmp.insert)
          tmp.dm[[k]] <- tmp.dm[[k]][order(as.numeric(row.names(tmp.dm[[k]]))), ]
        }
      }
      dm.trap[[s]] <- tmp.dm
      if (RSF) {
        if (any(!tcovnms %in% names(rsfDF[[s]]))){
          rsfMissing <- tcovnms[which(!tcovnms %in% names(rsfDF[[s]]))]
          for (miss in 1:length(rsfMissing)) {
            rsfDF[[s]][,rsfMissing[miss]] <- 0
          }
        }
        dm.rsf[[s]] <- my.model.matrix(mod2, rsfDF[[s]])[ , -1, drop=FALSE]
      }
    }
    if (any(paste0("session:", t.nms) %in% 
            attr(terms(model[[2]]), "term.labels"))) {
      t.nms.sess <- t.nms[which(paste0("session:", t.nms) %in% 
                                attr(terms(model[[2]]), "term.labels"))]
      tmpTsess <- rep(1:ns, each = length(t.nms.sess))
      tmpTcovs <- rep(t.nms.sess, ns)
      names.beta.trap1 <- paste("t.beta.", tmpTcovs, ".session", tmpTsess, sep = "")
      if (length(t.nms) > length(t.nms.sess)) {
        t.nms.nosess <- t.nms[!t.nms %in% t.nms.sess]
        names.beta.trap <- c(names.beta.trap1, paste("t.beta.", t.nms.nosess, sep = ""))
      }
      else {
        names.beta.trap <- names.beta.trap1
      }
      pars.beta.trap <- rep(0,length(names.beta.trap))
    }
    else {
      names.beta.trap <- paste("t.beta.", t.nms, sep = "")
      pars.beta.trap <- rep(0,length(names.beta.trap))
    }
  }

  
  #clean formatting to here

  if (length(dens.fx) > 0) {
        dcovnms <- colnames(ssDF[[1]])
        dCovMissing <- dens.fx[which(!dens.fx %in% dcovnms)]
        if (length(dCovMissing) > 0) {
            stop("I can't find these covariates in 'ssDF'", for (i in dCovMissing) print(i))
        }
    }
    if (DorN == "N") {
        if (length(dens.fx) > 0) {
            dIPP <- TRUE
            mod1 <- update(model[[1]], ~. - sex - session)
            for (s in 1:ns) {
                dm.den[[s]] <- model.matrix(mod1, ssDF[[s]])[,-1,drop=FALSE]
                if (s == 1)
                  d.nms <- colnames(dm.den[[s]])
            }
            if("Session" %in% all.vars(model[[1]])) {
                tmpDsess <- rep(1:ns, each = length(d.nms))
                tmpDcovs <- rep(d.nms, ns)
                names.beta.den <- paste("d.beta.", tmpDcovs,
                  ".session", tmpDsess, sep = "")
                pars.beta.den <- rep(0,length(names.beta.den))
                names.n0 <- paste0("n0.session", 1:ns)
                pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
            }
            if("session" %in% all.vars(model[[1]])) {
                names.beta.den <- paste0("d.beta.", d.nms, sep = "")
                pars.beta.den <- rep(0,length(names.beta.den))
                names.n0 <- paste("n0.session", 1:ns)
                pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
            }
            if(!("session" %in% all.vars(model[[1]])) & !("Session" %in%
                all.vars(model[[1]]))) {
                names.beta.den <- paste0("d.beta.", d.nms, sep = "")
                pars.beta.den <- rep(0,length(names.beta.den))
                names.n0 <- paste0("n0.")
                pars.n0 <- log(mean(unlist(lapply(scrFrame$caphist,nrow))))
            }
        }
        else {
            dHPP <- TRUE
            names.beta.den <- NULL
            pars.beta.den <- NULL
            for (s in 1:ns) {
                dm.den[[s]] <- my.model.matrix(~1, ssDF[[s]])
            }
            if ("Session" %in% all.vars(model[[1]])) {
                n0Session <- TRUE
                names.n0 <- paste0("n0.session", 1:ns)
                pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
            }
            if ("session" %in% all.vars(model[[1]])) {
                n0Session <- TRUE
                names.n0 <- paste0("n0.session", 1:ns)
                pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
            }
            if (!("session" %in% all.vars(model[[1]])) & !("Session" %in%
                all.vars(model[[1]]))) {
                names.n0 <- paste0("n0.")
                pars.n0 <- log(mean(unlist(lapply(scrFrame$caphist,nrow))))
            }
        }
    }
    if (DorN == "D") {
        mod1 <- update(model[[1]], ~. - sex)
        for (s in 1:ns) {
          dm.den[[s]] <- model.matrix(mod1, as.data.frame(ssDF[[s]]))
        }
        d.nms <- colnames(dm.den[[1]])
        names.beta.den <- paste("d.beta", d.nms, sep = ".")
        chx <- grep(fixed=TRUE,"Intercept", names.beta.den)
        if (length(chx) > 0)
            names.beta.den[chx] <- "d0.(Intercept)"
        pars.d0 <- log(mean((unlist(lapply(scrFrame$caphist,nrow)))/unlist(lapply(ssDF, nrow))))
        pars.beta.den <- c(pars.d0, rep(0.1, length(names.beta.den) - 1))
        pars.n0 <- NULL
        names.n0 <- NULL
    }

    if (distmet == "ecol" && length(allvars.dist) == 0) {
        message("You specified 'ecological distance' (distmet='ecol') but provided no\ncost surface.\n    Euclidean distance will be used.")
    }
    if (length(allvars.dist) > 0) {
        ccovnms <- colnames(costDF[[1]])
        cCovMissing <- allvars.dist[which(!allvars.dist %in% ccovnms)]
        if (length(cCovMissing) > 0) {
            stop("I cant find theses covariates in 'costDF'",
                for (i in cCovMissing) print(i))
        }
    }

    if(distmet == "ecol"){
      mod4 <- update(model[[4]], ~. - sex - session)
      for (s in 1:ns) {
        dm.cost[[s]] <- model.matrix(mod4, as.data.frame(costDF[[s]]))
      }
      c.nms <- colnames(dm.cost[[1]])
      names.dist <- paste("c.beta", c.nms, sep = ".")
      chx <- grep(fixed=TRUE,"Intercept", names.dist)
      if (length(chx) > 0)
        names.dist[chx] <- "c0.(Intercept)"
      pars.c0 <- 0.01
      pars.dist <- c(pars.c0, rep(0.01, length(names.dist) - 1))
    }
    

    if (!smallslow) {
        if (distmet == "euc") {
            for (s in 1:ns) {
                D[[s]] <- e2dist(scrFrame$traps[[s]][, c("X",
                  "Y")], ssDF[[s]][, c("X", "Y")])
            }
        }
    }
    if ("indCovs" %in% names(scrFrame)) {
        if ("sex" %in% names(scrFrame$indCovs[[1]])) {
            anySex <- TRUE
        }
    }
    tmp.p0.names <- "p0.(Intercept)"
    if (sum(var.p0.1, var.p0.2, var.p0.3, var.p0.4) == 0) {
        tmp.p0.names <- "p0.(Intercept)"
        pDot <- TRUE
    }
    if (var.p0.1 & !var.p0.2) {
        tmp.p0.names <- c("p0.(Intercept)", "p0.male")
        pJustsex <- TRUE
    }
    if (!var.p0.1 & var.p0.2) {
        if (ns > 1) {
            tmp.p0.names <- c("p0.(Intercept)", paste0("p0.session", 2:ns))
            pJustsesh <- TRUE
        }
        else {
            tmp.p0.names <- "p0.(Intercept)"
            pDot <- TRUE
        }
    }
    if (var.p0.1 & var.p0.2) {
        tmp.p0.names <- c("p0.(Intercept)", "p0.male")
        if (ns > 1) {
            tmp.p0.names <- c(tmp.p0.names, paste0("p0.session",
                2:ns))
            pJustsesh <- TRUE
        }
        else {
            tmp.p0.names <- tmp.p0.names
            pJustsex <- TRUE
        }
        pJustsex <- TRUE
    }
    if (var.p0.3) {
        tmp.p0.names <- c(tmp.p0.names, paste0("p0.t", 2:hiK))
        pTime <- TRUE
    }
    if (var.p0.4) {
        if (ns > 1) {
            tmp.p0.names <- c("p0.(Intercept)", paste0("p0.f.session", 2:ns),
                paste0("p0.m.session", 1:ns))
            pBothsexnsesh <- TRUE
        }
        else {
            tmp.p0.names <- c("p0.(Intercept)", "p0.male")
            pJustsex <- TRUE
        }
    }
    names.p0 <- tmp.p0.names
    pars.p0 <- rep(0, length(names.p0))
    
    st.p0 <- qlogis(0.1*(sum(sapply(scrFrame$caphist,sum))/
             sum(sapply(scrFrame$caphist,
                        function(x) prod(dim(x)[c(1,3)])))))
    pars.p0[1] <- st.p0
    if (any(var.p0.1, var.p0.1, var.sig.1) && !anySex)
        stop("Sex defined in a model but no sex data provided.")
    if (var.b.1 & !var.b.2 & !var.b.3 & !var.b.4) {
        pars.p0 <- c(pars.p0, 0)
        names.p0 <- c(names.p0, "p.behav")
        bDot <- TRUE
    }
    if (var.b.2 & !var.b.4) {
        pars.p0 <- c(pars.p0, 0, 0)
        names.p0 <- c(names.p0, "p.behav.f", "p.behav.m")
        bJustsex <- TRUE
    }
    if (var.b.3 & !var.b.4) {
        pars.p0 <- c(pars.p0, rep(0, ns))
        names.p0 <- c(names.p0, paste0("p.behav.session", 1:ns))
        bJustsesh <- TRUE
    }
    if (var.b.4) {
        pars.p0 <- c(pars.p0, rep(0, 2 * ns))
        names.p0 <- c(names.p0, paste0("p.behav.f.session", 1:ns),
            paste0("p.behav.m.session", 1:ns))
        bBothsexnsesh <- TRUE
    }

    #Setup: sigma
    pars.sig <- NULL
    names.sig <- NULL
    var.sig.2 <- "session" %in% allvars.sig
    var.sig.3 <- any(c("sex:session", "session:sex") %in% 
                     attr(terms(model[[3]]), "term.labels"))
    
    tmp.mm <- model.matrix(model[[3]], scrFrame$sigCovs)
    names.sig <- paste("sig.",colnames(tmp.mm),sep="")
    pars.sig <- rep(0.1, length(names.sig))
    pars.sig[1] <- log(0.5 * mmdm)

    if (anySex) {
        if (sexmod == "constant") {
            pars.sex <- 0
            names.sex <- "psi.constant"
        }
        if (sexmod == "session") {
            pars.sex <- rep(0, ns)
            names.sex <- paste("psi", 1:ns, sep = "")
        }
    }
    else {
        pars.sex <- NULL
        names.sex <- NULL
    }
    
    #create starting values objects (pv: values, pn: names) 
    pv <- round(c(pars.p0, pars.sig, pars.beta.trap, pars.beta.den, pars.dist, 
                  pars.n0, pars.sex), 2)
    pn <- c(names.p0, names.sig, names.beta.trap, names.beta.den, names.dist, 
            names.n0, names.sex)
    if (!is.null(start.vals)) {
      if (length(pv) == length(start.vals)) {
        pv <- start.vals
      }
      else {
        message(
          "The number of starting values provided doesnt match the \n\n 
           number of parameters in the model. Randomly generated values \n\n 
           are being used. Use getStarts = T to get correct length.")
      }
    }
    if (getStarts == TRUE) {
      oSCR.start <- list(parameters = pn, values = pv)
      return(oSCR.start)
    }
    
    #create the prevcap objects
    if (pBehave) {
      prevcap <- do.prevcap(scrFrame)
    }

    #Do the trimming
    get.trims <- do.trim(scrFrame, ssDF, trimS)
    trimR <- get.trims$trimR
    trimC <- get.trims$trimC
    
    msLL.nosex <- function(pv = pv, pn = pn, D = D, hiK = hiK, nG = nG, 
                           nK = nK, dm.den = dm.den, dm.trap = dm.trap) {
        
      alpha0 <- array(0, dim = c(ns, hiK, 2))
      tmpP <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.(Intercept)", names.p0)]]
        if (pDot & !pTime) {
            alpha0[, , ] <- tmpP
        }
        if (pDot & pTime) {
            tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1] <- tmpP + tmpT
            }
        }
        if (pJustsesh & !pTime) {
            tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",
                names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1] <- tmpP + tmpSS[s]
            }
        }
        if (pJustsesh & pTime) {
            tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
            tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",
                names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1] <- tmpP + tmpSS[s] + tmpT
            }
        }
        BRmat <- array(0, c(ns, hiK, 1))
        if (bDot) {
            BRmat[, , 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav", names.p0)]]
        }
        if (bJustsesh) {
            for (k in 1:hiK) {
                BRmat[, k, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.session", names.p0)]]
            }
        }
        alpha0[, , 2] <- alpha0[, , 1] + BRmat[, , 1]

        #sigma
        sig.beta <- pv[grep("sig.",pn)]
        alphsig <- model.matrix(model[[3]],scrFrame$sigCovs) %*% sig.beta
        alphsig <- 1/(2 * exp(alphsig)^2)
        
        if (trap.covs) {
            t.beta <- matrix(NA, ns, length(t.nms))
            if (any(paste0("session:", t.nms) %in% attr(terms(model[[2]]),
                "term.labels"))) {
                for (s in 1:ns) {
                  if (length(t.nms) > length(t.nms.sess)) {
                    t.beta[s, ] <- pv[pn %in% c(names.beta.trap[grep(fixed=TRUE,paste("session",
                      s, sep = ""), names.beta.trap)], names.beta.trap[as.vector(unlist(sapply(t.nms.nosess,
                      function(x) {
                        grep(fixed=TRUE,x, names.beta.trap)
                      })))])]
                  }
                  else {
                    t.beta[s, ] <- pv[pn %in% c(names.beta.trap[grep(fixed=TRUE,paste("session",
                      s, sep = ""), names.beta.trap)])]
                  }
                }
            }
            else {
                for (s in 1:ns) {
                  t.beta[s, ] <- pv[pn %in% names.beta.trap]
                }
            }
        }
        if (DorN == "N") {
            if (dIPP) {
                d.beta <- matrix(NA, ns, length(d.nms))
                if ("session" %in% all.vars(model[[1]])) {
                  for (s in 1:ns) {
                    d.beta[s, ] <- pv[pn %in% names.beta.den[grep(fixed=TRUE,paste("session",
                      s, sep = ""), names.beta.den)]]
                  }
                }
                else {
                  for (s in 1:ns) {
                    d.beta[s, ] <- pv[pn %in% names.beta.den]
                  }
                }
            }
        }
        else {
            d.beta <- pv[pn %in% names.beta.den]
        }
        if (distmet == "ecol") {
            dist.beta <- pv[pn %in% names.dist]
        }
        if (n0Session)
            n0 <- exp(pv[pn %in% names.n0])
        if (!n0Session)
            n0 <- rep(exp(pv[pn %in% names.n0]), ns)
        outLik <- 0
        if (predict) {
            preds <- list()
            lik.bits <- list()
            ss.bits <- list()
        }
        for (s in 1:length(scrFrame$caphist)) {
            Ys <- scrFrame$caphist[[s]]
            if (telem){ #check if telemetry exists
              Ytels <- YYtel[[s]]
              cap.tel <- scrFrame$telemetry$cap.tel[[s]]  #index of captured ind w/ collars
              lik.marg.tel <- rep(NA, nrow(Ytels))
            } 
            if (predict)
                preds[[s]] <- matrix(NA, nrow = nrow(Ys) + 1,
                  ncol = nrow(ssDF[[s]]))

            zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
            Ys <- abind(Ys, zeros, along = 1)

            if (distmet == "ecol") {
                cost <- exp(dm.cost[[s]] %*% dist.beta)
                costR <- rasterFromXYZ(cbind(costDF[[s]][, c(1,
                  2)], cost))
                if (is.null(PROJ)) {
                  projection(costR) <- "+proj=utm +zone=12 +datum=WGS84"
                }
                else {
                  projection(costR) <- PROJ
                }
                tr <- transition(costR, transitionFunction = function(x) (1/(mean(x))),
                  direction = directions)
                trLayer <- geoCorrection(tr, scl = F)
                D[[s]] <- costDistance(trLayer,
                                       as.matrix(scrFrame$traps[[s]][,c("X", "Y")]), 
                                       as.matrix(ssDF[[s]][, c("X","Y")]))
            }
            if (smallslow) {
                if (distmet == "euc") {
                  D <- list()
                  D[[s]] <- e2dist(scrFrame$traps[[s]][, c("X",
                    "Y")], ssDF[[s]][, c("X", "Y")])
                }
            }
            if (telem){
              # only Euclidean distance for telemetry fixes
              Drsf[[s]] <- e2dist(rsfDF[[s]][, c("X", "Y")], rsfDF[[s]][, c("X", "Y")])
            }
            lik.marg <- rep(NA, nrow(Ys))
            if (!is.null(trimS)) {
                pixels.prior <- rep(T, nG[s])
                pixels.post <- apply(D[[s]], 2, min) <= trimS
                pixels <- (pixels.prior & pixels.post)
                pixels <- ifelse(pixels, 1, 0)
            }
            else {
                pixels <- rep(1, nG[s])
            }
            if (DorN == "N") {
                if (dIPP) {
                  d.s <- exp(dm.den[[s]] %*% d.beta[s, ])
                  pi.s <- (d.s * pixels)/sum(d.s * pixels)
                }
                if (dHPP) {
                  pi.s <- pixels/(sum(pixels))
                }
            }
            else {
                d.s <- exp(dm.den[[s]] %*% d.beta)
                pi.s <- (d.s * pixels)/sum(d.s * pixels)
            }
            
            # some collared ind captured, so keep lik.cond for combining later
            if (telemetry == "dep"){
              lik.cond.tel <- matrix(0,nrow=length(cap.tel),ncol=nG[s])
            }
            
            Kern <- exp(-alphsig[s] * D[[s]]^theta)

            for (i in 1:nrow(Ys)) {
              if (plotit) {
                pp <- sum(trimR[[s]][[i]][[k]])
                plot(ssDF[[s]][, c("X", "Y")], pch = 16, col = "grey", cex = 0.5, asp = 1,
                     main = paste("Session:", s, " Individual: ", i, " traps: ", sum(pp), sep = " "))
                points(ssDF[[s]][trimC[[s]][[i]], c("X", "Y")], pch = 16, col = 2, cex = mycex)
                points(scrFrame$traps[[s]][,c("X", "Y")], pch = 15, col = "grey", cex = 1)
                points(scrFrame$traps[[s]][trimR[[s]][[i]][[k]], c("X")],
                       scrFrame$traps[[s]][trimR[[s]][[i]][[k]],c("Y")],
                       pch = 15, col = 1, cex = 1)
                pp <- apply(Ys[i,,],1,sum)>0
                points(scrFrame$traps[[s]][pp,c("X", "Y")], pch = 15, col = 3, cex = mycex)
              }


                 lik.cond <- numeric(nG[s])
#               if (multicatch)
#                  Pm <- matrix(0, sum(trimR[[s]][[i]][[k]]) + 1, sum(trimC[[s]][[i]]))
#                if (!multicatch)
#                  Pm <- matrix(0, sum(trimR[[s]][[i]][[k]]), sum(trimC[[s]][[i]]))

           for (k in 1:nK[s]) {
            
            if(!("removed" %in% names(scrFrame$indCovs[[s]]))){
              dead <- 1
            }else{
              if(i < nrow(Ys)){
                dead <- ifelse(k > scrFrame$indCovs[[s]]$removed[i],0,1)
              }else{
                dead <- 1
              }
            }
            
 
                   if (pBehave) {
                    a0 <- alpha0[s,k,1] * (1 - c(prevcap[[s]][i,,k])) +
                          alpha0[s,k,2] * c(prevcap[[s]][i,,k])
                  }
                  else {
                    a0 <- rep(alpha0[s, k, 1], nrow(D[[s]]))
                  }
                  if (trap.covs) {
                    a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
                  }
                  if (encmod == "B")
                    probcap <- c(dead * plogis(a0[trimR[[s]][[i]][[k]]])) *
                      Kern[trimR[[s]][[i]][[k]], trimC[[s]][[i]] ]
                  if (encmod %in% c("P","CLOG"))
                    probcap <- c(dead * exp(a0[trimR[[s]][[i]][[k]]])) *
                      Kern[trimR[[s]][[i]][[k]], trimC[[s]][[i]]]
                  if (!multicatch){
                    if (encmod == "B") {
                      probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                        trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                        1, probcap[1:length(probcap)], log = TRUE))
                    }
                    if (encmod == "P") {
                      probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,
                        trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                        probcap[1:length(probcap)], log = TRUE))
                    }
                    if (encmod == "CLOG") {
                      probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                        trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                        1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
                    }
                  }
                  else {
                    probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
                    probcap <- t(t(probcap)/colSums(probcap))
                    vvv <- rep(c(Ys[i, trimR[[s]][[i]][[k]], k], 1 -
                      any(Ys[i, trimR[[s]][[i]][[k]], k] > 0)), sum(trimC[[s]][[i]]))
                    vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv ==
                      1])
                    probcap[1:length(probcap)] <- vvv
                  }
            #      if (!is.null(scrFrame$trapOperation)) {
            #        probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]],  k]
            #      }
            #cat("dim probca: ", dim(probcap),fill=TRUE)
            #cat("length lik.cond: ", length(lik.cond),fill=TRUE)
            #cat("sum trimC: ", sum(trimC[[s]][[i]]), fill=TRUE)

    ####        lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + colSums(probcap)


          if(is.matrix(probcap))
               lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + colSums(probcap)
          if(!is.matrix(probcap))
               lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + probcap
#            if(i==nrow(Ys) & k==1){
#                b<<- probcap
#                return(0)
#                }
            #Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
                }  # end k loop
#                lik.cond <- numeric(nG[s])
#    lik.cond[trimC[[s]][[i]]] <- exp(colSums(Pm, na.rm = T))
                 
                 if (telemetry == "dep"){
                   if (i %in% cap.tel){
                     lik.cond.tel[match(i,cap.tel),] <- lik.cond
                   }
                 }
                 
                 
                lik.cond[trimC[[s]][[i]]] <- exp(lik.cond[trimC[[s]][[i]]])  ####colSums(Pm, na.rm = T))
              #  b<<- lik.cond
              #  return(0)
                lik.marg[i] <- sum(lik.cond * pi.s)

                if (predict)
                  preds[[s]][i, ] <- (lik.cond * pi.s)/lik.marg[i]
            }
            
            
            if(telem){
              
              if(RSF){
                rsf.lam0 <- dm.rsf[[s]] %*% c(t.beta[s,])
                rsf.lam0 <- array(rsf.lam0,dim=c(nrow(rsfDF[[s]]),nrow(rsfDF[[s]])))
              } else {
                rsf.lam0 <- 0
              }
              for (i in 1:nrow(Ytels)){
                
                probs <- t(exp(rsf.lam0 - alphsig[s] * Drsf[[s]]^theta))
                denom <- rowSums(probs)
                probs <- t(probs/denom)
                
                lik.marg.tel[i] <- sum( exp(Ytels[i,,drop=F] %*% log(probs)) * as.vector(pi.s) )
                #browser()
                if (telemetry == "dep"){
                  if (i <= length(cap.tel)){
                    # combine conditional likelihoods if some collared ind were captured
                    lik.cond.tot <- (Ytels[i,,drop=F] %*% log(probs)) + lik.cond.tel[i,]
                    #lik.cond.tot[trimC[[s]][[cap.tel[i]]]] <- exp(lik.cond.tot[trimC[[s]][[cap.tel[i]]]])
                    lik.cond.tot[lik.cond.tot != 0] <- exp(lik.cond.tot[lik.cond.tot != 0])
                    
                    # fix marginal likelihoods
                    lik.marg[cap.tel[i]] <- sum(lik.cond.tot * as.vector(pi.s)) 
                    lik.marg.tel[i] <- 1
                    
                    if (predict){
                      preds[[s]][cap.tel[i], ] <- lik.cond.tot * as.vector(pi.s) / lik.marg[cap.tel[i]]
                    }
                    
                  }
                }
                
              }
            }
            
            if (!predict) {
                if (DorN == "N") {
                  nv <- c(rep(1, length(lik.marg) - 1), n0[s])
                  part1 <- lgamma((nrow(Ys) - 1) + n0[s] + 1) -
                    lgamma(n0[s] + 1)
                  part2 <- sum(nv * log(lik.marg))
                }
                if (DorN == "D") {
                  nv <- c(rep(1, length(lik.marg) - 1), 1)
                  atheta <- 1 - lik.marg[nrow(Ys)]
                  nind <- nrow(Ys) - 1
                  part1 <- nind * log(sum(d.s * pixels)) - sum(d.s *
                    pixels) * atheta
                  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
                }
              
              if(telem){
                part3 <- sum(log(lik.marg.tel))
              } else {
                part3 <- 0
              }
              
              ll <- -1 * (part1 + part2 + part3)
              outLik <- outLik + ll
            }
            if (predict) {
                lik.bits[[s]] <- cbind(lik.mar = lik.marg)
                ss.bits[[s]] <- cbind(pi.s, d.s, lik.cond)
                colnames(ss.bits[[s]]) <- c("pi.s", "d.s", "lik.cond")
            }
        }
        if (!predict) {
            out <- outLik
            return(out)
        }
        if (predict) {
            return(list(preds = preds, lik.bits = lik.bits, ss.bits = ss.bits,
                ssDF = ssDF, data = scrFrame$caphist, traps = scrFrame$traps))
        }
    }
    msLL.sex <- function(pv, pn, scrFrame, D, Y, nG, nK, hiK, dm.den, dm.trap) {
        alpha0 <- array(0, c(ns, hiK, 2, 2))
        tmpP <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.(Intercept)", names.p0)]]
        if (pDot & !pTime) {
            alpha0[, , , 1] <- tmpP
        }
        if (pTime) {
            tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP + tmpT
                alpha0[s, , 2, 1] <- tmpP + tmpT
            }
        }
        if (pJustsex & !pTime & !pJustsesh & !pBothsexnsesh) {
            tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP
                alpha0[s, , 2, 1] <- tmpP + tmpS
            }
        }
        if (pJustsex & pTime & !pJustsesh & !pBothsexnsesh) {
            tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
            tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP + tmpT
                alpha0[s, , 2, 1] <- tmpP + tmpT + tmpS
            }
        }
        if (pJustsesh & !pTime & !pJustsex & !pBothsexnsesh) {
            tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",
                names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP + tmpSS[s]
                alpha0[s, , 2, 1] <- tmpP + tmpSS[s]
            }
        }
        if (pJustsesh & pTime & !pJustsex & !pBothsexnsesh) {
            tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
            tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",
                names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP + tmpT + tmpSS[s]
                alpha0[s, , 2, 1] <- tmpP + tmpT + tmpSS[s]
            }
        }
        if (pJustsesh & pJustsex & !pTime & !pBothsexnsesh) {
            tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
            tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",
                names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP + tmpSS[s]
                alpha0[s, , 2, 1] <- tmpP + tmpSS[s] + tmpS
            }
        }
        if (pJustsesh & pJustsex & pTime & !pBothsexnsesh) {
            tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
            tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
            tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",
                names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP + tmpT + tmpSS[s]
                alpha0[s, , 2, 1] <- tmpP + tmpT + tmpSS[s] +
                  tmpS
            }
        }
        if (pBothsexnsesh & !pTime) {
            tmpSSF <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.f.session",
                names.p0)]])
            tmpSSM <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.m.session", names.p0)]]
            for (k in 1:hiK) {
                alpha0[, k, 1, 1] <- tmpP + tmpSSF
                alpha0[, k, 2, 1] <- tmpP + tmpSSM
            }
        }
        if (pBothsexnsesh & pTime) {
            stop("model with time varying parameters AND a sex-session interaction not implemented")
            tmpSS <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.session",
                names.p0)]])
            tmpS <- pv[pn %in% names.p0[grep(fixed=TRUE,"p0.male", names.p0)]]
            tmpT <- c(0, pv[pn %in% names.p0[grep(fixed=TRUE,"p0.t", names.p0)]])
            for (s in 1:ns) {
                alpha0[s, , 1, 1] <- tmpP + tmpT + tmpSS[s]
                alpha0[s, , 2, 1] <- tmpP + tmpT + tmpSS[s] +
                  tmpS
            }
        }
        BRmat <- array(0, c(ns, hiK, 2, 1))
        if (bDot) {
            BRmat[, , , 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav",
                names.p0)]]
        }
        if (bJustsex) {
            BRmat[, , 1, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.f",
                names.p0)]]
            BRmat[, , 2, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.m",
                names.p0)]]
        }
        if (bJustsesh) {
            for (k in 1:hiK) {
                BRmat[, k, 1, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.session",
                  names.p0)]]
                BRmat[, k, 2, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.session",
                  names.p0)]]
            }
        }
        if (bBothsexnsesh) {
            for (k in 1:hiK) {
                BRmat[, k, 1, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.f.session",
                  names.p0)]]
                BRmat[, k, 2, 1] <- pv[pn %in% names.p0[grep(fixed=TRUE,"p.behav.m.session",
                  names.p0)]]
            }
        }
        alpha0[, , 1, 2] <- alpha0[, , 1, 1] + BRmat[, , 1, 1]
        alpha0[, , 2, 2] <- alpha0[, , 2, 1] + BRmat[, , 2, 1]

        #sigma
        sig.beta <- pv[grep("sig.",pn)]
        alphsig <- model.matrix(model[[3]],scrFrame$sigCovs) %*% sig.beta
        alphsig <- 1/(2 * exp(alphsig)^2)
        alphsig <- matrix(alphsig, ns, 2, byrow=FALSE)

        if (trap.covs) {
            t.beta <- matrix(NA, ns, length(t.nms))
            if (any(paste0("session:", t.nms) %in% attr(terms(model[[2]]),
                "term.labels"))) {
                for (s in 1:ns) {
                  if (length(t.nms) > length(t.nms.sess)) {
                    t.beta[s, ] <- pv[pn %in% c(names.beta.trap[grep(fixed=TRUE,paste("session",
                      s, sep = ""), names.beta.trap)], names.beta.trap[as.vector(unlist(sapply(t.nms.nosess,
                      function(x) {
                        grep(fixed=TRUE,x, names.beta.trap)
                      })))])]
                  }
                  else {
                    t.beta[s, ] <- pv[pn %in% c(names.beta.trap[grep(fixed=TRUE,paste("session",
                      s, sep = ""), names.beta.trap)])]
                  }
                }
            }
            else {
                for (s in 1:ns) {
                  t.beta[s, ] <- pv[pn %in% names.beta.trap]
                }
            }
        }
        if (DorN == "N") {
            if (dIPP) {
                d.beta <- matrix(NA, ns, length(d.nms))
                if ("session" %in% all.vars(model[[1]])) {
                  for (s in 1:ns) {
                    d.beta[s, ] <- pv[pn %in% names.beta.den[grep(fixed=TRUE,paste("session",
                      s, sep = ""), names.beta.den)]]
                  }
                }
                else {
                  for (s in 1:ns) {
                    d.beta[s, ] <- pv[pn %in% names.beta.den]
                  }
                }
            }
        }
        else {
            d.beta <- pv[pn %in% names.beta.den]
        }
        if (distmet == "ecol") {
            dist.beta <- pv[pn %in% names.dist]
        }
        if (n0Session)
            n0 <- exp(pv[pn %in% names.n0])
        if (!n0Session)
            n0 <- rep(exp(pv[pn %in% names.n0]), ns)
        if (sexmod == "constant")
            psi.sex <- rep(plogis(pv[grep(fixed=TRUE,"psi", pn)]), ns)
        if (sexmod == "session")
            psi.sex <- plogis(pv[grep(fixed=TRUE,"psi", pn)])
        outLik <- 0
        if (predict) {
            preds <- list()
            lik.bits <- list()
            ss.bits <- list()
        }
        for (s in 1:length(scrFrame$caphist)) {
            Ys <- scrFrame$caphist[[s]]
            if (telem){ #check if telemetry exists
              Ytels <- YYtel[[s]]
              sxtel <- scrFrame$telemetry$indCovs[[s]]$sex + 1
              cap.tel <- scrFrame$telemetry$cap.tel[[s]]  #index of captured ind w/ collars
              lik.marg.tel <- rep(NA, nrow(Ytels))
            } 
            if (predict)
                preds[[s]] <- matrix(NA, nrow = nrow(Ys) + 1,
                  ncol = nrow(ssDF[[s]]))

            zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
            Ys <- abind(Ys, zeros, along = 1)

            sx <- c(scrFrame$indCovs[[s]]$sex + 1, NA)
            
            if (distmet == "ecol") {
                cost <- exp(dm.cost[[s]] %*% dist.beta)
                costR <- rasterFromXYZ(cbind(costDF[[s]][, c(1,
                  2)], cost))
                if (is.null(PROJ)) {
                  projection(costR) <- "+proj=utm +zone=12 +datum=WGS84"
                }
                else {
                  projection(costR) <- PROJ
                }
                tr <- transition(costR, transitionFunction = function(x) (1/(mean(x))),
                  direction = directions)
                trLayer <- geoCorrection(tr, scl = F)
                D[[s]] <- costDistance(trLayer, as.matrix(scrFrame$traps[[s]][,
                  c("X", "Y")]), as.matrix(ssDF[[s]][, c("X",
                  "Y")]))
            }
            if (smallslow) {
                if (distmet == "euc") {
                  D <- list()
                  D[[s]] <- e2dist(scrFrame$traps[[s]][, c("X",
                    "Y")], ssDF[[s]][, c("X", "Y")])
                }
            }
            if (telem){
              # only Euclidean distance for telemetry fixes
              Drsf[[s]] <- e2dist(rsfDF[[s]][, c("X", "Y")], rsfDF[[s]][, c("X", "Y")])
            }
            lik.marg <- lik.marg1 <- lik.marg2 <- rep(NA, nrow(Ys))
            if (!is.null(trimS)) {
                pixels.prior <- rep(T, nG[s])
                pixels.post <- apply(D[[s]], 2, min) <= trimS
                pixels <- (pixels.prior & pixels.post)
                pixels <- ifelse(pixels, 1, 0)
            }
            else {
                pixels <- rep(1, nG[s])
            }
            if (DorN == "N") {
                if (dIPP) {
                  d.s <- exp(dm.den[[s]] %*% d.beta[s, ])
                  pi.s <- (d.s * pixels)/sum(d.s * pixels)
                }
                if (dHPP) {
                  pi.s <- pixels/(sum(pixels))
                }
            }
            else {
                d.s <- exp(dm.den[[s]] %*% d.beta)
                pi.s <- (d.s * pixels)/sum(d.s * pixels)
            }
            
            # some collared ind captured, so keep lik.cond for combining later
            if (telemetry == "dep"){
                lik.cond.tel <- matrix(0,nrow=length(cap.tel),ncol=nG[s])
            }
            
            for (i in 1:nrow(Ys)) {
                if (plotit) {
                  pp <- sum(trimR[[s]][[i]][[k]])
                  plot(ssDF[[s]][, c("X", "Y")], pch = 16, col = "grey", cex = 0.5, asp = 1,
                       main = paste("Session:", s, " Individual: ", i, " traps: ", sum(pp), sep = " "))
                  points(ssDF[[s]][trimC[[s]][[i]], c("X", "Y")], pch = 16, col = 2, cex = mycex)
                  points(scrFrame$traps[[s]][,c("X", "Y")], pch = 15, col = "grey", cex = 1)
                  points(scrFrame$traps[[s]][trimR[[s]][[i]][[k]], c("X")],
                         scrFrame$traps[[s]][trimR[[s]][[i]][[k]],c("Y")],
                         pch = 15, col = 1, cex = 1)
                  pp <- apply(Ys[i,,],1,sum)>0
                  points(scrFrame$traps[[s]][pp,c("X", "Y")], pch = 15, col = 3, cex = mycex)
                }
                if (!is.na(sx[i])) {


#                     if (multicatch)
#                  Pm <- Pm1 <- Pm2 <- matrix(0, sum(trimR[[s]][[i]][[k]]) +                     1, sum(trimC[[s]][[i]]))
#                if (!multicatch)
#                  Pm <- Pm1 <- Pm2 <- tmpPm <- matrix(0, sum(trimR[[s]][[i]][[k]]), sum(trimC[[s]][[i]]))

                  lik.cond <- numeric(nG[s])


                    for (k in 1:nK[s]) {
                      if(!("removed" %in% names(scrFrame$indCovs[[s]]))){
                        dead <- 1
                      }else{
                        if(i < nrow(Ys)){
                          dead <- ifelse(k > scrFrame$indCovs[[s]]$removed[i],0,1)
                        }else{
                          dead <- 1
                        }
                      }
                      
                      




                      if (pBehave) {
                      a0 <- alpha0[s, k, sx[i], 1] * (1 - c(prevcap[[s]][i,
                        , k])) + alpha0[s, k, sx[i], 2] * c(prevcap[[s]][i,
                        , k])
                    }
                    else {
                      a0 <- rep(alpha0[s, k, sx[i], 1], nrow(D[[s]]))
                    }
                    if (trap.covs) {
                      a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,
                        ]))
                    }
                    if (encmod == "B")
                      probcap <- c(dead*plogis(a0[trimR[[s]][[i]][[k]]])) *
                        exp(-alphsig[s, sx[i]] * D[[s]][trimR[[s]][[i]][[k]],
                          trimC[[s]][[i]]]^theta)
                    if (encmod %in% c("P","CLOG"))
                      probcap <- c(dead*exp(a0[trimR[[s]][[i]][[k]]])) *
                        exp(-alphsig[s, sx[i]] * D[[s]][trimR[[s]][[i]][[k]],
                          trimC[[s]][[i]]]^theta)
                    if (!multicatch) {
                      if (encmod == "B") {
                        probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          1, probcap[1:length(probcap)], log = TRUE))
                      }
                      if (encmod == "P") {
                        probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          probcap[1:length(probcap)], log = TRUE))
                      }
                      if (encmod == "CLOG") {
                        probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
                      }
                    }
                    else {
                      probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
                      probcap <- t(t(probcap)/colSums(probcap))
                      vvv <- rep(c(Ys[i, trimR[[s]][[i]][[k]], k],
                        1 - any(Ys[i, trimR[[s]][[i]][[k]], k] > 0)),
                        sum(trimC[[s]][[i]]))
                      vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv ==
                        1])
                      probcap[1:length(probcap)] <- vvv
                    }
#                    if (!is.null(scrFrame$trapOperation)) {
#                      probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]],k]
#                    }

          if(is.matrix(probcap))
               lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + colSums(probcap)
          if(!is.matrix(probcap))
               lik.cond[trimC[[s]][[i]]]<- lik.cond[trimC[[s]][[i]]] + probcap

               #####Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]

           }


           #                  lik.cond <- numeric(nG[s])
                  if (telemetry == "dep"){
                    if (i %in% cap.tel){
                      lik.cond.tel[match(i,cap.tel),] <- lik.cond
                    }
                  }
                  
                  lik.cond[trimC[[s]][[i]]] <- exp(lik.cond[trimC[[s]][[i]]])  ####colSums(Pm, na.rm = T))
######            lik.cond[trimC[[s]][[i]]] <- exp(colSums(Pm,na.rm = T))
                  tmpPsi <- (sx[i] == 1) * (1 - psi.sex[s]) +
                    (sx[i] == 2) * psi.sex[s]
                  lik.marg[i] <- sum(lik.cond * pi.s) * tmpPsi
                  if (predict) {
                    preds[[s]][i, ] <- lik.cond * pi.s * tmpPsi/lik.marg[i]
                  }
                }
                else {


#               if (multicatch)
#                  Pm <- Pm1 <- Pm2 <- matrix(0, sum(trimR[[s]][[i]][[k]]) +                     1, sum(trimC[[s]][[i]]))
#                if (!multicatch)
#                  Pm <- Pm1 <- Pm2 <- tmpPm <- matrix(0, sum(trimR[[s]][[i]][[k]]),                     sum(trimC[[s]][[i]]))

                  lik.cond1<- lik.cond2 <- numeric(nG[s])


                  for (k in 1:nK[s]) {




                      if (pBehave) {
                      a0.1 <- alpha0[s, k, 1, 1] * (1 - c(prevcap[[s]][i,
                        , k])) + alpha0[s, k, 1, 2] * c(prevcap[[s]][i,
                        , k])
                      a0.2 <- alpha0[s, k, 2, 1] * (1 - c(prevcap[[s]][i,
                        , k])) + alpha0[s, k, 2, 2] * c(prevcap[[s]][i,
                        , k])
                    }
                    else {
                      a0.1 <- rep(alpha0[s, k, 1, 1], nrow(D[[s]]))
                      a0.2 <- rep(alpha0[s, k, 2, 1], nrow(D[[s]]))
                    }
                    if (trap.covs) {
                      a0.1 <- a0.1 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
                      a0.2 <- a0.2 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
                    }
                    if (encmod == "B")
                      probcap <- c(plogis(a0.1[trimR[[s]][[i]][[k]]])) *
                        exp(-alphsig[s, 1] * D[[s]][trimR[[s]][[i]][[k]],
                          trimC[[s]][[i]]]^theta)
                    if (encmod %in% c("P","CLOG"))
                      probcap <- c(exp(a0.1[trimR[[s]][[i]][[k]]])) *
                        exp(-alphsig[s, 1] * D[[s]][trimR[[s]][[i]][[k]],
                          trimC[[s]][[i]]]^theta)
                    if (!multicatch) {
                      if (encmod == "B") {
                        probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          1, probcap[1:length(probcap)], log = TRUE))
                      }
                      if (encmod == "P") {
                        probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          probcap[1:length(probcap)], log = TRUE))
                      }
                      if (encmod == "CLOG") {
                        probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
                      }
                    }
                    else {
                      probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
                      probcap <- t(t(probcap)/colSums(probcap))
                      vvv <- rep(c(Ys[i, trimR[[s]][[i]][[k]], k],
                        1 - any(Ys[i, trimR[[s]][[i]][[k]], k] > 0)),
                        sum(trimC[[s]][[i]]))
                      vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv ==
                        1])
                      probcap[1:length(probcap)] <- vvv
                    }
#                    if (!is.null(scrFrame$trapOperation)) {
#                      probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]][[k]],  k]
#                    }


               #Pm1[1:length(Pm1)] <- Pm1[1:length(Pm1)] +  probcap[1:length(probcap)]
         if(is.matrix(probcap))
               lik.cond1[trimC[[s]][[i]]]<- lik.cond1[trimC[[s]][[i]]] + colSums(probcap)
          if(!is.matrix(probcap))
               lik.cond1[trimC[[s]][[i]]]<- lik.cond1[trimC[[s]][[i]]] + probcap

#               lik.cond1[trimC[[s]][[i]]]<- lik.cond1[trimC[[s]][[i]]] + colSums(probcap)


                    if (encmod == "B")
                      probcap <- c(plogis(a0.2[trimR[[s]][[i]][[k]]])) *
                        exp(-alphsig[s, 2] * D[[s]][trimR[[s]][[i]][[k]],
                          trimC[[s]][[i]]]^theta)
                    if (encmod %in% c("P","CLOG"))
                      probcap <- c(exp(a0.2[trimR[[s]][[i]][[k]]])) *
                        exp(-alphsig[s, 2] * D[[s]][trimR[[s]][[i]][[k]],
                          trimC[[s]][[i]]]^theta)
                    if (!multicatch) {
                      if (encmod == "B") {
                        probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          1, probcap[1:length(probcap)], log = TRUE))
                      }
                      if (encmod == "P") {
                        probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          probcap[1:length(probcap)], log = TRUE))
                      }
                      if (encmod == "CLOG") {
                        probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,
                          trimR[[s]][[i]][[k]], k], sum(trimC[[s]][[i]])),
                          1, 1-exp(-probcap[1:length(probcap)]), log = TRUE))
                      }
                    }
                    else {
                      probcap <- rbind(probcap, rep(1, sum(trimC[[s]][[i]])))
                      probcap <- t(t(probcap)/colSums(probcap))
                      vvv <- rep(c(Ys[i, trimR[[s]][[i]][[k]], k],
                        1 - any(Ys[i, trimR[[s]][[i]][[k]], k] > 0)),
                        sum(trimC[[s]][[i]]))
                      vvv[vvv == 1] <- log(probcap[1:length(probcap)][vvv ==
                        1])
                      probcap[1:length(probcap)] <- vvv
                    }
#                    if (!is.null(scrFrame$trapOperation)) {
#                      probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]],
#                        k]
#                    }

         if(is.matrix(probcap))
               lik.cond2[trimC[[s]][[i]]]<- lik.cond2[trimC[[s]][[i]]] + colSums(probcap)
          if(!is.matrix(probcap))
               lik.cond2[trimC[[s]][[i]]]<- lik.cond2[trimC[[s]][[i]]] + probcap

#               lik.cond2[trimC[[s]][[i]]]<- lik.cond2[trimC[[s]][[i]]] + colSums(probcap)

               ###Pm2[1:length(Pm2)] <- Pm2[1:length(Pm2)] +  probcap[1:length(probcap)]
                  }

                  ###lik.cond1 <- lik.cond2 <- numeric(nG[s])
             #     lik.cond1[trimC[[s]][[i]]] <- exp(colSums(Pm1, na.rm = T))
             #     lik.cond2[trimC[[s]][[i]]] <- exp(colSums(Pm2, na.rm = T))

                  lik.cond1[trimC[[s]][[i]]] <- exp(lik.cond1[trimC[[s]][[i]]])  ####colSums(Pm, na.rm = T))
                  lik.cond2[trimC[[s]][[i]]] <- exp(lik.cond2[trimC[[s]][[i]]])  ####colSums(Pm, na.rm = T))


                  lik.marg1[i] <- sum(lik.cond1 * pi.s)
                  lik.marg2[i] <- sum(lik.cond2 * pi.s)
                  lik.marg[i] <- lik.marg1[i] * (1 - psi.sex[s]) +
                    lik.marg2[i] * psi.sex[s]
                  if (predict) {
                    lik.cond <- (lik.cond1 * (1 - psi.sex[s]) +
                      lik.cond2 * psi.sex[s])
                    preds[[s]][i, ] <- (lik.cond * pi.s)/lik.marg[i]
                  }
                }
            }
            
            if(telem){

              if(RSF){
                rsf.lam0 <- dm.rsf[[s]] %*% c(t.beta[s,])
                rsf.lam0 <- array(rsf.lam0,dim=c(nrow(rsfDF[[s]]),nrow(rsfDF[[s]])))
              } else {
                rsf.lam0 <- 0
              }
              for (i in 1:nrow(Ytels)){
                
                probs <- t(exp(rsf.lam0 - alphsig[s, sxtel[i]] * Drsf[[s]]^theta))
                denom <- rowSums(probs)
                probs <- t(probs/denom)
                
                lik.marg.tel[i] <- sum( exp(Ytels[i,,drop=F] %*% log(probs)) * as.vector(pi.s) )
                #browser()
                if (telemetry == "dep"){
                  if (i <= length(cap.tel)){
                    # combine conditional likelihoods if some collared ind were captured
                    lik.cond.tot <- (Ytels[i,,drop=F] %*% log(probs)) + lik.cond.tel[i,]
                    #lik.cond.tot[trimC[[s]][[cap.tel[i]]]] <- exp(lik.cond.tot[trimC[[s]][[cap.tel[i]]]])
                    lik.cond.tot[is.na(lik.cond.tot)] <- 0
                    lik.cond.tot[lik.cond.tot != 0] <- exp(lik.cond.tot[lik.cond.tot != 0])
                  
                    tmpPsi <- (sx[cap.tel[i]] == 1) * (1 - psi.sex[s]) + (sx[cap.tel[i]] == 2) * psi.sex[s]
                    # fix marginal likelihoods
                    lik.marg[cap.tel[i]] <- sum(lik.cond.tot * as.vector(pi.s)) * tmpPsi
                    lik.marg.tel[i] <- 1
                  
                    if (predict){
                      preds[[s]][cap.tel[i], ] <- lik.cond.tot * as.vector(pi.s) * tmpPsi/lik.marg[cap.tel[i]]
                    }
                  
                  }
                }
                              
              }
            }

            if (!predict) {
                if (DorN == "N") {
                  nv <- c(rep(1, length(lik.marg) - 1), n0[s])
                  part1 <- lgamma((nrow(Ys) - 1) + n0[s] + 1) -
                    lgamma(n0[s] + 1)
                  part2 <- sum(nv * log(lik.marg))
                }
                if (DorN == "D") {
                  nv <- c(rep(1, length(lik.marg) - 1), 1)
                  atheta <- 1 - lik.marg[nrow(Ys)]
                  nind <- nrow(Ys) - 1
                  part1 <- nind * log(sum(d.s * pixels)) - sum(d.s *
                    pixels) * atheta
                  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
              }
              if(telem){
                part3 <- sum(log(lik.marg.tel))
              } else {
                part3 <- 0
              }
              
              ll <- -1 * (part1 + part2 + part3)
              outLik <- outLik + ll
            }
            if (predict) {
                lik.bits[[s]] <- cbind(lik.mar = lik.marg)
                ss.bits[[s]] <- cbind(pi.s, d.s, lik.cond)
                colnames(ss.bits[[s]]) <- c("pi.s", "d.s", "lik.cond")
            }
        }
        if (!predict) {
            out <- outLik
            return(out)
        }
        if (predict) {
            return(list(preds = preds, ss.bits = ss.bits, lik.bits = lik.bits,
                        ssDF = ssDF, data = scrFrame$caphist, traps = scrFrame$traps))
        }

    }   # end likelihood

    if(getStarts == FALSE) {
      if (!predict) {
        message("Fitting model: D", paste(model)[1], 
                            ", p0", paste(model)[2], 
                         ", sigma", paste(model)[3],
                           ", asu", paste(model)[4], sep = " ")
        if (!anySex) {
          if (telem){
            message("Telemetry: ",telemetry)
          }
          message("Using ll function 'msLL.nosex' \nHold on tight!")
          message(Sys.time())
          message(paste(pn, " ", sep = " | "))
          message(" ")
          myfit <- suppressWarnings(
                     nlm(msLL.nosex, p = pv, pn = pn, D = D, nG = nG, nK = nK,
                         hiK = hiK, dm.den = dm.den, dm.trap = dm.trap,
                         hessian = hessian, print.level = print.level, 
                         iterlim = 200))
        }
        else {
          if (telem){
            message("Telemetry: ",telemetry)
          }
          message("Using ll function 'msLL.sex' \nHold on tight!")
          message(Sys.time())
          message(paste(pn, " ", sep = " | "))
          message(" ")
          myfit <- suppressWarnings(
                     nlm(msLL.sex, scrFrame=scrFrame, p = pv, pn = pn, D = D, nG = nG, nK = nK,
                         hiK = hiK, dm.den = dm.den, dm.trap = dm.trap,
                         hessian = hessian, print.level = print.level, 
                         iterlim = 200))
        }
        links <- rep(NA, length(pn))
        pars <- myfit$estimate
        if (encmod == "B") {
          links[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- "(logit)"
        }
        else {
          links[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- "(log)"
        }
        links[grep(fixed=TRUE,"sig.(Intercept)", pn)] <- "(log)"
        links[grep(fixed=TRUE,"n0.", pn)] <- "(log)"
        links[grep(fixed=TRUE,"d0.(Intercept)", pn)] <- "(log)"
        links[grep(fixed=TRUE,"c0.(Intercept)", pn)] <- "(log)"
        links[grep(fixed=TRUE,"psi", pn)] <- "(logit)"
        links[grep(fixed=TRUE,"beta", pn)] <- "(Identity)"
        trans.mle <- rep(0, length(pv))
        if (encmod == "B") {
          trans.mle[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- plogis(pars[grep(fixed=TRUE,"p0.(Intercept)", pn)])
        }
        else {
          trans.mle[grep(fixed=TRUE,"p0.(Intercept)", pn)] <- exp(pars[grep(fixed=TRUE,"p0.(Intercept)", pn)])
        }
        trans.mle[grep(fixed=TRUE,"sig.(Intercept)", pn)] <- exp(pars[grep(fixed=TRUE,"sig.(Intercept)", pn)])
        trans.mle[grep(fixed=TRUE,"n0.", pn)] <- exp(pars[grep(fixed=TRUE,"n0.", pn)])
        trans.mle[grep(fixed=TRUE,"d0.(Intercept)", pn)] <- exp(pars[grep(fixed=TRUE,"d0.(Intercept)", pn)])
        trans.mle[grep(fixed=TRUE,"c0.(Intercept)", pn)] <- exp(pars[grep(fixed=TRUE,"c0.(Intercept)", pn)])
        trans.mle[grep(fixed=TRUE,"psi", pn)] <- plogis(pars[grep(fixed=TRUE,"psi", pn)])
        trans.mle[grep(fixed=TRUE,"beta", pn)] <- pars[grep(fixed=TRUE,"beta", pn)]
        if (pBehave) {
          links[grep(fixed=TRUE,"pBehav", pn)] <- "(Identity)"
          trans.mle[grep(fixed=TRUE,"pBehav", pn)] <- pars[grep(fixed=TRUE,"pBehav", pn)]
        }
        std.err <- rep(rep(NA, length(pv)))
        trans.se <- rep(NA, length(pv))
        if("hessian" %in% names(myfit)) {
          if(sum(myfit$hessian) != 0){
          #Need a check for this error and return mles and a warning
          #Error in solve.default(myfit$hessian) :
          #Lapack routine dgesv: system is exactly singular: U[1,1] = 0
          std.err <- sqrt(diag(solve(myfit$hessian)))
        }
        else {
          warning("Something went wrong! Try better starting values.")
        }
        }
        outStats <- data.frame(parameters = pn, link = links, mle = myfit$estimate, 
                               std.er = std.err, mle.tr = trans.mle, se.tr = trans.se)
        VcV <- NULL
        if(DorN == "N") {
        ## write some code to return per session density 
        }
        else{
        ## write some code to return per session abundance 
            
        }
        endtime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
        output <- list(call = cl, model=model.call,rawOutput = myfit, 
                       outStats = outStats, 
                       coef.mle = data.frame(param = pn, mle = myfit$estimate),
                       Area = areaS, nObs = unlist(lapply(scrFrame$caphist,nrow)),
                       mmdm = mmdm, nll = myfit$minimum, 
                       AIC = 2 * myfit$minimum + 2 * length(myfit$estimate),
                       started = starttime, ended = endtime,
                       proctime = (proc.time() - ptm)[3]/60, scrFrame = scrFrame,
                       ssDF = ssDF, costDF = costDF, rsfDF = rsfDF)
        class(output) <- "oSCR.fit"
        return(output)
      }
      if (predict) {
          message("Predicting model: D", paste(model)[1], 
                                 ", p0", paste(model)[2], 
                              ", sigma", paste(model)[3],
                                ", asu", paste(model)[4], sep = " ")
          if (!anySex) {
            if (telem){
              message("Telemetry: ",telemetry)
            }
            message("Using ll function 'msLL.nosex' \nHold on tight!")
            message(Sys.time())
            message(paste(pn, " ", sep = " | "))
            message(" ")
            myfit <- msLL.nosex(p = start.vals, pn = pn, D = D, hiK = hiK, 
                                nG = nG, nK = nK, dm.den = dm.den, 
                                dm.trap = dm.trap)
          }
          else {
            if (telem){
              message("Telemetry: ",telemetry)
            }
            message("Using ll function 'msLL.sex' \nHold on tight!")
            message(Sys.time())
            message(paste(pn, " ", sep = " | "))
            message(" ")
            myfit <- msLL.sex(scrFrame=scrFrame, p = start.vals, pn = pn, D = D, nK = nK, nG = nG, 
                              hiK = hiK, dm.den = dm.den, dm.trap = dm.trap)
          }
          return(myfit)
        }
    }
}
