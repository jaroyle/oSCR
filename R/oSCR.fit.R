oSCR.fit <-
function(scrFrame, model = list(D~1, p0~1, sig~1, asu~1), ssDF = NULL, costDF = NULL,
         distmet=c("euc","user","ecol")[1], sexmod = c('constant','session')[1],
         encmod = c("B","P")[1], DorN = c('D','N')[1], directions = 8, Dmat = NULL,
         trimS = NULL, start.vals = NULL, PROJ = NULL, pxArea = 1, plotit = F,
         mycex = 0.5, tester = F, pl = 0, nlmgradtol = 1e-6, nlmstepmax = 10,
         predict=FALSE, smallslow = FALSE, multicatch=FALSE,hessian=T, print.level = 0,
         getStarts = FALSE){
##NOTES: 'session' = n0 different
##       'Session' = betas differ too!
##       Reommend trimming the State space to trim value!


################################################################################
#                            Setting things up!                                 #
################################################################################

# this is a function that make design matrices that are fully 'contrast' parameterizations
my.model.matrix <- function(form,data){
  mdm <-suppressWarnings(model.matrix(form, data, contrasts.arg =
  lapply(data.frame(data[,sapply(data.frame(data), is.factor)]), contrasts, contrasts = FALSE)))
  return(mdm)
}
# here I want to find a max obserevd movement for trimming:
#  if(is.null(trimS)){
#    max.dist <- 0
    max.dist <- NULL
    for(i in 1:length(scrFrame$caphist)){
     for(j in 1:nrow(scrFrame$caphist[[i]])){
       where <- apply(scrFrame$caphist[[i]][j,,],1,sum)>0
       max.dist <- c(max.dist,max(0,dist(scrFrame$traps[[i]][where,c("X","Y")]),na.rm=T))
     }
    }
    mmdm <- mean(max.dist[max.dist>0],na.rm=T)
################################################################################
# Some setting and checks
#

  ptm <- proc.time() # Start
  starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y") # Date stamp for start
  cl <- match.call(expand.dots = TRUE) # The call given to the function

  if(!require(abind)) stop("need to install package 'abind'")
  if(!require(Formula)) stop("need to load package 'Formula'")
  if(distmet=="ecol"){
   if(!require(raster)) stop("need to install package 'raster'")
   if(!require(gdistance)) stop("need to install package 'gdistance'")
  }
  if(!inherits(scrFrame,"scrFrame")){
    stop("Data must be of class 'scrFrame'")
  }
  if(encmod=="B" & max(unlist(lapply(scrFrame$caphist,max)))>1){
    stop("Data in caphist must be Binary")
  }

  if(distmet=="ecol"){
   if(is.null(PROJ)){
     message("Projection not provided, using default: '+proj=utm +zone=12 +datum=WGS84'")
   }
  }

  if(!is.null(ssDF) & length(ssDF)!=length(scrFrame$caphist))
    stop("A 'state space' object must be provided for EACH session.")
#need to provide a warning NOT an error!
#  if(length(scrFrame$traps)!=length(scrFrame$caphist))
#    warning("A 'traps' object must be provided for EACH session.")

  if(multicatch){
   for(s in 1:length(scrFrame$caphist)){
     captures<- apply(scrFrame$caphist[[s]],c(1,3),sum)  # Number of captures per individual and occ
     if(any(captures>1))
      stop("error: multicatch system cannot have > 1 capture.")
   }
  }
## ADD A CHECK FOR POISSON vs. BINOMIAL DATA RELATIVE TO SELECTED ENCOUNTER MODEL!
  maxY <- unlist(lapply(scrFrame$caphist,max))
  if(any(maxY > 1) & encmod == "B")
      stop("caphist must be binary when using the Binomial encounter model")

  if(all(maxY == 1) & encmod == "P")
      stop("caphist looks binary but Poisson encounter model is selected")

  pars.p0 <- NULL ; names.p0 <- NULL
  pars.sig <- NULL ; names.sig <- NULL
  pars.beta.trap <- NULL ; names.beta.trap <- NULL
  pars.beta.den <- NULL ; names.beta.den <- NULL
  pars.n0 <- NULL ; names.n0 <- NULL
  pars.beta.den <- NULL ; names.beta.den <- NULL
  pars.dist <- NULL ; names.dist <- NULL
  pars.dist <- NULL ; names.dist <- NULL
  singleS <- NULL ; singleT <- NULL ; singleG <- NULL

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

  if(length(model)==3){model[[4]] <- formula(~1)}
  for(i in 1:4){
    model[[i]] <- update.formula(model[[i]],NULL~.)
  }


################################################################################
# Make a state-space if one doesnt exist (remove - building this is important!)

  # make a ssDF if one doesnt exist! USE 'make.basic.S' function!
  if(is.null(ssDF)){
    message("Generating a state space based on traps")
    dHPP <- TRUE
    ssDF <- make.ssDF(scrFrame, buffer,res)
  }

  ns <- length(scrFrame$caphist)
  nt <- length(scrFrame$traps)
  nK <- unlist(lapply(scrFrame$caphist,function(x) dim(x)[3]))
  hiK <- max(nK)
  nG <- unlist(lapply(ssDF,nrow))

  nnn <- all(unlist(lapply(ssDF,function(x){"session" %in% names(x)})))
  areaS <- NULL
  if("session" %in% all.vars(model[[1]]) & (!nnn)){
   for(s in 1:ns){
    ssDF[[s]]$session <- factor(rep(s,nrow(ssDF[[s]])),levels=1:ns)
   }
  }


################################################################################
# Some general settings
#

  allvars.D <- all.vars(model[[1]])
  dens.fx <- allvars.D[!allvars.D %in% c("D","session","Session")]

  allvars.T <- all.vars(model[[2]])
  trap.fx <- allvars.T[!allvars.T %in% c("p0","session","Session","sex","t","T","b")]

  allvars.p0a <- all.vars(model[[2]])
  allvars.p0 <- allvars.p0a[!allvars.p0a=="p0"]
  allvars.siga <- all.vars(model[[3]])
  allvars.sig <- allvars.siga[!allvars.siga=="sig"]
  var.p0.1 <- "sex" %in% allvars.p0
  var.p0.2 <- "session" %in% allvars.p0
  var.p0.3 <- "t" %in% allvars.p0
  var.p0.4 <- any(c("sex:session", "session:sex") %in% attr(terms(model[[2]]),"term.labels"))
  var.sig.1 <- "sex" %in% allvars.sig
  var.sig.2 <- "session" %in% allvars.sig
  var.sig.3 <- any(c("sex:session", "session:sex") %in% attr(terms(model[[3]]),"term.labels"))
  var.b.1 <- "b" %in% attr(terms(model[[2]]),"term.labels")
  var.b.2 <- any(c("b:sex", "sex:b") %in% attr(terms(model[[2]]),"term.labels"))
  var.b.3 <- any(c("b:session", "session:b") %in% attr(terms(model[[2]]),"term.labels"))
  var.b.4 <- any(c("b:session:sex", "b:sex:session", "sex:session:b", "sex:b:session",
                   "session:b:sex","session:sex:b") %in% attr(terms(model[[2]]),"term.labels"))
  pBehave <- any(c(var.b.1,var.b.2,var.b.3,var.b.4))

  allvars.dist <- all.vars(model[[4]])
  allvars.dist <- allvars.dist[!allvars.dist=="asu"]

  for(s in 1:ns){
   if(!is.null(trimS)){
     pixels.prior <- rep(T,nG[s])
     pixels.post <- apply(e2dist(scrFrame$traps[[s]][,c("X","Y")],
                                 ssDF[[s]][,c("X","Y")]),2,min)<=trimS
     pixels <- (pixels.prior & pixels.post)
     pixels <- ifelse(pixels,1,0)
   }else{
     pixels <- rep(1,nG[s])
   }
   areaS <- c(areaS,sum(pixels) * pxArea)
  }


################################################################################
# design matrix for the trap covariates

  if(length(trap.fx)>0){
  #check all covariates are present:
    trap.covs <- TRUE
    tcovnms <- colnames(scrFrame$trapCovs[[1]][[1]]) # covs MUST be the same!
    tCovMissing <- trap.fx[which(!trap.fx %in% tcovnms)]

   if(length(tCovMissing)>0){
     stop("I cant find theses covariates in 'scrFrame$trapCovs'",
          for(i in tCovMissing)print(i))
   }

  #make the design matrix
    mod2 <- update(model[[2]],~. - sex - session - Session - t - b - 1)
   if(any(c("session","Session") %in% allvars.T)) tSession <- TRUE
    for(s in 1:ns){
      tmp.dm <- list()
    for(k in 1:nK[s]){
      tmp.dm[[k]] <- my.model.matrix(mod2,scrFrame$trapCovs[[s]][[k]])
     if(s==1 && k==1) t.nms <- colnames(tmp.dm[[k]])
     if(nrow(tmp.dm[[k]])!=nrow(scrFrame$trapCovs[[s]][[k]])){#deal with NAs!
         mis <-setdiff(rownames(scrFrame$trapCovs[[s]][[k]]),
                       rownames(my.model.matrix(mod2,scrFrame$trapCovs[[s]][[k]])))
         tmp.insert <- matrix(NA,length(mis),ncol(tmp.dm[[k]]))
         row.names(tmp.insert) <- mis
         tmp.dm[[k]] <- rbind(tmp.dm[[k]],tmp.insert)
         tmp.dm[[k]] <- tmp.dm[[k]][order(as.numeric(row.names(tmp.dm[[k]]))),]
     }#end if
    }#end k
     dm.trap[[s]] <- tmp.dm
   }#end s

   #set up the parameters to be estimated
   if("Session" %in% all.vars(model[[2]])){
     tmpTsess <- rep(1:ns,each=length(t.nms))
     tmpTcovs <- rep(t.nms,ns)
     names.beta.trap <- paste("t.beta.",tmpTcovs,".sess",tmpTsess,sep="")
     pars.beta.trap <- rnorm(length(names.beta.trap))
   }else{
     names.beta.trap <- paste("t.beta.",t.nms,sep="")
     pars.beta.trap <- rnorm(length(names.beta.trap))
   }
  }

################################################################################
# design matrix for the density - can have session specififcty

  # check all covariates are present:
  if(length(dens.fx)>0){
    dcovnms <- colnames(ssDF[[1]]) # covs MUST be the same!
    dCovMissing <- dens.fx[which(!dens.fx %in% dcovnms)]
   if(length(dCovMissing)>0){
     stop("I cant find theses covariates in 'ssDF'",
          for(i in dCovMissing)print(i))
   }
  }

  #make the design matrix
  if(DorN=="N"){
   if(length(dens.fx)>0){
     dIPP <- TRUE
     mod1 <- update(model[[1]],~. - sex - Session - session - 1)
    for(s in 1:ns){
      dm.den[[s]] <- my.model.matrix(mod1,ssDF[[s]])
      if(s==1) d.nms <- colnames(dm.den[[s]])
    }
    if("Session" %in% all.vars(model[[1]])){
      tmpDsess <- rep(1:ns,each=length(d.nms))
      tmpDcovs <- rep(d.nms,ns)
      names.beta.den <- paste("d.beta.",tmpDcovs,".sess",tmpDsess,sep="")
      pars.beta.den <- rnorm(length(names.beta.den))
      names.n0 <- paste0("n0.sess",1:ns)
      pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
    }
    if("session" %in% all.vars(model[[1]])){
      names.beta.den <- paste0("d.beta.",d.nms,sep="")
      pars.beta.den <- rnorm(length(names.beta.den))
      names.n0 <- paste("n0.sess",1:ns)
      pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
    }
    if(!("session" %in% all.vars(model[[1]])) &
       !("Session" %in% all.vars(model[[1]]))){
      names.beta.den <- paste0("d.beta.",d.nms,sep="")
      pars.beta.den <- rnorm(length(names.beta.den))
      names.n0 <- paste0("n0.")
      pars.n0 <- log(mean(unlist(lapply(scrFrame$caphist,nrow))))
    }
   }else{
     dHPP <- TRUE
     names.beta.den <- NULL
     pars.beta.den <- NULL
    for(s in 1:ns){
      dm.den[[s]] <- my.model.matrix(~1,ssDF[[s]])
    }
    if("Session" %in% all.vars(model[[1]])){
      n0Session <- TRUE
      names.n0 <- paste0("n0.sess",1:ns)
      pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
    }
    if("session" %in% all.vars(model[[1]])){
      n0Session <- TRUE
      names.n0 <- paste0("n0.sess",1:ns)
      pars.n0 <- log(unlist(lapply(scrFrame$caphist,nrow)))
    }
    if(!("session" %in% all.vars(model[[1]])) &
       !("Session" %in% all.vars(model[[1]]))){
      names.n0 <- paste0("n0.")
      pars.n0 <- log(mean(unlist(lapply(scrFrame$caphist,nrow))))
    }
   }
  }

  if(DorN=="D"){
    mod1 <- update(model[[1]],~. - sex)
   for(s in 1:ns){
     dm.den[[s]] <- model.matrix(mod1,as.data.frame(ssDF[[s]]))
   }
    d.nms <- colnames(dm.den[[1]])
    names.beta.den <- paste("d.beta",d.nms,sep=".")
    chx <- grep("Intercept",names.beta.den)
    if(length(chx)>0)
      names.beta.den[chx] <- "d0."
    pars.d0 <- log(mean((unlist(lapply(scrFrame$caphist,nrow)))/unlist(lapply(ssDF,nrow))))
    pars.beta.den <- c(pars.d0,rep(0.1,length(names.beta.den)-1))
    pars.n0 <- NULL
    names.n0 <- NULL
  }


################################################################################
# make appropriate distance matrices need to extent to ecological distance
#

# right now the cost is the same across sessions!
  if(distmet=="ecol" && length(allvars.dist)==0){
    message("You specified 'ecological distance' (distmet='ecol') but provided no\ncost surface.
    Euclidean distance will be used.")
  }
  if(length(allvars.dist)>0){
    ccovnms <- colnames(costDF[[1]]) # covs MUST be the same!
    cCovMissing <- allvars.dist[which(!allvars.dist %in% ccovnms)]
   if(length(cCovMissing)>0){
     stop("I cant find theses covariates in 'costDF'",
          for(i in cCovMissing) print(i))
   }
  }

  mod4 <- update(model[[4]],~. - sex - session - 1)

  if(distmet == "ecol"){
  # need an error check here to see if cost coords match S
   for(s in 1:ns){
     dm.cost[[s]] <- my.model.matrix(mod4,costDF[[s]])
     names.dist <- paste("c.beta.",allvars.dist,sep="")
     pars.dist <- rep(0,length(names.dist))
   }
  }
  if(!smallslow){
   if(distmet == "euc"){# add a 'else' for providing your own!
    for(s in 1:ns){
      D[[s]] <- e2dist(scrFrame$traps[[s]][,c("X","Y")], ssDF[[s]][,c("X","Y")])
    }
   }
  }


################################################################################
# p0: define sex and/or session specific parameters
#  var.p0.1 <- "sex" only
#  var.p0.2 <- "session" only
#  var.p0.3 <- "t" only
#  var.p0.4 <- "sex:session" / "session:sex"

  if("indCovs" %in% names(scrFrame)){
   if("sex" %in% names(scrFrame$indCovs[[1]])){
     anySex <- TRUE
   }
  }

  tmp.p0.names <- "p0.int"

  if(sum(var.p0.1,var.p0.2,var.p0.3,var.p0.4)==0){
    tmp.p0.names <- "p0.int"
    pDot <- TRUE
  }

  if(var.p0.1 & !var.p0.2){
    tmp.p0.names <- c("p0.int","p0.male")
    pJustsex <- TRUE
  }

  if(!var.p0.1 & var.p0.2){
    if(ns>1){
      tmp.p0.names <- c("p0.int",paste0("p0.sess",2:ns))
      pJustsesh <- TRUE
    }else{
      tmp.p0.names <- "p0.int"
      pDot <- TRUE
    }
  }

  if(var.p0.1 & var.p0.2){
    tmp.p0.names <- c("p0.int","p0.male")
    if(ns>1){
      tmp.p0.names <- c(tmp.p0.names,paste0("p0.sess",2:ns))
      pJustsesh <- TRUE
    }else{
      tmp.p0.names <- tmp.p0.names
      pJustsex <- TRUE
    }
    pJustsex <- TRUE
  }

  if(var.p0.3){
    tmp.p0.names <- c(tmp.p0.names,paste0("p0.t",2:hiK))
    pTime <- TRUE
  }
  if(var.p0.4){
    if(ns>1){
      tmp.p0.names <- c("p0.int",paste0("p0.f.sess",2:ns),paste0("p0.m.sess",1:ns))
      pBothsexnsesh <- TRUE
    }else{
      tmp.p0.names <- c("p0.int","p0.male")
      pJustsex <- TRUE
    }
  }

  names.p0 <- tmp.p0.names
  pars.p0 <- rep(0,length(names.p0))
  pars.p0[1] <- -1.5 #strting value for p

  if(any(var.p0.1, var.p0.1, var.sig.1) && !anySex)
   stop("Sex defined in a model but no sex data provided.")

  #experimental!
  # var.b.1 <- "b" %in% attr(terms(model[[2]]),"term.labels")
  # var.b.2 <- c("b:sex", "sex:b") %in% attr(terms(model[[2]]),"term.labels"))
  # var.b.3 <- c("b:session", "session:b") %in% attr(terms(model[[2]]),"term.labels"))
  # var.b.4 <- c("b:session:sex", "b:sex:session", "sex:session:b", "sex:b:session",
  #               "session:b:sex","session:sex:b") %in% attr(terms(model[[2]]),"term.labels"))

  if(var.b.1 & !var.b.2 & !var.b.3 & !var.b.4){
    pars.p0 <- c(pars.p0,0)
    names.p0 <- c(names.p0,"p.behav")
    bDot <- TRUE
  }

  if(var.b.2 & !var.b.4){
    pars.p0 <- c(pars.p0,0,0)
    names.p0 <- c(names.p0,"p.behav.f","p.behav.m")
    bJustsex <- TRUE
  }

  if(var.b.3 & !var.b.4){
    pars.p0 <- c(pars.p0,rep(0,ns))
    names.p0 <- c(names.p0,paste0("p.behav.sess",1:ns))
    bJustsesh <- TRUE
  }

  if(var.b.4){
    pars.p0 <- c(pars.p0,rep(0,2*ns))
    names.p0 <- c(names.p0, paste0("p.behav.f.sess",1:ns),paste0("p.behav.m.sess",1:ns))
    bBothsexnsesh <- TRUE
  }


################################################################################
# define sex and/or session specific parameters
# NB: sig = 1/(2 * sigma^2) | sigma = sqrt(1/2*sig)
#
# var.sig.1 <- "sex" %in% allvars.sig
# var.sig.2 <- "session" %in% allvars.sig
# var.sig.3 <- any(c("sex:session", "session:sex") %in% attr(terms(model[[3]]),"term.labels"))

  tmp.sig.names <- "sig.int"

  if(sum(var.sig.1,var.sig.2,var.sig.3)==0){
    aDot <- TRUE
  }

  if(var.sig.2 & !var.sig.3){
    if(ns>1){
      tmp.sig.names <- c(tmp.sig.names,paste0("sig.sess.",2:ns))
      aJustsesh <- TRUE
    }else{
      aDot <- TRUE
    }
  }

  if(var.sig.1 & !var.sig.3){
    tmp.sig.names <- c(tmp.sig.names,"sig.male")
    aJustsex <- TRUE
  }

  if(var.sig.3){
    if(ns>1){
      tmp.sig.names <- c(tmp.sig.names,paste0("sig.f.sess",2:ns),paste0("sig.m.sess",1:ns))
      aBothsexnsesh <- TRUE
    }else{
      aJustsex <- TRUE
    }
  }

  names.sig <- tmp.sig.names
  pars.sig <- rep(0.1,length(names.sig))
  pars.sig[1] <- log(0.5 * mmdm)

################################################################################
# some objects to fit models
#

  if(scrFrame$type=='scr'){
    YY <- scrFrame$caphist
  }
  #need to add conversion function here from scrbook
  if(scrFrame$type=='secr'){
    YY <- 'secr2scr'
  }

  if(anySex){
   if(sexmod=='constant'){
     pars.sex <- 0
     names.sex <- "psi.constant"
   }
   if(sexmod=='session'){
     pars.sex <- rep(0,ns)
     names.sex <- paste("psi",1:ns,sep="")
   }
  }else{
    pars.sex <- NULL
    names.sex <- NULL
  }
  pv <- round(c(pars.p0,pars.sig,pars.beta.trap,pars.beta.den,pars.dist,pars.n0,pars.sex),2)
  pn <- c(names.p0,names.sig,names.beta.trap,names.beta.den,names.dist,names.n0,names.sex)
  if(!is.null(start.vals)){
   if(length(pv)==length(start.vals)){
     pv <- start.vals
   }else{
   message("The number of starting values provided doesnt match the \n
           number of parameters in the model. Randomly generated values \n
           are being used. Use getStarts = T to get correct length.")
   }
  }

################################################################################
#     Likelihood functions - these can probably be dumped somewhere else?      #
################################################################################


if(pBehave){
    prevcap<- list()

    for(s in 1:length(YY)){
      Ys <- YY[[s]]

     prevcap[[s]] <- array(0, dim=c(dim(Ys)[1],dim(Ys)[2],dim(Ys)[3]))
       first<- matrix(0,dim(Ys)[1],dim(Ys)[2])
      for(i in 1:dim(Ys)[1]){
       for(j in 1:dim(Ys)[2]){
        if(sum(Ys[i,j,])>0){
          first[i,j]<- min((1:(dim(Ys)[3]))[Ys[i,j,]>0] )
          prevcap[[s]][i,j,1:first[i,j]]<- 0
        if(first[i,j] < dim(Ys)[3])
          prevcap[[s]][i,j,(first[i,j]+1):(dim(Ys)[3])] <- 1
        }
       }
      }
       zeros <- array(0,c(1,dim(prevcap[[s]])[2],dim(prevcap[[s]])[3]))
       prevcap[[s]] <- abind(prevcap[[s]],zeros,along=1)
   }

   }


   ## All the trim stuff is done here
   ## Note: the "trap operation" should be done here too, REMOVE TRAPS THAT ARE NOT OPERATIONAL _OR_ far away

trimR<- trimC<- list()

for(s in 1:length(YY)){
    Ys<- YY[[s]]
# multicatch block 1 [looks like not needed]
     if(!multicatch){
       zeros <- array(0,c(1,dim(Ys)[2],dim(Ys)[3]))
       Ys <- abind(Ys,zeros,along=1)
     }
     if(multicatch){
       zeros <- array(0,c(1, dim(Ys)[2], dim(Ys)[3]))
       Ys <- abind(Ys,zeros,along=1)
     }


    trimR[[s]]<- list()
    trimC[[s]]<- list()
    for(i in 1:nrow(Ys)){

     if(is.null(trimS)){
        pp <- rep(T,ncol(Ys))
        trimC[[s]][[i]] <- rep(T,nG[s])
        trimR[[s]][[i]] <- pp
      }else{
       if(i<nrow(Ys)){
         pp <- apply(Ys[i,,],1,sum)>0
         trimR[[s]][[i]] <- apply(
                   rbind(
                    rep(trimS+2,nrow(scrFrame$traps[[s]])),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          scrFrame$traps[[s]][,c("X","Y")])),2,min)<=(2*trimS)
         trimC[[s]][[i]] <- apply(
                   rbind(
                    rep(trimS+2,nG[s]),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          ssDF[[s]][,c("X","Y")])),2,min,na.rm=T)<=trimS
       }else{
         pp <- rep(T,ncol(Ys))
         trimC[[s]][[i]] <- apply(
                   rbind(
                    rep(trimS+2,nG[s]),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          ssDF[[s]][,c("X","Y")])),2,min,na.rm=T)<=trimS
         #trimC <- rep(T,nG[[s]])
         trimR[[s]][[i]] <- pp
       }
   }

   # multicatch block 2
      if(multicatch)#need all traps!
       trimR[[s]][[i]] <- rep(T,length(trimR[[s]][[i]]))

}
}



#
#
#
#
#
#
# likelihood function definition for no sex model
################################################################################
# NO sex
#
#
#
  msLL.nosex <- function(pv=pv, pn=pn, YY=YY, D=D, hiK=hiK, nG=nG, nK=nK, dm.den=dm.den, dm.trap=dm.trap){

    alpha0 <- array(0,dim=c(ns,hiK,2))#[session,k,BR]
    tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]

    if(pDot & !pTime){
      alpha0[,,] <- tmpP
    }
    if(pDot & pTime){
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1] <- tmpP + tmpT
      }
    }

    if(pJustsesh & !pTime){
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1] <- tmpP + tmpSS[s]
      }
    }

    if(pJustsesh & pTime){
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1] <- tmpP + tmpSS[s] + tmpT
      }
    }

    #add BR
    BRmat <- array(0, c(ns, hiK, 1))#[sess,K,sex,BR]

    if(bDot){
      BRmat[,,1] <- pv[pn%in%names.p0[grep("p.behav",names.p0)]]
    }

    if(bJustsesh){
      for(k in 1:hiK){
        BRmat[,k,1] <- pv[pn%in%names.p0[grep("p.behav.sess",names.p0)]]
        BRmat[,k,1] <- pv[pn%in%names.p0[grep("p.behav.sess",names.p0)]]
      }
    }

    alpha0[,,2] <- alpha0[,,1] + BRmat[,,1]


    #sig
    alphsig <- numeric(ns)
    if(aDot){
      tmpA <- pv[pn%in%names.sig[grep("sig.int",names.sig)]]
      for(s in 1:ns){
        alphsig[s] <- tmpA
      }
    }
    if(aJustsesh){
      tmpA <- pv[pn%in%names.sig[grep("sig.int",names.sig)]]
      tmpSS <- c(0,pv[pn%in%names.sig[grep("sig.sess",names.sig)]])
      for(s in 1:ns){
        alphsig[s] <- tmpA + tmpSS[s]
      }
    }

    alphsig <- 1/(2*exp(alphsig)^2)

    #trap betas
    if(trap.covs){
      t.beta <- matrix(NA,ns,length(t.nms))
     if("Session" %in% all.vars(model[[2]])){
      for(s in 1:ns){
        t.beta[s,] <- pv[pn%in%names.beta.trap[grep(paste("sess",s,sep=""),names.beta.trap)]]
      }
     }else{
      for(s in 1:ns){
        t.beta[s,] <- pv[pn%in%names.beta.trap]
      }
     }
    }

  #density betas
    if(DorN=="N"){
     if(dIPP){
       d.beta <- matrix(NA,ns,length(d.nms))
      if("Session" %in% all.vars(model[[1]])){
       for(s in 1:ns){
         d.beta[s,] <- pv[pn%in%names.beta.den[grep(paste("sess",s,sep=""),names.beta.den)]]
       }
      }else{
       for(s in 1:ns){
        d.beta[s,] <- pv[pn%in%names.beta.den]
       }
      }
     }
    }else{
      d.beta <- pv[pn%in%names.beta.den]
    }

  #the distance measurements
    if(distmet=="ecol"){
      dist.beta <- pv[pn%in%names.dist]
    }

  #n0
    if(n0Session) n0 <- exp(pv[pn%in%names.n0])
    if(!n0Session) n0 <- rep(exp(pv[pn%in%names.n0]),ns)

    outLik <- 0
    
    if(predict)
       preds<- list()
    
  # calculate likelihood
    for(s in 1:length(YY)){
      Ys <- YY[[s]]
      if(predict)
         preds[[s]]<- matrix(NA,nrow=nrow(Ys)+1,ncol=nrow(ssDF[[s]]))

# multicatch block 1 [looks like not needed]
     if(!multicatch){
       zeros <- array(0,c(1,dim(Ys)[2],dim(Ys)[3]))
       Ys <- abind(Ys,zeros,along=1)
     }
     if(multicatch){
       zeros <- array(0,c(1, dim(Ys)[2], dim(Ys)[3]))
       Ys <- abind(Ys,zeros,along=1)
     }

     if(distmet=="ecol"){
       cost <- exp(dm.cost[[s]] %*% exp(dist.beta)) #could exponentiate the beta?
       costR <- rasterFromXYZ(cbind(costDF[[s]][,c(1,2)],cost))
       if(is.null(PROJ)){
         projection(costR) <- '+proj=utm +zone=12 +datum=WGS84'
       }else{
         projection(costR) <- PROJ
       }
       tr <- transition(costR, transitionFunction=function(x) (1/(mean(x))), direction = directions)
       trLayer <- geoCorrection(tr, scl = F)
       D[[s]] <- costDistance(trLayer,
                              as.matrix(scrFrame$traps[[s]][,c("X","Y")]),
                              as.matrix(ssDF[[s]][,c(c("X","Y"))]))
     }

     if(smallslow){
      if(distmet == "euc"){# add a 'else' for providing your own!
        D <- list()
        D[[s]] <- e2dist(scrFrame$traps[[s]][,c("X","Y")], ssDF[[s]][,c("X","Y")])
      }
     }

     lik.marg <- rep(NA,nrow(Ys))

     #rescale the pi.s when using trimS
     if(!is.null(trimS)){
       pixels.prior <- rep(T,nG[s])
       pixels.post <- apply(D[[s]],2,min)<=trimS
       pixels <- (pixels.prior & pixels.post)
       pixels <- ifelse(pixels,1,0)
     }else{
       pixels <- rep(1,nG[s])
     }

     if(DorN=="N"){
      if(dIPP){
        d.s <- exp(dm.den[[s]] %*% d.beta[s,])
        pi.s <- (d.s*pixels)/sum(d.s*pixels)
      }
      if(dHPP){
        pi.s <- pixels/(sum(pixels))
      }
     }else{
       d.s <- exp(dm.den[[s]] %*% d.beta)
       pi.s <- (d.s*pixels)/sum(d.s*pixels)
     }

     Kern <- exp(-alphsig[s] * D[[s]]^2)

     for(i in 1:nrow(Ys)){

       #################################
       #visualize the local evaluations
       if(plotit){
         plot(ssDF[[s]][,c("X","Y")],pch=16,col="grey",cex=0.5,asp=1,
         main=paste("Session:",s," Individual: ",i," traps: ",sum(pp),sep=" "))
         points(ssDF[[s]][trimC[[s]][[i]],c("X","Y")],pch=16,col=2,cex=mycex)
         points(ssDF[[s]][trimC[[s]][[i]],],pch=16,col=2,cex=mycex)
         points(scrFrame$traps[[s]][trimR[[s]][[i]],c("X","Y")],pch=3,col=4,cex=mycex,lwd=mycex)
         points(scrFrame$traps[[s]][pp,c("X")],scrFrame$traps[[s]][pp,c("Y")],pch=16,col=3,cex=1.5)
       }
       #################################
# multicatch block 3
      if(multicatch)
        Pm <-  matrix(0,sum(trimR[[s]][[i]])+1,sum(trimC[[s]][[i]]))
      if(!multicatch)
        Pm <-  matrix(0,sum(trimR[[s]][[i]]),sum(trimC[[s]][[i]]))
      
      for(k in 1:nK[s]){
       if(pBehave){
         a0 <- alpha0[s,k,1] * (1-c(prevcap[[s]][i,,k])) + alpha0[s,k,2] * c(prevcap[[s]][i,,k])
       }else{
         a0 <- rep(alpha0[s,k,1],nrow(D[[s]]))
       }
       if(trap.covs){
         a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
       }
       if(encmod=="B")
         probcap <- c(plogis(a0[trimR[[s]][[i]]])) * Kern[trimR[[s]][[i]],trimC[[s]][[i]]]
         #probcap <- c(plogis(a0[trimR[[s]][[i]]])) * exp(-alphsig[s] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)

       if(encmod=="P")
         probcap <- c(plogis(a0[trimR[[s]][[i]]])) * Kern[trimR[[s]][[i]],trimC[[s]][[i]]]
         #probcap <- c(exp(a0[trimR[[s]][[i]]])) * exp(-alphsig[s] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)
       
# multicatch block 4
        if(!multicatch){
         if(encmod=="B"){
           probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])),1, 
                                           probcap[1:length(Pm)],log = TRUE))}
         if(encmod=="P"){
           probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])), 
                                           probcap[1:length(Pm)],log = TRUE))}

       }else{
        probcap<- rbind(probcap,rep(1, sum(trimC[[s]][[i]]) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        ###probcap<-   exp(probcap)/(1+sum(exp(probcap)))   #NOTE: need a trap mask multiplied here
        # probcap[1:length(probcap)] <- c(dbinom(rep(   c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k]) ),sum(trimC))  ,    1,
        #                                   probcap[1:length(Pm)],log = TRUE))
        vvv<-  rep(   c(Ys[i,trimR[[s]][[i]],k],1-any(Ys[i,trimR[[s]][[i]],k]>0) ),sum(trimC[[s]][[i]]))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
       }
       if(!is.null(scrFrame$trapOperation)){
         probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]],k]
       }
        Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
      }
       lik.cond <- numeric(nG[s])
#if(!is.matrix(Pm)) <--- WHAT OS THIS?
       lik.cond[trimC[[s]][[i]]] <- exp(colSums(Pm,na.rm=T))
       lik.marg[i] <- sum(lik.cond * pi.s)
       if(predict)
         preds[[s]][i,]<- (lik.cond*pi.s)/lik.marg[i]
   }

     ###Liklihood:
     if(!predict){
     ##Binomial:
      if(DorN == "N"){
        nv <- c(rep(1, length(lik.marg) - 1), n0[s])
        part1 <- lgamma((nrow(Ys)-1) + n0[s] + 1) - lgamma(n0[s] + 1)
        part2 <- sum(nv * log(lik.marg))
      }
     ##Poisson:
      if(DorN == "D"){
        nv <- c(rep(1, length(lik.marg) - 1), 1)
        atheta <- 1 - lik.marg[nrow(Ys)]
        nind <- nrow(Ys) - 1
        part1 <- nind * log(sum(d.s*pixels)) - sum(d.s*pixels) * atheta
        part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
      }
      ll <- -1 * (part1 + part2)
      outLik <- outLik + ll
     }


}  # end i loop here
if(!predict){
      out <- outLik
      return(out)
    }
    if(predict){
       return(list(preds=preds, pi.s=pi.s, ssDF=ssDF, data=YY, traps=scrFrame$traps, d.s=d.s, lik.marg=lik.marg,
                    lik.cond=lik.cond)  )
    }
  }




################################################################################
# WITH sex
#

# need to add sex specific trap level covariate coefficients.
msLL.sex <- function(pv, pn, YY, D, Y, nG, nK, hiK, dm.den, dm.trap) {

    alpha0 <- array(0, c(ns, hiK, 2, 2))#[sess,K,sex,BR]
    tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]

    if(pDot & !pTime){
      alpha0[,,,1] <- tmpP
    }

    if(pTime){
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1,1] <- tmpP + tmpT
        alpha0[s,,2,1] <- tmpP + tmpT
      }
    }

    if(pJustsex & !pTime & !pJustsesh & !pBothsexnsesh){
      tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
      for(s in 1:ns){
        alpha0[s,,1,1] <- tmpP
        alpha0[s,,2,1] <- tmpP + tmpS
      }
    }

    if(pJustsex & pTime & !pJustsesh & !pBothsexnsesh){
       tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
       tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
       for(s in 1:ns){
         alpha0[s,,1,1] <- tmpP + tmpT
         alpha0[s,,2,1] <- tmpP + tmpT + tmpS
       }
    }

    if(pJustsesh & !pTime & !pJustsex & !pBothsexnsesh){
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1,1] <- tmpP + tmpSS[s]
        alpha0[s,,2,1] <- tmpP + tmpSS[s]
      }
    }

    if(pJustsesh & pTime & !pJustsex & !pBothsexnsesh){
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1,1] <- tmpP + tmpT + tmpSS[s]
        alpha0[s,,2,1] <- tmpP + tmpT + tmpSS[s]
      }
    }

    if(pJustsesh & pJustsex & !pTime & !pBothsexnsesh){
      tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1,1] <- tmpP + tmpSS[s]
        alpha0[s,,2,1] <- tmpP + tmpSS[s] + tmpS
      }
    }

    if(pJustsesh & pJustsex & pTime & !pBothsexnsesh){
      tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1,1] <- tmpP + tmpT + tmpSS[s]
        alpha0[s,,2,1] <- tmpP + tmpT + tmpSS[s] + tmpS
      }
    }

    if(pBothsexnsesh & !pTime){
      tmpSSF <- c(0,pv[pn%in%names.p0[grep("p0.f.sess",names.p0)]])
      tmpSSM <- pv[pn%in%names.p0[grep("p0.m.sess",names.p0)]]
      for(k in 1:hiK){
        alpha0[,k,1,1] <- tmpP + tmpSSF
        alpha0[,k,2,1] <- tmpP + tmpSSM
      }
    }

    if(pBothsexnsesh & pTime){
      stop("model with time varying parameters AND a sex-session interaction not implemented")
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      for(s in 1:ns){
        alpha0[s,,1,1] <- tmpP + tmpT + tmpSS[s]
        alpha0[s,,2,1] <- tmpP + tmpT + tmpSS[s] + tmpS
      }
    }
#add BR
    BRmat <- array(0, c(ns, hiK, 2, 1))#[sess,K,sex,BR]

    if(bDot){
      BRmat[,,,1] <- pv[pn%in%names.p0[grep("p.behav",names.p0)]]
    }

    if(bJustsex){
      BRmat[,,1,1] <- pv[pn%in%names.p0[grep("p.behav.f",names.p0)]]
      BRmat[,,2,1] <- pv[pn%in%names.p0[grep("p.behav.m",names.p0)]]
    }
    if(bJustsesh){
      for(k in 1:hiK){
        BRmat[,k,1,1] <- pv[pn%in%names.p0[grep("p.behav.sess",names.p0)]]
        BRmat[,k,2,1] <- pv[pn%in%names.p0[grep("p.behav.sess",names.p0)]]
      }
    }

    if(bBothsexnsesh){
      for(k in 1:hiK){
        BRmat[,k,1,1] <- pv[pn%in%names.p0[grep("p.behav.f.sess",names.p0)]]
        BRmat[,k,2,1] <- pv[pn%in%names.p0[grep("p.behav.m.sess",names.p0)]]
      }
    }

    alpha0[,,1,2] <- alpha0[,,1,1] + BRmat[,,1,1]
    alpha0[,,2,2] <- alpha0[,,2,1] + BRmat[,,2,1]


    #sig
    alphsig <- matrix(0,ns,2)
    tmpA <- pv[pn%in%names.sig[grep("sig.int",names.sig)]]
    if(aDot){
      alphsig[,] <- tmpA
    }

    if(aJustsex & !aJustsesh){
      tmpSex <- pv[pn%in%names.sig[grep("sig.male",names.sig)]]
      alphsig[,1] <- tmpA
      alphsig[,2] <- tmpA + tmpSex
    }

    if(aJustsesh & !aJustsex){
      tmpSS <- c(0,pv[pn%in%names.sig[grep("sig.sess",names.sig)]])
      for(s in 1:ns){
       alphsig[s,1] <- tmpA + tmpSS[s]
       alphsig[s,2] <- tmpA + tmpSS[s]
     }
    }

    if(aJustsesh & aJustsex){
      tmpSS <- c(0,pv[pn%in%names.sig[grep("sig.sess",names.sig)]])
      tmpSex <- pv[pn%in%names.sig[grep("sig.male",names.sig)]]
      for(s in 1:ns){
        alphsig[s,1] <- tmpA + tmpSS[s]
        alphsig[s,2] <- tmpA + tmpSS[s] + tmpSex
      }
    }

    if(aBothsexnsesh){
      tmpSF <- c(0,pv[pn%in%names.sig[grep("sig.f.sess",names.sig)]])
      tmpSM <- pv[pn%in%names.sig[grep("sig.m.sess",names.sig)]]
      alphsig[,1] <- tmpA + tmpSF
      alphsig[,2] <- tmpA + tmpSM
    }
    alphsig <- 1/(2*exp(alphsig)^2)

  #trap betas  (no sex effect on covariates 0 could add 'Sex' for the interaction?)
    if(trap.covs){
      t.beta <- matrix(NA,ns,length(names.beta.trap))
     if("Session" %in% all.vars(model[[2]])){
      for(s in 1:ns){
        t.beta[s,] <- pv[pn%in%names.beta.trap[grep(paste("sess",s,sep=""),names.beta.trap)]]
      }
     }else{
      for(s in 1:ns){
        t.beta[s,] <- pv[pn%in%names.beta.trap]
      }
     }
    }

  #density betas
    if(DorN=="N"){
     if(dIPP){
       d.beta <- matrix(NA,ns,length(d.nms))
      if("Session" %in% all.vars(model[[1]])){
       for(s in 1:ns){
         d.beta[s,] <- pv[pn%in%names.beta.den[grep(paste("sess",s,sep=""),names.beta.den)]]
       }
      }else{
       for(s in 1:ns){
        d.beta[s,] <- pv[pn%in%names.beta.den]
       }
      }
     }
    }else{
      d.beta <- pv[pn%in%names.beta.den]
    }

   #the distance measurements
    if(distmet=="ecol"){
      dist.beta <- pv[pn%in%names.dist]
    }

   #n0
    if(n0Session) n0 <- exp(pv[pn%in%names.n0])
    if(!n0Session) n0 <- rep(exp(pv[pn%in%names.n0]),ns)

   #psi
    if(sexmod=='constant') psi.sex <- rep(plogis(pv[grep("psi",pn)]),ns)
    if(sexmod=='session') psi.sex <- plogis(pv[grep("psi",pn)])

    outLik <- 0
    
    if(predict)
       preds<- list()

 # calculate likelihood
    for(s in 1:length(YY)){
     Ys <- YY[[s]]
     if(predict)
       preds[[s]]<- matrix(NA,nrow=nrow(Ys)+1,ncol=nrow(ssDF[[s]]))

#   andy comment this out 2/11 2016
#     if(pBehave){
#       prevcap <- array(0, dim=c(dim(Ys)[1],dim(Ys)[2],dim(Ys)[3]))
#       first<- matrix(0,dim(Ys)[1],dim(Ys)[2])
#      for(i in 1:dim(Ys)[1]){
#       for(j in 1:dim(Ys)[2]){
#        if(sum(Ys[i,j,])>0){
#          first[i,j]<- min((1:(dim(Ys)[3]))[Ys[i,j,]>0] )
#          prevcap[i,j,1:first[i,j]]<- 0
#        if(first[i,j] < dim(Ys)[3])
#          prevcap[i,j,(first[i,j]+1):(dim(Ys)[3])] <- 1
#        }
#       }
#      }
#       zeros <- array(0,c(1,dim(prevcap)[2],dim(prevcap)[3]))
#       prevcap <- abind(prevcap,zeros,along=1)
#     }

# multicatch block 1 [not needed]
     if(!multicatch){
       zeros <- array(0,c(1,dim(Ys)[2],dim(Ys)[3]))
       Ys <- abind(Ys,zeros,along=1)
     }
     if(multicatch){
       zeros <- array(0,c(1, dim(Ys)[2], dim(Ys)[3]))
       Ys <- abind(Ys,zeros,along=1)
     }
     sx <- c(scrFrame$indCovs[[s]]$sex+1,NA)
     #1=female, 2=male

     if(distmet=="ecol"){
       cost <- exp(dm.cost[[s]] %*% exp(dist.beta)) #should exponentiate the beta?
       costR <- rasterFromXYZ(cbind(costDF[[s]][,c(1,2)],cost))
       if(is.null(PROJ)){
         projection(costR) <- '+proj=utm +zone=12 +datum=WGS84'
       }else{
         projection(costR) <- PROJ
       }
       tr <- transition(costR, transitionFunction=function(x) (1/(mean(x))),
                        direction = directions)
       trLayer <- geoCorrection(tr, scl = F)
       D[[s]] <- costDistance(trLayer,as.matrix(scrFrame$traps[[s]][,c("X","Y")]),as.matrix(ssDF[[s]][,c("X","Y")]))
     }
     if(smallslow){
      if(distmet == "euc"){# add a 'else' for providing your own!
        D <- list()
        D[[s]] <- e2dist(scrFrame$traps[[s]][,c("X","Y")], ssDF[[s]][,c("X","Y")])
      }
     }

      lik.marg <- lik.marg1 <- lik.marg2 <- rep(NA,nrow(Ys))
     #rescale the pi.s when using trimS
     if(!is.null(trimS)){
       pixels.prior <- rep(T,nG[s])
       pixels.post <- apply(D[[s]],2,min)<=trimS
       pixels <- (pixels.prior & pixels.post)
       pixels <- ifelse(pixels,1,0)
     }else{
       pixels <- rep(1,nG[s])
     }

     if(DorN=="N"){
      if(dIPP){
        d.s <- exp(dm.den[[s]] %*% d.beta[s,])
        pis <- (d.s*pixels)/sum(d.s*pixels)
      }
      if(dHPP){
        pi.s <- pixels/(sum(pixels))
      }
     }else{
       d.s <- exp(dm.den[[s]] %*% d.beta)
       pi.s <- (d.s*pixels)/sum(d.s*pixels)
     }

     for(i in 1:nrow(Ys)){

       #################################
       #visualize the local evaluations
       if(plotit){
         plot(ssDF[[s]][,c("X","Y")],pch=16,col="grey",cex=0.5,asp=1,
         main=paste("Session:",s," Individual: ",i," traps: ",sum(pp),sep=" "))
         points(ssDF[[s]][trimC[[s]][[i]],c("X","Y")],pch=16,col=2,cex=mycex)
         points(ssDF[[s]][trimC[[s]][[i]],],pch=16,col=2,cex=mycex)
         points(scrFrame$traps[[s]][trimR[[s]][[i]],c("X","Y")],pch=3,col=4,cex=mycex,lwd=mycex)
         points(scrFrame$traps[[s]][pp,c("X")],scrFrame$traps[[s]][pp,c("Y")],pch=16,col=3,cex=1.5)
       }
       #################################
# multicatch block 3
      if(multicatch)
       Pm <- Pm1 <- Pm2 <- matrix(0,sum(trimR[[s]][[i]])+1,sum(trimC[[s]][[i]]))
      if(!multicatch)
       Pm <- Pm1 <- Pm2 <- tmpPm <- matrix(0,sum(trimR[[s]][[i]]),sum(trimC[[s]][[i]]))

      if(!is.na(sx[i])){
       for(k in 1:nK[s]){
        if(pBehave){
          a0 <- alpha0[s,k,sx[i],1] * (1-c(prevcap[[s]][i,,k])) + alpha0[s,k,sx[i],2] * c(prevcap[[s]][i,,k])
        }else{
          a0 <- rep(alpha0[s,k,sx[i],1],nrow(D[[s]]))
        }
        if(trap.covs){
          a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
        }
        if(encmod=="B")
          probcap <- c(plogis(a0[trimR[[s]][[i]]])) * exp(-alphsig[s,sx[i]] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)

        if(encmod=="P")
          probcap <- c(exp(a0[trimR[[s]][[i]]])) * exp(-alphsig[s,sx[i]] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)

## Multicatch block 4
        if(!multicatch){
         if(encmod=="B"){
          probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])),1,
                                          probcap[1:length(Pm)],log = TRUE))}
         if(encmod=="P"){
          probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])),
                                          probcap[1:length(Pm)],log = TRUE))}
        }else{
        probcap<- rbind(probcap,rep(1, sum(trimC[[s]][[i]]) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        ###probcap<-   exp(probcap)/(1+sum(exp(probcap)))   #NOTE: need a trap mask multiplied here
        # probcap[1:length(probcap)] <- c(dbinom(rep(   c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k]) ),sum(trimC))  ,    1,
        #                                   probcap[1:length(Pm)],log = TRUE))
        vvv<-  rep(   c(Ys[i,trimR[[s]][[i]],k],1-any(Ys[i,trimR[[s]][[i]],k]>0) ),sum(trimC[[s]][[i]]))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
        }
        if(!is.null(scrFrame$trapOperation)){
          probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]],k]
        }
         Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
       }
        lik.cond <- numeric(nG[s])
        lik.cond[trimC[[s]][[i]]] <- exp(colSums(Pm,na.rm=T))
        tmpPsi <- (sx[i]==1) * (1-psi.sex[s]) + (sx[i]==2) * psi.sex[s]
        lik.marg[i] <- sum(lik.cond * pi.s) * tmpPsi
        if(predict){
           preds[[s]][i,]<- lik.cond*pi.s/lik.marg[i]
        }
      }else{
       for(k in 1:nK[s]){
        if(pBehave){
          a0.1 <- alpha0[s,k,1,1] * (1-c(prevcap[[s]][i,,k])) + alpha0[s,k,1,2] * c(prevcap[[s]][i,,k])
          a0.2 <- alpha0[s,k,2,1] * (1-c(prevcap[[s]][i,,k])) + alpha0[s,k,2,2] * c(prevcap[[s]][i,,k])
        }else{
          a0.1 <- rep(alpha0[s,k,1,1],nrow(D[[s]]))
          a0.2 <- rep(alpha0[s,k,2,1],nrow(D[[s]]))
        }
        if(trap.covs){
          a0.1 <- a0.1 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
          a0.2 <- a0.2 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
        }

# multicatch block 4 repeated
        #mixture #1
        if(encmod=="B")
          probcap <- c(plogis(a0.1[trimR[[s]][[i]]])) * exp(-alphsig[s,1] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)
        if(encmod=="P")
          probcap <- c(exp(a0.1[trimR[[s]][[i]]])) * exp(-alphsig[s,1] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)
         if(!multicatch){
          if(encmod=="B"){
           probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])),1,
                                           probcap[1:length(probcap)],log = TRUE))}
          if(encmod=="P"){
           probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])),
                                           probcap[1:length(probcap)],log = TRUE))}

         }else{
        probcap<- rbind(probcap,rep(1, sum(trimC[[s]][[i]]) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        vvv<-  rep(   c(Ys[i,trimR[[s]][[i]],k],1-any(Ys[i,trimR[[s]][[i]],k]>0) ),sum(trimC[[s]][[i]]))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
        # In Chris' code
        #   probcap<- rbind(probcap,rep(1, sum(trimC)))
        #   probcap<- probcap/colSums(probcap)
        #   probcap[1:length(probcap)] <- c(dbinom(rep(c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k])),sum(trimC)),1,
        #                                   probcap[1:length(probcap)],log = TRUE))
         }
         if(!is.null(scrFrame$trapOperation)){
           probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]],k]
         }
         Pm1[1:length(Pm1)] <- Pm1[1:length(Pm1)] + probcap[1:length(probcap)]

         #mixture #1
         if(encmod=="B")
           probcap <- c(plogis(a0.2[trimR[[s]][[i]]])) * exp(-alphsig[s,2] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)
         if(encmod=="P")
           probcap <- c(exp(a0.2[trimR[[s]][[i]]])) * exp(-alphsig[s,2] * D[[s]][trimR[[s]][[i]],trimC[[s]][[i]]]^2)

         if(!multicatch){
          if(encmod=="B"){
           probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])),1,
                                           probcap[1:length(probcap)],log = TRUE))}
          if(encmod=="P"){
           probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR[[s]][[i]],k],sum(trimC[[s]][[i]])),
                                           probcap[1:length(probcap)],log = TRUE))}

         }else{
        probcap<- rbind(probcap,rep(1, sum(trimC[[s]][[i]]) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        vvv<-  rep(   c(Ys[i,trimR[[s]][[i]],k],1-any(Ys[i,trimR[[s]][[i]],k]>0) ),sum(trimC[[s]][[i]]))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
       #    probcap<- rbind(probcap,rep(1, sum(trimC)))
       #    probcap<- probcap/colSums(probcap)
       #    probcap[1:length(probcap)] <- c(dbinom(rep(c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k])),sum(trimC)),1,
       #                                    probcap[1:length(probcap)],log = TRUE))
         }
         if(!is.null(scrFrame$trapOperation)){
           probcap <- probcap * scrFrame$trapOperation[[s]][trimR[[s]][[i]],k]
         }
         Pm2[1:length(Pm2)] <- Pm2[1:length(Pm2)] + probcap[1:length(probcap)]
       }
        lik.cond1 <- lik.cond2 <- numeric(nG[s])
        lik.cond1[trimC[[s]][[i]]] <- exp(colSums(Pm1,na.rm=T))
        lik.cond2[trimC[[s]][[i]]] <- exp(colSums(Pm2,na.rm=T))
        lik.marg1[i] <- sum(lik.cond1 * pi.s)
        lik.marg2[i] <- sum(lik.cond2 * pi.s)
        lik.marg[i]<- lik.marg1[i] * (1-psi.sex[s]) + lik.marg2[i] * psi.sex[s]
        if(predict){
          lik.cond <- (lik.cond1 * (1-psi.sex[s]) + lik.cond2 * psi.sex[s]) 
          preds[[s]][i,]<- (lik.cond*pi.s)/lik.marg[i]
        }
       }
      }  # end loop over n.individuals ... i index
 
     ###Liklihood:
     if(!predict){
     ##Binomial:
      if(DorN == "N"){
        nv <- c(rep(1, length(lik.marg) - 1), n0[s])
        part1 <- lgamma((nrow(Ys)-1) + n0[s] + 1) - lgamma(n0[s] + 1)
        part2 <- sum(nv * log(lik.marg))
      }
     ##Poisson:
      if(DorN == "D"){
        nv <- c(rep(1, length(lik.marg) - 1), 1)
        atheta <- 1 - lik.marg[nrow(Ys)]
        nind <- nrow(Ys) - 1
        part1 <- nind * log(sum(d.s*pixels)) - sum(d.s*pixels) * atheta
        part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
      }
      ll <- -1 * (part1 + part2)
      outLik <- outLik + ll
  }
   }  # end loop over sessions
    if(!predict){
      out <- outLik
      return(out)
    }
    if(predict){ 
        return(list(preds=preds, pi.s=pi.s, ssDF=ssDF, data=YY, traps=scrFrame$traps, d.s=d.s, lik.marg=lik.marg,
                    lik.cond=lik.cond)  )
     }  # end predict
  } # end msLL.sex likelihood


################################################################################
#                 Choosing and fitting the appropriate model                   #
################################################################################
## Fitting functions:
##  - msLL.nosex
##  - msLL.sex
  if(getStarts==TRUE){
    oSCR.start <- list("parameters" = pn, "values" = pv)
    return(oSCR.start)
  }else{

  if(!predict){
  message("Fitting model: D",paste(model)[1],", p0",paste(model)[2],", sigma",
           paste(model)[3],", asu",paste(model)[4], sep=" ")

  if(!anySex){
     message("Using ll function 'msLL.nosex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- suppressWarnings(
               nlm(msLL.nosex,p=pv,pn=pn,YY=YY,D=D,nG=nG,nK=nK, hiK=hiK, dm.den=dm.den,
                   dm.trap=dm.trap,hessian=T, print.level = print.level))
  }else{
     message("Using ll function 'msLL.sex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- suppressWarnings(
               nlm(msLL.sex,p=pv,pn=pn,YY=YY,D=D,nG=nG,nK=nK, hiK=hiK,
                   dm.den=dm.den,dm.trap=dm.trap,hessian=T,
                   print.level = print.level) )
  }


################################################################################
#                   Post processing of the outputs etc...                      #
################################################################################
  links <- rep(NA,length(pn))
  pars <- myfit$estimate
  links[grep("p0.int",pn)] <- "(logit)"
  links[grep("sig.int",pn)] <- "(log)"
  links[grep("n0.",pn)] <- "(log)"
  links[grep("d0.",pn)] <- "(log)"
  links[grep("psi",pn)] <- "(logit)"
  links[grep("beta",pn)] <- "(Identity)"
  trans.mle <- rep(0,length(pv))
  trans.mle[grep("p0.int",pn)] <- plogis(pars[grep("p0.int",pn)])
  trans.mle[grep("sig.int",pn)] <- exp(pars[grep("sig.int",pn)])
  trans.mle[grep("n0.",pn)] <- exp(pars[grep("n0.",pn)])
  trans.mle[grep("d0.",pn)] <- exp(pars[grep("d0.",pn)])
  trans.mle[grep("psi",pn)] <- plogis(pars[grep("psi",pn)])
  trans.mle[grep("beta",pn)] <- pars[grep("beta",pn)]
  if(pBehave){
    links[grep("pBehav",pn)] <- "(Identity)"
    trans.mle[grep("pBehav",pn)] <- pars[grep("pBehav",pn)]
  }
  #see my hessian!
  sese <- rep(rep(NA,length(pv)))
  trans.se <- rep(NA,length(pv))
  if("hessian"%in%names(myfit)){
   if(sum(myfit$hessian)!=0){
     std.err <- sqrt(diag(solve(myfit$hessian)))
   }
  }
  outStats <- data.frame( parameters=pn,
                          link = links,
                          mle=round(myfit$estimate,3),
                          std.er = std.err,
                          mle.tr = round(trans.mle,3),
                          se.tr = trans.se)
  VcV <- NULL
  if(DorN=="N"){
    ED <- (exp(pars[grep("n0.",pn)]) + unlist(lapply(scrFrame$caphist,nrow)))/areaS
  }else{
    ED <- NULL
  }
  endtime <- format(Sys.time(), "%H:%M:%S %d %b %Y") # Date stamp for start
  
  output <- list(call = cl,
                 #scrFrame = scrFrame,
                 #G = NULL,
                 #start = pv,
                 #link = c("n0 (log)", "p0 (logit)", "sig (log)"),
                 rawOutput = myfit,
                 outStats = outStats,
                 coef.mle = data.frame(param = pn, mle = myfit$estimate),
                 Area = areaS,
                 ED = ED,
                 nObs = unlist(lapply(scrFrame$caphist,nrow)),
                 mmdm = mmdm,
                 #VcV = VcV,
                 nll = myfit$minimum,
                 AIC = 2*myfit$minimum + 2* length(myfit$estimate),
                 started = starttime,
                 ended = endtime,
                 proctime = (proc.time() - ptm)[3])
   class(output) <- "oSCR.fit"
   return(output)
  }

  if(predict){
  message("Predicting model: D",paste(model)[1],", p0",paste(model)[2],", sigma",
           paste(model)[3],", cost",paste(model)[4], sep=" ")
  
# msLL.nosex <- function(pv=pv, pn=pn, YY=YY, D=D, hiK=hiK, nG=nG, nK=nK, dm.den=dm.den, dm.trap=dm.trap){

  if(!anySex){
     message("Using ll function 'msLL.nosex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- msLL.nosex(p=start.vals,pn=pn,YY=YY,D=D, 
                  hiK=hiK, nG=nG,nK=nK, dm.den=dm.den,dm.trap=dm.trap)
 
 }else{
     message("Using ll function 'msLL.sex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- msLL.sex(p=start.vals,pn=pn,YY=YY,D=D,nK=nK,
              nG=nG,hiK=hiK, dm.den=dm.den,dm.trap=dm.trap)
  }
  return(myfit)
  }
 }
}
