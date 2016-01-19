oSCR.fit <-
function(scrFrame, model = list(D~1, p0~1, a1~1, path~1), ssDF = NULL, costDF = NULL,
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
#                            Setting thing up!                                 #
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
    #trimS <- 6*max.dist
#  }
   mmdm <- mean(max.dist,na.rm=T)
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
  pars.a1 <- NULL ; names.a1 <- NULL
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
    ssDF <- list()
   for(i in 1:length(scrFrame$traps)){ # make a state space for every set of traps
    dd <- as.matrix(dist(scrFrame$traps[[i]]))
    dd[dd==0] <- NA
    nnd <- mean(apply(dd,1,min,na.rm=T))
    ssDF[[i]] <- expand.grid(seq(min(scrFrame$traps[[i]][,1])-4*nnd,
                                 max(scrFrame$traps[[i]][,1])+4*nnd,nnd/4),
                             seq(min(scrFrame$traps[[i]][,2])-4*nnd,
                                 max(scrFrame$traps[[i]][,2])+4*nnd,nnd/4))
   }
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
  allvars.a1a <- all.vars(model[[3]])
  allvars.a1 <- allvars.a1a[!allvars.a1a=="a1"]
  var.p0.1 <- "sex" %in% allvars.p0
  var.p0.2 <- "session" %in% allvars.p0
  var.p0.3 <- "t" %in% allvars.p0
  var.a1.1 <- "sex" %in% allvars.a1
  var.a1.2 <- "session" %in% allvars.a1
  pBehave <- "b" %in% all.vars(model[[2]])

  allvars.dist <- all.vars(model[[4]])
  allvars.dist <- allvars.dist[!allvars.dist=="path"]

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
#  var.p0.1 <- "sex" %in% allvars.p0
#  var.p0.2 <- "session" %in% allvars.p0
#

  if("indCovs" %in% names(scrFrame)){
   if("sex" %in% names(scrFrame$indCovs[[1]])){
     anySex <- TRUE
   }
  }

  tmp.a0.name1 <- "p0.int"
  tmp.a0.name2 <- ifelse(var.p0.1,"p0.male",NA)
  if(var.p0.2){
   if(ns>1){
     tmp.a0.name3 <- paste0("p0.sess",2:ns)
   }
  }else{
    tmp.a0.name3 <- NA
  }
  if(var.p0.3){
    tmp.a0.name4 <- paste0("p0.t",2:hiK)
    pTime <- TRUE
  }else{
    tmp.a0.name4 <- NA
  }
  names.p0 <- c(tmp.a0.name1,tmp.a0.name2,tmp.a0.name3,tmp.a0.name4)
  names.p0 <- names.p0[!is.na(names.p0)]
  pars.p0 <- rep(0,length(names.p0))
  pars.p0[1] <- qlogis(0.05)

  if(var.p0.1 && var.p0.2){
    pBothsexnsesh <- TRUE
  }else{
   if(var.p0.1){
     pJustsex <- TRUE
   }else{
    if(var.p0.2){
      pJustsesh <- TRUE
    }else{
      pDot <- TRUE
    }
   }
  }


#  if(var.p0.1 && var.p0.2){
#    pars.p0 <- rnorm(ns*2,qlogis(0.1),0.2)#"p.ss"
#    tmpPsex <- rep(c(1,2),ns)
#    tmpPsess <- rep(1:ns,each=2)
#    names.p0 <- paste("p0.sex",tmpPsex,"session",tmpPsess,sep="")
#    pBothsexnsesh <- TRUE
#  }else{
#   if(var.p0.1){
#     pars.p0 <- rnorm(2,qlogis(0.1),0.2)#"p.sex"
#     names.p0 <- c("p0.sex1","p0.sex2")
#     pJustsex <- TRUE
#   }else{
#    if(var.p0.2){
#      pars.p0 <- rnorm(ns,qlogis(0.1),0.2)#"p.ses"
#      tmpPsess <- 1:ns
#      names.p0 <-  paste("p0.session",tmpPsess,sep="")
#      pJustsesh <- TRUE
#    }else{
#      pars.p0 <- rnorm(1,qlogis(0.1),0.2)#"p."
#      names.p0 <- c("p0.")
#      pDot <- TRUE
#    }
#   }
#  }
#
  if(any(var.p0.1, var.a1.1) && !anySex)
   stop("Sex defined in a model but no sex data provided.")


################################################################################
# p0 can vary by occasion!
#

### To do:
###  - relax k_g = K
###  - make the 'T' = trend work
#  var.p0.t <- "t" %in% allvars.p0
#  var.p0.T <- "T" %in% allvars.p0
#  if(var.p0.t){
#    pars.p0 <- rep(pars.p0,hiK)
#    names.p0 <- paste(rep(names.p0,each=hiK),"..t",1:hiK,sep="")
#    pTime <- TRUE
#  }
#  if(var.p0.T){
#    pars.p0 <- c(pars.p0,0)
#    names.p0 <- c(names.p0,"T.trend") ## This need working up
#  }
  if(pBehave){
    pars.p0 <- c(pars.p0,0)
    names.p0 <- c(names.p0,"p.behav")
  }


################################################################################
# define sex and/or session specific parameters
# NB: a1 = 1/(2 * sigma^2) | sigma = sqrt(1/2*a1)
#

######### CS: attempt to fix sex:session thing
  #var.a1 = sex var.a2 = session
  tmp.a1.name1 <- "a1.int"
  tmp.a1.name2 <- ifelse(var.a1.1,"a1.male",NA)
  if(var.a1.2){
   if(ns>1){
     tmp.a1.name3 <- paste0("a1.sess",2:ns)
   }
  }else{
    tmp.a1.name3 <- NA
  }
  names.a1 <- c(tmp.a1.name1,tmp.a1.name2,tmp.a1.name3)
  names.a1 <- names.a1[!is.na(names.a1)]
  pars.a1 <- rep(0,length(names.a1))
  pars.a1[1] <- log(1/(mmdm^2))
  
  if(var.a1.1 && var.a1.2){
    aBothsexnsesh <- TRUE
  }else{
   if(var.a1.1){
     aJustsex <- TRUE
   }else{
    if(var.a1.2){
      aJustsesh <- TRUE
    }else{
      aDot <- TRUE
    }
   }
  }
  
############################################### 
#  if(var.a1.1 && var.a1.2){
#    pars.a1 <- rnorm(ns*2,0,0.2)#"a1.ss"
#    tmpAsex <- rep(c(1,2),ns)
#    tmpAsess <- rep(1:ns,each=2)
#    names.a1 <- paste("a1.sex",tmpAsex,"session",tmpAsess,sep="")
#    aBothsexnsesh <- TRUE
#  }else{
#   if(var.a1.1){
#     pars.a1 <- rnorm(2,0,0.2)#"a1.sex"
#     names.a1 <- c("a1.sex1","a1.sex2")
#     aJustsex <- TRUE
#   }else{
#    if(var.a1.2){
#      pars.a1 <- rnorm(ns,0.1,0.2)#"a1.sess"
#      tmpAsess <- 1:ns
#      names.a1 <- paste("a1.session",tmpAsess,sep="")
#      aJustsesh <- TRUE
#    }else{
#      cnames.a1 <- c("a1.")
#      aDot <- TRUE
#    }
#   }
#  }
#


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
  pv <- round(c(pars.p0,pars.a1,pars.beta.trap,pars.beta.den,pars.dist,pars.n0,pars.sex),2)
  pn <- c(names.p0,names.a1,names.beta.trap,names.beta.den,names.dist,names.n0,names.sex)
  if(!is.null(start.vals)){
   if(length(pv)==length(start.vals)){
     pv <- start.vals
   }else{
   message("The starting values provided are not the required length, using generated values")
   }
  }

################################################################################
#     Likelihood functions - these can probably be dumped somewhere else?      #
################################################################################


################################################################################
# NO sex  
#

  msLL.nosex <- function(pv=pv, pn=pn, YY=YY, D=D, hiK=hiK, nG=nG, nK=nK, dm.den=dm.den, dm.trap=dm.trap){

  #p0
    alpha0 <- array(NA,dim=c(ns,hiK,2))
    if(pDot & !pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpPB <- ifelse(pBehave,pv[pn%in%names.p0[grep("p.behav",names.p0)]],0)
      alpha0[,,1] <- tmpP
      alpha0[,,2] <- alpha0[,,1] + tmpPB
    }
    if(pDot & pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      tmpPB <- ifelse(pBehave,pv[pn%in%names.p0[grep("p.behav",names.p0)]],0)
     for(s in 1:ns){
      alpha0[s,,1] <- tmpP + tmpT
      alpha0[s,,2] <- alpha0[s,,1] + tmpPB
     }
    }
    if(pJustsesh & !pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      tmpPB <- ifelse(pBehave,pv[pn%in%names.p0[grep("p.behav",names.p0)]],0)
     for(s in 1:ns){
       alpha0[s,,1] <- tmpP + tmpSS[s]
       alpha0[s,,2] <- alpha0[s,,1] + tmpPB
     }
    }
    if(pJustsesh & pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      tmpPB <- ifelse(pBehave,pv[pn%in%names.p0[grep("p.behav",names.p0)]],0)
     for(s in 1:ns){
       alpha0[s,,1] <- tmpP + tmpSS[s] + tmpT
       alpha0[s,,2] <- alpha0[s,,1] + tmpPB
     }
    }

  #a1
    alpha1 <- numeric(ns)
    if(aDot){
      tmpA <- pv[pn%in%names.a1[grep("a1.int",names.a1)]]
      for(s in 1:ns){
        alpha1[s] <- exp(tmpA)
      }
    }
    if(aJustsesh){
      tmpA <- pv[pn%in%names.a1[grep("a1.int",names.a1)]]
      tmpSS <- c(0,pv[pn%in%names.a1[grep("a1.sess",names.a1)]])
      for(s in 1:ns){
        alpha1[s] <- exp(tmpA + tmpSS[s])
      }
    }

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

  # calculate likelihood
    for(s in 1:length(YY)){
      Ys <- YY[[s]]
     if(pBehave){
       prevcap <- array(0, dim=c(dim(Ys)[1],dim(Ys)[2],dim(Ys)[3]))
       first<- matrix(0,dim(Ys)[1],dim(Ys)[2])
      for(i in 1:dim(Ys)[1]){
       for(j in 1:dim(Ys)[2]){
        if(sum(Ys[i,j,])>0){
          first[i,j]<- min((1:(dim(Ys)[3]))[Ys[i,j,]>0] )
          prevcap[i,j,1:first[i,j]]<- 0
        if(first[i,j] < dim(Ys)[3])
          prevcap[i,j,(first[i,j]+1):(dim(Ys)[3])] <- 1
        }
       }
      }
       zeros <- array(0,c(1,dim(prevcap)[2],dim(prevcap)[3]))
       prevcap <- abind(prevcap,zeros,along=1)
     }

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

     for(i in 1:nrow(Ys)){
      if(is.null(trimS)){
        pp <- rep(T,ncol(Ys))
        trimC <- rep(T,nG[s])
        trimR <- pp
      }else{
       if(i<nrow(Ys)){
         pp <- apply(Ys[i,,],1,sum)>0
         trimR <- apply(
                   rbind(
                    rep(trimS+2,nrow(scrFrame$traps[[s]])),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          scrFrame$traps[[s]][,c("X","Y")])),2,min)<=(2*trimS)
         trimC <- apply(
                   rbind(
                    rep(trimS+2,nG[s]),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          ssDF[[s]][,c("X","Y")])),2,min,na.rm=T)<=trimS
       }else{
         pp <- rep(T,ncol(Ys))
         trimC <- apply(
                   rbind(
                    rep(trimS+2,nG[s]),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          ssDF[[s]][,c("X","Y")])),2,min,na.rm=T)<=trimS
         #trimC <- rep(T,nG[[s]])
         trimR <- pp
       }
      }
# multicatch block 2
      if(multicatch)#need all traps!
       trimR <- rep(T,length(trimR))

       #################################
       #visualize the local evaluations
       if(plotit){
         plot(ssDF[[s]][,c("X","Y")],pch=16,col="grey",cex=0.5,asp=1,
         main=paste("Session:",s," Individual: ",i," traps: ",sum(pp),sep=" "))
         points(ssDF[[s]][trimC,c("X","Y")],pch=16,col=2,cex=mycex)
         points(ssDF[[s]][trimC,],pch=16,col=2,cex=mycex)
         points(scrFrame$traps[[s]][trimR,c("X","Y")],pch=3,col=4,cex=mycex,lwd=mycex,)
         points(scrFrame$traps[[s]][pp,c("X")],scrFrame$traps[[s]][pp,c("Y")],pch=16,col=3,cex=1.5)
       }
       #################################
# multicatch block 3
       if(multicatch)
        Pm <-  matrix(0,sum(trimR)+1,sum(trimC))
      if(!multicatch)
        Pm <-  matrix(0,sum(trimR),sum(trimC))
      for(k in 1:nK[s]){
       if(pBehave){
         a0 <- alpha0[s,k,1] * (1-c(prevcap[i,,k])) + alpha0[s,k,2] * c(prevcap[i,,k])
       }else{
         a0 <- rep(alpha0[s,k,1],nrow(D[[s]]))
       }
       if(trap.covs){
         a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
       }
       if(encmod=="B")
         probcap <- c(plogis(a0[trimR])) * exp(-alpha1[s] * D[[s]][trimR,trimC]^2)
       if(encmod=="P")
         probcap <- c(exp(a0[trimR])) * exp(-alpha1[s] * D[[s]][trimR,trimC]^2)

# multicatch block 4
        if(!multicatch){
         if(encmod=="B"){
           probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR,k],sum(trimC)),1,
                                           probcap[1:length(Pm)],log = TRUE))}
         if(encmod=="P"){
           probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR,k],sum(trimC)),
                                           probcap[1:length(Pm)],log = TRUE))}
                                           
       }else{
        probcap<- rbind(probcap,rep(1, sum(trimC) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        ###probcap<-   exp(probcap)/(1+sum(exp(probcap)))   #NOTE: need a trap mask multiplied here
        # probcap[1:length(probcap)] <- c(dbinom(rep(   c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k]) ),sum(trimC))  ,    1,
        #                                   probcap[1:length(Pm)],log = TRUE))
        vvv<-  rep(   c(Ys[i,trimR,k],1-any(Ys[i,trimR,k]>0) ),sum(trimC))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
       }
       if(!is.null(scrFrame$trapOperation)){
         probcap <- probcap * scrFrame$trapOperation[[s]][trimR,k]
       }
        Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
      }
       lik.cond <- numeric(nG[s])
if(!is.matrix(Pm)) browser()
       lik.cond[trimC] <- exp(colSums(Pm,na.rm=T))
       lik.marg[i] <- sum(lik.cond * pi.s)
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
     if(predict){ #needs to be ammended - currently WRONG!
       tmp.post <- matrix(NA,nG[s],nrow(Ys))
      for(i in 1:nrow(Ys)){
        Pm <-  matrix(0,length(trimR),length(trimC))
       for(k in 1:dim(Ys)[3]){
        probcap <- c(plogis(a0)) * exp(-alpha1[s] * D[[s]]^2)
        probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,,k],nG[s]),1,probcap[1:length(Pm)],log = TRUE))
        Pm[1:length(tmpPm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
       }
        lik.cond <- exp(colSums(Pm))
        tmp.post[,i]<- (lik.cond*(1/nG[s]))/lik.marg[i]
      }
      posterior[[s]] <- cbind(ssDF[[s]][,c("X","Y")],tmp.post)
     }
    }
    if(!predict){
      out <- outLik
      return(out)
    }
    if(predict){
      return(posterior)
    }
  }


################################################################################
# WITH sex
#

# need to add sex specific trap level covariate coefficients.
msLL.sex <- function(pv, pn, YY, D, Y, nG, nK, hiK, dm.den, dm.trap) {

    alpha0 <- array(NA, c(ns, hiK, 2, 2))
    tmpPB <- ifelse(pBehave,pv[pn%in%names.p0[grep("p.behav",names.p0)]],0)
    if(pDot & !pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP
       alpha0[s,,2,1] <- tmpP
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB
     }
    }
    if(pDot & pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP + tmpT
       alpha0[s,,2,1] <- tmpP + tmpT
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB
     }
    }
    if(pJustsex & !pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpS <- c(0,pv[pn%in%names.p0[grep("p0.male",names.p0)]])
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,1] <- tmpP + tmpS
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB + tmpS
     }
    }
    if(pJustsex & pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP + tmpT
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,1] <- tmpP + tmpT + tmpS
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB
     }
    }
    if(pJustsesh & !pTime){ # no sex here but still sex likelihood!
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP + tmpSS[s]
       alpha0[s,,2,1] <- alpha0[s,,1,1]
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB
     }
    }
    if(pJustsesh & pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP + tmpT + tmpSS[s]
       alpha0[s,,2,1] <- alpha0[s,,1,1]
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB
     }
    }
    if(pBothsexnsesh & !pTime){ # no sex here but still sex likelihood!
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP + tmpSS[s] 
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,1] <- tmpP + tmpSS[s] + tmpS
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB
     }
    }
    if(pBothsexnsesh & pTime){
      tmpP <- pv[pn%in%names.p0[grep("p0.int",names.p0)]]
      tmpSS <- c(0,pv[pn%in%names.p0[grep("p0.sess",names.p0)]])
      tmpS <- pv[pn%in%names.p0[grep("p0.male",names.p0)]]
      tmpT <- c(0,pv[pn%in%names.p0[grep("p0.t",names.p0)]])
     for(s in 1:ns){
       alpha0[s,,1,1] <- tmpP + tmpT + tmpSS[s]
       alpha0[s,,1,2] <- alpha0[s,,1,1] + tmpPB
       alpha0[s,,2,1] <- tmpP + tmpT + tmpSS[s] + tmpS
       alpha0[s,,2,2] <- alpha0[s,,2,1] + tmpPB
     }
    }

  #a1
    alpha1 <- matrix(NA,ns,2)
    if(aDot){
      tmpA <- pv[pn%in%names.a1[grep("a1.int",names.a1)]]
      alpha1[] <- tmpA
    }
    if(aJustsex){
#      alpha1[,1] <- exp(pv[pn%in%names.a1[grep("sex1",names.a1)]])
#      alpha1[,2] <- exp(pv[pn%in%names.a1[grep("sex2",names.a1)]])
      tmpA <- pv[pn%in%names.a1[grep("a1.int",names.a1)]]
      tmpS <- pv[pn%in%names.a1[grep("a1.male",names.a1)]]
      alpha1[,1] <- tmpA
      alpha1[,2] <- alpha1[,1] + tmpS
    }
    if(aJustsesh){ # no sex here but still sex likelihood!
      tmpA <- pv[pn%in%names.a1[grep("a1.int",names.a1)]]
      tmpSS <- c(0,pv[pn%in%names.a1[grep("a1.sess",names.a1)]])
     for(s in 1:ns){
       alpha1[s,1] <- tmpA + tmpSS[s]
       alpha1[s,2] <- alpha1[s,1]
     }
    }
    if(aBothsexnsesh){
      tmpA <- pv[pn%in%names.a1[grep("a1.int",names.a1)]]
      tmpS <- pv[pn%in%names.a1[grep("a1.male",names.a1)]]
      tmpSS <- c(0,pv[pn%in%names.a1[grep("a1.sess",names.a1)]])
     for(s in 1:ns){
#      alpha1[s,1] <- exp(pv[pn%in%names.a1[grep(paste("sex1session",s,sep=""),names.a1)]])
#      alpha1[s,2] <- exp(pv[pn%in%names.a1[grep(paste("sex2session",s,sep=""),names.a1)]])
      alpha1[s,1] <- tmpA + tmpSS[s]
      alpha1[s,2] <- alpha1[s,1] + tmpS
     }
    }
    alpha1 <- exp(alpha1)

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

 # calculate likelihood
    for(s in 1:length(YY)){
     Ys <- YY[[s]]
     if(predict)
      tmp.post <- matrix(NA,nG[s],nrow(Ys))

     if(pBehave){
       prevcap <- array(0, dim=c(dim(Ys)[1],dim(Ys)[2],dim(Ys)[3]))
       first<- matrix(0,dim(Ys)[1],dim(Ys)[2])
      for(i in 1:dim(Ys)[1]){
       for(j in 1:dim(Ys)[2]){
        if(sum(Ys[i,j,])>0){
          first[i,j]<- min((1:(dim(Ys)[3]))[Ys[i,j,]>0] )
          prevcap[i,j,1:first[i,j]]<- 0
        if(first[i,j] < dim(Ys)[3])
          prevcap[i,j,(first[i,j]+1):(dim(Ys)[3])] <- 1
        }
       }
      }
       zeros <- array(0,c(1,dim(prevcap)[2],dim(prevcap)[3]))
       prevcap <- abind(prevcap,zeros,along=1)
     }
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
           if(is.null(trimS)){
        pp <- rep(T,ncol(Ys))
        trimC <- rep(T,nG[s])
        trimR <- pp
      }else{
       if(i<nrow(Ys)){
         pp <- apply(Ys[i,,],1,sum)>0
         trimR <- apply(
                   rbind(
                    rep(trimS+2,nrow(scrFrame$traps[[s]])),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          scrFrame$traps[[s]][,c("X","Y")])),2,min)<=(2*trimS)
         trimC <- apply(
                   rbind(
                    rep(trimS+2,nG[s]),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          ssDF[[s]][,c("X","Y")])),2,min,na.rm=T)<=trimS
       }else{
         pp <- rep(T,ncol(Ys))
         trimC <- apply(
                   rbind(
                    rep(trimS+2,nG[s]),
                     e2dist(matrix(unlist(scrFrame$traps[[s]][pp,c("X","Y")]),sum(pp),2),
                                          ssDF[[s]][,c("X","Y")])),2,min,na.rm=T)<=trimS
         #trimC <- rep(T,nG[[s]])
         trimR <- pp
       }
      }
# multicatch block 2
      if(multicatch)#need all traps!
       trimR <- rep(T,length(trimR))

       #################################
       #visualize the local evaluations
       if(plotit){
         plot(ssDF[[s]][,c("X","Y")],pch=16,col="grey",cex=0.5,asp=1,
         main=paste("Session:",s," Individual: ",i," traps: ",sum(pp),sep=" "))
         points(ssDF[[s]][trimC,c("X","Y")],pch=16,col=2,cex=mycex)
         points(ssDF[[s]][trimC,],pch=16,col=2,cex=mycex)
         points(scrFrame$traps[[s]][trimR,c("X","Y")],pch=3,col=4,cex=mycex,lwd=mycex,)
         points(scrFrame$traps[[s]][pp,c("X")],scrFrame$traps[[s]][pp,c("Y")],pch=16,col=3,cex=1.5)
       }
       #################################
# multicatch block 3
      if(multicatch)
       Pm <- Pm1 <- Pm2 <- matrix(0,sum(trimR)+1,sum(trimC))
      if(!multicatch)
       Pm <- Pm1 <- Pm2 <- tmpPm <- matrix(0,sum(trimR),sum(trimC))

      if(!is.na(sx[i])){
       for(k in 1:nK[s]){
        if(pBehave){
          a0 <- alpha0[s,k,sx[i],1] * (1-c(prevcap[i,,k])) + alpha0[s,k,sx[i],2] * c(prevcap[i,,k])
        }else{
          a0 <- rep(alpha0[s,k,sx[i],1],nrow(D[[s]]))
        }
        if(trap.covs){
          a0 <- a0 + (dm.trap[[s]][[k]] %*% c(t.beta[s,]))
        }
        if(encmod=="B")
          probcap <- c(plogis(a0[trimR])) * exp(-alpha1[s,sx[i]] * D[[s]][trimR,trimC]^2)
        if(encmod=="P")
          probcap <- c(exp(a0[trimR])) * exp(-alpha1[s,sx[i]] * D[[s]][trimR,trimC]^2)

## Multicatch block 4
        if(!multicatch){
         if(encmod=="B"){
          probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR,k],sum(trimC)),1,
                                          probcap[1:length(Pm)],log = TRUE))}
         if(encmod=="P"){
          probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR,k],sum(trimC)),
                                          probcap[1:length(Pm)],log = TRUE))}
        }else{
        probcap<- rbind(probcap,rep(1, sum(trimC) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        ###probcap<-   exp(probcap)/(1+sum(exp(probcap)))   #NOTE: need a trap mask multiplied here
        # probcap[1:length(probcap)] <- c(dbinom(rep(   c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k]) ),sum(trimC))  ,    1,
        #                                   probcap[1:length(Pm)],log = TRUE))
        vvv<-  rep(   c(Ys[i,trimR,k],1-any(Ys[i,trimR,k]>0) ),sum(trimC))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
        }
        if(!is.null(scrFrame$trapOperation)){
          probcap <- probcap * scrFrame$trapOperation[[s]][trimR,k]
        }
         Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
       }
        lik.cond <- numeric(nG[s])
        lik.cond[trimC] <- exp(colSums(Pm,na.rm=T))
        tmpPsi <- (sx[i]==1) * (1-psi.sex[s]) + (sx[i]==2) * psi.sex[s]
        lik.marg[i] <- sum(lik.cond * pi.s) * tmpPsi
        if(predict){
          tmp.post[,i] <- (lik.cond*pi.s)/lik.marg[i]
        }
      }else{
       for(k in 1:nK[s]){
        if(pBehave){
          a0.1 <- alpha0[s,k,1,1] * (1-c(prevcap[i,,k])) + alpha0[s,k,1,2] * c(prevcap[i,,k])
          a0.2 <- alpha0[s,k,2,1] * (1-c(prevcap[i,,k])) + alpha0[s,k,2,2] * c(prevcap[i,,k])
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
          probcap <- c(plogis(a0.1[trimR])) * exp(-alpha1[s,1] * D[[s]][trimR,trimC]^2)
        if(encmod=="P")
          probcap <- c(exp(a0.1[trimR])) * exp(-alpha1[s,1] * D[[s]][trimR,trimC]^2)
         if(!multicatch){
          if(encmod=="B"){
           probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR,k],sum(trimC)),1,
                                           probcap[1:length(probcap)],log = TRUE))}
          if(encmod=="P"){
           probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR,k],sum(trimC)),
                                           probcap[1:length(probcap)],log = TRUE))}
                                           
         }else{
        probcap<- rbind(probcap,rep(1, sum(trimC) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        vvv<-  rep(   c(Ys[i,trimR,k],1-any(Ys[i,trimR,k]>0) ),sum(trimC))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
        # In Chris' code
        #   probcap<- rbind(probcap,rep(1, sum(trimC)))
        #   probcap<- probcap/colSums(probcap)
        #   probcap[1:length(probcap)] <- c(dbinom(rep(c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k])),sum(trimC)),1,
        #                                   probcap[1:length(probcap)],log = TRUE))
         }
         if(!is.null(scrFrame$trapOperation)){
           probcap <- probcap * scrFrame$trapOperation[[s]][trimR,k]
         }
         Pm1[1:length(Pm1)] <- Pm1[1:length(Pm1)] + probcap[1:length(probcap)]

         #mixture #1
         if(encmod=="B")
           probcap <- c(plogis(a0.2[trimR])) * exp(-alpha1[s,2] * D[[s]][trimR,trimC]^2)
         if(encmod=="P")
           probcap <- c(exp(a0.2[trimR])) * exp(-alpha1[s,2] * D[[s]][trimR,trimC]^2)

         if(!multicatch){
          if(encmod=="B"){
           probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,trimR,k],sum(trimC)),1,
                                           probcap[1:length(probcap)],log = TRUE))}
          if(encmod=="P"){
           probcap[1:length(probcap)] <- c(dpois(rep(Ys[i,trimR,k],sum(trimC)),
                                           probcap[1:length(probcap)],log = TRUE))}
                                           
         }else{
        probcap<- rbind(probcap,rep(1, sum(trimC) ) )
        probcap<- t(t(probcap)/colSums(probcap))
        vvv<-  rep(   c(Ys[i,trimR,k],1-any(Ys[i,trimR,k]>0) ),sum(trimC))      # I think this and below is computing the multinomial likelihood
        vvv[vvv==1]<-    log( probcap[1:length(Pm)][vvv==1] )
        probcap[1:length(Pm)]<- vvv
       #    probcap<- rbind(probcap,rep(1, sum(trimC)))
       #    probcap<- probcap/colSums(probcap)
       #    probcap[1:length(probcap)] <- c(dbinom(rep(c(Ys[i,trimR,k],1-sum(Ys[i,trimR,k])),sum(trimC)),1,
       #                                    probcap[1:length(probcap)],log = TRUE))
         }
         if(!is.null(scrFrame$trapOperation)){
           probcap <- probcap * scrFrame$trapOperation[[s]][trimR,k]
         }
         Pm2[1:length(Pm2)] <- Pm2[1:length(Pm2)] + probcap[1:length(probcap)]
       }
        lik.cond1 <- lik.cond2 <- numeric(nG[s])
        lik.cond1[trimC] <- exp(colSums(Pm1,na.rm=T))
        lik.cond2[trimC] <- exp(colSums(Pm2,na.rm=T))
        lik.marg1[i] <- sum(lik.cond1 * pi.s)
        lik.marg2[i] <- sum(lik.cond2 * pi.s)
        lik.marg[i]<- lik.marg1[i] * (1-psi.sex[s]) + lik.marg2[i] * psi.sex[s]
       }
      }
      if(predict)
       posterior[[s]] <- cbind(ssDF[[s]][,c("X","Y")],tmp.post)

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
     if(predict){ #needs to be ammended - currently WRONG!
       tmp.post <- matrix(NA,nG[s],nrow(Ys))
      for(i in 1:nrow(Ys)){
        Pm <-  matrix(0,length(trimR),length(trimC))
       for(k in 1:dim(Ys)[3]){
         probcap <- c(plogis(a0)) * exp(-alpha1[s] * D[[s]]^2)
         probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i,,k],nG[s]),1,probcap[1:length(Pm)],log = TRUE))
         Pm[1:length(tmpPm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
       }
        lik.cond <- exp(colSums(Pm))
        tmp.post[,i]<- (lik.cond*(1/nG[s]))/lik.marg[i]
      }
      posterior[[s]] <- cbind(ssDF[[s]][,c("X","Y")],tmp.post)
     }
   }
    if(!predict){
      out <- outLik
      return(out)
    }
    if(predict){
      return(posterior)
    }
  }


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
           paste(model)[3],", cost",paste(model)[4], sep=" ")

  if(!anySex){
     message("Using ll function 'msLL.nosex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- nlm(msLL.nosex,p=pv,pn=pn,YY=YY,D=D,nG=nG,nK=nK, hiK=hiK, dm.den=dm.den,
                  dm.trap=dm.trap,hessian=T, print.level = print.level)
  }else{
     message("Using ll function 'msLL.sex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- nlm(msLL.sex,p=pv,pn=pn,YY=YY,D=D,nG=nG,nK=nK, hiK=hiK,
                  dm.den=dm.den,dm.trap=dm.trap,hessian=T,
                  print.level = print.level)
  }


  
################################################################################
#                   Post processing of the outputs etc...                      #
################################################################################

  links <- rep(NA,length(pn))
  pars <- myfit$estimate
  links[grep("p0.int",pn)] <- "(logit)"
  links[grep("a1.int",pn)] <- "(log)"
  links[grep("n0.",pn)] <- "(log)"
  links[grep("d0.",pn)] <- "(exp)"
  links[grep("psi",pn)] <- "(logit)"
  links[grep("beta",pn)] <- "(Identity)"
  trans.mle <- rep(0,length(pv))
  trans.mle[grep("p0.int",pn)] <- plogis(pars[grep("p0.int",pn)])
  trans.mle[grep("a1.int",pn)] <- exp(pars[grep("a1.int",pn)])
  trans.mle[grep("n0.",pn)] <- exp(pars[grep("n0.",pn)])
  trans.mle[grep("d0.",pn)] <- exp(pars[grep("d0.",pn)])
  trans.mle[grep("psi",pn)] <- plogis(pars[grep("psi",pn)])
  trans.mle[grep("beta",pn)] <- pars[grep("beta",pn)]
  if(pBehave){
    links[grep("pBehav",pn)] <- "(Identity)"
    trans.mle[grep("pBehav",pn)] <- pars[grep("pBehav",pn)]
  }
  #see my hessian!
    if("hessian"%in%names(myfit)){
#    sese <- sqrt(diag(solve(myfit$hessian)))
    # add 95%CI's here
    #  trans.se[which(pn %in% "p0.")] <- plogis(pars[which(pn %in% "p0.")])
    #  trans.se[which(pn %in% "a1.")] <- exp(pars[which(pn %in% "a1.")])
    #  trans.se[which(pn %in% "n0")] <- exp(pars[which(pn %in% "n0")])
    trans.se <- rep(NA,length(pv))
  }else{
    sese <- rep(rep(NA,length(pv)))
    trans.se <- rep(NA,length(pv))
  }
  outStats <- data.frame( parameters=pn,
                          link = links,
                          mle=round(myfit$estimate,3),
                          #se = sese,
                          mle.tr = round(trans.mle,3),
                          se.tr = trans.se)
  VcV <- NULL
  idx <- outStats[,"parameters"] == "a1.int"
  sigma <- sqrt(1/(2*exp(outStats[idx,"mle"])))
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
                 #link = c("n0 (log)", "p0 (logit)", "a1 (log)"),
                 rawOutput = myfit,
                 outStats = outStats,
                 Area = areaS,
                 ED = ED,
                 nObs = unlist(lapply(scrFrame$caphist,nrow)),
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

  if(!anySex){
     message("Using ll function 'msLL.nosex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- msLL.nosex(p=start.vals,pn=pn,YY=YY,D=D,
              nG=nG,K=K,hiK=hiK, dm.den=dm.den,dm.trap=dm.trap)
  }else{
     message("Using ll function 'msLL.sex' \nHold on tight!")
     message(paste(pn," ",sep=" | "))
     message(" ")
     myfit <- msLL.sex(p=start.vals,pn=pn,YY=YY,D=D,
              nG=nG,K=K, hiK=hiK, dm.den=dm.den,dm.trap=dm.trap)
  }
  return(myfit)
  }
 }
}


