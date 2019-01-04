# process telemetry data into raster cell frequencies and also compute
# "sbar" -- individual activity center estimates
telemetry.processor <- function(ssDF, teldata){
  library(FNN)
  nsess<- length(teldata)
  sbar<- list()
  dlst<- list()
  nfreq<- list()
  for(s in 1:nsess){
    s.grid <- as.vector(get.knnx(ssDF[[s]][,c("X","Y")],teldata[[s]][,c("X","Y")],1)$nn.index)
    sbar[[s]]<- aggregate(teldata[[s]][,2:3],list(teldata[[s]][,1]),mean)
    tmp<-table(s.grid, teldata[[s]][,1])
    tmp2<- matrix(0,nrow=nrow(ssDF[[s]]), ncol=nrow(sbar[[s]]) )
    tmp2[as.numeric(rownames(tmp)),]<- tmp
    nfreq[[s]]<- t(tmp2)
    dlst[[s]]<-t( e2dist(sbar[[s]][,2:3], ssDF[[s]][,c("X","Y")])  )
    #lambda<- exp(-(1/(2*sigma*sigma))*d*d)
  }
  list("dlst"=dlst, "nfreq" = nfreq,  "sbar"=sbar)
}


# Distance between traps and activity centers
# Faster than using loops
e2dist <- function (x, y){
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(y)) y <- as.matrix(y)
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


#image.scale <-
#function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks",
#    "ranges"))
#{
#    # sort out the location
#    n <- length(col)
#    usr <- par("usr")
#    mx <- mean(usr[1:2]); my <- mean(usr[3:4])
#    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
#    if (missing(x))
#        x <- mx + 1.05*dx/2	# default x to right of image
#    else if (is.list(x)) {
#        if (length(x$x) == 2)
#          size <- c(diff(x$x), -diff(x$y)/n)
#        y <- x$y[1]
#        x <- x$x[1]
#    } else x <- x[1]
#    if (is.null(size))
#        if (is.null(y)) {
#          size <- 0.618*dy/n	# default size, golden ratio
#          y <- my + 0.618*dy/2	# default y to give centred scale
#        } else size <- (y-my)*2/n
#    if (length(size)==1)
#        size <- rep(size, 2)	# default square boxes
#    if (is.null(y))
#        y <- my + n*size[2]/2
#    # draw the image scale
#    i <- seq(along = col)
#    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2],
#        col = rev(col), xpd = TRUE)
#    # sort out the labels
#    rng <- range(z, na.rm = TRUE)
#    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
#    bks <- formatC(bks, format="f", digits=digits)
#    labels <- match.arg(labels)
#    if (labels == "breaks")
#        ypts <- y - c(0, i) * size[2]
#    else {
#        bks <- paste(bks[-1], bks[-(n+1)], sep = " - ")
#        ypts <- y - (i - 0.5) * size[2]
#    }
#    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj =
#        ifelse(size[1]>0, 0, 1), xpd = TRUE)
#}
#

################################################################################
## Make a previous capture object for use with a behavioral response model
do.prevcap <- function(scrFrame){
  prevcap <- list()
  for (s in 1:length(scrFrame$caphist)) {
    Ys <- scrFrame$caphist[[s]]
    prevcap[[s]] <- array(0, dim = c(dim(Ys)[1], dim(Ys)[2], dim(Ys)[3]))
    first <- matrix(0, dim(Ys)[1], dim(Ys)[2])
    for (i in 1:dim(Ys)[1]) {
      for (j in 1:dim(Ys)[2]) {
        if (sum(Ys[i, j, ]) > 0) {
          first[i, j] <- min((1:(dim(Ys)[3]))[Ys[i, j, ] > 0])
          prevcap[[s]][i, j, 1:first[i, j]] <- 0
          if (first[i, j] < dim(Ys)[3])
            prevcap[[s]][i, j, (first[i, j] + 1):(dim(Ys)[3])] <- 1
        }
      }  
    }
    zeros <- array(0, c(1, dim(prevcap[[s]])[2], dim(prevcap[[s]])[3]))
    prevcap[[s]] <- abind(prevcap[[s]], zeros, along = 1)
  }
  return(prevcap)
}

################################################################################
## Trim to do local evaluations - identify traps and s's within trim of capture

do.trim <- function(scrFrame, ssDF, trimS){
  trimR <- trimC <- list()
  Nocc <- scrFrame$occasions
  
  for (s in 1:length(scrFrame$caphist)) {#each session
    Ys <- scrFrame$caphist[[s]]
    Xs <- scrFrame$traps[[s]]
    Os <- scrFrame$trapOperation[[s]]
    Ss <- ssDF[[s]]
    zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
    Ys <- abind::abind(Ys, zeros, along = 1)
    trimR[[s]] <- list()
    trimC[[s]] <- list()
    
    for (i in 1:nrow(Ys)) { #each individual
      trimR[[s]][[i]] <- list()
      trimC[[s]][[i]] <- list()
      
      #if no trim: global evaluation
      if (is.null(trimS)) {
        pp <- rep(T, ncol(Ys))
        trimC[[s]][[i]] <- rep(T, nrow(ssDF[[s]]))
        
        for (k in 1:Nocc[s]) {#each occasion
          trimR[[s]][[i]][[k]] <- pp
        }
      }
      #if trim: local evaluation
      else {
        if (i < nrow(Ys)) {
          if (dim(Ys)[3]>1) {
            pp <- apply(Ys[i, , ], 1, sum) > 0
          }
          else {
            pp <- Ys[i, , 1] > 0
          }
          for (k in 1:Nocc[s]) {
            tmp.a <- rep(trimS*3 + 2, nrow(Xs))
            tmp.b <- e2dist(Xs[pp, c("X", "Y"), drop = FALSE], 
                            Xs[ , c("X", "Y")])
            if (!is.null(Os)) {
              tmp.ls <- (apply(rbind(tmp.a, tmp.b), 2, min) <= (2 * trimS)) & 
                        (Os[,k] == 1)
            }
            else {
              tmp.ls <- (apply(rbind(tmp.a, tmp.b), 2, min) <= (2 * trimS))
            }
            trimR[[s]][[i]][[k]] <- tmp.ls 
          }
          tmp.a <- rep(trimS + 2, nrow(Ss))
          tmp.b <- e2dist(Xs[pp, c("X", "Y"), drop = FALSE], 
                          Ss[, c("X","Y")]) 
          tmp.ls <- apply(rbind(tmp.a, tmp.b),2, min, na.rm = T) <= trimS
          trimC[[s]][[i]] <- tmp.ls
        }#end i
        else {
          pp <- rep(T, ncol(Ys))
          for (k in 1:Nocc[s]) {
            if (!is.null(Os)) {
              trimR[[s]][[i]][[k]] <- pp & (Os[,k] == 1) 
            }
            else {
              trimR[[s]][[i]][[k]] <- pp
            }
          }
          tmp.a <- rep(trimS + 2, nrow(Ss))
          tmp.b <- e2dist(Xs[pp, c("X", "Y"), drop = FALSE],
                          Ss[, c("X","Y")])
          tmp.ls <- apply(rbind(tmp.a, tmp.b), 2, min, na.rm = T) <= trimS
          trimC[[s]][[i]] <- tmp.ls
        }
      }
    }
  }
  trim.obj <- list(trimR = trimR, trimC = trimC)
  return(trim.obj)
}
################################################################################


################################################################################
## PRINTING
print.oSCR.fitList <- function(x){
  df.out <- data.frame(model = names(x),
                       logL = round(unlist(lapply(x, function(y) y$rawOutput$minimum)),2),
                       K = unlist(lapply(x, function(y) length(y$rawOutput$estimate))),
                       AIC = round(unlist(lapply(x, function(y) y$AIC)),2))
  df.out$dAIC <- round(df.out$AIC - min(df.out$AIC),2)
  rownames(df.out) <- NULL
  print(df.out)
}


print.oSCR.modSel <- function(x){
  if(class(x) == "oSCR.modSel"){
    cat(" AIC model table:",fill=TRUE)
    cat("",fill=TRUE)
    #    tmp.df1 <- data.frame(model = x[[1]][,1],round(x[[1]][,-1],2))
    #    colnames(tmp.df1) <- colnames(x$aic.tab)
    print(x[[1]], digits=2)
    #    cat("",fill=TRUE)
    #    cat("",fill=TRUE)
    #    cat(" Table of coefficients:",fill=TRUE)
    #    cat("",fill=TRUE)
    #    tmp.df2 <- data.frame(model = x[[2]][,1],round(x[[2]][,-1],2))
    #    colnames(tmp.df2) <- colnames(x$coef.tab)
    #    print(tmp.df2)
  }else{
    print(x)
  }
}




################################################################################
# functions required for sub routines in get.real()

get.err <- function(x, p, vcv, nms, j, err=1, d.factor=1,ftype){
  require("car")
  if(err==1) out <- deltaMethod(object=p, vcov.=vcv, g=make.expr(x,nms,d.factor,ftype))     
  if(err==2) out <- deltaMethod(object=p, vcov.=vcv, g=make.expr.M(x,nms,j,d.factor,ftype))
  if(err==3) out <- deltaMethod(object=p, vcov.=vcv, g=make.expr.F(x,nms,j,d.factor,ftype))
  if(err==4) out <- deltaMethod(object=p, vcov.=vcv, g=make.expr.p(x,nms,ftype))
  if(err==5) out <- deltaMethod(object=p, vcov.=vcv, g=make.expr.sig(x,nms,ftype))
  rownames(out) <- NULL
  return(out)
}

#####################################################
## making the expressions to pass to deltaMethod

  #####################################
  #No sex:
  make.expr <- function(x,nms,d.factor,ftype){
    if(length(nms)==1){
      if(ftype=="r")  expr <- paste0("exp(d0) * ",d.factor)
      if(ftype=="lp") expr <- paste0("d0 + log(",d.factor,")")
    }else{  
      if(ftype=="r")  expr <- c("exp(d0", paste0(" + ",nms[-1],"*",x[-1]),")*",d.factor)
      if(ftype=="lp") expr <- c(paste0("d0*",x[1]), paste0(" + ",nms[-1],"*",x[-1]),"+log(",d.factor,")")
    }
    expr <- paste(expr,collapse="")
    return(expr)
  }

  #####################################
  #Sex (mixture):
  
  #male
  make.expr.M <- function(x,nms,j,d.factor,ftype){
    if(length(nms)==1){
      if(ftype=="r")  expr <- paste0("(exp(psi.",j,")/(1+exp(psi.",j,"))) * exp(d0) * ",d.factor)
      if(ftype=="lp") expr <- paste0("log((exp(psi.",j,")/(1+exp(psi.",j,"))) * ",d.factor,") + d0")
    }else{
      if(ftype=="r"){
        expr <- paste(c(paste0("(exp(psi.",j,")/(1+exp(psi.",j,"))) * exp(d0*",x[1]), 
                        paste0(" + ",nms[-1],"*",x[-1]),
                        paste0(") * ", d.factor)), collapse = "")
      }
      if(ftype=="lp"){
        expr <- paste(c(paste0("log((exp(psi.",j,")/(1+exp(psi.",j,"))) * ",d.factor,") + d0*",x[1]), 
                        paste0(" + ",nms[-1],"*",x[-1])), collapse = "")
      }
    }
    return(expr)
  }
  #female
  make.expr.F <- function(x,nms,j,d.factor,ftype){
    if(length(nms)==1){
      if(ftype=="r")  expr <- paste0("(1-exp(psi.",j,")/(1+exp(psi.",j,"))) * exp(d0) * ",d.factor)
      if(ftype=="lp") expr <- paste0("log((1-exp(psi.",j,")/(1+exp(psi.",j,"))) * ",d.factor,") + d0")
    }else{
      if(ftype=="r"){
        expr <- paste(c(paste0("(1-exp(psi.",j,")/(1+exp(psi.",j,"))) * exp(d0*",x[1]), 
                        paste0(" + ",nms[-1],"*",x[-1]),
                        paste0(") * ",d.factor)), collapse = "")
      }
      if(ftype=="lp"){
        expr <- paste(c(paste0("log((1-exp(psi.",j,")/(1+exp(psi.",j,"))) * ",d.factor,") + d0*",x[1]), 
                        paste0(" + ",nms[-1],"*",x[-1])), collapse = "")
      } 
    }
    return(expr)
  }

  #####################################
  #Detection:
  make.expr.p <- function(x,nms,ftype){
    if(length(nms)==1){
      if(ftype=="r")  expr <- "exp(p0)/(1+exp(p0))"
      if(ftype=="lp") expr <- "p0"
    }else{
      if(ftype=="r"){
        expr <- paste(c(paste0("exp(p0*",x[1]),
                        paste0(" + ",nms[-1],"*",x[-1]),
                        paste0(") / (1+exp(p0*",x[1]),
                        paste0(" + ",nms[-1],"*",x[-1]),"))"), collapse = "")
      }
      if(ftype=="lp"){
        expr <- paste(c(paste0("p0*",x[1]),
                        paste0(" + ",nms[-1],"*",x[-1])), collapse = "")
      }
    }
    return(expr)
  }
  
  #####################################
  #Sigma
  make.expr.sig <- function(x,nms,ftype){
    if(length(nms)==1){
      expr <- "exp(sig0)"
      if(ftype=="r")  expr <- "exp(sig0)"
      if(ftype=="lp") expr <- "sig0"
    }else{
      expr <- c("exp(sig0", paste0("+",nms[-1],"*",x[-1]),") / (1+exp(p0", paste0("+",nms[-1],"*",x[-1]),"))")
      if(ftype=="r"){
        expr <- paste(c(paste0("exp(sig0*",x[1]),
                        paste0(" + ",nms[-1],"*",x[-1]),")"), collapse = "")
      }
      if(ftype=="lp"){
        expr <- paste(c(paste0("sig0*",x[1]),
                        paste0(" + ",nms[-1],"*",x[-1])), collapse = "")
      }
    }
    return(expr)
  }
###############################################################
  
  
  
#####################################################
## Handling offsets

is.offset <- function(f){
  off <- list()
  off[["offset"]] <- "offset" %in% names(attributes(terms.formula(f)))
  off[["pos"]] <- attributes(terms.formula(f))$offset 
  off[["name"]] <- ifelse(is.null(off$pos),NA,all.vars(f)[off$pos])
  return(off)
}
  