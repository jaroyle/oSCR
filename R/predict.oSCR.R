predict.oSCR <- function (scr.fit, scrFrame = NULL, ssDF = NULL, costDF = NULL, 
                          rsfDF = NULL, override.trim=FALSE) {
  mles <- scr.fit$rawOutput$estimate
  call <- scr.fit$call
  oSCR.fit2 <- oSCR.fit
  call.fix <- names(call)[!names(call) %in% c("", "scrFrame", "ssDF", "costDF", "rsfDF")]
  if(override.trim)
    call.fix <- call.fix[call.fix!="trimS"]
  if (length(call.fix>0)) {
    for (i in 1:length(call.fix)) {
      formals(oSCR.fit2)[[call.fix[i]]] <- call[[call.fix[i]]]
    }
  }
  if(is.null(scrFrame)) {
    sf <- scr.fit$scrFrame
  }
  else {
    sf <- scrFrame
  }
  if(is.null(ssDF)) {
    ss <- scr.fit$ssDF
  }
  else {
    ss <- ssDF
  }
  if(is.null(costDF)) {
    cs <- scr.fit$costDF
  }
  else {
    cs <- costDF
  }
  if(is.null(rsfDF)) {
    rs <- scr.fit$rsfDF
  }
  else {
    rs <- ssDF
  }
  
  out <- oSCR.fit2(scrFrame = sf, ssDF = ss, costDF = cs, 
                   rsfDF = rs, start.vals = mles, predict = TRUE)
  nsess <- length(out$preds)
  r <- list()
  total <- list()
  pbar<- list()
  for (s in 1:nsess) {
    nguys <- dim(out$preds[[s]])[1]
    Nhat <- sum(out$ss.bits[[s]][, "d.s"])
    n0 <- out$ss.bits[[s]][, "d.s"] * out$ss.bits[[s]][,"lik.cond"]
    tmp <- sp::SpatialPoints(out$ssDF[[s]][,c("X","Y")])
    tmp <- try(sp::points2grid(tmp))
    if( class(tmp) == "try-error") {
      cat("Cannot rasterize state-space",fill=TRUE)
      pbar[[s]]<- cbind(out$ssDF[[s]][, c("X","Y")], 
                        pbar = (1 - out$ss.bits[[s]][, "lik.cond"]))
    }
    else {
      pbar[[s]] <- raster::rasterFromXYZ(
                     cbind(out$ssDF[[s]][, c("X", "Y")], 
                           pbar = (1 - out$ss.bits[[s]][, "lik.cond"])))
    }
    cat("Nhat: ", sum(n0) + nguys - 1, fill = TRUE)
    out$preds[[s]][nguys, ] <- n0
    total[[s]] <- apply(out$preds[[s]], 2, sum)
    cat("sum of predictions: ", sum(total[[s]]), fill = TRUE)
    if (class(tmp) == "try-error") {
      r[[s]] <- cbind(out$ssDF[[s]][, c("X", "Y")], total[[s]])
    }
    else {
      r[[s]] <- rasterFromXYZ(cbind(out$ssDF[[s]][, c("X", "Y")], total[[s]]))
    }
  }
  return(list(r = r, ssN = total, preds = out$preds, ssDF = out$ssDF, pbar=pbar))
}
