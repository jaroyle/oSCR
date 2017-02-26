predict.oSCR <-
function (scrFrame, scr.fit, ssDF, costDF = NULL, rsfDF = NULL, override.trim=FALSE)
{
    mles <- scr.fit$rawOutput$estimate
    call <- scr.fit$call
    oSCR.fit2 <- oSCR.fit
    call.fix <- names(call)[!names(call) %in% c("", "scrFrame",
        "ssDF", "costDF", "rsfDF")]

    if(override.trim)
       call.fix<- call.fix[call.fix!="trimS"]

    for (i in 1:length(call.fix)) {
        formals(oSCR.fit2)[[call.fix[i]]] <- call[[call.fix[i]]]
    }
    out <- oSCR.fit2(scrFrame = scrFrame, ssDF = ssDF, costDF = costDF, rsfDF = rsfDF,
        start.vals = mles, predict = TRUE)
    nsess <- length(out$preds)
    r <- list()
    total <- list()
    pbar<- list()
    for (s in 1:nsess) {
        nguys <- dim(out$preds[[s]])[1]
        Nhat <- sum(out$ss.bits[[s]][, "d.s"])
        library(raster)
        n0 <- out$ss.bits[[s]][, "d.s"] * out$ss.bits[[s]][,
        "lik.cond"]



tmp<- SpatialPoints(out$ssDF[[s]][,c("X","Y")])
tmp<- try(points2grid(tmp))
if( class(tmp) == "try-error"){
    cat("Cannot rasterize state-space",fill=TRUE)
    pbar[[s]]<- cbind(out$ssDF[[s]][, c("X",
            "Y")], pbar = (1 - out$ss.bits[[s]][, "lik.cond"]))
}else{
        pbar[[s]] <- rasterFromXYZ(cbind(out$ssDF[[s]][, c("X",
            "Y")], pbar = (1 - out$ss.bits[[s]][, "lik.cond"])))
     }


        cat("Nhat: ", sum(n0) + nguys - 1, fill = TRUE)
        out$preds[[s]][nguys, ] <- n0
        total[[s]] <- apply(out$preds[[s]], 2, sum)
        cat("sum of predictions: ", sum(total[[s]]), fill = TRUE)
        if( class(tmp) == "try-error"){
r[[s]]<-   cbind(out$ssDF[[s]][, c("X",             "Y")], total[[s]])
 }else{
        r[[s]] <- rasterFromXYZ(cbind(out$ssDF[[s]][, c("X",
                                                        "Y")], total[[s]]))
        }
        r[[s]] <- rasterFromXYZ(cbind(out$ssDF[[s]][, c("X",
            "Y")], total[[s]]))
    }
    return(list(r = r, ssN = total, preds = out$preds, ssDF = out$ssDF, pbar=pbar))
}
