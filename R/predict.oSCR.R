predict.oSCR <-    function (scrFrame, scr.fit, ssDF, costDF = NULL)
{
    mles <- scr.fit$rawOutput$estimate
    call <- scr.fit$call
    oSCR.fit2 <- oSCR.fit
    call.fix <- names(call)[!names(call) %in% c("", "scrFrame",
        "ssDF", "costDF")]
    for (i in 1:length(call.fix)) {
        formals(oSCR.fit2)[[call.fix[i]]] <- call[[call.fix[i]]]
    }
    out <- oSCR.fit2(scrFrame = scrFrame, ssDF = ssDF,
        costDF = costDF, start.vals = mles, predict = TRUE)
    nsess<-length(out$preds)

    r<- list()
    total<-list()
    
    for(s in 1:nsess){
    nguys<- dim(out$preds[[s]])[1]
    Nhat<- sum(out$ss.bits[[s]][,"d.s"])
    library(raster)
# number of individuals not detected at each pixel
# this is wrong, should involve lik.condl ... ie Pr(y=0 | s)
n0<- out$ss.bits[[s]][,"d.s"]*out$ss.bits[[s]][,"lik.cond"]
cat("Nhat: ", sum(n0) + nguys  - 1, fill=TRUE)
# probability that a guy at pixel s is NOT captured
##plot(rasterFromXYZ(cbind(out1$ssDF[[1]],out1$lik.cond)))
out$preds[[s]][nguys,]<- n0
# Does not include the "n0" weighting yet
total[[s]]<- apply(out$preds[[s]],2,sum)
# check
cat("sum of predictions: ", sum(total[[s]]) , fill=TRUE)
r[[s]]<- rasterFromXYZ(cbind(out$ssDF[[s]][,c("X","Y")],total[[s]]))

}


return(list(r=r, ssN=total, preds=out$preds, ssDF=out$ssDF))
}
