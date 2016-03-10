predict.oSCR <-
function(scrFrame, ssDF, scr.fit){
 
  mles <- scr.fit$rawOutput$estimates
  mod <- scr.fit$call
  mc<- mod$multicatch
  tr<- mod$trimS
  prediction <- oSCR.fit(scrFrame=scrFrame, ssDF=ssDF, model=mod, start.vals=mles, trimS=tr, 
  multicatch=scr.fit$call$mc,
  predict=TRUE)
  return(prediction)
 
}
