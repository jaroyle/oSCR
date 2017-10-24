oSCR.parfit <- function(mods, ncores=2,
                        scrFrame, ssDF, encmod = c("B", "P", "CLOG")[1], multicatch = FALSE, theta = 2, 
                        trimS = NULL, DorN = c("D", "N")[1], sexmod = c("constant", "session")[1], 
                        costDF = NULL, distmet = c("euc", "user", "ecol")[1], directions = 8, 
                        PROJ = NULL, rsfDF = NULL, RSF = FALSE, telemetry = c("none","ind","dep")[1],
                        se = TRUE, predict = FALSE, start.vals = NULL, getStarts = FALSE, pxArea = 1, 
                        plotit = F, mycex = 1, nlmgradtol = 1e-06, nlmstepmax = 10, smallslow = FALSE, 
                        print.level = 0){
  library(doParallel)
  wd <- getwd()
  par.wrap <- function(model.number, mods,
                       scrFrame, ssDF, encmod, multicatch, theta, trimS, DorN, sexmod, 
                       costDF, distmet, directions, PROJ, rsfDF, RSF, telemetry, se, predict, start.vals, 
                       getStarts, pxArea, plotit, mycex, nlmgradtol, nlmstepmax, smallslow, print.level){
    mod <- list(as.formula(paste(mods[[1]][[model.number]])), #density
                as.formula(paste(mods[[2]][[model.number]])), #detection
                as.formula(paste(mods[[3]][[model.number]])), #sigma
                as.formula(paste(mods[[4]][[model.number]]))) #asu
    
    fm <- oSCR.fit(model = mod, scrFrame, ssDF, encmod, multicatch, theta, trimS, DorN, sexmod, 
                   costDF, distmet, directions, PROJ, rsfDF, RSF, telemetry, se, predict, start.vals, 
                   getStarts, pxArea, plotit, mycex, nlmgradtol, nlmstepmax, smallslow, print.level)
    return(fm)
  }
  nmods <- nrow(mods)
  cl <- makeCluster(ncores)  
  registerDoParallel(cl)
  out <- foreach(i = 1:nmods) %dopar% {
    library(oSCR)
    tmp <- par.wrap(i, mods, scrFrame, ssDF, encmod, multicatch, theta, 
                    trimS, DorN, sexmod, costDF, distmet, directions, PROJ, 
                    rsfDF, RSF, telemetry, se, predict, start.vals, getStarts, 
                    pxArea, plotit, mycex, nlmgradtol, nlmstepmax, smallslow, 
                    print.level)
    save("tmp",file=paste0(paste0(wd,"/model_",i,".RData")))
    return(tmp)
  }
  stopCluster(cl)
  save("out",file="models.RData")
  out<- fitList.oSCR(out,rename=TRUE)
  return(out)
}
