connectivity.surface <- function(scr.fit, scrFrame, ssDF, costDF=NULL, PROJ = NULL,
                                 directions=8){
  
  if(is.null(PROJ))
    PROJ <- "+proj=utm +zone=12 +datum=WGS84"
  dwc = list()
  Ec = list()
  dm.cost = list()
  D = list()
  prY = list()
  pred = predict.oSCR(scrFrame,scr.fit,ssDF,costDF)
  asu.mod = formula(scr.fit$model[[4]])
  for (s in 1:length(costDF)) {
    dm.cost[[s]] = model.matrix(asu.mod, as.data.frame(costDF[[s]]))
  }
  c.nms = colnames(dm.cost[[1]])
  names.dist = paste("c.beta", c.nms, sep = ".")
  chx = grep(fixed=TRUE,"Intercept", names.dist)
  if (length(chx) > 0)
    names.dist[chx] <- "c0.(Intercept)"
  tmp.par = data.frame(parameters = names.dist)
  hat.a2 = merge(tmp.par, scr.fit$outStats[,c("parameters","mle")],by.x=T)[,2]

  
  sig.mod = formula(scr.fit$model[[3]])
  if("sex" %in% all.vars(sig.mod)){
    tmp.sig = c(scr.fit$outStats$mle[grep(fixed=T,"sig.(Intercept)", scr.fit$outStats$parameters)],
                scr.fit$outStats$mle[grep(fixed=T,"sig.(Intercept)",scr.fit$outStats$parameters)]+
                scr.fit$outStats$mle[grep(fixed=T,"sig.male",scr.fit$outStats$parameters)])  
    a1 = 1/(2*exp(tmp.sig)^2)
    psi = scr.fit$outStats$mle[grep("psi",scr.fit$outStats$parameters)]         
  } 
  #need to add session varying sigma
  #need to add session varying sex
  
  for(i in 1:length(ssDF)){
    
    cost = exp(dm.cost[[i]] %*% c(hat.a2))
    costR = rasterFromXYZ(cbind(costDF[[i]][, c(1,2)], cost))
    projection(costR) <- PROJ
    tr = transition(costR, 
                     transitionFunction = function(x) (1/(mean(x))),
                     direction = directions)
    trLayer = geoCorrection(tr, scl = F)
    D[[i]] = costDistance(trLayer, 
                          as.matrix(ssDF[[i]][,c("X", "Y")]), 
                          as.matrix(ssDF[[i]][, c("X","Y")]))

    # 1. density weighted connectivity
    dwc[[i]] = exp(-a1[1] * D[[i]]^2) %*% (pred$ssN[[i]] * (1 - plogis(psi))) +
               exp(-a1[2] * D[[i]]^2) %*% (pred$ssN[[i]] * plogis(psi))
      
    # 2. expected connectivity for this need to change PREDICT
    Ec[[1]] = NULL
  }
  return(list(dwc,Ec))
}
