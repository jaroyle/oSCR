make.mods <-
function(density,detection,sigma,cost=c(~1) ){
 
  # aa <- paste0("~",density)
  # bb <- paste0("~",detection)
  # cc <- paste0("~",sigma)
  # dd <- paste0("~",cost)
  aa<- density
  bb<- detection
  cc<- sigma
  dd<- cost
   mm <- expand.grid(dd,cc,bb,aa)[,c(4,3,2,1)]
  return(mm)
 
}
