flatten <-
function(y3d){
  if(any(y3d[!is.na(y3d) & y3d>1])) return("error: some encounters > 1")
  y2d<- apply(y3d,c(1,3),sum)
  y2d[y2d>1 & !is.na(y2d)]<- 1
  y2d
}
