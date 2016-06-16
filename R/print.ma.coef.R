print.ma.coef <- function(ma){
  class(ma) <- "data.frame"
  nums <- vapply(ma, is.numeric, FUN.VALUE = logical(1))
  ma[,nums] <- round(ma[,nums], digits = 2)
  print(ma)
  cat("", fill=TRUE)
  cat("* Denotes 'shrinkage' values.", fill=TRUE)
}
