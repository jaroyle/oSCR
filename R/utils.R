



# Distance between traps and activity centers
# Faster than using loops
e2dist <- function (x, y)
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}






image.scale <-
function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks", 
    "ranges"))
{
    # sort out the location
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2]); my <- mean(usr[3:4])
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    if (missing(x))
        x <- mx + 1.05*dx/2	# default x to right of image
    else if (is.list(x)) {
        if (length(x$x) == 2) 
          size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    } else x <- x[1]
    if (is.null(size))
        if (is.null(y)) {
          size <- 0.618*dy/n	# default size, golden ratio
          y <- my + 0.618*dy/2	# default y to give centred scale
        } else size <- (y-my)*2/n
    if (length(size)==1)
        size <- rep(size, 2)	# default square boxes
    if (is.null(y))
        y <- my + n*size[2]/2
    # draw the image scale
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2], 
        col = rev(col), xpd = TRUE)
    # sort out the labels
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format="f", digits=digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
        ypts <- y - c(0, i) * size[2]
    else {
        bks <- paste(bks[-1], bks[-(n+1)], sep = " - ")
        ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj =
        ifelse(size[1]>0, 0, 1), xpd = TRUE) 
}





#fitList - a function for creating a list of oSCR model output objects
#          formattted as an 'oSCR.fitList' for convenient post processeing


# x: a simple list of oSCR.fit objects
# rename: remame objects based on the fitted model

fitList <- function(x,rename=F){

  if(rename==FALSE & is.null(names(x)))
   fl.names <- 1:length(x)

  if(rename==FALSE & !is.null(names(x)))
   fl.names <- names(x)

  if(rename==TRUE){
    n1 <- paste(lapply(x,function(z) x$call$model[2]))
    n1 <- ifelse(n1 %in% "NULL",".",n1)
    n2 <- paste(lapply(x,function(z) x$call$model[3]))
    n2 <- ifelse(n2 %in% "NULL",".",n2)
    n3 <- paste(lapply(x,function(z) x$call$model[4]))
    n3 <- ifelse(n3 %in% "NULL",".",n3)
    n4 <- paste(lapply(x,function(z) x$call$model[5]))
    n4 <- ifelse(n4 %in% "NULL",".",n4)

    fl.names <- paste("D(",n1,") ","p(",n2,") ",
                      "sig(",n3,") ","asu(",n4,")", sep="")
  }
  names(x) <- fl.names
  class(x) <- "oSCR.fitList"
  return(x)
}





#print.oSCR.fit - prints a summary table of an oSCR.fitList object
# x: an oSCR.fitList

print.oSCR.fitList <- function(x){
  if(class(x) == "oSCR.fitList"){

    df.out <- data.frame(model = names(x),
                         logL = unlist(lapply(x, function(y) y$rawOutput$minimum)),
                         K = unlist(lapply(x, function(y) length(y$rawOutput$estimate))),
                         AIC = unlist(lapply(x, function(y) y$AIC)))

    df.out$dAIC <- df.out$AIC - min(df.out$AIC)
  }else{print(x)}

  print(df.out,digits=2)
}





#modSel - a function for generating an ordered (by deltaAIC) model selection
#         table from an 'oSCR.fitList'

# x: an 'oSCR.fitList object

modSel <- function(x){

  if(class(x)!="oSCR.fitList")
    stop("Object must be an oSCR.fitList")

  df.out <- data.frame(model = names(x),
                       logL = unlist(lapply(x, function(y) y$rawOutput$minimum)),
                       K = unlist(lapply(x, function(y) length(y$rawOutput$estimate))),
                       AIC = unlist(lapply(x, function(y) y$AIC)))

  df.out$dAIC <- df.out$AIC - min(df.out$AIC)
  (df.out <- df.out[order(df.out$dAIC),])
  return(df.out)
}



