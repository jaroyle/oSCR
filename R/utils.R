



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


################################################################################
## PRINTING

#print.oSCR.fit - prints a summary table of an oSCR.fitList object
# x: an oSCR.fitList

## This should be change to print interesting summaries from the model list
## currently does just the same as modSel

print.oSCR.fitList <- function(x){

    df.out <- data.frame(model = names(x),
                         logL = round(unlist(lapply(x, function(y) y$rawOutput$minimum)),2),
                         K = unlist(lapply(x, function(y) length(y$rawOutput$estimate))),
                         AIC = round(unlist(lapply(x, function(y) y$AIC)),2))

    df.out$dAIC <- round(df.out$AIC - min(df.out$AIC),2)
    rownames(df.out) <- NULL
    print(df.out)
}






print.oSCR.modSel <- function(x){
  if(class(x) == "oSCR.modSel"){
cat(" AIC model table:",fill=TRUE)
cat("",fill=TRUE)
tmp.df1 <- data.frame(model = x[[1]][,1],round(x[[1]][,-1],2))
colnames(tmp.df1) <- colnames(x$aic.tab)
print(tmp.df1)
cat("",fill=TRUE)
cat("",fill=TRUE)
cat(" Table of coefficients:",fill=TRUE)
cat("",fill=TRUE)
tmp.df2 <- data.frame(model = x[[2]][,1],round(x[[2]][,-1],2))
colnames(tmp.df2) <- colnames(x$coef.tab)
print(tmp.df2)
  }else{
print(x)
}
}


