\name{oSCR}
\alias{oSCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oSCR(scrFrame, model = list(D ~ 1, p0 ~ 1, a1 ~ 1, path ~ 1), ssDF = NULL, costDF = NULL, distmet = c("euc", "user", "ecol")[1], sexmod = c("constant", "session")[1], DorN = c("D", "N")[1], directions = 8, Dmat = NULL, trimS = NULL, res = NULL, smallslow = FALSE, predict = FALSE, start.vals = NULL, PROJ = NULL, plotit = F, mycex = 0.5, tester = F, pl = 0, nlmgradtol = 1e-06, nlmstepmax = 10, multicatch = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scrFrame}{
%%     ~~Describe \code{scrFrame} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{ssDF}{
%%     ~~Describe \code{ssDF} here~~
}
  \item{costDF}{
%%     ~~Describe \code{costDF} here~~
}
  \item{distmet}{
%%     ~~Describe \code{distmet} here~~
}
  \item{sexmod}{
%%     ~~Describe \code{sexmod} here~~
}
  \item{DorN}{
%%     ~~Describe \code{DorN} here~~
}
  \item{directions}{
%%     ~~Describe \code{directions} here~~
}
  \item{Dmat}{
%%     ~~Describe \code{Dmat} here~~
}
  \item{trimS}{
%%     ~~Describe \code{trimS} here~~
}
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{smallslow}{
%%     ~~Describe \code{smallslow} here~~
}
  \item{predict}{
%%     ~~Describe \code{predict} here~~
}
  \item{start.vals}{
%%     ~~Describe \code{start.vals} here~~
}
  \item{PROJ}{
%%     ~~Describe \code{PROJ} here~~
}
  \item{plotit}{
%%     ~~Describe \code{plotit} here~~
}
  \item{mycex}{
%%     ~~Describe \code{mycex} here~~
}
  \item{tester}{
%%     ~~Describe \code{tester} here~~
}
  \item{pl}{
%%     ~~Describe \code{pl} here~~
}
  \item{nlmgradtol}{
%%     ~~Describe \code{nlmgradtol} here~~
}
  \item{nlmstepmax}{
%%     ~~Describe \code{nlmstepmax} here~~
}
  \item{multicatch}{
%%     ~~Describe \code{multicatch} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (scrFrame, model = list(D ~ 1, p0 ~ 1, a1 ~ 1, path ~ 
    1), ssDF = NULL, costDF = NULL, distmet = c("euc", "user", 
    "ecol")[1], sexmod = c("constant", "session")[1], DorN = c("D", 
    "N")[1], directions = 8, Dmat = NULL, trimS = NULL, res = NULL, 
    smallslow = FALSE, predict = FALSE, start.vals = NULL, PROJ = NULL, 
    plotit = F, mycex = 0.5, tester = F, pl = 0, nlmgradtol = 1e-06, 
    nlmstepmax = 10, multicatch = FALSE) 
{
    if (multicatch) {
        for (s in 1:length(scrFrame$caphist)) {
            captures <- apply(scrFrame$caphist[[s]], c(1, 3), 
                sum)
            if (any(captures > 1)) 
                return("error: multicatch system cannot have > 1 capture", 
                  fill = TRUE)
        }
    }
    my.model.matrix <- function(form, data) {
        mdm <- suppressWarnings(model.matrix(form, data, contrasts.arg = lapply(data.frame(data[, 
            sapply(data.frame(data), is.factor)]), contrasts, 
            contrasts = FALSE)))
        return(mdm)
    }
    if (is.null(trimS)) {
        max.dist <- 0
        for (i in 1:length(scrFrame$caphist)) {
            for (j in 1:nrow(scrFrame$caphist[[i]])) {
                where <- apply(scrFrame$caphist[[i]][j, , ], 
                  1, sum) > 0
                max.dist <- max(max.dist, max(0, dist(scrFrame$traps[[i]][where, 
                  ])))
            }
        }
    }
    ptm <- proc.time()
    starttime <- format(Sys.time(), "\%H:\%M:\%S \%d \%b \%Y")
    cl <- match.call(expand.dots = TRUE)
    if (!require(abind)) 
        stop("need to install package 'abind'")
    if (!require(Formula)) 
        stop("need to load package 'Formula'")
    if (distmet == "ecol") {
        if (!require(raster)) 
            stop("need to install package 'raster'")
        if (!require(gdistance)) 
            stop("need to install package 'gdistance'")
    }
    if (!inherits(scrFrame, "scrFrame")) {
        stop("Data must be of class 'scrFrame'")
    }
    if (max(unlist(lapply(scrFrame$caphist, max))) > 1) {
        stop("Data in caphist must be Binary")
    }
    if (distmet == "ecol") {
        if (is.null(PROJ)) {
            message("Projection not provided, using default: '+proj=utm +zone=12 +datum=WGS84'")
        }
    }
    pars.p0 <- NULL
    names.p0 <- NULL
    pars.a1 <- NULL
    names.a1 <- NULL
    pars.beta.trap <- NULL
    names.beta.trap <- NULL
    pars.beta.den <- NULL
    names.beta.den <- NULL
    pars.d0 <- NULL
    names.d0 <- NULL
    pars.beta.den <- NULL
    names.beta.den <- NULL
    pars.dist <- NULL
    names.dist <- NULL
    pars.dist <- NULL
    names.dist <- NULL
    singleS <- NULL
    singleT <- NULL
    singleG <- NULL
    sess.ss <- NULL
    sess.ss.nG <- NULL
    sess.ssArea <- NULL
    sess.trapCovs <- list()
    D <- list()
    YY <- list()
    dm.den <- list()
    tmp.dm <- list()
    dm.trap <- list()
    dm.cost <- list()
    posterior <- list()
    dHPP <- FALSE
    dIPP <- FALSE
    dSession <- FALSE
    trap.covs <- FALSE
    pDot <- FALSE
    pTime <- FALSE
    pJustsex <- FALSE
    pJustsesh <- FALSE
    pBothsexnsesh <- FALSE
    anySex <- FALSE
    aDot <- FALSE
    aJustsex <- FALSE
    aJustsesh <- FALSE
    aBothsexnsesh <- FALSE
    if (length(model) == 3) {
        model[[4]] <- formula(~1)
    }
    for (i in 1:4) {
        model[[i]] <- update.formula(model[[i]], NULL ~ .)
    }
    K <- ifelse(is.list(scrFrame$nocc), scrFrame$nocc[[1]], scrFrame$nocc)
    nsess.trap <- length(scrFrame$traps)
    allvars.D <- all.vars(model[[1]])
    dens.fx <- allvars.D[!allvars.D \%in\% c("D")]
    allvars.T <- all.vars(model[[2]])
    trap.fx <- allvars.T[!allvars.T \%in\% c("p0", "session", "Session", 
        "sex", "t", "T")]
    allvars.p0a <- all.vars(model[[2]])
    allvars.p0 <- allvars.p0a[!allvars.p0a == "p0"]
    allvars.a1a <- all.vars(model[[3]])
    allvars.a1 <- allvars.a1a[!allvars.a1a == "a1"]
    var.p0.1 <- "sex" \%in\% allvars.p0
    var.p0.2 <- "session" \%in\% allvars.p0
    var.a1.1 <- "sex" \%in\% allvars.a1
    var.a1.2 <- "session" \%in\% allvars.a1
    allvars.dist <- all.vars(model[[4]])
    allvars.dist <- allvars.dist[!allvars.dist == "path"]
    if (is.null(ssDF)) {
        message("Generating a state space based on traps")
        dHPP <- TRUE
        ssDF <- list()
        for (i in 1:length(scrFrame$traps)) {
            dd <- as.matrix(dist(scrFrame$traps[[i]]))
            dd[dd == 0] <- NA
            nnd <- mean(apply(dd, 1, min, na.rm = T))
            ssDF[[i]] <- expand.grid(seq(min(scrFrame$traps[[i]][, 
                1]) - 4 * nnd, max(scrFrame$traps[[i]][, 1]) + 
                4 * nnd, nnd/4), seq(min(scrFrame$traps[[i]][, 
                2]) - 4 * nnd, max(scrFrame$traps[[i]][, 2]) + 
                4 * nnd, nnd/4))
        }
    }
    ns <- length(scrFrame$caphist)
    ng <- length(ssDF)
    nt <- length(scrFrame$traps)
    if (ns == 1) {
        if (ng > 1) 
            message(">1 'state-space' provided, using the first ONLY.")
        if (nt > 1) 
            message(">1 'traps' object provided, using the first ONLY.")
        singleS <- TRUE
        singleT <- TRUE
        singleG <- TRUE
        sess.ss.nG <- nrow(ssDF[[1]])
        sess.ssArea <- ifelse(is.null(res), 1 * sess.ss.nG, res * 
            sess.ss.nG)
        nnn <- names(ssDF[[1]])
    }
    if (ns > 1) {
        singleS <- FALSE
        if (ng == 1 && nt == 1) {
            message("Same state space will be used for all sessions.")
            message("Same traps will be used for all sessions.")
            singleG <- TRUE
            singleT <- TRUE
            sess.ss.nG <- unlist(lapply(ssDF, nrow))
            sess.ssArea <- ifelse(is.null(res), 1 * sess.ss.nG, 
                res * sess.ss.nG)
        }
        if (ng == ns && nt == ns) {
            singleG <- FALSE
            singleT <- FALSE
            nnn <- names(ssDF[[1]])
            for (s in 2:ns) {
                if (!identical(names(ssDF[[s]]), nnn)) 
                  stop("ssDF must have the same covariate names and be in the same order")
            }
            sess.ss.nG <- unlist(lapply(ssDF, nrow))
            sess.ssArea <- ifelse(is.null(res), 1 * sess.ss.nG, 
                res * sess.ss.nG)
        }
    }
    if (nt != ng) 
        stop("Number of trap arrays '!=' number of state-spaces")
    if ("session" \%in\% all.vars(model[[1]]) & !("session" \%in\% 
        nnn)) {
        for (g in 1:ns) {
            ssDF[[g]]$session <- factor(rep(g, nrow(ssDF[[g]])), 
                levels = 1:ns)
        }
    }
    mod2 <- update(model[[2]], ~. - sex - session - Session - 
        t - 1)
    if (length(trap.fx) > 0) {
        trap.covs <- TRUE
        tcovnms <- colnames(scrFrame$trapCovs[[1]][[1]])
        tCovMissing <- trap.fx[!trap.fx \%in\% tcovnms]
        print(tCovMissing)
        if (length(tCovMissing) > 0) {
            stop("I cant find theses covariates in 'scrFrame$trapCovs'", 
                for (i in tCovMissing) print(i))
        }
        if (any(c("session", "Session") \%in\% allvars.T)) 
            tSession <- TRUE
        if (singleT) {
            for (g in 1:ns) {
                for (k in 1:(dim(scrFrame$caphist[[g]])[3])) {
                  if (length(scrFrame$trapCovs[[1]]) > 1) {
                    tmp.dm[[k]] <- my.model.matrix(mod2, scrFrame$trapCovs[[1]][[k]])
                    if (g == 1 && k == 1) 
                      t.nms <- colnames(tmp.dm[[k]])
                    if (nrow(tmp.dm[[k]]) != nrow(scrFrame$trapCovs[[1]][[k]])) {
                      mis <- setdiff(rownames(scrFrame$trapCovs[[1]][[k]]), 
                        rownames(my.model.matrix(mod2, scrFrame$trapCovs[[1]][[k]])))
                      tmp.insert <- matrix(NA, length(mis), ncol(tmp.dm[[k]]))
                      row.names(tmp.insert) <- mis
                      tmp.dm[[k]] <- rbind(tmp.dm[[k]], tmp.insert)
                      tmp.dm[[k]] <- tmp.dm[[k]][order(as.numeric(row.names(tmp.dm[[k]]))), 
                        ]
                    }
                  }
                  else {
                    tmp.dm[[k]] <- my.model.matrix(mod2, scrFrame$trapCovs[[1]][[1]])
                    if (g == 1 && k == 1) 
                      t.nms <- colnames(tmp.dm[[k]])
                    if (nrow(tmp.dm[[k]]) != nrow(scrFrame$trapCovs[[1]][[1]])) {
                      mis <- setdiff(rownames(scrFrame$trapCovs[[1]][[1]]), 
                        rownames(my.model.matrix(mod2, scrFrame$trapCovs[[1]][[1]])))
                      tmp.insert <- matrix(NA, length(mis), ncol(tmp.dm[[k]]))
                      row.names(tmp.insert) <- mis
                      tmp.dm[[k]] <- rbind(tmp.dm[[k]], tmp.insert)
                      tmp.dm[[k]] <- tmp.dm[[k]][order(as.numeric(row.names(tmp.dm[[k]]))), 
                        ]
                    }
                  }
                }
                dm.trap[[g]] <- tmp.dm
            }
        }
        else {
            for (g in 1:ns) {
                for (k in 1:(dim(scrFrame$caphist[[g]])[3])) {
                  tmp.dm[[k]] <- my.model.matrix(mod2, scrFrame$trapCovs[[g]][[k]])
                  if (g == 1 && k == 1) 
                    t.nms <- colnames(tmp.dm[[k]])
                  if (nrow(tmp.dm[[k]]) != nrow(scrFrame$trapCovs[[g]][[k]])) {
                    cat("case b", fill = TRUE)
                    mis <- setdiff(rownames(scrFrame$trapCovs[[g]][[k]]), 
                      rownames(my.model.matrix(mod2, scrFrame$trapCovs[[g]][[k]])))
                    tmp.insert <- matrix(NA, length(mis), ncol(tmp.dm[[k]]))
                    row.names(tmp.insert) <- mis
                    tmp.dm[[k]] <- rbind(tmp.dm[[k]], tmp.insert)
                    tmp.dm[[k]] <- tmp.dm[[k]][order(as.numeric(row.names(tmp.dm[[k]]))), 
                      ]
                    if (!is.matrix(tmp.dm[[k]])) 
                      tmp.dm[[k]] <- matrix(tmp.dm[[k]], ncol = 1)
                  }
                }
                dm.trap[[g]] <- tmp.dm
            }
        }
        if ("Session" \%in\% all.vars(model[[2]])) {
            tmpTsess <- rep(1:ns, each = length(t.nms))
            tmpTcovs <- rep(t.nms, ns)
            names.beta.trap <- paste("t.beta.", tmpTcovs, ".sess", 
                tmpTsess, sep = "")
            pars.beta.trap <- rnorm(length(names.beta.trap))
        }
        else {
            names.beta.trap <- paste("t.beta.", t.nms, sep = "")
            pars.beta.trap <- rnorm(length(names.beta.trap))
        }
    }
    if (DorN == "D") {
        if (length(dens.fx) > 0) {
            dcovnms <- colnames(ssDF[[1]])
            dCovMissing <- dens.fx[which(!dens.fx \%in\% dcovnms)]
            if (length(dCovMissing) > 0) {
                stop("I cant find theses covariates in 'ssDF'", 
                  for (i in dCovMissing) print(i))
            }
        }
        if (("session" \%in\% all.vars(model[[1]]))) 
            mod1 <- update(model[[1]], ~. - sex - 1)
        if (!("session" \%in\% all.vars(model[[1]]))) 
            mod1 <- update(model[[1]], ~. - sex)
        for (g in 1:ns) {
            dm.den[[g]] <- model.matrix(mod1, ssDF[[g]])
        }
        d.nms <- colnames(dm.den[[1]])
        names.beta.den <- paste("d.beta", d.nms, sep = ".")
        pars.d0 <- log(mean((unlist(lapply(scrFrame$caphist, 
            nrow)))/unlist(lapply(ssDF, nrow))))
        pars.beta.den <- c(pars.d0, rep(0.1, length(names.beta.den) - 
            1))
    }
    if (distmet == "ecol" && length(allvars.dist) == 0) {
        message("You specified 'ecological distance' (distmet='ecol') but provided no\ncost surface.\n    Euclidean distance will be used.")
    }
    if (length(allvars.dist) > 0) {
        ccovnms <- colnames(costDF[[1]])
        cCovMissing <- allvars.dist[which(!allvars.dist \%in\% 
            ccovnms)]
        if (length(cCovMissing) > 0) {
            stop("I cant find theses covariates in 'costDF'", 
                for (i in cCovMissing) print(i))
        }
    }
    mod4 <- update(model[[4]], ~. - sex - session - 1)
    if (distmet == "ecol") {
        for (g in 1:ns) {
            dm.cost[[g]] <- my.model.matrix(mod4, costDF[[g]])
            names.dist <- paste("c.beta.", allvars.dist, sep = "")
            pars.dist <- rep(0, length(names.dist))
        }
    }
    if (!smallslow) {
        if (distmet == "euc") {
            for (g in 1:ns) {
                D[[g]] <- e2dist(scrFrame$traps[[g]][, c("X", 
                  "Y")], ssDF[[g]][, c("X", "Y")])
            }
        }
    }
    if ("indCovs" \%in\% names(scrFrame)) {
        if ("sex" \%in\% names(scrFrame$indCovs[[1]])) {
            anySex <- TRUE
        }
    }
    if (var.p0.1 && var.p0.2) {
        pars.p0 <- rnorm(ns * 2, qlogis(0.1), 0.2)
        tmpPsex <- rep(c(1, 2), ns)
        tmpPsess <- rep(1:ns, each = 2)
        names.p0 <- paste("p0.sex", tmpPsex, "session", tmpPsess, 
            sep = "")
        pBothsexnsesh <- TRUE
    }
    else {
        if (var.p0.1) {
            pars.p0 <- rnorm(2, qlogis(0.1), 0.2)
            names.p0 <- c("p0.sex1", "p0.sex2")
            pJustsex <- TRUE
        }
        else {
            if (var.p0.2) {
                pars.p0 <- rnorm(ns, qlogis(0.1), 0.2)
                tmpPsess <- 1:ns
                names.p0 <- paste("p0.session", tmpPsess, sep = "")
                pJustsesh <- TRUE
            }
            else {
                pars.p0 <- rnorm(1, qlogis(0.1), 0.2)
                names.p0 <- c("p0.")
                pDot <- TRUE
            }
        }
    }
    if (any(var.p0.1, var.a1.1) && !anySex) 
        stop("Sex defined in a model but no sex data provided.")
    var.p0.t <- "t" \%in\% allvars.p0
    var.p0.T <- "T" \%in\% allvars.p0
    if (var.p0.t) {
        pars.p0 <- rep(pars.p0, K)
        names.p0 <- paste(rep(names.p0, each = K), "..t", 1:K, 
            sep = "")
        pTime <- TRUE
    }
    if (var.p0.T) {
        pars.p0 <- c(pars.p0, 0)
        names.p0 <- c(names.p0, "T.trend")
    }
    if (var.a1.1 && var.a1.2) {
        pars.a1 <- rnorm(ns * 2, 0, 0.2)
        tmpAsex <- rep(c(1, 2), ns)
        tmpAsess <- rep(1:ns, each = 2)
        names.a1 <- paste("a1.sex", tmpAsex, "session", tmpAsess, 
            sep = "")
        aBothsexnsesh <- TRUE
    }
    else {
        if (var.a1.1) {
            pars.a1 <- rnorm(2, 0, 0.2)
            names.a1 <- c("a1.sex1", "a1.sex2")
            aJustsex <- TRUE
        }
        else {
            if (var.a1.2) {
                pars.a1 <- rnorm(ns, 0.1, 0.2)
                tmpAsess <- 1:ns
                names.a1 <- paste("a1.session", tmpAsess, sep = "")
                aJustsesh <- TRUE
            }
            else {
                pars.a1 <- rnorm(1, 0, 0.2)
                names.a1 <- c("a1.")
                aDot <- TRUE
            }
        }
    }
    if (scrFrame$type == "scr") {
        YY <- scrFrame$caphist
    }
    if (scrFrame$type == "secr") {
        YY <- "secr2scr"
    }
    nG <- sess.ss.nG
    if (anySex) {
        if (sexmod == "constant") {
            pars.sex <- 0
            names.sex <- "psi.constant"
        }
        if (sexmod == "session") {
            pars.sex <- rep(0, ns)
            names.sex <- paste("psi", 1:ns, sep = "")
        }
    }
    else {
        pars.sex <- NULL
        names.sex <- NULL
    }
    pv <- round(c(pars.p0, pars.a1, pars.beta.trap, pars.beta.den, 
        pars.dist, pars.sex), 2)
    pn <- c(names.p0, names.a1, names.beta.trap, names.beta.den, 
        names.dist, names.sex)
    if (!is.null(start.vals)) {
        if (length(pv) == length(start.vals)) {
            pv <- start.vals
        }
        else {
            message("The starting values provided are not the required length, using generated values")
        }
    }
    msLL.nosex <- function(pv = pv, pn = pn, YY = YY, nG = nG, 
        K = K, dm.den = dm.den, dm.trap = dm.trap) {
        alpha0 <- matrix(NA, ns, K)
        if (pDot & !pTime) {
            tmpP <- pv[pn \%in\% names.p0[grep("p0.", names.p0)]]
            alpha0[] <- tmpP
        }
        if (pDot & pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep("p0.", names.p0)]]
                alpha0[s, ] <- tmpP
            }
        }
        if (pJustsesh & !pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep(paste("session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, ] <- tmpP
            }
        }
        if (pJustsesh & pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep(paste("session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, ] <- tmpP
            }
        }
        if (aDot) {
            alpha1 <- exp(rep(pv[pn \%in\% names.a1], ns))
        }
        if (aJustsesh) {
            alpha1 <- exp(pv[pn \%in\% names.a1])
        }
        if (trap.covs) {
            t.beta <- matrix(NA, ns, length(t.nms))
            if ("Session" \%in\% all.vars(model[[2]])) {
                for (s in 1:ns) {
                  t.beta[s, ] <- pv[pn \%in\% names.beta.trap[grep(paste("sess", 
                    s, sep = ""), names.beta.trap)]]
                }
            }
            else {
                for (s in 1:ns) {
                  t.beta[s, ] <- pv[pn \%in\% names.beta.trap]
                }
            }
        }
        if (distmet == "ecol") {
            dist.beta <- pv[pn \%in\% names.dist]
        }
        d.beta <- pv[pn \%in\% names.beta.den]
        if (DorN == "N") 
            d0 <- ifelse(dSession, exp(pv[pn \%in\% names.d0]), 
                rep(exp(pv[pn \%in\% names.d0]), ns))
        outLik <- 0
        for (s in 1:length(YY)) {
            Ys <- YY[[s]]
            if (!multicatch) {
                zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
                Ys <- abind(Ys, zeros, along = 1)
            }
            if (multicatch) {
                zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
                Ys <- abind(Ys, zeros, along = 1)
            }
            if (predict) 
                tmp.post <- matrix(NA, nG[s], nrow(Ys))
            if (distmet == "ecol") {
                cost <- exp(dm.cost[[s]] \%*\% dist.beta)
                costR <- rasterFromXYZ(cbind(costDF[[s]][, c(1, 
                  2)], cost))
                if (is.null(PROJ)) {
                  projection(costR) <- "+proj=utm +zone=12 +datum=WGS84"
                }
                else {
                  projection(costR) <- PROJ
                }
                tr <- transition(costR, transitionFunction = function(x) (1/(mean(x))), 
                  direction = directions)
                trLayer <- geoCorrection(tr, scl = F)
                D <- list()
                D[[s]] <- costDistance(trLayer, as.matrix(scrFrame$traps[[s]][, 
                  c("X", "Y")]), as.matrix(ssDF[[s]][, c("X", 
                  "Y")]))
            }
            if (smallslow) {
                if (distmet == "euc") {
                  D <- list()
                  D[[s]] <- e2dist(scrFrame$traps[[s]][, c("X", 
                    "Y")], ssDF[[s]][, c("X", "Y")])
                }
            }
            tmpD <- e2dist(scrFrame$traps[[s]][, c("X", "Y")], 
                ssDF[[s]][, c("X", "Y")])
            if (DorN == "D") {
                d.s <- exp(dm.den[[s]] \%*\% d.beta)
                pi.s <- d.s/sum(d.s)
            }
            if (DorN == "N") {
                d.s <- exp(dm.den[[s]] \%*\% d.beta)
                pi.s <- d.s/sum(d.s)
            }
            lik.marg <- rep(NA, nrow(Ys))
            trimS <- ifelse(is.null(trimS), max(tmpD), trimS)
            for (i in 1:nrow(Ys)) {
                if (i < nrow(Ys)) {
                  pp <- apply(Ys[i, , ], 1, sum) > 0
                }
                else {
                  pp <- apply(Ys[i, , ], 1, sum) >= 0
                }
                inflate <- (1 - pp) * 2 * trimS
                trimC <- apply(tmpD + inflate, 2, min, na.rm = T) < 
                  trimS
                inflate <- (1 - trimC) * 2 * trimS
                trimR <- apply(t(t(tmpD) + inflate), 1, min, 
                  na.rm = T) < trimS
                Pm <- matrix(0, sum(trimR), sum(trimC))
                if (multicatch) {
                  Pm <- matrix(0, sum(trimR) + 1, sum(trimC))
                }
                if (plotit) {
                  plot(scrFrame$traps[[s]][trimR, c("X", "Y")], 
                    pch = 3, col = 4, asp = 1, main = paste("Session:", 
                      s, " Individual: ", i, " traps: ", sum(pp), 
                      sep = " "))
                  points(ssDF[[s]][trimC, ], pch = 16, col = 2, 
                    cex = mycex)
                  points(scrFrame$traps[[s]][pp, c("X")], scrFrame$traps[[s]][pp, 
                    c("Y")], pch = 16, col = 3, cex = 1.5)
                }
                kk <- ifelse(length(dim(Ys)) == 2, 1, dim(Ys)[3])
                for (k in 1:kk) {
                  if (!trap.covs) {
                    a0 <- rep(alpha0[s, k], nrow(D[[s]]))
                  }
                  else {
                    a0 <- alpha0[s, k] + (dm.trap[[s]][[k]] \%*\% 
                      c(t.beta[s, ]))
                  }
                  probcap <- c(plogis(a0[trimR])) * exp(-alpha1[s] * 
                    D[[s]][trimR, trimC]^2)
                  if (!multicatch) {
                    probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i, 
                      trimR, k], sum(trimC)), 1, probcap[1:length(Pm)], 
                      log = TRUE))
                  }
                  if (multicatch) {
                    probcap <- rbind(probcap, rep(1, sum(trimC)))
                    probcap <- probcap/colSums(probcap)
                    probcap[1:length(probcap)] <- c(dbinom(rep(c(Ys[i, 
                      trimR, k], 1 - sum(Ys[i, trimR, k])), sum(trimC)), 
                      1, probcap[1:length(Pm)], log = TRUE))
                  }
                  if (!is.null(scrFrame$trapOperation)) {
                    probcap <- probcap * scrFrame$trapOperation[[s]][trimR, 
                      k]
                  }
                  Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
                }
                lik.cond <- numeric(nG[s])
                if (!is.matrix(Pm)) 
                  browser()
                lik.cond[trimC] <- exp(colSums(Pm, na.rm = T))
                lik.marg[i] <- sum(lik.cond * pi.s)
                if (predict) {
                  tmp.post[, i] <- (lik.cond * pi.s)/lik.marg[i]
                }
            }
            if (!predict) {
                if (DorN == "N") {
                  nv <- c(rep(1, length(lik.marg) - 1), d0[s])
                  part1 <- lgamma((nrow(Ys) - 1) + d0[s] + 1) - 
                    lgamma(d0[s] + 1)
                  part2 <- sum(nv * log(lik.marg))
                }
                if (DorN == "D") {
                  nv <- c(rep(1, length(lik.marg) - 1), 1)
                  atheta <- 1 - lik.marg[nrow(Ys)]
                  nind <- nrow(Ys) - 1
                  part1 <- nind * log(sum(d.s)) - sum(d.s) * 
                    atheta
                  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
                }
                ll <- -1 * (part1 + part2)
                outLik <- outLik + ll
            }
            if (predict) {
                posterior[[s]] <- cbind(ssDF[[s]][, c("X", "Y")], 
                  tmp.post)
            }
        }
        if (!predict) {
            out <- outLik
            Pm <- probcap <- Ys <- probcap <- NULL
            return(out)
        }
        if (predict) {
            return(posterior)
        }
    }
    msLL.sex <- function(pv, pn, YY, nG, K, dm.den, dm.trap) {
        alpha0 <- array(NA, c(ns, K, 2))
        if (pDot & !pTime) {
            tmpP <- pv[pn \%in\% names.p0[grep("p0.", names.p0)]]
            alpha0[] <- tmpP
        }
        if (pDot & pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep("p0.", names.p0)]]
                alpha0[s, , 1] <- tmpP
                alpha0[s, , 2] <- tmpP
            }
        }
        if (pJustsex & !pTime) {
            tmpP <- pv[pn \%in\% names.p0[grep("sex1", names.p0)]]
            alpha0[, , 1] <- tmpP
            tmpP <- pv[pn \%in\% names.p0[grep("sex2", names.p0)]]
            alpha0[, , 2] <- tmpP
        }
        if (pJustsex & pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep("sex1", names.p0)]]
                alpha0[s, , 1] <- tmpP
                tmpP <- pv[pn \%in\% names.p0[grep("sex2", names.p0)]]
                alpha0[s, , 2] <- tmpP
            }
        }
        if (pJustsesh & !pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep(paste("session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, , ] <- tmpP
            }
        }
        if (pJustsesh & pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep(paste("session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, , 1] <- tmpP
                alpha0[s, , 2] <- tmpP
            }
        }
        if (pBothsexnsesh & !pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep(paste("sex1session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, , 1] <- tmpP
                tmpP <- pv[pn \%in\% names.p0[grep(paste("sex2session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, , 2] <- tmpP
            }
        }
        if (pBothsexnsesh & pTime) {
            for (s in 1:ns) {
                tmpP <- pv[pn \%in\% names.p0[grep(paste("sex1session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, , 1] <- tmpP
                tmpP <- pv[pn \%in\% names.p0[grep(paste("sex2session", 
                  s, sep = ""), names.p0)]]
                alpha0[s, , 2] <- tmpP
            }
        }
        alpha1 <- matrix(NA, length(YY), 2)
        if (aDot) {
            alpha1[] <- exp(pv[pn \%in\% names.a1])
        }
        if (aJustsex) {
            alpha1[, 1] <- exp(pv[pn \%in\% names.a1[grep("sex1", 
                names.a1)]])
            alpha1[, 2] <- exp(pv[pn \%in\% names.a1[grep("sex2", 
                names.a1)]])
        }
        if (aJustsesh) {
            for (s in 1:ns) {
                alpha1[s, ] <- exp(pv[pn \%in\% names.a1[grep(paste("session", 
                  s, sep = ""), names.a1)]])
            }
        }
        if (aBothsexnsesh) {
            for (s in 1:ns) {
                alpha1[s, 1] <- exp(pv[pn \%in\% names.a1[grep(paste("sex1session", 
                  s, sep = ""), names.a1)]])
                alpha1[s, 2] <- exp(pv[pn \%in\% names.a1[grep(paste("sex2session", 
                  s, sep = ""), names.a1)]])
            }
        }
        if (trap.covs) {
            t.beta <- matrix(NA, ns, length(names.beta.trap))
            if ("Session" \%in\% all.vars(model[[2]])) {
                for (s in 1:ns) {
                  t.beta[s, ] <- pv[pn \%in\% names.beta.trap[grep(paste("sess", 
                    s, sep = ""), names.beta.trap)]]
                }
            }
            else {
                for (s in 1:ns) {
                  t.beta[s, ] <- pv[pn \%in\% names.beta.trap]
                }
            }
        }
        if (distmet == "ecol") {
            dist.beta <- pv[pn \%in\% names.dist]
        }
        d.beta <- pv[pn \%in\% names.beta.den]
        if (DorN == "N") 
            d0 <- ifelse(dSession, exp(pv[pn \%in\% names.d0]), 
                rep(exp(pv[pn \%in\% names.d0]), ns))
        if (sexmod == "constant") 
            psi.sex <- rep(plogis(pv[grep("psi", pn)]), ns)
        if (sexmod == "session") 
            psi.sex <- plogis(pv[grep("psi", pn)])
        outLik <- 0
        for (s in 1:length(YY)) {
            Ys <- YY[[s]]
            zeros <- array(0, c(1, dim(Ys)[2], dim(Ys)[3]))
            Ys <- abind(Ys, zeros, along = 1)
            sx <- c(scrFrame$indCovs[[s]]$sex + 1, NA)
            if (predict) 
                tmp.post <- matrix(NA, nG[s], nrow(Ys))
            if (distmet == "ecol") {
                cost <- exp(dm.cost[[s]] \%*\% dist.beta)
                costR <- rasterFromXYZ(cbind(costDF[[s]][, c(1, 
                  2)], cost))
                if (is.null(PROJ)) {
                  projection(costR) <- "+proj=utm +zone=12 +datum=WGS84"
                }
                else {
                  projection(costR) <- PROJ
                }
                tr <- transition(costR, transitionFunction = function(x) (1/(mean(x))), 
                  direction = directions)
                trLayer <- geoCorrection(tr, scl = F)
                D <- list()
                D[[s]] <- costDistance(trLayer, as.matrix(scrFrame$traps[[s]][, 
                  c("X", "Y")]), as.matrix(ssDF[[s]][, c("X", 
                  "Y")]))
            }
            if (smallslow) {
                if (distmet == "euc") {
                  D <- list()
                  D[[s]] <- e2dist(scrFrame$traps[[s]][, c("X", 
                    "Y")], ssDF[[s]][, c("X", "Y")])
                }
            }
            tmpD <- e2dist(scrFrame$traps[[s]][, c("X", "Y")], 
                ssDF[[s]][, c("X", "Y")])
            if (DorN == "D") {
                d.s <- exp(dm.den[[s]] \%*\% d.beta)
                pi.s <- d.s/sum(d.s)
            }
            if (DorN == "N") {
                d.s <- exp(dm.den[[s]] \%*\% d.beta)
                pi.s <- d.s/sum(d.s)
            }
            lik.marg <- lik.marg1 <- lik.marg2 <- rep(NA, nrow(Ys))
            trimS <- ifelse(is.null(trimS), max(tmpD), trimS)
            for (i in 1:nrow(Ys)) {
                if (i < nrow(Ys)) {
                  pp <- apply(Ys[i, , ], 1, sum) > 0
                }
                else {
                  pp <- apply(Ys[i, , ], 1, sum) >= 0
                }
                inflate <- (1 - pp) * 2 * trimS
                trimC <- apply(tmpD + inflate, 2, min, na.rm = T) < 
                  trimS
                inflate <- (1 - trimC) * 2 * trimS
                trimR <- apply(t(t(tmpD) + inflate), 1, min, 
                  na.rm = T) < trimS
                Pm <- Pm1 <- Pm2 <- tmpPm <- NULL
                if (plotit) {
                  plot(scrFrame$traps[[s]][trimR, c("X", "Y")], 
                    pch = 3, col = 4, asp = 1, main = paste("Session:", 
                      s, " Individual: ", i, " traps: ", sum(pp), 
                      sep = " "))
                  points(ssDF[[s]][trimC, ], pch = 16, col = 2, 
                    cex = mycex)
                  points(scrFrame$traps[[s]][pp, c("X")], scrFrame$traps[[s]][pp, 
                    c("Y")], pch = 16, col = 3, cex = 1.5)
                }
                if (!is.na(sx[i])) {
                  Pm <- matrix(0, sum(trimR), sum(trimC))
                  kk <- ifelse(length(dim(Ys)) == 2, 1, dim(Ys)[3])
                  for (k in 1:kk) {
                    if (!trap.covs) {
                      a0 <- rep(alpha0[s, k, sx[i]], nrow(D[[s]]))
                    }
                    else {
                      a0 <- alpha0[s, k, sx[i]] + (dm.trap[[s]][[k]] \%*\% 
                        c(t.beta[s, ]))
                    }
                    probcap <- c(plogis(a0[trimR])) * exp(-alpha1[s, 
                      sx[i]] * D[[s]][trimR, trimC]^2)
                    probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i, 
                      trimR, k], sum(trimC)), 1, probcap[1:length(probcap)], 
                      log = TRUE))
                    if (!is.null(scrFrame$trapOperation)) {
                      probcap <- probcap * scrFrame$trapOperation[[s]][trimR, 
                        k]
                    }
                    Pm[1:length(Pm)] <- Pm[1:length(Pm)] + probcap[1:length(probcap)]
                  }
                  lik.cond <- numeric(nG[s])
                  lik.cond[trimC] <- exp(colSums(Pm, na.rm = T))
                  tmpPsi <- (sx[i] == 1) * (1 - psi.sex[s]) + 
                    (sx[i] == 2) * psi.sex[s]
                  lik.marg[i] <- sum(lik.cond * pi.s) * tmpPsi
                  if (predict) {
                    tmp.post[, i] <- (lik.cond * pi.s)/lik.marg[i]
                  }
                }
                else {
                  Pm1 <- Pm2 <- matrix(0, sum(trimR), sum(trimC))
                  kk <- ifelse(length(dim(Ys)) == 2, 1, dim(Ys)[3])
                  for (k in 1:kk) {
                    if (!trap.covs) {
                      a0.1 <- rep(alpha0[s, k, 1], length(trimR))
                      a0.2 <- rep(alpha0[s, k, 2], length(trimR))
                    }
                    else {
                      a0.1 <- alpha0[s, k, 1] + (dm.trap[[s]][[k]] \%*\% 
                        c(t.beta[s, ]))
                      a0.2 <- alpha0[s, k, 2] + (dm.trap[[s]][[k]] \%*\% 
                        c(t.beta[s, ]))
                    }
                    probcap <- c(plogis(a0.1[trimR])) * exp(-alpha1[s, 
                      1] * D[[s]][trimR, trimC]^2)
                    probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i, 
                      trimR, k], sum(trimC)), 1, probcap[1:length(probcap)], 
                      log = TRUE))
                    if (!is.null(scrFrame$trapOperation)) {
                      probcap <- probcap * scrFrame$trapOperation[[s]][trimR, 
                        k]
                    }
                    Pm1[1:length(Pm1)] <- Pm1[1:length(Pm1)] + 
                      probcap[1:length(probcap)]
                    probcap <- c(plogis(a0.2[trimR])) * exp(-alpha1[s, 
                      2] * D[[s]][trimR, trimC]^2)
                    probcap[1:length(probcap)] <- c(dbinom(rep(Ys[i, 
                      trimR, k], sum(trimC)), 1, probcap[1:length(probcap)], 
                      log = TRUE))
                    if (!is.null(scrFrame$trapOperation)) {
                      probcap <- probcap * scrFrame$trapOperation[[s]][trimR, 
                        k]
                    }
                    Pm2[1:length(Pm2)] <- Pm2[1:length(Pm2)] + 
                      probcap[1:length(probcap)]
                  }
                  lik.cond <- lik.cond1 <- lik.cond2 <- numeric(nG[s])
                  lik.cond1[trimC] <- exp(colSums(Pm1, na.rm = T))
                  lik.cond2[trimC] <- exp(colSums(Pm2, na.rm = T))
                  lik.cond <- lik.cond1 * (1 - psi.sex[s]) + 
                    lik.cond2 * psi.sex[s]
                  lik.marg1[i] <- sum(lik.cond1 * pi.s)
                  lik.marg2[i] <- sum(lik.cond2 * pi.s)
                  lik.marg[i] <- lik.marg1[i] * (1 - psi.sex[s]) + 
                    lik.marg2[i] * psi.sex[s]
                  if (predict) {
                    tmp.post[, i] <- (lik.cond * pi.s)/lik.marg[i]
                  }
                }
            }
            if (!predict) {
                if (DorN == "N") {
                  nv <- c(rep(1, length(lik.marg) - 1), d0[s])
                  part1 <- lgamma((nrow(Ys) - 1) + d0[s] + 1) - 
                    lgamma(d0[s] + 1)
                  part2 <- sum(nv * log(lik.marg))
                }
                if (DorN == "D") {
                  nv <- c(rep(1, length(lik.marg) - 1), 1)
                  atheta <- 1 - lik.marg[nrow(Ys)]
                  nind <- nrow(Ys) - 1
                  part1 <- nind * log(sum(d.s)) - sum(d.s) * 
                    atheta
                  part2 <- sum(nv[1:nind] * log(lik.marg[1:nind]))
                }
                ll <- -1 * (part1 + part2)
                outLik <- outLik + ll
            }
            if (predict) {
                posterior[[s]] <- cbind(ssDF[[s]][, c("X", "Y")], 
                  tmp.post)
            }
        }
        if (!predict) {
            out <- outLik
            Pm <- probcap <- Ys <- probcap <- NULL
            return(out)
        }
        if (predict) {
            return(posterior)
        }
    }
    if (!predict) {
        message("Fitting model: D", paste(model)[1], ", p0", 
            paste(model)[2], ", sigma", paste(model)[3], ", cost", 
            paste(model)[4], sep = " ")
        if (!anySex) {
            message("Using ll function 'msLL.nosex' \nHold on tight!")
            message(paste(pn, " ", sep = " | "))
            message(" ")
            myfit <- nlm(msLL.nosex, p = pv, pn = pn, YY = YY, 
                nG = sess.ss.nG, K = K, dm.den = dm.den, dm.trap = dm.trap, 
                hessian = T, stepmax = 10)
        }
        else {
            message("Using ll function 'msLL.sex' \nHold on tight!")
            message(paste(pn, " ", sep = " | "))
            message(" ")
            myfit <- nlm(msLL.sex, p = pv, pn = pn, YY = YY, 
                nG = sess.ss.nG, K = K, dm.den = dm.den, dm.trap = dm.trap, 
                hessian = T, stepmax = 10)
        }
        links <- rep(NA, length(pn))
        pars <- myfit$estimate
        links[grep("p0.", pn)] <- "(logit)"
        links[grep("a1.", pn)] <- "(log)"
        links[grep("psi", pn)] <- "(logit)"
        links[grep("beta", pn)] <- "(Identity)"
        mle.se <- NULL
        trans.mle <- rep(0, length(pv))
        trans.mle[grep("p0.", pn)] <- plogis(pars[grep("p0.", 
            pn)])
        trans.mle[grep("a1.", pn)] <- exp(pars[grep("a1.", pn)])
        trans.mle[grep("psi", pn)] <- plogis(pars[grep("psi", 
            pn)])
        trans.mle[grep("d.beta", pn)] <- pars[grep("d.beta", 
            pn)]
        if ("hessian" \%in\% names(myfit)) {
            trans.se <- rep(NA, length(pv))
        }
        else {
            sese <- rep(rep(NA, length(pv)))
            trans.se <- rep(NA, length(pv))
        }
        outStats <- data.frame(parameters = pn, link = links, 
            mle = round(myfit$estimate, 3), se = round(myfit$estimate, 
                3) * 0, mle.tr = round(trans.mle, 3), se.tr = trans.se)
        VcV <- NULL
        endtime <- format(Sys.time(), "\%H:\%M:\%S \%d \%b \%Y")
        idx <- outStats[, "parameters"] == "a1."
        sigma <- sqrt(1/(2 * exp(outStats[idx, "mle"])))
        output <- list(call = cl, rawOutput = myfit, outStats = outStats, 
            sigma = sigma, multicatch = multicatch, Area = sess.ssArea, 
            nll = myfit$minimum, AIC = 2 * myfit$minimum + 2 * 
                length(myfit$estimate), started = starttime, 
            ended = endtime, proctime = (proc.time() - ptm)[3])
        class(output) <- "oSCR.fit"
        return(output)
    }
    if (predict) {
        message("Predicting model: D", paste(model)[1], ", p0", 
            paste(model)[2], ", sigma", paste(model)[3], ", cost", 
            paste(model)[4], sep = " ")
        if (!anySex) {
            message("Using ll function 'msLL.nosex' \nHold on tight!")
            message(paste(pn, " ", sep = " | "))
            message(" ")
            myfit <- msLL.nosex(p = start.vals, pn = pn, YY = YY, 
                nG = sess.ss.nG, K = K, dm.den = dm.den, dm.trap = dm.trap)
        }
        else {
            message("Using ll function 'msLL.sex' \nHold on tight!")
            message(paste(pn, " ", sep = " | "))
            message(" ")
            myfit <- msLL.sex(p = start.vals, pn = pn, YY = YY, 
                nG = sess.ss.nG, K = K, dm.den = dm.den, dm.trap = dm.trap)
        }
        return(myfit)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
