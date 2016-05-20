secr2oscr <-
function (edf, sess.col = 1, id.col = 2, occ.col = 3, trap.col = 4, 
    sex.col = NULL, tdf = NULL, K = NULL, ntraps = NULL, remove.zeros = FALSE, 
    remove.extracaps = FALSE, sex.nacode = NULL, tdf.sep = "/") 
{
    if (!is.null(sex.col) & is.null(sex.nacode)) {
        ux <- length(unique(edf[, sex.col]))
        if (ux > 2) {
            cat("error: more than 2 sex codes, no sex.nacode specified", 
                fill = TRUE)
            return(NULL)
        }
    }
    if (!is.null(sex.col) & !is.null(sex.nacode)) {
        edf[edf[, sex.col] %in% sex.nacode, sex.col] <- NA
        ux <- length(unique(edf[, sex.col][!is.na(edf[, sex.col])]))
        if (ux > 2) {
            cat("error: more than 2 sex codes, no sex.nacode specified", 
                fill = TRUE)
            return(NULL)
        }
    }
    Xsess <- edf[, sess.col]
    if (!is.numeric(Xsess)) {
        Xsess <- as.numeric(as.factor(as.character(Xsess)))
    }
    Xid <- edf[, id.col]
    if (!is.numeric(Xid)) {
        Xid <- as.numeric(as.factor(as.character(Xid)))
    }
    Xocc <- edf[, occ.col]
    if (!is.numeric(Xocc)) {
        Xocc <- as.numeric(as.factor(as.character(Xocc)))
    }
    Xtrap <- edf[, trap.col]
      trap.names <- Xtrap
    if (!is.numeric(Xtrap)) {
        Xtrap <- as.numeric(as.factor(as.character(Xtrap)))
    }
    out <- cbind(Xsess, Xid, Xocc, Xtrap)
    colnames(out) <- c("session", "individual", "occasion", "trap")
    if (!is.null(sex.col)) {
        Xsex <- edf[, sex.col]
        if (!is.numeric(Xsex)) {
            Xsex <- as.numeric(as.factor(as.character(Xsex))) - 
                1
        }
        out <- cbind(out, sex = Xsex)
    }
    trap.names <- trap.names[order(out[, "individual"])]
    out <- out[order(out[, "individual"]), ]
    nsess <- max(out[, "session"])
    caphist <- list()
    nn <- list()
    usex <- list()
    for (s in 1:nsess) {
        xx <- out[out[, "session"] == s, ]
        nind <- max(xx[, "individual"])
        if (is.null(K[s])) {
            Ks <- max(xx[, "occasion"])
        }
        else {
            Ks <- K[s]
        }
        if (is.null(ntraps[s])) {
            ntrapss <- max(xx[, "trap"])
        }
        else {
            ntrapss <- ntraps[s]
        }
        y3d <- array(0, c(nind, ntrapss, Ks))
        for (obs in 1:nrow(xx)) {
            y3d[xx[obs, "individual"], xx[obs, "trap"], xx[obs, 
                "occasion"]] <- y3d[xx[obs, "individual"], xx[obs, 
                "trap"], xx[obs, "occasion"]] + 1
        }
        caphist[[s]] <- y3d
        nn[[s]] <- apply(y3d, c(1), sum)
        if (remove.zeros) 
            caphist[[s]] <- caphist[[s]][nn[[s]] > 0, , ]
        if (remove.extracaps) 
            caphist[[s]][caphist[[s]] > 1] <- 1
        nn[[s]] <- apply(caphist[[s]], c(1), sum)
        if (any(nn[[s]] == 0)) 
            cat("Some individuals in session", s, " have 0 captures", 
                fill = TRUE)
        if (max(caphist[[s]] > 1)) 
            cat("Some individuals in session", s, " captured > 1 time in a trap/occasion", 
                fill = TRUE)
        if(!is.null(sex.col)) usex[[s]] <- xx[!duplicated(xx[, c("individual", "sex")]),             c("individual", "sex")]
    }
    if (is.null(K)) {
        K <- NULL
        for (s in 1:length(caphist)) {
            K[s] <- dim(caphist[[s]])[3]
        }
    }
    traplocs <- NULL
    trapopp <- NULL
    trapcovs <- NULL
    if (!is.null(tdf)) {
        if (is.list(tdf)) 
            ntdf <- length(tdf)
        if (!is.list(tdf)) 
            ntdf <- 1
        if (ntdf != nsess) {
            if (length(K)>1 & var(K) != 0) {
                cat("Error: variable trap operation period indicated but not provided", 
                  fill = TRUE)
                return(NULL)
            }
            cat("Warning: # tdf files not equal to number of sessions. Assuming trap coords constant across sessions", 
                fill = TRUE)
            tdfx <- list()
            for (s in 1:nsess) {
                tdfx[[s]] <- tdf[[1]]
            }
            tdf <- tdfx
        }
        if (is.list(tdf)) {
            all.tcnames <- NULL
            traplocs <- list(NULL)
            trapopp <- list(NULL)
            trapcovs <- list(NULL)
            for (s in 1:length(tdf)) {
                allnames <- tdf[[s]][, 1]
                if (any(is.na(match(trap.names[out[, 1] == s], 
                  allnames)))) {
                  cat("some trap names in EDF not in TDF", fill = TRUE)
                  return(NULL)
                }
                traplocs[[s]] <- as.matrix(tdf[[s]][, 2:3])
                colnames(traplocs[[s]])<-c("X","Y")
                xx <- tdf[[s]][, 4:ncol(tdf[[s]])]
                is.trapcovs <- any(xx[1, ] == tdf.sep)
                if (is.trapcovs) {
                  xx.check <- (1:ncol(xx))[xx[1, ] == tdf.sep]
                  tc.nams <- dimnames(xx)[[2]][(xx.check + 1):ncol(xx)]
                  trapcovs[[s]] <- as.matrix(xx[, (xx.check + 
                    1):ncol(xx)])
                  colnames(trapcovs[[s]]) <- tc.nams
                  trapopp[[s]] <- as.matrix(xx[, 1:(xx.check - 
                    1)])
                  all.tcnames <- c(all.tcnames, tc.nams)
                }
                else {
                  trapopp[[s]] <- as.matrix(xx[, 1:ncol(xx)])
                }
            }
        }
        else {
            allnames <- tdf[, 1]
            if (any(is.na(match(trap.names, allnames)))) {
                cat("some trap names in EDF not in TDF", fill = TRUE)
                return(NULL)
            }
            all.tcnames <- NULL
            traplocs <- list(NULL)
            trapopp <- list(NULL)
            trapcovs <- list(NULL)
            traplocs[[1]] <- as.matrix(tdf[, 2:3])
            xx <- tdf[, 4:ncol(tdf)]
            is.trapcovs <- any(xx[1, ] == tdf.sep)
            if (is.trapcovs) {
                xx.check <- (1:ncol(xx))[xx[1, ] == tdf.sep]
                tc.nams <- dimnames(xx)[[2]][(xx.check + 1):ncol(xx)]
                trapcovs[[1]] <- as.matrix(xx[, (xx.check + 1):ncol(xx)])
                colnames(trapcovs[[1]]) <- tc.nams
                trapopp[[1]] <- as.matrix(xx[, 1:(xx.check - 
                  1)])
                all.tcnames <- c(all.tcnames, tc.nams)
            }
            else {
                trapopp[[1]] <- as.matrix(xx[, 1:ncol(xx)])
            }
        }
        for (s in 1:length(trapopp)) {
            if (any(is.na(trapopp))) {
                return("Error: missing values not allowed in trap operation")
            }
        }
    }
    if(!is.null(sex.col)) sex.oscr = list(data.frame(sex = usex[[1]][, 2]), data.frame(sex = usex[[2]][,         2]))
    else  sex.oscr = NULL
    
    list(edf = out, y3d = caphist, sex = sex.oscr, traplocs = traplocs, 
        trapopp = trapopp, trapcovs = trapcovs)
}
