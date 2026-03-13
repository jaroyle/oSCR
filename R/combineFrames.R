#' Combine Multiple Single-Session scrFrame Objects into a Multi-Session scrFrame
#'
#' This function is the inverse of \code{subFrame}: it takes a list of
#' \code{scrFrame} objects (typically single-session) and combines them into a
#' single multi-session \code{scrFrame} suitable for use with \code{oSCR.fit}.
#'
#' @param sfList A list of objects of class \code{scrFrame}, one per session.
#'   Each element should ideally be a single-session scrFrame, but multi-session
#'   frames are also supported (their sessions are appended in order).
#'
#' @return An object of class \code{scrFrame} containing all sessions from the
#'   input list, with the following components:
#'   \item{caphist}{List of n x J x K capture history arrays, one per session.}
#'   \item{traps}{List of trap coordinate data frames (columns X and Y), one per session.}
#'   \item{indCovs}{List of individual covariate data frames, one per session,
#'     or NULL if absent in all input frames.}
#'   \item{trapCovs}{List of lists of occasion-level trap covariate data frames,
#'     one per session, or NULL if absent in all input frames.}
#'   \item{trapOperation}{List of J x K binary trap operation matrices, one per
#'     session, or NULL if absent in all input frames.}
#'   \item{telemetry}{Combined telemetry object (fixfreq, indCovs, cap.tel),
#'     or NULL if absent in all input frames.}
#'   \item{nocc}{Maximum number of occasions across all sessions.}
#'   \item{type}{Taken from the first scrFrame in the list.}
#'
#' @details
#' Components that are NULL in some frames but present in others are padded with
#' NULL for the missing sessions, matching the behaviour of \code{make.scrFrame}.
#' A warning is issued in this case so you are aware of the mismatch.
#'
#' Telemetry objects are combined by appending the per-session elements of
#' \code{fixfreq}, \code{indCovs}, and \code{cap.tel} across all input frames.
#'
#' @author Based on oSCR package conventions by Sutherland, Royle & Linden.
#'
#' @seealso \code{\link{subFrame}}, \code{\link{make.scrFrame}}, \code{\link{oSCR.fit}}
#'
#' @examples
#' ## Simulate two single-session scrFrames and combine them
#' # library(oSCR)
#' # d1 <- sim.SCR.ms(N=50, K=5, alpha0=-1, sigma=0.5, nSession=1)
#' # d2 <- sim.SCR.ms(N=60, K=5, alpha0=-1, sigma=0.5, nSession=1)
#' # sf1 <- make.scrFrame(caphist=d1$y3d, traps=d1$traplocs)
#' # sf2 <- make.scrFrame(caphist=d2$y3d, traps=d2$traplocs)
#' # sf.combined <- combineFrames(list(sf1, sf2))
#'
#' ## Round-trip check: split and recombine
#' # data(rbs)
#' # sf1 <- subFrame(rbs$scrFrame, 1)
#' # sf2 <- subFrame(rbs$scrFrame, 2)
#' # sf.back <- combineFrames(list(sf1, sf2))

combineFrames <- function(sfList) {

  # --- Input checks ---
  if (!is.list(sfList) || length(sfList) < 1) {
    stop("'sfList' must be a non-empty list of scrFrame objects.")
  }
  if (!all(sapply(sfList, inherits, "scrFrame"))) {
    stop("All elements of 'sfList' must be of class 'scrFrame'.")
  }

  # Components that are per-session lists (excluding scalars nocc, type)
  session.components <- c("caphist", "traps", "indCovs", "trapCovs",
                          "trapOperation", "telemetry")

  # Check which components are present (non-NULL) in at least one frame
  present <- sapply(session.components, function(comp) {
    any(sapply(sfList, function(sf) !is.null(sf[[comp]])))
  })

  # Warn if a component is present in some frames but not others
  partial <- sapply(session.components, function(comp) {
    vals <- sapply(sfList, function(sf) !is.null(sf[[comp]]))
    any(vals) && !all(vals)
  })
  if (any(partial)) {
    warning(
      "The following components are present in some scrFrames but not others; ",
      "missing sessions will be padded with NULL:\n  ",
      paste(session.components[partial], collapse = ", ")
    )
  }

  # --- Helper: flatten one component across all frames into a single list ---
  # Each scrFrame stores per-session data as a list; we concatenate them.
  flatten_component <- function(comp) {
    if (!present[comp]) return(NULL)
    out <- list()
    for (sf in sfList) {
      if (is.null(sf[[comp]])) {
        # Pad with one NULL per session in this frame
        n_sess <- length(sf$caphist)
        for (s in seq_len(n_sess)) out <- c(out, list(NULL))
      } else {
        out <- c(out, sf[[comp]])
      }
    }
    out
  }

  # --- Helper: combine telemetry (itself a named list of per-session lists) ---
  combine_telemetry <- function() {
    if (!present["telemetry"]) return(NULL)
    tel.out <- list(fixfreq = list(), indCovs = list(), cap.tel = list())
    has.cap.tel <- FALSE
    for (sf in sfList) {
      if (is.null(sf$telemetry)) {
        n_sess <- length(sf$caphist)
        for (s in seq_len(n_sess)) {
          tel.out$fixfreq  <- c(tel.out$fixfreq,  list(NULL))
          tel.out$indCovs  <- c(tel.out$indCovs,  list(NULL))
          tel.out$cap.tel  <- c(tel.out$cap.tel,  list(NULL))
        }
      } else {
        tel.out$fixfreq <- c(tel.out$fixfreq, sf$telemetry$fixfreq)
        tel.out$indCovs <- c(tel.out$indCovs, sf$telemetry$indCovs)
        if (!is.null(sf$telemetry$cap.tel)) {
          tel.out$cap.tel <- c(tel.out$cap.tel, sf$telemetry$cap.tel)
          has.cap.tel <- TRUE
        } else {
          n_sess <- length(sf$caphist)
          for (s in seq_len(n_sess)) {
            tel.out$cap.tel <- c(tel.out$cap.tel, list(NULL))
          }
        }
      }
    }
    if (!has.cap.tel) tel.out$cap.tel <- NULL
    tel.out
  }

  # --- Assemble the combined scrFrame ---
  sf.combined <- list()
  sf.combined$caphist       <- flatten_component("caphist")
  sf.combined$traps         <- flatten_component("traps")
  sf.combined$indCovs       <- flatten_component("indCovs")
  sf.combined$trapCovs      <- flatten_component("trapCovs")
  sf.combined$trapOperation <- flatten_component("trapOperation")
  sf.combined$telemetry     <- combine_telemetry()

  # nocc: maximum occasions across all sessions
  sf.combined$nocc <- max(sapply(sf.combined$caphist, function(x) dim(x)[3]))

  # occasions: per-session occasion counts, required by do.trim (called when trimS is set)
  sf.combined$occasions <- unlist(lapply(sf.combined$caphist, function(x) dim(x)[3]))

  # mmdm: per-session mean maximum distance moved — concatenate across all frames.
  # If present in some frames but not others, warn and pad missing sessions with NA.
  has.mmdm <- sapply(sfList, function(sf) !is.null(sf$mmdm))
  if (any(has.mmdm)) {
    if (!all(has.mmdm)) {
      warning("'mmdm' is present in some scrFrames but not others; missing sessions padded with NA.")
    }
    sf.combined$mmdm <- unlist(lapply(sfList, function(sf) {
      if (is.null(sf$mmdm)) rep(NA_real_, length(sf$caphist)) else sf$mmdm
    }))
  }

  # type: inherit from first frame (should be "scr" for all)
  sf.combined$type <- sfList[[1]]$type

  class(sf.combined) <- "scrFrame"
  return(sf.combined)
}
