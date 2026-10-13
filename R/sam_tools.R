#-------------------------------------------------------------------------------
# Wrappers around stockassessment functions that work on the sam fit stored
# in an FLSAM object (see samfit) and return data.frames
#-------------------------------------------------------------------------------

#- Convert a samforecast table into a long data.frame
.samforecast2df <- function(x){
  tab   <- attr(x,"tab")
  cols  <- do.call(rbind,strsplit(colnames(tab),":"))
  stats <- unique(cols[,2])
  est   <- stats[!stats %in% c("low","high")][1]
  res   <- do.call(rbind,lapply(unique(cols[,1]),function(q){
             data.frame(year=as.numeric(rownames(tab)),quantity=q,
                        value=tab[,paste(q,est,sep=":")],
                        lbnd=tab[,paste(q,"low",sep=":")],
                        ubnd=tab[,paste(q,"high",sep=":")],
                        row.names=NULL)}))
  attr(res,"samforecast") <- x
  return(res)
}

#' Short-term forecasts from an FLSAM assessment
#'
#' \code{FLSAM.forecast} runs the stochastic forecast of
#' \code{\link[stockassessment]{forecast}} and \code{FLSAM.modelforecast} the
#' model-based forecast of \code{\link[stockassessment]{modelforecast}} on the
#' sam fit stored in an \code{\link{FLSAM-class}} object.
#'
#' @param sam An \code{FLSAM} object with a stored fit (see \code{\link{samfit}}).
#' @param ... Arguments passed to \code{stockassessment::forecast} or
#'   \code{stockassessment::modelforecast}, e.g. \code{fscale}, \code{fval},
#'   \code{catchval}, \code{ave.years}, \code{rec.years}, \code{nosim}.
#' @return A \code{data.frame} with columns \code{year}, \code{quantity}
#'   (fbar, rec, ssb, catch), \code{value} (median or most likely trajectory),
#'   \code{lbnd} and \code{ubnd}. The original \code{samforecast} object is
#'   stored in the attribute \code{"samforecast"} and can be plotted with
#'   \code{plot(attr(x, "samforecast"))}.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)
#' fc  <- FLSAM.forecast(sam, fscale=c(1,1,1), ave.years=2009:2011, rec.years=2000:2011)
#' fc
#' }
#' @rdname FLSAM.forecast
#' @export
FLSAM.forecast <- function(sam,...){
  .samforecast2df(stockassessment::forecast(samfit(sam),...))
}

#' @rdname FLSAM.forecast
#' @export
FLSAM.modelforecast <- function(sam,...){
  .samforecast2df(stockassessment::modelforecast(samfit(sam),...))
}

#' Reference points from an FLSAM assessment
#'
#' Estimates deterministic reference points with
#' \code{\link[stockassessment]{deterministicReferencepoints}} on the sam fit
#' stored in an \code{\link{FLSAM-class}} object. MSY-based reference points
#' require a stock-recruitment model (\code{srr} other than 0); per-recruit
#' reference points (e.g. \code{"Max"}, \code{"0.1dYPR"} for F0.1,
#' \code{"0.35SPR"}) do not.
#'
#' @param sam An \code{FLSAM} object with a stored fit (see \code{\link{samfit}}).
#' @param referencepoints Character vector of reference points, see
#'   \code{stockassessment::deterministicReferencepoints}.
#' @param ... Further arguments passed on, e.g. \code{nYears}, \code{aveYears},
#'   \code{selYears}.
#' @return A \code{data.frame} with columns \code{refpoint}, \code{quantity},
#'   \code{value}, \code{lbnd}, \code{ubnd}. The original object is stored in
#'   the attribute \code{"sam_referencepoints"}.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)
#' FLSAM.refpts(sam, c("Max","0.1dYPR","0.35SPR"))
#' }
#' @export
FLSAM.refpts <- function(sam,referencepoints,...){
  rp  <- stockassessment::deterministicReferencepoints(samfit(sam),referencepoints,...)
  res <- do.call(rbind,lapply(names(rp$tables),function(q){
           tab <- rp$tables[[q]]
           if(is.null(dim(tab)) || ncol(tab)<3) return(NULL)
           data.frame(refpoint=rownames(tab),quantity=q,value=tab[,1],
                      lbnd=tab[,2],ubnd=tab[,3],row.names=NULL)}))
  attr(res,"sam_referencepoints") <- rp
  return(res)
}

#' Yield and spawning biomass per recruit
#'
#' Calculates yield and SSB per recruit curves and the per-recruit reference
#' points (Fmax, F0.1 and F at a given SPR) with
#' \code{\link[stockassessment]{ypr}}.
#'
#' @param sam An \code{FLSAM} object with a stored fit (see \code{\link{samfit}}).
#' @param ... Arguments passed to \code{stockassessment::ypr}, e.g.
#'   \code{Flimit}, \code{Fdelta}, \code{aveYears}, \code{sprProp}.
#' @return A \code{data.frame} with columns \code{fbar}, \code{ssb} and
#'   \code{yield} (per recruit); the reference points are in the attribute
#'   \code{"refpts"}.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)
#' y <- FLSAM.ypr(sam)
#' attr(y, "refpts")
#' }
#' @export
FLSAM.ypr <- function(sam,...){
  y   <- stockassessment::ypr(samfit(sam),...)
  res <- data.frame(fbar=y$fbar,ssb=y$ssb,yield=y$yield)
  attr(res,"refpts") <- unlist(y[c("fmax","f01","fsprProp","sprProp")])
  return(res)
}

#' Catch by fleet
#'
#' Estimated catch (in weight) per catch fleet with confidence intervals, from
#' \code{\link[stockassessment]{catchbyfleettable}}.
#'
#' @param sam An \code{FLSAM} object with a stored fit (see \code{\link{samfit}}).
#' @param obs.show Also return the observed catch.
#' @return A \code{data.frame} with columns \code{year}, \code{fleet},
#'   \code{value}, \code{lbnd}, \code{ubnd} (and \code{obs}).
#' @export
FLSAM.catchbyfleet <- function(sam,obs.show=FALSE){
  tab    <- stockassessment::catchbyfleettable(samfit(sam),obs.show=obs.show)
  fleets <- names(sam@control@fleets)[sam@control@fleets==0]
  nf     <- ncol(tab)/ifelse(obs.show,4,3)
  if(nf!=length(fleets)) fleets <- sub("^Catch\\((.*)\\)$","\\1",colnames(tab)[seq_len(nf)])
  #- columns are blocks of estimate, low, high (and observed) by fleet
  res    <- do.call(rbind,lapply(seq_along(fleets),function(i){
              cols <- (0:3)*nf+i
              df   <- data.frame(year=as.numeric(rownames(tab)),fleet=fleets[i],
                                 value=tab[,cols[1]],lbnd=tab[,cols[2]],ubnd=tab[,cols[3]],row.names=NULL)
              if(obs.show) df$obs <- tab[,cols[4]]
              df}))
  return(res)
}

#' Text description of the assessment model
#'
#' Returns the description of the model configuration as written by
#' \code{\link[stockassessment]{modelDescription}}.
#'
#' @param sam An \code{FLSAM} object with a stored fit (see \code{\link{samfit}}).
#' @return A character string (printed with \code{cat}, returned invisibly).
#' @export
model.description <- function(sam){
  fit <- samfit(sam)
  #- modelDescription pastes per fleet when maxAgePlusGroup is fleet specific: use the stock value
  fit$conf$maxAgePlusGroup <- fit$conf$maxAgePlusGroup[1]
  txt <- paste(utils::capture.output(stockassessment::modelDescription(fit)),collapse="\n")
  txt <- sub("\nNULL$","",txt)
  cat(txt,"\n")
  return(invisible(txt))
}
