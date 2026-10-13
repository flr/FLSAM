setOldClass("sam")
setClassUnion("samOrNULL", c("sam","NULL"))

#' FLSAM class
#'
#' Result of an \code{\link{FLSAM}} stock assessment.
#'
#' @slot name,desc,range As in \code{FLComp}: name, description and age/year/fbar range.
#' @slot control The \code{\link{FLSAM.control-class}} object used.
#' @slot nopar Number of estimated (fixed-effect) parameters.
#' @slot n.states Number of states (log(N) ages plus log(F) states).
#' @slot states The F-state coupling matrix (see \code{FLSAM.control}).
#' @slot components Estimated proportions of partial-survey components (if any).
#' @slot stock.wt,catch.wt,mat,m Estimated biological processes when the
#'   corresponding process models are switched on in the control object.
#' @slot nlogl Negative log-likelihood at the optimum.
#' @slot vcov Covariance matrix of the fixed-effect parameters.
#' @slot rescov Covariance matrix of all reported (derived) quantities.
#' @slot obscov List of observation covariance matrices for fleets with
#'   correlated observations.
#' @slot params \code{data.frame} with \code{name}, \code{value} and
#'   \code{std.dev} of all parameters, random effects and derived quantities.
#' @slot stock.n,harvest Estimated numbers and fishing mortality at age.
#' @slot residuals \code{data.frame} with \code{year}, \code{fleet}, \code{age},
#'   \code{log.obs}, \code{log.mdl} and the one-observation-ahead residual
#'   \code{std.res}.
#' @slot info Version and run information.
#' @slot fit The underlying \code{stockassessment} fit (class \code{sam}), or
#'   \code{NULL}. It gives access to all functions of the stockassessment
#'   package, see \code{\link{samfit}}.
#'
#' @seealso \code{\link{FLSAM}}, \code{\link{samfit}}, \code{\link{accessors}}
#' @rdname FLSAM-class
#' @aliases sam-class samOrNULL-class
#' @exportClass FLSAM
setClass("FLSAM",
         representation(
           "FLComp",
           control    = "FLSAM.control",
           nopar      = "integer",
           n.states   = "integer",
           states     = "matrix",
           components = "matrix",
           stock.wt   = "FLQuant",
           catch.wt   = "FLQuant",
           mat        = "FLQuant",
           m          = "FLQuant",
           nlogl      = "numeric",
           vcov       = "matrix",
           rescov     = "matrix",
           obscov     = "list",
           params     = "data.frame",
           stock.n    = "FLQuant",
           harvest    = "FLQuant",
           residuals  = "data.frame",
           info       = "matrix",
           fit        = "samOrNULL"),
         prototype=prototype(fit=NULL),
         validity=function(object){
           # Everything is fine
           return(TRUE)}
)

#' Run a SAM stock assessment from FLR objects
#'
#' Converts the FLR data and configuration into SAM format
#' (\code{\link{FLSAM2SAM}}, \code{\link{ctrl2conf}}), fits the model with
#' \code{\link[stockassessment]{sam.fit}} and converts the results back to an
#' \code{\link{FLSAM-class}} object (\code{\link{SAM2FLR}}).
#'
#' @param stcks An \code{FLStock}, or an \code{FLStocks} with elements
#'   \code{residual} and optionally \code{sum}.
#' @param tun An \code{FLIndices} object with the survey indices. The
#'   \code{index.var} slots are passed on as observation weights.
#' @param ctrl An \code{\link{FLSAM.control-class}} object.
#' @param catch.vars Optional \code{FLQuants} with weights for the catch
#'   observations (elements \code{residual} and optionally \code{sum}).
#' @param return.fit If \code{TRUE}, return the raw \code{sam} fit instead of
#'   an \code{FLSAM} object.
#' @param starting.values Optional earlier assessment whose estimates are used
#'   as starting values: an \code{FLSAM} object (with or without stored fit,
#'   e.g. last year's assessment or a retrospective peel) or a \code{sam} fit.
#'   Fixed-effect parameters are copied when their number matches the current
#'   configuration; the states (\code{logN}, \code{logF}) are matched by year,
#'   and years that are not in the earlier assessment get the value of the
#'   nearest year. Parameters that do not match keep their default starting
#'   value (with a message). See \code{\link{startingValues}}.
#' @param set.pars Optional \code{data.frame} with columns \code{par},
#'   \code{number} and \code{value} to fix parameters at a given value
#'   (on the SAM, mostly log, scale). \code{number} is the zero-based index
#'   within the parameter vector (i.e. the coupling number). \code{par} is either
#'   an FLSAM name (see \code{\link{setParsTable}}) or any fixed-effect
#'   parameter name of SAM (see \code{names(stockassessment::defpar(...))}),
#'   e.g. \code{"logSdLogTotalObs"} or \code{"logPhiSW"}.
#' @param keep.fit Store the \code{sam} fit in the \code{fit} slot (default
#'   \code{TRUE}); set to \code{FALSE} to save memory.
#' @param ... Further arguments passed to \code{\link[stockassessment]{sam.fit}},
#'   e.g. \code{newtonsteps}, \code{rel.tol}, \code{lower}, \code{upper},
#'   \code{rm.unidentified}, \code{penalizeSpline}, \code{silent}.
#' @return An \code{\link{FLSAM-class}} object (or a \code{sam} fit when
#'   \code{return.fit=TRUE}).
#' @seealso \code{\link{FLSAM.control}}, \code{\link{retro}}, \code{\link{looi}},
#'   \code{\link{samfit}}, \code{\link{setParsTable}}
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)
#' plot(sam)
#' # fix the observation variance of coupling 1 and the F correlation parameter
#' sp  <- data.frame(par=c("obs.vars","cor.F"), number=c(1,0), value=c(log(0.3),0.5))
#' sam.fixed <- FLSAM(NSH, NSH.tun, NSH.ctrl, set.pars=sp)
#' lr.test(sam.fixed, sam)
#' }
#' @export
FLSAM <-function(stcks,tun,ctrl,catch.vars=NULL,return.fit=F,starting.values=NULL,set.pars=NULL,keep.fit=TRUE,...){
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(is(stcks, "FLStock")) stcks <- FLStocks(residual=stcks)
  ctrl   <- .upgradeControl(ctrl)
  dataS  <- FLSAM2SAM(stcks,tun,ctrl@sumFleets,catch.vars,.recapture(ctrl))
  confS  <- ctrl2conf(ctrl,dataS)
  parS   <- defpar(dataS,confS)
  if(!is.na(ctrl@hockeyStickCurve[1]) & ctrl@srr==63)
    parS$rec_pars[3] <- log(ctrl@hockeyStickCurve[1])
  if(!is.null(starting.values))
    parS <- startingValues(parS,dataS,starting.values)

  if(!is.null(set.pars)){
    sp    <- .setPars(parS,set.pars)
    parS  <- sp$par
    fit   <- sam.fit(dataS,confS,parS,map=sp$map,sim.condRE=ctrl@simulate,...)
  } else {
    fit   <- sam.fit(dataS,confS,parS,sim.condRE=ctrl@simulate,...)
  }

  #- Check if you want variance and covariance info
  #- Get FLSAM object back
  if(return.fit){
    res <- fit
  } else {
    res <- try(SAM2FLR(fit,ctrl,keep.fit=keep.fit))
  }

  return(res)
}

#' Starting values from an earlier assessment
#'
#' Fills a SAM parameter list (as returned by
#' \code{\link[stockassessment]{defpar}}) with the estimates of an earlier
#' assessment, which may cover different years (e.g. one year less of data).
#' Used by \code{FLSAM(..., starting.values=)}.
#'
#' Fixed-effect parameters are copied when their length matches. The random
#' effects \code{logN} and \code{logF} are matched by year: years of the new
#' data that are not in the earlier assessment get the value of the nearest
#' earlier year (or the first year). Everything else keeps its default.
#'
#' @param par A SAM parameter list, e.g. \code{stockassessment::defpar(data, conf)}.
#' @param data The SAM data list of the new assessment (for its years).
#' @param start An \code{FLSAM} object (with or without stored fit) or a \code{sam} fit.
#' @param verbose Report the parameters that keep their default value.
#' @return The SAM parameter list with starting values.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' #- last year's assessment as starting values for this year's assessment
#' sam.last <- FLSAM(window(NSH,end=2010),
#'                   FLIndices(lapply(NSH.tun,window,end=2011)),NSH.ctrl)
#' sam.now  <- FLSAM(NSH,NSH.tun,NSH.ctrl,starting.values=sam.last)
#' }
#' @export
startingValues <- function(par,data,start,verbose=TRUE){
  #- Earlier estimates as list + their years
  if(is(start,"sam")){
    pl <- start$pl; yrsOld <- start$data$years; fromParams <- FALSE
  } else if(is(start,"FLSAM")){
    if(.hasSlot(start,"fit") && !is.null(start@fit)){
      pl <- start@fit$pl; yrsOld <- start@fit$data$years; fromParams <- FALSE
    } else {
      #- no stored fit (e.g. saved with an older version): use the params slot
      p  <- start@params
      pl <- lapply(split(p$value,factor(p$name,levels=unique(p$name))),identity)
      yrsOld <- start@range["minyear"]:start@range["maxyear"]; fromParams <- TRUE
    }
  } else {
    stop("starting.values must be an FLSAM object or a sam fit")
  }
  yrsNew <- data$years
  kept   <- character()
  for(n in intersect(names(par),names(pl))){
    new <- par[[n]]; old <- pl[[n]]
    if(length(new)==0) next
    if(n %in% c("logN","logF")){
      #- states x years: match by year
      nr <- nrow(new)
      if(is.null(dim(old))){
        if(length(old) != nr*length(yrsOld)){ kept <- c(kept,n); next }
        old <- matrix(old,nrow=nr)
      }
      if(nrow(old)!=nr || ncol(old)!=length(yrsOld)){ kept <- c(kept,n); next }
      idx <- pmax(findInterval(yrsNew,yrsOld),1)       # same year, else nearest earlier year
      par[[n]][] <- old[,idx,drop=FALSE]
    } else {
      #- the params slot holds some parameters twice (as parameter and as derived quantity)
      if(fromParams && length(old)==2*length(new)) old <- old[seq_along(new)]
      if(length(old)==length(new)){
        par[[n]][] <- old
      } else {
        kept <- c(kept,n)
      }
    }
  }
  if(verbose && length(kept)>0)
    message("starting.values: default starting values kept for ",paste(kept,collapse=", "))
  return(par)}

#' Parameters that can be fixed with set.pars
#'
#' Returns the translation table between the FLSAM names that can be used in
#' the \code{par} column of \code{set.pars} (see \code{\link{FLSAM}}) and the
#' SAM parameter names. Any other fixed-effect SAM parameter name can be used
#' directly.
#'
#' @return A character matrix with columns \code{FLSAM} and \code{SAM}.
#' @examples
#' setParsTable()
#' @export
setParsTable <- function(){
  matrix(c("logN.vars","logSdLogN",
           "logP.vars","logSdLogP",
           "catchabilities","logFpar",
           "power.law.exps","logQpow",
           "f.vars","logSdLogFsta",
           "obs.vars","logSdLogObs",
           "cor.obs","transfIRARdist",
           "cor.F","itrans_rho",
           "srr","rec_pars",
           "scalePars","logScale",
           "predVarObsLink","predVarObs",
           "XtraSd","logXtraSd"),ncol=2,byrow=T,dimnames=list(NULL,c("FLSAM","SAM")))
}

#- Translate set.pars into SAM names (column SAM added, unknown parameters dropped with a warning)
.setParsSAM <- function(par,set.pars){
  tbl       <- setParsTable()
  random    <- c("logN","logF","logP","logSW","logCW","logitMO","logNM","missing")
  set.pars$SAM <- ifelse(set.pars$par %in% tbl[,"FLSAM"],tbl[match(set.pars$par,tbl[,"FLSAM"]),"SAM"],as.character(set.pars$par))
  ok        <- set.pars$SAM %in% setdiff(names(par),random)
  ok[ok]    <- (set.pars$number[ok]+1) <= sapply(par[set.pars$SAM[ok]],length)
  if(any(!ok))
    warning(paste(set.pars$par[!ok],"cannot be set"))
  return(set.pars[ok,,drop=FALSE])}

#- Fix parameters: returns updated starting values and TMB map
.setPars <- function(par,set.pars){
  set.pars  <- .setParsSAM(par,set.pars)
  map       <- list()
  for(iPar in unique(set.pars$SAM)){
    idx           <- set.pars$number[set.pars$SAM==iPar]+1
    par[[iPar]][idx] <- set.pars$value[set.pars$SAM==iPar]
    mp            <- seq_along(par[[iPar]])
    mp[idx]       <- NA
    map[[iPar]]   <- factor(mp)
  }
  return(list(par=par,map=map))}

#' Internal fitting helpers
#'
#' Copies of internal \pkg{stockassessment} helpers used by
#' \code{sam.fitfast}, a stripped-down version of
#' \code{\link[stockassessment]{sam.fit}} without the variance calculation.
#'
#' @param parameters A SAM parameter list.
#' @param dat,data A SAM data list.
#' @param conf A SAM configuration list.
#' @param lower,upper Bounds.
#' @param sim.condRE Simulate conditional on random effects.
#' @param ... Passed to \code{MakeADFun}.
#' @keywords internal
#' @rdname sam.fitfast
#' @export
getLowerBounds<-function(parameters){
    list(sigmaObsParUS=rep(-10,length(parameters$sigmaObsParUS)))
}

#' @rdname sam.fitfast
#' @export
getUpperBounds<-function(parameters){
    list(sigmaObsParUS=rep(10,length(parameters$sigmaObsParUS)))
}

#' @rdname sam.fitfast
#' @export
clean.void.catches<-function(dat, conf){
  cfidx <- which(dat$fleetTypes==0)
  aidx <- unique(dat$aux[dat$aux[,2]%in%cfidx,3]-conf$minAge+1)
  faidx <- as.matrix(expand.grid(cfidx, aidx))
  faidx <- faidx[which(conf$keyLogFsta[faidx]== -1),,drop=FALSE]
  rmidx <- paste0(dat$aux[,2],"x",dat$aux[,3]-conf$minAge+1) %in%  paste0(faidx[,1],"x",faidx[,2])
  dat$aux <- dat$aux[!rmidx,]
  dat$logobs <- dat$logobs[!rmidx]
  dat$weight <- dat$weight[!rmidx]
  dat$nobs<-sum(!rmidx)
  dat$minAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=min))
  dat$maxAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=max))
  newyear<-min(as.numeric(dat$aux[,"year"])):max(as.numeric(dat$aux[,"year"]))
  newfleet<-min(as.numeric(dat$aux[,"fleet"])):max(as.numeric(dat$aux[,"fleet"]))
  mmfun<-function(f,y, ff){idx<-which(dat$aux[,"year"]==y & dat$aux[,"fleet"]==f); ifelse(length(idx)==0, NA, ff(idx)-1)}
  dat$idx1<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=min)
  dat$idx2<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=max)
  dat
}


#' @rdname sam.fitfast
#' @export
sam.fitfast <- function(data, conf, parameters, lower=getLowerBounds(parameters), upper=getUpperBounds(parameters), sim.condRE=TRUE, ...){
  data <- clean.void.catches(data,conf)
  tmball <- c(data, conf, simFlag=as.numeric(sim.condRE))
  if(is.null(tmball$resFlag)){tmball$resFlag <- 0}
  nmissing <- sum(is.na(data$logobs))
  parameters$missing <- numeric(nmissing)
  if(length(conf$keyVarLogP)>1){
    ran <- c("logN", "logF","logPS", "missing")
    obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",...)
  } else {
    ran <- c("logN", "logF", "missing")
    obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",...)
  }

  lower2<-rep(-Inf,length(obj$par))
  upper2<-rep(Inf,length(obj$par))
  for(nn in names(lower)) lower2[names(obj$par)==nn]=lower[[nn]]
  for(nn in names(upper)) upper2[names(obj$par)==nn]=upper[[nn]]

  opt <- try(nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=0, eval.max=1000, iter.max=500),lower=lower2,upper=upper2))
  rep <- obj$report()
#  sdrep <- TMB::sdreport(obj,opt$par)
#
#  # Last two states
#  idx <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
#  sdrep$estY <- sdrep$value[idx]
#  sdrep$covY <- sdrep$cov[idx,idx]
#
#  idx <- c(which(names(sdrep$value)=="beforeLastLogN"),which(names(sdrep$value)=="beforeLastLogF"))
#  sdrep$estYm1 <- sdrep$value[idx]
#  sdrep$covYm1 <- sdrep$cov[idx,idx]
#
#  pl <- as.list(sdrep,"Est")
#  plsd <- as.list(sdrep,"Std")
#
#  #sdrep$cov<-NULL # save memory
#
#  fit <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=conf, opt=opt, obj=obj, rep=rep, low=lower2, hig=upper2)

  if(!is(opt, "try-error")){
    sdrep <- list(par.fixed=length(opt$par),cov.fixed=matrix(NA),cov=matrix(NA))
    pl    <- c(lapply(split(opt$par,names(opt$par)),unname),list(logPS=rep$logP),list(logN=rep$logN),list(logF=log(rep$totF[-which(duplicated(conf$keyLogFsta[1,])),])))
  } else {
    sdrep <- NA
    pl    <- NA
  }
  fit <- list(sdrep=sdrep, pl=pl, plsd=NA, data=data, conf=conf, opt=opt, obj=obj, rep=rep, low=lower2, hig=upper2)
return(fit)}


# class
#' @rdname FLSAMs
#' @exportClass FLSAMs
setClass("FLSAMs", contains="FLlst")

# constructor
#' @rdname FLSAMs
setGeneric("FLSAMs", function(object, ...){
	standardGeneric("FLSAMs")
	}
)

#' @rdname FLSAMs
setMethod("FLSAMs", signature(object="ANY"), function(object, ...){
	lst1 <- list(...)
	nlst <- length(lst1)
	lst <- list()
	length(lst) <- nlst + 1
	lst[[1]] <- object
	lst[-1] <- lst1
	FLSAMs(lst)
})

#' @rdname FLSAMs
setMethod("FLSAMs", signature(object="FLSAM"), function(object, ...) {
    lst <- c(object, list(...))
    FLSAMs(lst)
})

#' @rdname FLSAMs
setMethod("FLSAMs", "missing", function(...){
	if(missing(...)){
		new("FLSAMs")
	} else {
		lst <- list(...)
		FLSAMs(lst)
	}
})

#' @rdname FLSAMs
setMethod("FLSAMs", "list", function(object){
        #Check list contents are FLSAM objects
        for(i in seq(object)) {
          if(!is(object[[i]],"FLSAM")) {
            stop(sprintf("Object %i in list is not an FLSAM",i)) }} 
        #Setup unique names
        list.names <- names(object)
        obj.names <- sapply(object,name)
        outnames <- if(is.null(list.names)) {obj.names} else {list.names}
        if(any(table(outnames) !=1) | any(outnames=="") |is.null(outnames)) outnames <- as.character(seq(object))
	new("FLSAMs", .Data=object,names=outnames)
})

#' @rdname FLSAMs
setMethod("FLSAMs", "FLSAMs", function(object){
	return(object)
}) # }}}
