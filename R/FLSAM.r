setClass("FLSAM",
         representation(
           "FLComp",
           control    = "FLSAM.control",
           nohess     = "logical",
           nopar      = "integer",
           n.states   = "integer",
           states     = "matrix",
           nlogl      = "numeric",
           vcov       = "matrix",
           rescov     = "matrix",
           obscov     = "list",
           params     = "data.frame",
           stock.n    = "FLQuant",
           harvest    = "FLQuant",
           residuals  = "data.frame",
           info       = "matrix"),
         prototype=prototype(),
         validity=function(object){
           # Everything is fine
           return(TRUE)}
)



FLSAM <-function(stck,tun,ctrl,batch.mode=FALSE,pin.sam=NULL,map=NULL,bounds=NULL,return.fit=FALSE){
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  #General Setup
  inputSAM      <- new("FLSAMinput")
  inputSAM@stck <- stck
  inputSAM@tun  <- tun
  inputSAM@ctrl <- ctrl

  #- Check validity of objects
  #if(any(c(validObject(stck),validObject(tun),validObject(ctrl),validObject(inputSAM))==F)) stop("Validity check failed")
  
  #Write output files
  fls2run     <- FLR2SAM(stck,tun,ctrl,pin.sam,map)
  data        <- fls2run[["data"]]
  parameters  <- fls2run[["parameters"]]
  map         <- fls2run[["map"]]
  
  #- Get the data in place to run the simulation
  nmissing    <- sum(is.na(data$logobs))
  parameters$missing <- numeric(nmissing)

  tmball      <- c(data)
  ran         <- c("logN", "logF", "missing")
  tmb.stem    <- .get.stem(ctrl)

  #- Find the location of the dll
  if(length(ctrl@sam.binary)!=0) {
    tmb.exec <- ctrl@sam.binary
    if(!file.exists(tmb.exec)) {stop(sprintf("Cannot find specified sam executable (binary) file: %s",tmb.exec))}
  } else if (.Platform$OS.type=="unix") {
    tmb.exec <- file.path(system.file("bin", "linux", package="FLSAM",
                 mustWork=TRUE), tmb.stem)
  } else if (.Platform$OS.type == "windows") {
    tmb.exec <- file.path(system.file("bin", "windows", package="FLSAM", mustWork=TRUE),
    sprintf(tmb.stem))
  } else {
    stop(sprintf("Platform type, %s, is not currently supported.",R.version$os))
  }
  
  #- Load and run the model
  dyn.load(dynlib(tmb.exec))
  if(is.null(map))
    obj     <- MakeADFun(tmball, parameters, random=ran, DLL="sam")
  if(!is.null(map))
    obj     <- MakeADFun(tmball, parameters, random=ran, DLL="sam",map=map)

  lower2    <- rep(-Inf,length(obj$par))
  upper2    <- rep(Inf,length(obj$par))
  if(!is.null(bounds)){
    for(nn in rownames(bounds)) lower2[names(obj$par)==nn]=bounds[nn,"lower"]
    for(nn in rownames(bounds)) upper2[names(obj$par)==nn]=bounds[nn,"upper"]
  }

  #- Fit the model
  opt       <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower2,upper=upper2)
  if(opt$convergence==0) print("Assessment converged")
  if(opt$convergence!=0) warning("Assessment did not converge")
  
  # Take a few extra newton steps
  for(i in seq_len(3)) {
    g <- as.numeric( obj$gr(opt$par) )
    h <- optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  }

  #- Check if you want variance and covariance info
  if(ctrl@nohess){
    rep           <- obj$report()
    fit           <- list(data=data,conf=parameters,opt=opt,obj=obj,rep=rep)
  }
  if(!ctrl@nohess){
    rep           <- obj$report()
    sdrep         <- sdreport(obj,opt$par)

    # Last two states
    idx           <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
    sdrep$estY    <- sdrep$value[idx]
    sdrep$covY    <- sdrep$cov[idx,idx]

    idx           <- c(which(names(sdrep$value)=="beforeLastLogN"),which(names(sdrep$value)=="beforeLastLogF"))
    sdrep$estYm1  <- sdrep$value[idx]
    sdrep$covYm1  <- sdrep$cov[idx,idx]
    
    pl            <- as.list(sdrep,"Est")
    plsd          <- as.list(sdrep,"Std")
    fit           <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=parameters, opt=opt, obj=obj, rep=rep)
  }

  #- Get FLSAM object back
  if(return.fit){
    res <- fit
  } else {
    res <- SAM2FLR(fit,ctrl)
  }

  return(res)
}


#---------------------------------------------------
# Validity checks
#---------------------------------------------------
setClass("FLSAMinput",
    representation(
      "list",
      stck     = "FLStock",
      ctrl     = "FLSAM.control",
      tun      = "FLIndices"),
  	prototype=prototype(),
  	validity=function(object){
                	# Everything is fine
                	return(TRUE)}
)

setValidity("FLSAMinput",
  function(object){

  #- Pull structure apart
  stck <- object@stck
  ctrl <- object@ctrl
  tun  <- object@tun
  
  #- Check ranges between stck and ctrl
  idxstock <- grep("maxyear",names(range(stck))); idxctrl <- grep("maxyear",names(ctrl@range))
  if(any(apply(cbind(range(stck)[-idxstock],ctrl@range[-idxctrl]),1,duplicated)[2,]==F)) stop("Range in stock object different from range in control object")
  if(ctrl@range[idxctrl] > (range(stck)[idxstock]+1)) stop("Year range in control larger than in stock+1. Only 1 year ahead tuning series are implemented")

  #- Check if number of fleets specified matches stck and tun
  if(length(ctrl@fleets) != (length(dimnames(stck@catch)$area)+length(tun))) stop("Number of fleets specified does not match stock and tun object")
  if(all(ctrl@fleets %in% c(0,3))) stop("Assessment cannot (currently) be run with only catch and SSB index, other survey fleet needs to be specified too.
                                         Contact the developers of SAM and FLSAM if you require this addition")
  
  #- Check if age dimensions match in stck ctrl and tun
  rng <- abs(range(do.call(rbind,lapply(tun,function(x){cbind(range(x)["min"],range(x)["max"])})),na.rm=T))
  if(rng[1] < range(stck)["min"] | rng[2] > range(stck)["max"]) stop("Age range in tuning series outside age range stock")
  
  #- Check if year dimensions match in stck ctrl and tun
  rng <- range(do.call(rbind,lapply(tun,function(x){cbind(range(x)["minyear"],range(x)["maxyear"])})),na.rm=T)
  if(rng[2] > (range(stck)["maxyear"]+1)) stop("Year range in tuning series larger than in stock+1. Only 1 year ahead tuning series are implemented")
  if(rng[1] < range(stck)["minyear"]) stop("Year range in tuning series smaller than in stock")

  #- Check if iter dimension is not larger than 1
  if(dims(stck)$iter > 1 | any(unlist(lapply(tun,function(x)return(dims(x)$iter)))>1)) stop("Only 1 iteration is allowed")
  }
)

# class
setClass("FLSAMs", contains="FLlst")

# constructor
setGeneric("FLSAMs", function(object, ...){
	standardGeneric("FLSAMs")
	}
)

setMethod("FLSAMs", signature(object="ANY"), function(object, ...){
	lst1 <- list(...)
	nlst <- length(lst1)
	lst <- list()
	length(lst) <- nlst + 1
	lst[[1]] <- object
	lst[-1] <- lst1
	FLSAMs(lst)
})

setMethod("FLSAMs", signature(object="FLSAM"), function(object, ...) {
    lst <- c(object, list(...))
    FLSAMs(lst)
})

setMethod("FLSAMs", "missing", function(...){
	if(missing(...)){
		new("FLSAMs")
	} else {
		lst <- list(...)
		FLSAMs(lst)
	}
})

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

setMethod("FLSAMs", "FLSAMs", function(object){
	return(object)
}) # }}}
