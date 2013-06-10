setClass("FLSAM",
  	representation(
      "FLComp",
      control  = "FLSAM.control",
      nopar    = "integer",
      nlogl    = "numeric",
      vcov     = "matrix",
      params   = "data.frame",
      stock.n  = "FLQuant",
      harvest  = "FLQuant",
      residuals = "data.frame",
      info     = "matrix"),
  	prototype=prototype(),
  	validity=function(object){
                	# Everything is fine
                	return(TRUE)}
)


FLSAM <-function(stck,tun,ctrl,run.dir=tempdir(),batch.mode=FALSE,pin.sam=NULL) {
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  #General Setup
  inputSAM      <- new("FLSAMinput")
  inputSAM@stck <- stck
  inputSAM@tun  <- tun
  inputSAM@ctrl <- ctrl

  #- Check validity of objects
  if(any(c(validObject(stck),validObject(tun),validObject(ctrl),validObject(inputSAM))==F)) stop("Validity check failed")
  
  #Write output files
  FLR2SAM(stck,tun,ctrl,run.dir,pin.sam)

  #Run SAM
  if(is.null(pin.sam)==TRUE){
    rtn <- runSAM(ctrl, run.dir)
  } else {
    rtn <- runSAM(ctrl,run.dir,use.pin=TRUE)
  }
  if(rtn!=0) {
    if(batch.mode) {
      return(NULL)
    } else {
      stop(sprintf("An error occurred while running ADMB. Return code %s.",rtn))
    }
  }

  #Read the results back in
  res <- SAM2FLR(ctrl,run.dir)

  return(res)
}

#---------------------------------------------------
# We're ready! Run the executable
#---------------------------------------------------
runSAM <- function(ctrl,run.dir=tempdir(),use.pin=FALSE){
  admb.stem <- .get.admb.stem(ctrl) 
  admb.args <-  "-nr 2 -noinit -iprint 5"
  if(ctrl@nohess) {admb.args <- paste(admb.args,"-nohess")}
  if(use.pin) {
     admb.args <- paste(admb.args,"-usepin")
  } else {
     #Pin files are a contaigon if you don't want them!
     unlink(file.path(run.dir,sprintf("%s.pin",admb.stem)))
  }

  #If binary is specified, use it. Otherwise copy it from
  #the package distribution
  if(length(ctrl@sam.binary)!=0) {  
    admb.exec <- ctrl@sam.binary
    if(!file.exists(admb.exec)) {stop(sprintf("Cannot find specified sam executable (binary) file: %s",admb.exec))}
  } else if (.Platform$OS.type=="unix") {
    admb.exec <- file.path(system.file("bin", "linux", package="FLSAM",
                   mustWork=TRUE), admb.stem)
  } else if (.Platform$OS.type == "windows") {
    admb.exec <- file.path(system.file("bin", "windows", package="FLSAM", mustWork=TRUE),
      sprintf("%s.exe",admb.stem))
  } else {
    stop(sprintf("Platform type, %s, is not currently supported.",R.version$os))
  }
  file.copy(admb.exec, run.dir)

  #Run!
  cmd <- sprintf("./%s %s" ,basename(admb.exec),admb.args)
  olddir <- setwd(run.dir)
  rtn <- system(cmd)
  setwd(olddir)
  return(rtn)
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
  rng <- range(do.call(rbind,lapply(tun,function(x){cbind(range(x)["min"],range(x)["max"])})),na.rm=T)
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
