setClass("FLSAM",
  	representation(
      "FLComp",
      control  = "FLSAM.control",
      nopar    = "integer",
      nlogl    = "numeric",
      params   = "data.frame",
      stock.n  = "FLQuant",
      harvest  = "FLQuant",
      residuals = "data.frame"),
  	prototype=prototype(),
  	validity=function(object){
                	# Everything is fine
                	return(TRUE)}
)


FLSAM <-function(stck,tun,ctrl,run.dir=tempdir(),batch.mode=FALSE) {
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  #General Setup
  admb.stem <- "sam" 
  
  inputSAM      <- new("FLSAMinput")
  inputSAM@stck <- stck
  inputSAM@tun  <- tun
  inputSAM@ctrl <- ctrl

  #- Check validity of objects
  if(any(c(validObject(stck),validObject(tun),validObject(ctrl),validObject(inputSAM))==F)) stop("Validity check failed")
  
  #Write output files
  FLR2SAM(stck,tun,ctrl,run.dir)

  #---------------------------------------------------
  # We're ready! Run the executable
  #---------------------------------------------------
  admb.args <-  "-nr 2 -noinit -iprint 1"
  #Platform specific issues
  if (.Platform$OS.type=="unix") {
    admb.exec <- file.path(system.file("bin", "linux", package="FLSAM",
                   mustWork=TRUE), admb.stem)
    file.copy(admb.exec, run.dir)
  } else if (.Platform$OS.type == "windows") {
    admb.exec <- file.path(system.file("bin", "windows", package="FLSAM", mustWork=TRUE),
      sprintf("%s.exe",admb.stem))
    file.copy(admb.exec, run.dir)
  } else {
    stop(sprintf("Platform type, %s, is not currently supported.",R.version$os))
  }

  #Run!
  cmd <- sprintf("./%s -nr 2 -noinit -iprint 1" ,basename(admb.exec))
  olddir <- setwd(run.dir)
  rtn <- system(cmd)
  setwd(olddir)
  if(rtn!=0) {
    if(batch.mode) {
      return(NULL)
    } else {
      stop(sprintf("An error occurred while running ADMB. Return code %s.",rtn))
    }
  }

  #---------------------------------------------------
  # Now read the results from the assessment
  #---------------------------------------------------
  res <- SAM2FLR(ctrl,run.dir)

  return(res)
}

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
  if(any(apply(cbind(range(stck),ctrl@range),1,duplicated)[2,]==F)) stop("Range in stock object different from range in control object")

  #- Check if number of fleets specified matches stck and tun
  if(length(ctrl@fleets) != (length(dimnames(stck@catch)$area)+length(tun))) stop("Number of fleets specified does not match stock and tun object")
  
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
