
setMethod("update", signature(object="FLSAM"),
  function(object,stck,tun,run.dir=tempdir(),batch.mode=FALSE)
  {

  #- Check validity of objects
  ctrl <- object@control
  if(!validObject(object))    stop("FLSAM object invalid")
  if(!validObject(stck))      stop("FLStock object invalid")
  if(!validObject(tun))       stop("FLIndices object invalid")


#  inputSAM      <- new("FLSAMinput")
#  inputSAM@stck <- stck
#  inputSAM@tun  <- tun
#  inputSAM@ctrl <- ctrl
#  if(!all(c(validObject(stck),validObject(tun),validObject(object),validObject(inputSAM)))) {
#        stop("Validity checks on input objects failed") }


  #-Test if ranges of FLStock are smaller than of FLSAM
  if(range(stck)["minyear"] < range(object)["minyear"] |
     range(stck)["maxyear"] > range(object)["maxyear"])
      stop("Year ranges of FLStock exceed those of the FLSAM stck,
      only previously estimated parameters can be used")

  #-Define starting and final year
  strt  <- range(stck)["minyear"]
  nd    <- range(stck)["maxyear"]

  #-Check for ranges of FLIndices stck
  if(all(unlist(lapply(tun,function(x){return(range(x)["maxyear"])})) > nd))
    warning("All indices longer than maxyear in FLStock,
    is the FLIndices stck trimmed to the right years?")
  if(any(unlist(lapply(tun,function(x){return(range(x)["maxyear"])})) > nd)){
    params2pin(object,strt,nd+1,save.dir=run.dir)
  } else {
    params2pin(object,strt,nd,  save.dir=run.dir)
  }
  
  #-Update FLSAM.control stck by copying directly from the oject. Modifications
  # may be necessary for revised year ranges 
  ctrl       <- object@control
  ctrl@range["maxyear"] <- max(sapply(tun,range)["maxyear",],stck@range["maxyear"])
  ctrl@range["minyear"] <- min(sapply(tun,range)["minyear",],stck@range["minyear"])
  
  #-Run the assessment, "by hand" taking advantage of the usepin argument
  FLR2SAM(stck,tun,ctrl,run.dir)
  rtn <- runSAM(ctrl, run.dir, use.pin=TRUE)
  if(rtn!=0) {
    if(batch.mode) {
      return(NULL)
    } else {
      stop(sprintf("An error occurred while running ADMB. Return code %s.",rtn))
    }
  }
  res <- SAM2FLR(ctrl,run.dir)

  return(res)
  }
)


params2pin <- function(object,start=NULL,end=NULL,save.dir=tempdir()){

  #---------------------------------------------------
  # Book keeping
  #---------------------------------------------------
  if(!is(object,"FLSAM")) stop("Supplied object must be an FLSAM")
  run.time  <- Sys.time()
  strt      <- ifelse(is.null(start)==T,range(object)["minyear"],start)
  nd        <- ifelse(is.null(end)==T  ,range(object)["maxyear"],end)

  #---------------------------------------------------
  # Truncate data if necessary
  #---------------------------------------------------

  #-Extract stock.n and harvest values from SAM object
  n         <- log(window(object@stock.n,strt,nd))
  f         <- log(window(object@harvest[!duplicated(object@control@states[names(which(object@control@fleets==0)),]),],strt,nd))

  #-Define how many parameters are estimated
  n.states  <- length(unique(object@control@states[names(which(object@control@fleets==0)),]))
  ages      <- dims(object)$age

  idxNs     <- c(mapply(seq,from=seq(1,length(n)+length(f),ages+n.states),
                            to  =seq(1,length(n)+length(f),ages+n.states)+ages-1,
                            by  =1))
  idxFs     <- c(mapply(seq,from=seq(1,length(n)+length(f),ages+n.states)+ages,
                            to  =seq(1,length(n)+length(f),ages+n.states)+n.states+ages-1,
                            by  =1))
                            
  #-Create vector with combined n and f values
  u         <- numeric(length(n)+length(f))
  u[idxNs]  <- n
  u[idxFs]  <- f
  U         <- data.frame(name="U",value=u,std.dev=NA)

  #-Merge U's with other parameters (even if truncated)
  par.list  <-  c("logFpar","logSdLogFsta","logSdLogN","logSdLogObs","rec_loga",
                  "rec_logb","rho","logScaleSSB","logPowSSB","logSdSSB","U")
  pars      <- object@params
  uidx      <- which(pars$name == "U")
  parsTop   <- pars[1:(uidx[1]-1),]
  pars      <- rbind(parsTop,U,if(nrow(pars)>max(uidx)){pars[(rev(uidx)[1]+1):nrow(pars),]})

  #---------------------------------------------------
  # Create pin file
  #---------------------------------------------------
  pin.file  <- file.path(save.dir,"sam.pin")
  unlink(pin.file)    #Get rid of anything that is already there
  for(iPar in par.list){
    #Get pars to write
    par.vals <- pars[pars$name==iPar,"value"] 
    if(length(par.vals)==0) par.vals <- 0 
    #Write to file
    cat("#",iPar,": \n ",file=pin.file,append=TRUE,sep="")
    cat(par.vals," \n ",file=pin.file,append=TRUE)
  }
invisible(NULL)}





