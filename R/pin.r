
setMethod("update", signature(object="FLSAM"),
  function(object,stck,tun,run.dir=tempdir(),batch.mode=FALSE)
  {

  #- Check validity of objects
  ctrl <- object@control
  inputSAM      <- new("FLSAMinput")
  inputSAM@stck <- stck
  inputSAM@tun  <- tun
  inputSAM@ctrl <- ctrl
  if(!all(c(validObject(stck),validObject(tun),validObject(object),validObject(inputSAM)))) {
        stop("Validity checks on input objects failed") }

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
  
  #-Update FLSAM.control stck
  #ctrl        <- FLSAM.control(stck,tun)
  #for(iSlot in c("states","logN.vars","catchabilities","power.law.exps",
  #               "f.vars","obs.vars","srr","cor.F","nohess","timeout"))
  # {slot(ctrl,iSlot) <- slot(object@control,iSlot) }

  
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


params2pin <- function(FLSAM,start=NULL,end=NULL,save.dir=tempdir()){

  #---------------------------------------------------
  # Book keeping
  #---------------------------------------------------
  run.time  <- Sys.time()
  strt      <- ifelse(is.null(start)==T,range(FLSAM)["minyear"],start)
  nd        <- ifelse(is.null(end)==T  ,range(FLSAM)["maxyear"],end)

  #---------------------------------------------------
  # Truncate data if necessary
  #---------------------------------------------------

  #-Extract stock.n and harvest values from SAM object
  n         <- log(window(FLSAM@stock.n,strt,nd))
  f         <- log(window(FLSAM@harvest[!duplicated(FLSAM@control@states[names(which(FLSAM@control@fleets==0)),]),],strt,nd))

  #-Define how many parameters are estimated
  n.states  <- length(unique(FLSAM@control@states[names(which(FLSAM@control@fleets==0)),]))
  ages      <- dims(FLSAM)$age

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
  pars      <- FLSAM@params
  uidx      <- which(pars$name == "U")
  parsTop   <- pars[1:(uidx[1]-1),]
  pars      <- rbind(parsTop,U,if(nrow(pars)>max(uidx)){pars[(rev(uidx)[1]+1):nrow(pars),]})

  parnames  <- ac(unique(pars$name))
  parvals   <- pars$value

  #---------------------------------------------------
  # Create pin file
  #---------------------------------------------------
  pin.file  <- file.path(save.dir,"sam.pin")
  for(iPar in parnames){
    cat("#",iPar,": \n ",file=pin.file,append=TRUE,sep="")
    cat(parvals[which(pars$name %in% iPar)]," \n ",file=pin.file,append=TRUE)
  }
invisible(NULL)}





