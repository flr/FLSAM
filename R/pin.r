
setMethod("update", signature(object="FLStock"),
  function(object,y,z,run.dir=tempdir())
  {

  if(!validObject(object)) stop("FLStock object invalid")
  if(!validObject(y))      stop("FLSAM object invalid")
  if(!validObject(z))      stop("FLIndices object invalid")

  #-Test if ranges of FLStock are smaller than of FLSAM
  if(range(object)["minyear"] < range(y)["minyear"] |
     range(object)["maxyear"] > range(y)["maxyear"])
      stop("Year ranges of FLStock exceed those of the FLSAM object,
      only previously estimated parameters can be used")

  #-Define starting and final year
  strt  <- range(object)["minyear"]
  nd    <- range(object)["maxyear"]

  #-Check for ranges of FLIndices object
  if(all(unlist(lapply(z,function(x){return(range(x)["maxyear"])})) > nd))
    warning("All indices longer than maxyear in FLStock,
    is the FLIndices object trimmed to the right years?")
  if(any(unlist(lapply(z,function(x){return(range(x)["maxyear"])})) > nd)){
    params2pin(  y,strt,nd+1,save.dir=run.dir)
  } else {
      params2pin(y,strt,nd,  save.dir=run.dir)
    }
  
  #-Update FLSAM.control object
  ctrl        <- FLSAM.control(object,z)
  for(iSlot in c("states","logN.vars","catchabilities","power.law.exps",
                 "f.vars","obs.vars","srr","cor.F","nohess","timeout"))
    slot(ctrl,iSlot) <- slot(y@control,iSlot)

  
  #-Update the assessment
  ynew        <- FLSAM(object,z,ctrl,run.dir,batch.mode=T,pin.dir=run.dir)

  return(ynew)
  }
)




params2pin <- function(FLSAM,start=NULL,end=NULL,save.dir=tempdir()){

  #---------------------------------------------------
  # Book keeping
  #---------------------------------------------------
  run.time  <- Sys.time()
  strt      <- ifelse(is.null(start)==T,range(FLSAM)["minyear"],start)
  nd        <- ifelse(is.null(end)==T  ,range(FLSAM)["maxyear"],end)

  .file.header <- function(fname,ftime) {
     cat("# Auto generated file\n", file=fname)
     cat(sprintf("# Datetime : %s\n\n",ftime),file=fname,append=TRUE)
  }
  
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





