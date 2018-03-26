
setMethod("update", signature(object="FLSAM"),
  function(object,stck,tun,ctrl.new)
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

  if(object@control@scaleNoYears!=0){
    object@control@scaleNoYears <- length(object@control@scaleYears[which(object@control@scaleYears %in% strt:nd)])
    object@control@scaleYears   <- object@control@scaleYears[which(object@control@scaleYears %in% strt:nd)]
    object@control@scalePars    <- object@control@scalePars[ac(object@control@scaleYears),]
    object@control              <- update(object@control)
  }

  #-Update FLSAM.control stck by copying directly from the oject. Modifications
  # may be necessary for revised year ranges 
  ctrl@range["maxyear"] <- max(sapply(tun,range)["maxyear",],stck@range["maxyear"])
  ctrl@range["minyear"] <- min(sapply(tun,range)["minyear",],stck@range["minyear"])

  #-Run the assessment, "by hand" taking advantage of the usepin argument
  res <- FLSAM(stck,tun,ctrl.new,pin.sam=object)

  return(res)
  }
)



.update.params <- function(ctrl.old,ctrl.new,iSlot,parameters,sam,convDF){

  #- Select parameter conversion
  conv      <- convDF[which(convDF$FLRname == iSlot),]

  #- Check if no gaps exist in parameter values
  parVals   <- sort(na.omit(unique(c(slot(ctrl.new,iSlot)))))
  if(any(diff(parVals)>1)) stop("Run 'update' on control object first")
  if(iSlot != "srr"){
    #- Check if parameter binding has changed
    for(i in 0:(length(sort(na.omit(unique(c(slot(ctrl.new,iSlot))))))-1))
      parameters[[ac(conv$SAMname)]][i+1] <- ifelse(identical(which(slot(ctrl.old,iSlot) == i,arr.ind=T),which(slot(ctrl.new,iSlot) == i,arr.ind=T)),
                                                subset(params(sam),name==conv$SAMname)$value[i+1],
                                                  conv$default)
  }
  return(parameters)}


#- take default parameter values and overwrite with previous run converged parameter estimates
updateStart <- function(parameters,startVals){
  parametersUp <- parameters
  if("missing" %in% names(startVals))
    startVals <- startVals[-grep("missing",names(startVals))]
  for(iName in names(startVals)){
    if(identical(dim(parametersUp[[iName]]),dim(startVals[[iName]]))){
      parametersUp[[iName]] <- startVals[[iName]]
    } else {
      if(iName %in% c("logF","logN","logPS")){
        if((dim(parametersUp[[iName]])[2]-dim(startVals[[iName]])[2])==-1)
          parametersUp[[iName]] <- startVals[[iName]][,-ncol(startVals[[iName]])]
        if((dim(parametersUp[[iName]])[2]-dim(startVals[[iName]])[2])==1)
          parametersUp[[iName]] <- cbind(startVals[[iName]],startVals[[iName]][,ncol(startVals[[iName]])])
      }
      warnings(paste("Could not update according to starting value",iName,"as number of parameters specified does not match"))
    }
  }
  
  for(iName in names(startVals)){
    if(length(parametersUp[[iName]][which(parametersUp[[iName]] < -5)])>0)
      parametersUp[[iName]][which(parametersUp[[iName]] < -5)] <- parameters[[iName]][which(parametersUp[[iName]] < -5)]
  }

  return(parametersUp)}
  
#- Create defpar object from sam
FLSAM2defpar <- function(sam){
  ret<-list()
  ret$logFpar=numeric(max(sam@control@catchabilities,na.rm=T)+1)-5
  ret$logQpow=numeric(max(sam@control@power.law.exps,na.rm=T)+1)
  ret$logSdLogFsta=numeric(max(sam@control@f.vars)+1)-.7
  ret$logSdLogN=numeric(max(sam@control@logN.vars,na.rm=T)+1)-.35
  ret$logSdLogP=if(length(sam@control@logP.vars)>0){numeric(max(sam@control@logP.vars,na.rm=T)+1)-.7} else { numeric(0)}
  ret$logSdLogObs=numeric(max(sam@control@obs.vars,na.rm=T)+1)-.35
  ret$logSdLogTotalObs=numeric(sum(sam@control@cor.obs.Flag %in% c("ALN")))
  ret$transfIRARdist=if(all(is.na(sam@control@cor.obs)))numeric(0) else numeric(max(sam@control@cor.obs,na.rm=TRUE)+1)+0.05
  nbyfleet = (sam@control@cor.obs.Flag=="US")*(sam@control@range["max"]-sam@control@range["min"]+1-(sam@control@cor.obs.Flag=="ALN"))
  ret$sigmaObsParUS=numeric(sum(nbyfleet*(nbyfleet-1)/2,na.rm=TRUE))
  ret$rec_loga=if(sam@control@srr==0){numeric(0)}else{numeric(1)}
  ret$rec_logb=if(sam@control@srr==0){numeric(0)}else{numeric(1)}
  ret$itrans_rho=unlist(lapply(as.list(sam@control@cor.F),function(x){if(x==0){ ret <- numeric()} else { ret <- numeric(1)+.5}; return(ret)}))
  ret$rhop = if(length(sam@control@logP.vars)>0){0.5}else{numeric(0)}
  ret$logScale=if(sam@control@scaleNoYears==0){numeric(0)}else{numeric(max(sam@control@scalePars,na.rm=T)+1)}
  ret$logitReleaseSurvival=if(any(sam@control@fleets==5)){stop("Fleet type 5 FLSAM2parameters not implemented")
                           }else{numeric(0)}
  ret$logitRecapturePhi=if(any(sam@control@fleets==5)){stop("Fleet type 5 FLSAM2parameters not implemented")
                        }else{numeric(0)}
  ret$logAlphaSCB = if(length(sam@control@logP.vars)>0){sam@params[grep("logAlphaSCB",sam@params$name),"value"]}else{numeric(0)}
  ret$logF=matrix(0, nrow=max(sam@control@states,na.rm=T)+1,ncol=sam@control@range["maxyear"]-sam@control@range["minyear"]+1)
  ret$logN=matrix(0, nrow=sam@control@range["max"]-sam@control@range["min"]+1, ncol=sam@control@range["maxyear"]-sam@control@range["minyear"]+1)
  if(any(sam@control@fleets==6)){
    idxPart <- which(sam@control@fleets==6)
    ret$logPS=matrix(0, nrow=length(sam@control@logP.vars), ncol=ncol(sam@components))
  } else {
    ret$logPS=numeric()
  }
}

#- Move all estimated parameters from sam to a par object
FLSAM2par <- function(sam){
  ret               <- list()
  parNames          <- unique(sam@params$name)
  if("logFpar" %in% parNames)
    ret$logFpar       <- subset(sam@params,name=="logFpar")$value
  if("logQpow" %in% parNames)
    ret$logQpow       <- subset(sam@params,name=="logQpow")$value
  if("logSdLogFsta" %in% parNames)
    ret$logSdLogFsta  <- subset(sam@params,name=="logSdLogFsta")$value
  if("logSdLogN" %in% parNames)
    ret$logSdLogN     <- subset(sam@params,name=="logSdLogN")$value
  if("logSdLogP" %in% parNames)
    ret$logSdLogP     <- subset(sam@params,name=="logSdLogP")$value
  if("logSdLogObs" %in% parNames)
    ret$logSdLogObs   <- subset(sam@params,name=="logSdLogObs")$value
  if("logSdLogTotalObs" %in% parNames)
    ret$logSdLogTotalObs  <- subset(sam@params,name=="logSdLogTotalObs")$value
  if("transfIRARdist" %in% parNames)
    ret$transfIRARdist    <- subset(sam@params,name=="transfIRARdist")$value
  if("rec_loga" %in% parNames)
    ret$rec_loga      <- subset(sam@params,name=="rec_loga")$value
  if("rec_logb" %in% parNames)
    ret$rec_logb      <- subset(sam@params,name=="rec_logb")$value
  if("itrans_rho" %in% parNames)
    ret$itrans_rho    <- subset(sam@params,name=="itrans_rho")$value
  if("rhop" %in% parNames)
    ret$rhop          <- subset(sam@params,name=="rhop")$value
  if("logScale" %in% parNames)
    ret$logScale      <- subset(sam@params,name=="logScale")$value
  if("logAlphaSCB" %in% parNames)
    ret$logAlphaSCB   <- sam@params[grep("logAlphaSCB",sam@params$name),"value"]
  if("logF" %in% parNames)
    ret$logF          <- matrix(subset(sam@params,name=="logF")$value, nrow=max(sam@control@states,na.rm=T)+1,ncol=sam@control@range["maxyear"]-sam@control@range["minyear"]+1)
  if("logN" %in% parNames)
    ret$logN          <- matrix(subset(sam@params,name=="logN")$value, nrow=sam@control@range["max"]-sam@control@range["min"]+1, ncol=sam@control@range["maxyear"]-sam@control@range["minyear"]+1)
  if("logPS" %in% parNames)
    ret$logPS         <- matrix(subset(sam@params,name=="logPS")$value, nrow=length(sam@control@logP.vars), ncol=ncol(sam@components))

return(ret)}