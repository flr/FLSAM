
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


createPin <- function(sam,ctrl.new,...){

  if(!is(sam,"FLSAM")) stop("Supplied object must be an FLSAM")

  #- Generate empty list of parameters
  parameters <- list(logFpar = numeric(), logQpow = numeric(), logSdLogFsta = numeric(), logSdLogN = numeric(),
                     logSdLogObs = numeric(), logSdLogTotalObs = numeric(), transfIRARdist = numeric(), sigmaObsParUS = numeric(),
                     rec_loga = numeric(), rec_logb = numeric(), itrans_rho = numeric(), logScale = numeric(), logitReleaseSurvival = numeric(),
                     logitRecapturePhi = numeric(), logF = numeric(), logN = numeric())

  #- Create conversion between SAM and FLR for parameters
  convDF     <- data.frame(SAMname=names(parameters),
                           FLRname=c("catchabilities","power.law.exps","f.vars",
                                     "logN.vars","obs.vars",NA,"cor.obs",NA,"srr",
                                     "srr",NA,"scalePars",NA,NA,NA,NA),
                           default=c(-5,1,-0.7,-0.35,-0.35,0,0.05,0,0,0,0.5,0,0,0,0,0),
                           stringsAsFactors=F)

  #- Create new template
  if(length(unique(na.omit(c(ctrl.new@catchabilities)))) > 0)
    parameters$logFpar        <- rep(-5,length(unique(na.omit(c(ctrl.new@catchabilities)))))
  if(length(unique(na.omit(c(ctrl.new@power.law.exps)))) > 0)
    parameters$logQpow        <- rep(1,unique(na.omit(c(ctrl.new@power.law.exps))))
  if(length(unique(na.omit(c(ctrl.new@f.vars)))) > 0)
    parameters$logSdLogFsta   <- rep(-0.7,length(unique(na.omit(c(ctrl.new@f.vars)))))
  if(length(unique(na.omit(c(ctrl.new@logN.vars)))) > 0)
    parameters$logSdLogN      <- rep(-0.35,length(unique(na.omit(c(ctrl.new@logN.vars)))))
  if(length(unique(na.omit(c(ctrl.new@obs.vars)))) > 0)
    parameters$logSdLogObs    <- rep(-0.35,length(unique(na.omit(c(ctrl.new@obs.vars)))))
  if(any(ctrl.new@likFlag %in% "ALN"))
    parameters$logSdLogTotalObs <- 0
  if(length(unique(na.omit(c(ctrl.new@cor.obs))))>0)
    parameters$transfIRARdist <- rep(0.05,length(unique(na.omit(c(ctrl.new@cor.obs)))))
  if(any(ctrl.new@cor.obs.Flag %in% "US")){
    idx                       <- which(ctrl.new@cor.obs.Flag == "US")
    len                       <- length(na.omit(c((ctrl.new@cor.obs.Flag[idx]))))
    parameters$sigmaObsParUS  <- rep(0,len*(len-1)/2)
  }
  if(ctrl.new@srr > 0){
    parameters$rec_loga       <- 1
    parameters$rec_logb       <- 1
  }
  if(ctrl.new@cor.F)
    parameters$itrans_rho     <- 0.5
  if(length(unique(na.omit(c(ctrl.new@scalePars))))>0)
    parameters$logScale       <- rep(0,length(unique(na.omit(c(ctrl.new@scalePars)))))

  #- Fill template with old values
  for(iPar in names(parameters)){
    iSlot <- convDF[which(convDF$SAMname==iPar),"FLRname"]
    if(is.na(iSlot)==F){
      if(any(is.na(slot(ctrl.new,iSlot))==F)){
        parameters <- .update.params(sam@control,ctrl.new,iSlot,parameters,sam,convDF)
      }
    }
  }

  #- Fill template with old Fs and Ns
  strt      <- ctrl.new@range["minyear"]
  nd        <- ctrl.new@range["maxyear"]
  strtage   <- ctrl.new@range["min"]
  ndage     <- ctrl.new@range["max"]

  #- Check if year range has changed
  if(strt != range(sam)["minyear"] | nd != range(sam)["maxyear"]){
    idx     <- names(which(ctrl.new@fleets == 0))
    if(strtage == range(sam)["min"] & ndage == range(sam)["max"]){
      parameters$logF <- log(window(sam@harvest,start=strt,end=nd)[!duplicated(sam@control@states[idx,]),,drop=T])
      parameters$logN <- log(window(sam@stock.n,start=strt,end=nd)[                                     ,,drop=T])
    } else {
        idxpar        <- ctrl.new@states[idx,!duplicated(ctrl.new@states[idx,])]
      if(length(unique(ctrl.new@states[idx,]))>length(unique(sam@control@states[idx,])))
        idxpar        <- c(idxpar,rep(idxpar[length(idxpar)],length(unique(ctrl.new@states[idx,]))-length(unique(sam@control@states[idx,]))))
      parameters$logF <- log(window(trim(sam@harvest,age=strtage:ndage),start=strt,end=nd)[idxpar+1,,drop=T])
      parameters$logN <- log(window(trim(sam@stock.n,age=strtage:ndage),start=strt,end=nd)[        ,,drop=T])
      parameters$logF[which(is.na(parameters$logF))] <- 0
      parameters$logN[which(is.na(parameters$logN))] <- 0
    }
  } else {
      idx             <- names(which(ctrl.new@fleets == 0))
      parameters$logF <- log(sam@harvest[!duplicated(sam@control@states[idx,]),,drop=T])
      parameters$logN <- log(sam@stock.n[                                     ,,drop=T])
    }

  return(parameters)
}

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


