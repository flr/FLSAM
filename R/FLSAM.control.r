setClass("FLSAM.control",
 representation(
    name            ="character",
    desc            ="character",
    range	    ="numeric",   ## min and max age represented internally in the assessment
    fleets          ="numeric",  ## fleet meta data
    plus.group      ="logical",   ## we model the maximum age as a plus group or not?
    states          ="matrix",   ## matrix describing the coupling of the states
    logN.vars       ="vector",   ## vector coupling the logN variances
    logP.vars       ="vector",   ## vector coupling the logP variances
    catchabilities  ="matrix",   ## matrix coupling the catachability parameters
    power.law.exps  ="matrix",   ## matrix coupling the power law expoonents
    f.vars          ="matrix",   ## matrix of fishing mortality couplings
    obs.vars        ="matrix",   ## matrix coupling the observation variances
    srr             ="integer",   ## stock recruitment relationship
    scaleNoYears    ="integer",  ## number of years for which to apply a catch multiplier
    scaleYears      ="vector",   ## actual years for which scaling applies
    scalePars       ="matrix",   ## parameter binding over years and ages
    cor.F           ="vector",   ## should we estimate the Fs a being correlated?
    cor.obs         ="matrix",    ## parameter binding for correlated observations
    cor.obs.Flag    ="factor",    ## Type of observation correlation
    biomassTreat    ="vector",    ## way biomass indices are treated
    timeout         ="numeric", ## number of seconds before model timesout
    likFlag         ="factor",  ## Likelihood flags
    fixVarToWeight  ="logical",  ## Fixing variance to weight
    simulate        ="logical",  ## Simulate data
    residuals       ="logical",
    sumFleets       ="vector"),    ## Define which catch fleets should be summed over a certain period
  prototype=prototype(
    range           =as.numeric(1),   ## minimum age represented internally in the assessment
    plus.group      =as.logical(TRUE),   ## model the maximum age as a plus group?
    states          =as.matrix(0),   ## matrix describing the coupling of the states
    logN.vars       =as.vector(0),   ## vector of coupling of the logN variables
    logP.vars       =as.vector(0),   ## vector of coupling of the logP variables
    catchabilities  =as.matrix(0),   ## matrix characterising the coupling of the parameters
    power.law.exps  =as.matrix(0),   ## matrix of survey catchabilities
    f.vars          =as.matrix(0),   ## matrix of fishing mortality couplings
    obs.vars        =as.matrix(0),   ## matrix coupling the observation variances
    cor.F           =as.vector(0),   ## Don't estimate by default
    cor.obs         =as.matrix(0),   ## matrix of correlation in observations
    cor.obs.Flag    =factor(NA,levels=c("ID","AR","US")), ## Identity, auto-regressive or user defined
    biomassTreat    =as.vector(0),   ## vector of biomassTreatment
    srr             =as.integer(0),  ## stock recruitment relationship
    scaleNoYears    =as.integer(0),  ## number of years for which to apply a catch multiplier
    scaleYears      =as.vector(0),  ## actual years for which scaling applies
    scalePars       =as.matrix(NA),  ## parameter binding over years and ages
    timeout         =3600,           ## defaults to one hour
    likFlag         =factor(NA,levels=c("LN","ALN")),
    fixVarToWeight  =as.logical(FALSE),
    simulate        =as.logical(FALSE),
    residuals       =as.logical(TRUE),
    sumFleets       =as.vector(0)),    ## Define which catch fleets should be summed over a certain period
  validity=function(object){
               	# Everything is fine
               	return(TRUE)}
)



FLSAM.control <- function(stcks,tun,sumFleets=vector(),catch.vars=NULL,scaleYears=NULL) {
  #Default constructor
  #Create object
  if(length(stcks)>1 & length(sumFleets)==0) stop("sumFleets must be specified when stcks are provided")
  if(class(stcks)=="FLStock") stcks <- FLStocks(residual=stcks)
  ctrl <- new("FLSAM.control")

  ctrlSAM <- defcon(FLSAM2SAM(stcks,tun,sumFleets,catch.vars))
  if(any(!unlist(lapply(tun,function(x){return(x@type)})) %in% c("con","number","biomass","partial")))
    stop("index type not recognized. Use con, number or biomass")

  #Populate metadata
  ctrl@name  <- stcks[["residual"]]@name
  ctrl@desc  <- stcks[["residual"]]@desc
  ctrl@range <- stcks[["residual"]]@range
  
  ctrl@range["maxyear"] <- max(sapply(stcks,function(x) max(x@range["maxyear"])),
                            max(sapply(tun,function(x) max(x@range[c("maxyear")])))) 
  ctrl@range["minyear"] <- min(sapply(stcks,function(x) min(x@range["minyear"])),
                            min(sapply(tun,function(x) min(x@range[c("minyear")]))))
  ctrl@plus.group <- ifelse(is.na(stcks[["residual"]]@range["plusgroup"]==stcks[["residual"]]@range["max"])==F,stcks[["residual"]]@range["plusgroup"]==stcks[["residual"]]@range["max"],F)
  ctrl@fleets <-factor(unlist(lapply(tun,function(x){return(x@type)})),levels=c("con","number","biomass","partial"))
  ctrl@fleets <- c(na.omit(c(rep(0,dims(stcks[["residual"]])$area),as.numeric(ctrl@fleets),ifelse(length(sumFleets)==0,numeric(),7))))
  if(length(ctrl@fleets==4)>0)
    ctrl@fleets[ctrl@fleets==4] <- 6
  fleet.names <- c(na.omit(c(paste("catch",dimnames(stcks[["residual"]]@catch.n)$area),names(tun),ifelse(length(sumFleets)==0,character(),"sumFleet"))))
  names(ctrl@fleets) <- fleet.names

  #Setup coupling structures - default is for each parameter to be fitted independently
  default.coupling <- matrix(as.integer(NA),nrow=dims(stcks[["residual"]])$area+length(tun)+ifelse(length(sumFleets)==0,0,1),ncol=dims(stcks[["residual"]])$age,
                        dimnames=list(fleet=fleet.names,age=dimnames(stcks[["residual"]]@catch.n)$age))
  ctrl@states           <- default.coupling
  ctrl@catchabilities   <- default.coupling
  ctrl@power.law.exps   <- default.coupling
  ctrl@f.vars           <- default.coupling
  ctrl@obs.vars         <- default.coupling
  ctrl@cor.obs          <- default.coupling[,-1]
  colnames(ctrl@cor.obs)<- apply(cbind(dimnames(stcks[["residual"]]@catch.n)$age[-length(dimnames(stcks[["residual"]]@catch.n)$age)],dimnames(stcks[["residual"]]@catch.n)$age[-1]),1,paste,collapse="-")
  ctrl@cor.obs.Flag     <- factor(rep("ID",length(fleet.names)),levels=c("ID","AR","US"))
  ctrl@biomassTreat     <- rep(-1,length(fleet.names))
  ctrl@biomassTreat[which(ctrl@fleets %in% c(3,4))] <- 0
  ctrl@likFlag          <- factor(rep("LN",length(fleet.names)),levels=c("LN","ALN"))    #1=LN, 2=ALN
  ctrl@cor.F            <- rep(0,dims(stcks[["residual"]])$area)
  
  #Other variables
  ctrl@logN.vars        <- default.coupling[1,]
  ctrl@srr              <- as.integer(0)
  ctrl@logP.vars        <- NA
  if(length(ctrl@fleets==6)>0)
    ctrl@logP.vars      <- 1:(length(which(ctrl@fleets==6))-1)

  #Scaling variables
  if(length(scaleYears)>0){
    ctrl@scaleNoYears     <- length(scaleYears)
    ctrl@scaleYears       <- sort(scaleYears)
    ctrl@scalePars        <- matrix(as.integer(NA),nrow=ctrl@scaleNoYears,ncol=dims(stcks[["residual"]])$age,
                              dimnames=list(years=scaleYears,age=dimnames(stcks[["residual"]]@catch.n)$age))
  } else {
    ctrl@scaleNoYears     <- as.integer(0)
    ctrl@scaleYears       <- NA
    ctrl@scalePars        <- matrix(as.integer(NA),nrow=ctrl@scaleNoYears,ncol=dims(stcks[["residual"]])$age,
                              dimnames=list(years=scaleYears,age=dimnames(stcks[["residual"]]@catch.n)$age))
  }

  #Handle full coupling case
  ctrl@states[]           <- ctrlSAM$keyLogFsta
  ctrl@logN.vars[]        <- ctrlSAM$keyVarLogN
  ctrl@logP.vars          <- ctrlSAM$keyVarLogP
  ctrl@catchabilities[]   <- ctrlSAM$keyLogFpar
  ctrl@power.law.exps[]   <- ctrlSAM$keyQpow
  ctrl@f.vars[]           <- ctrlSAM$keyVarF
  ctrl@obs.vars[]         <- ctrlSAM$keyVarObs
  ctrl@srr                <- as.integer(ctrlSAM$stockRecruitmentModelCode)
  ctrl@cor.F[]            <- ctrlSAM$corFlag
  ctrl@cor.obs[]          <- ctrlSAM$keyCorObs
  ctrl@cor.obs.Flag[]     <- ctrlSAM$obsCorStruct
  ctrl@biomassTreat       <- ctrlSAM$keyBiomassTreat
  ctrl@likFlag[]          <- ctrlSAM$obsLikelihoodFlag
  ctrl@fixVarToWeight     <- ifelse(ctrlSAM$fixVarToWeight==0,FALSE,TRUE)
  ctrl@sumFleets          <- sumFleets

  ctrl <- update(ctrl)

  #Finished!
  return(ctrl)
}

#Function for dropping ages or surveys from a FLSAM.control
if (!isGeneric("drop.from.control")) {
    setGeneric("drop.from.control",function(object,...) standardGeneric("drop.from.control"))
}
setMethod("drop.from.control",signature(object="FLSAM.control"),
  function(object,fleets="missing",ages="missing") {
    #Drop the fleets first
    if(!missing("fleets")){
      whichFleet      <- which(names(object@fleets) %in% fleets)
      if(any(object@fleets[whichFleet]==6)){
        fleetsPart    <- names(object@fleets[which(object@fleets == 6)])
        object@logP.vars <- numeric()
        fleets        <- unique(c(fleets,fleetsPart))
      }
      for(slt.name in slotNames(object)){
        slt <- slot(object,slt.name)
        if(class(slt)=="matrix"){
          remaining.fleets <- setdiff(rownames(slt),fleets)
          slot(object,slt.name) <- slt[remaining.fleets,,drop=FALSE]
        }
      }
      whichFleet      <- which(names(object@fleets) %in% fleets)
      object@fleets   <- object@fleets[setdiff(names(object@fleets),fleets)]
      object@cor.obs.Flag <- object@cor.obs.Flag[-whichFleet]
      object@biomassTreat <- object@biomassTreat[-whichFleet]
      object@likFlag      <- object@likFlag[-whichFleet]
    }
    #Then drop the ages 
    if(!missing("ages")) {
      for(slt.name in slotNames(object)) {
        slt <- slot(object,slt.name)
        if(class(slt)=="matrix") {
          remaining.ages <- setdiff(colnames(slt),ages)
          if(slt.name == "cor.obs"){
            mincol <- unique(unlist(sapply(ages,grep,colnames(object@cor.obs))))
            slot(object,slt.name) <- slt[,-mincol]
          } else {
              slot(object,slt.name) <- slt[,remaining.ages,drop=FALSE]
            }
        }
      }
      object@logN.vars   <- object@logN.vars[setdiff(names(object@logN.vars),ages)]
    }
    #Finally, do an update
    object <- update(object)
    return(object)
})

setMethod("update", signature(object="FLSAM.control"),
  function(object){

  for(iSlt in slotNames(object)){
    if(class(slot(object,iSlt))=="matrix"){
      isNA <- which(slot(object,iSlt)==-1)
      slot(object,iSlt)[isNA] <- NA
      slot(object,iSlt)[]  <-  as.numeric(factor(slot(object,iSlt)))-1
      slot(object,iSlt)[isNA] <- -1
    }
    if(iSlt %in% c("logN.vars","scalePars","logP.vars")){
      isNA <- which(slot(object,iSlt)==-1)
      slot(object,iSlt)[]  <- as.numeric(factor(slot(object,iSlt)))-1
      slot(object,iSlt)[isNA] <- -1
    }
  }
  return(object)}
)

setValidity("FLSAM.control",
  function(object){
    #-1 All slots populated
    if(object@plus.group) if(is.na(object@range["plusgroup"])==T){stop("Specify plusgroup in object@range")}
    if(any(is.na(object@range[-3])==T)) stop("Specify values in object@range")
    if(any(is.na(object@logN.vars)==T)) stop("Not all ages in logN.vars specified")

    #-2 Settings within range
    if(object@srr < 0 | object@srr > 2) stop("SRR type outisde possible setting range")
    if(range(object@fleets)[1] < 0 | range(object@fleets)[2] > 4) stop(paste("Fleet type",names(object@fleets)[which(object@fleets < 0 | object@fleets > 4)],"outside possible setting range"))

    #-3 Dimension and dimnames check
    dmns <- list()
    dms  <- list()
    for(iSlt in slotNames(object)){
      if(class(slot(object,iSlt)) == "matrix"){
        dmns[[iSlt]]  <- dimnames(slot(object,iSlt))
        dms[[iSlt]]   <- dim(slot(object,iSlt))
      }
    }
    
    #check sumFleets year ranges
    if(any(apply(matrix(unlist(dmns),length(unlist(dmns[[1]])),length(dmns)),1,duplicated)[2,]==F)) stop("Dimnames are not the same for parameter sections")
    if(any(apply(matrix(unlist(dms), length(unlist(dms[[1]])), length(dms)), 1,duplicated)[2,]==F)) stop("Dimensions are not the same for parameter sections")
  }
)

