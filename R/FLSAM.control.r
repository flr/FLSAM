setClass("FLSAM.control",
 representation(
    name            ="character",
    desc            ="character",
    range	    ="numeric",   ## min and max age represented internally in the assessment
    fleets          ="numeric",  ## fleet meta data
    plus.group      ="logical",   ## we model the maximum age as a plus group or not?
    states          ="matrix",   ## matrix describing the coupling of the states
    logN.vars       ="vector",   ## vector coupling the logN variances
    catchabilities  ="matrix",   ## matrix coupling the catachability parameters
    power.law.exps  ="matrix",   ## matrix coupling the power law expoonents
    f.vars          ="matrix",   ## matrix of fishing mortality couplings
    obs.vars        ="matrix",   ## matrix coupling the observation variances
    obs.weight      ="matrix",    ## weight of observations
    srr             ="integer",   ## stock recruitment relationship
    scaleNoYears    ="integer",  ## number of years for which to apply a catch multiplier
    scaleYears      ="vector",   ## actual years for which scaling applies
    scalePars       ="matrix",   ## parameter binding over years and ages
    cor.F           ="vector",   ## should we estimate the Fs a being correlated?
    cor.obs         ="matrix",    ## parameter binding for correlated observations
    cor.obs.Flag    ="factor",    ## Type of observation correlation
    biomassTreat    ="vector",    ## way biomass indices are treated
    nohess          ="logical",   ## should the hessian be estimated?
    timeout         ="numeric", ## number of seconds before model timesout 
    likFlag         ="factor",  ## Likelihood flags
    fixVarToWeight  ="logical",  ## Fixing variance to weight
    simulate        ="logical",  ## Simulate data
    residuals       ="logical",
    sam.binary      ="character"),  ## path to sam binary
  prototype=prototype(
    range           =as.numeric(1),   ## minimum age represented internally in the assessment
    plus.group      =as.logical(TRUE),   ## model the maximum age as a plus group?
    states          =as.matrix(0),   ## matrix describing the coupling of the states
    logN.vars       =as.vector(0),   ## vector of coupling of the logN variables
    catchabilities  =as.matrix(0),   ## matrix characterising the coupling of the parameters
    power.law.exps  =as.matrix(0),   ## matrix of survey catchabilities
    f.vars          =as.matrix(0),   ## matrix of fishing mortality couplings
    obs.vars        =as.matrix(0),   ## matrix coupling the observation variances
    obs.weight      =as.matrix(NA),  ## weight set to NA by default
    cor.F           =as.vector(0),   ## Don't estimate by default
    cor.obs         =as.matrix(0),   ## matrix of correlation in observations
    cor.obs.Flag    =factor(NA,levels=c("ID","AR","US")), ## Identity, auto-regressive or user defined
    biomassTreat    =as.vector(0),   ## vector of biomassTreatment
    srr             =as.integer(0),  ## stock recruitment relationship
    scaleNoYears    =as.integer(0),  ## number of years for which to apply a catch multiplier
    scaleYears      =as.vector(0),  ## actual years for which scaling applies
    scalePars       =as.matrix(NA),  ## parameter binding over years and ages
    nohess          =as.logical(FALSE),          ## Estimate hessian by default 
    timeout         =3600,           ## defaults to one hour
    likFlag         =factor(NA,levels=c("LN","ALN")),
    fixVarToWeight  =as.logical(FALSE),
    simulate        =as.logical(FALSE),
    residuals       =as.logical(TRUE)),
  validity=function(object){
               	# Everything is fine
               	return(TRUE)}
)

FLSAM.control <- function(stck,tun,default="full",scaleYears=integer(),sam.binary="missing") {
  #Default constructor
  #Create object
  ctrl <- new("FLSAM.control")
  if(!missing("sam.binary")) ctrl@sam.binary <- sam.binary

  if(any(!sapply(tun,type) %in% c("con","number","biomass")))
    stop("index type not recognized. Use con, number or biomass")

  #Populate metadata
  ctrl@name  <- stck@name
  ctrl@desc  <- stck@desc
  ctrl@range <- stck@range
  ctrl@range["maxyear"] <- max(stck@range["maxyear"],
                            max(sapply(tun,function(x) max(x@range[c("maxyear")])))) 
  ctrl@plus.group <- ifelse(is.na(stck@range["plusgroup"]==stck@range["max"])==F,stck@range["plusgroup"]==stck@range["max"],F)
  ctrl@fleets <-factor(sapply(tun,type),levels=c("con","number","biomass"))
  ctrl@fleets <- c(0,as.numeric(ctrl@fleets))
  fleet.names <- c("catch",names(tun))
  names(ctrl@fleets) <- fleet.names

  #Setup coupling structures - default is for each parameter to be fitted independently
  default.coupling <- matrix(as.integer(NA),nrow=dims(stck)$area+length(tun),ncol=dims(stck)$age,
                        dimnames=list(fleet=fleet.names,age=dimnames(stck@catch.n)$age))
  ctrl@states           <- default.coupling
  ctrl@catchabilities   <- default.coupling
  ctrl@power.law.exps   <- default.coupling
  ctrl@f.vars           <- default.coupling
  ctrl@obs.vars         <- default.coupling
  ctrl@obs.weight       <- default.coupling
  ctrl@cor.obs          <- default.coupling[,-1]
  colnames(ctrl@cor.obs)<- apply(cbind(dimnames(stck@catch.n)$age[-length(dimnames(stck@catch.n)$age)],dimnames(stck@catch.n)$age[-1]),1,paste,collapse="-")
  ctrl@cor.obs.Flag     <- factor(rep("ID",length(fleet.names)),levels=c("ID","AR","US"))
  ctrl@biomassTreat     <- rep(-1,length(fleet.names))
  ctrl@biomassTreat[which(ctrl@fleets %in% c(3,4))] <- 0
  ctrl@likFlag          <- factor(rep("LN",length(fleet.names)),levels=c("LN","ALN"))    #1=LN, 2=ALN
  
  #Other variables
  ctrl@logN.vars        <- default.coupling[1,]
  ctrl@srr              <- as.integer(0)
  
  #Scaling variables
  if(length(scaleYears)>0){
    ctrl@scaleNoYears     <- length(scaleYears)
    ctrl@scaleYears       <- sort(scaleYears)
    ctrl@scalePars        <- matrix(as.integer(NA),nrow=ctrl@scaleNoYears,ncol=dims(stck)$age,
                              dimnames=list(years=scaleYears,age=dimnames(stck@catch.n)$age))
  } else {
    ctrl@scaleNoYears     <- as.integer(0)
    ctrl@scaleYears       <- NA
    ctrl@scalePars        <- matrix(as.integer(NA),nrow=ctrl@scaleNoYears,ncol=dims(stck)$age,
                              dimnames=list(years=scaleYears,age=dimnames(stck@catch.n)$age))
  }

  #Handle full coupling case
  if(default=="full"){
    idx                 <- lapply(tun,function(x){return(if(x@type == "biomass"){ 1
                                                         }else{ seq(as.numeric(range(x)["min"]),as.numeric(range(x)["max"]),1)})})
    idx[[names(ctrl@fleets)[which(ctrl@fleets==0)]]] <- seq(as.numeric(range(stck)["min"]),as.numeric(range(stck)["max"]),1)

    full.coupling       <- matrix(as.integer(NA),nrow=1+length(tun),ncol=dims(stck)$age,
                            dimnames=list(fleet=fleet.names,age=dimnames(stck@catch.n)$age))
    for(iFlt in names(ctrl@fleets[ctrl@fleets!=3]))
      full.coupling[iFlt,ac(idx[[iFlt]])] <- (sum(is.finite(full.coupling),na.rm=T)+1):(sum(is.finite(full.coupling),na.rm=T)+length(idx[[iFlt]]))-1
    for(iFlt in names(ctrl@fleets[ctrl@fleets==3]))
      full.coupling[iFlt,idx[[iFlt]]] <- (sum(is.finite(full.coupling),na.rm=T)+1):(sum(is.finite(full.coupling),na.rm=T)+length(idx[[iFlt]]))-1

    ctrl@states["catch",] <- 0:(length(ctrl@states["catch",])-1)                 #Fixed selection pattern
    ctrl@catchabilities   <- full.coupling   
    ctrl@catchabilities["catch",] <- NA        #Catchabilities don't make sense for "catch"
    ctrl@power.law.exps   <- default.coupling  #Don't fit power laws by default
    ctrl@obs.vars         <- full.coupling
    #Random walks default to a single variance
    ctrl@logN.vars[]      <- 0:(length(ctrl@logN.vars)-1)
    ctrl@f.vars["catch",] <- 0:(length(ctrl@f.vars["catch",])-1)
    if(length(scaleYears)>0)
      ctrl@scalePars[]    <- 0:(length(c(ctrl@scalePars))-1)
    ctrl@cor.obs[]        <- an(af(full.coupling[,-1]))-1
  }
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
    if(!missing("fleets")) {
      for(slt.name in slotNames(object)) {
        slt <- slot(object,slt.name)
        if(class(slt)=="matrix") {
            remaining.fleets <- setdiff(rownames(slt),fleets)
            slot(object,slt.name) <- slt[remaining.fleets,,drop=FALSE]}}
      whichFleet      <- which(names(object@fleets) %in% fleets)
      object@fleets   <- object@fleets[setdiff(names(object@fleets),fleets)]
      object@cor.obs.Flag <- object@cor.obs.Flag[-whichFleet]
      object@biomassTreat <- object@biomassTreat[-whichFleet]
      object@likFlag      <- object@likFlag[-whichFleet]}
    #Then drop the ages 
    if(!missing("ages")) {
      for(slt.name in slotNames(object)) {
        slt <- slot(object,slt.name)
        if(class(slt)=="matrix") {
          remaining.ages <- setdiff(colnames(slt),ages)
          if(slt == "cor.obs"){
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
      slot(object,iSlt)[]  <-  as.numeric(factor(slot(object,iSlt)))-1
    }
    if(iSlt %in% c("logN.vars","scalePars"))
      slot(object,iSlt)[]  <- as.numeric(factor(slot(object,iSlt)))-1
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
    if(any(apply(matrix(unlist(dmns),length(unlist(dmns[[1]])),length(dmns)),1,duplicated)[2,]==F)) stop("Dimnames are not the same for parameter sections")
    if(any(apply(matrix(unlist(dms), length(unlist(dms[[1]])), length(dms)), 1,duplicated)[2,]==F)) stop("Dimensions are not the same for parameter sections")
  }
)

