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
    srr             ="integer"),    ## stock recruitment relationship
  prototype=prototype(
    range           =as.numeric(1),   ## minimum age represented internally in the assessment
    plus.group      =as.logical(TRUE),   ## model the maximum age as a plus group?
    states          =as.matrix(0),   ## matrix describing the coupling of the states
    logN.vars       =as.vector(0),   ## vector of coupling of the logN variables
    catchabilities  =as.matrix(0),   ## matrix characterising the coupling of the parameters
    power.law.exps  =as.matrix(0),   ## matrix of survey catchabilities
    f.vars          =as.matrix(0),   ## matrix of fishing mortality couplings
    obs.vars        =as.matrix(0),   ## matrix coupling the observation variances
    srr             =as.integer(0)),    ## stock recruitment relationship
  validity=function(object){
               	# Everything is fine
               	return(TRUE)}
)

FLSAM.control <- function(stck,tun,default="full") {
  #Default constructor
  #Create object
  ctrl <- new("FLSAM.control")

  #Populate metadata
  ctrl@range <- stck@range
  ctrl@range["maxyear"] <- max(stck@range["maxyear"],
                            max(sapply(tun,function(x) max(x@range[c("maxyear")])))) 
  ctrl@plus.group <- stck@range["plusgroup"]==stck@range["max"]
  ctrl@fleets <-factor(sapply(tun,type),levels=c("con","number","biomass"))
  ctrl@fleets <- c(0,as.numeric(ctrl@fleets))
  names(ctrl@fleets) <- c("catch",sapply(tun,name))

  #Setup coupling structures - default is for each parameter to be fitted independently
  default.coupling <- matrix(as.integer(NA),nrow=1+length(tun),ncol=dims(stck)$age,
                        dimnames=list(fleet=c("catch",names(tun)),age=dimnames(stck@catch.n)$age))
  ctrl@states           <- default.coupling
  ctrl@catchabilities   <- default.coupling
  ctrl@power.law.exps   <- default.coupling
  ctrl@f.vars           <- default.coupling
  ctrl@obs.vars         <- default.coupling

  #Other variables
  ctrl@logN.vars            <- default.coupling[1,]
  ctrl@srr <- as.integer(0)

  #Handle full coupling case
  if(default=="full"){
    idx                 <- lapply(tun,function(x){return(if(x@type == "biomass"){ NA
                                                         }else{ seq(as.numeric(range(x)["min"]),as.numeric(range(x)["max"]),1)})})
    idx[[names(ctrl@fleets)[ctrl@fleets==0]]] <- seq(as.numeric(range(stck)["min"]),as.numeric(range(stck)["max"]),1)

    full.coupling       <- matrix(as.integer(NA),nrow=1+length(tun),ncol=dims(stck)$age,
                            dimnames=list(fleet=c("catch",names(tun)),age=dimnames(stck@catch.n)$age))
    for(iFlt in names(ctrl@fleets[ctrl@fleets!=3]))
      full.coupling[iFlt,ac(idx[[iFlt]])] <- (sum(is.finite(full.coupling),na.rm=T)+1):(sum(is.finite(full.coupling),na.rm=T)+length(idx[[iFlt]]))

    ctrl@states["catch",] <- 1                 #Fixed selection pattern
    ctrl@catchabilities   <- full.coupling   
    ctrl@catchabilities["catch",] <- NA        #Catchabilities don't make sense for "catch"
    ctrl@power.law.exps   <- default.coupling  #Don't fit power laws by default
    ctrl@obs.vars         <- full.coupling
    #Random walks default to a single variance
    ctrl@logN.vars[]      <- 1
    ctrl@f.vars["catch",] <- 1
    ctrl <- update(ctrl)
  }


  #Finished!
  return(ctrl)
}

#Function for dropping ages or surveys from a FLSAM.control
if (!isGeneric("drop")) {
    setGeneric("drop",function(object,...) standardGeneric("drop"))
}
setMethod("drop",signature(object="FLSAM.control"),
  function(object,fleets="missing",ages="missing") {
    #Drop the fleets first
    if(!missing("fleets")) {
      for(slt.name in slotNames(object)) {
        slt <- slot(object,slt.name)
        if(class(slt)=="matrix") {
            remaining.fleets <- setdiff(rownames(slt),fleets)
            slot(object,slt.name) <- slt[remaining.fleets,,drop=FALSE]}}
      object@fleets   <- object@fleets[setdiff(names(object@fleets),fleets)]}
    #Then drop the ages 
    if(!missing("ages")) {
      for(slt.name in slotNames(object)) {
        slt <- slot(object,slt.name)
        if(class(slt)=="matrix") {
            remaining.ages <- setdiff(colnames(slt),ages)
            slot(object,slt.name) <- slt[,remaining.ages,drop=FALSE]}}
      object@logN.vars   <- object@logN.vars[setdiff(names(object@logN.vars),ages)]}

    #Finally, do an update
    object <- update(object)
    return(object)
})

setMethod("update", signature(object="FLSAM.control"),
  function(object){

  for(iSlt in slotNames(object)){
    if(class(slot(object,iSlt))=="matrix"){
      slot(object,iSlt)[]  <-  as.numeric(factor(slot(object,iSlt)))
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
    if(any(apply(matrix(unlist(dmns),length(unlist(dmns[[1]])),length(dmns)),1,duplicated)[2,]==F)) stop("Dimnames are not the same for parameter sections")
    if(any(apply(matrix(unlist(dms), length(unlist(dms[[1]])), length(dms)), 1,duplicated)[2,]==F)) stop("Dimensions are not the same for parameter sections")
  }
)

