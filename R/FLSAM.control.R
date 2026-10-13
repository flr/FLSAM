#' FLSAM.control class
#'
#' Configuration object for an \code{\link{FLSAM}} assessment. It holds the
#' FLR representation of the SAM configuration list (see
#' \code{\link[stockassessment]{defcon}} and \code{\link[stockassessment]{sam.fit}}).
#' Most slots are coupling ("binding") matrices: rows are fleets, columns are
#' ages. Equal integers indicate that one parameter is shared, consecutive
#' integers starting at zero must be used, and \code{-1} (or \code{NA}) means
#' the parameter is not used. Run \code{update()} on the object after editing
#' couplings to renumber them consecutively.
#'
#' @slot name,desc Name and description.
#' @slot range Named vector with \code{min}, \code{max}, \code{plusgroup},
#'   \code{minyear}, \code{maxyear}, \code{minfbar}, \code{maxfbar}.
#' @slot fleets Named vector with the SAM fleet type of each fleet: 0 catch,
#'   2 age-structured survey, 3 biomass/SSB index, 6 partial survey (components),
#'   7 sum-of-fleets catch.
#' @slot plus.group Per fleet: is the oldest age a plus group (1) or not (0)
#'   (SAM: \code{maxAgePlusGroup}).
#' @slot states Coupling of the fishing mortality states (SAM: \code{keyLogFsta}).
#' @slot logN.vars Coupling of the log(N) process variances by age (SAM: \code{keyVarLogN}).
#' @slot logP.vars Coupling of the log(P) process variances of partial surveys (SAM: \code{keyVarLogP}).
#' @slot catchabilities Coupling of survey catchabilities (SAM: \code{keyLogFpar}).
#' @slot power.law.exps Coupling of density-dependent catchability power
#'   parameters (SAM: \code{keyQpow}). Also allowed for SSB and FSB indices.
#' @slot f.vars Coupling of the log(F) process variances (SAM: \code{keyVarF}).
#' @slot obs.vars Coupling of observation variances (SAM: \code{keyVarObs}).
#' @slot srr Stock-recruitment code (SAM: \code{stockRecruitmentModelCode}):
#'   0 random walk, 1 Ricker, 2 Beverton-Holt, 3 piece-wise constant,
#'   61 hockey stick, 62 AR(1), 63 bent hyperbola, 64/65 power functions,
#'   66 Shepherd, 67 Deriso, 68 Saila-Lorda, 69 sigmoidal Beverton-Holt,
#'   90-93 splines, and the depensatory variants (2xx, 4xx).
#' @slot scaleNoYears,scaleYears,scalePars Catch scaling: number of years,
#'   the years, and the year x age coupling of the scale parameters
#'   (SAM: \code{noScaledYears}, \code{keyScaledYears}, \code{keyParScaledYA}).
#' @slot cor.F Correlation structure of the F increments per catch fleet
#'   (SAM: \code{corFlag}): 0 independent, 1 compound symmetry, 2 AR(1),
#'   3 separable AR(1).
#' @slot cor.obs Coupling of the AR(1) observation correlation parameters
#'   (SAM: \code{keyCorObs}).
#' @slot cor.obs.Flag Observation covariance structure per fleet: \code{"ID"},
#'   \code{"AR"} or \code{"US"} (SAM: \code{obsCorStruct}).
#' @slot biomassTreat Treatment of biomass indices (SAM: \code{keyBiomassTreat}):
#'   -1 not a biomass index, 0 SSB index, 1 catch index, 2 FSB index,
#'   3 total catch, 4 total landings, 5 TSB index, 6 TSN index, 10 Fbar index.
#' @slot timeout Not used anymore (kept for backward compatibility).
#' @slot likFlag Observation likelihood per fleet: \code{"LN"} log-normal or
#'   \code{"ALN"} additive logistic normal (SAM: \code{obsLikelihoodFlag}).
#' @slot fixVarToWeight Logical per fleet: if weights are supplied
#'   (\code{catch.vars} or \code{index.var}) are they relative weights
#'   (\code{FALSE}) or fixed variances (\code{TRUE}). A single value is
#'   recycled over fleets (SAM: \code{fixVarToWeight}).
#' @slot fracMixF,fracMixN,fracMixObs Fraction of a t(3) distribution mixed
#'   into the logF increments, the logN increments (by age) and the
#'   observations (by fleet).
#' @slot constRecBreaks Break years for srr 3, or log-SSB knots for spline
#'   recruitment models.
#' @slot predVarObsLink Coupling of the prediction-variance link parameters
#'   (fleet x age).
#' @slot stockWeightModel,stockWeightMean,stockWeightObsVar Stock weight
#'   process model switch and couplings of its mean and observation variance.
#' @slot catchWeightModel,catchWeightMean,catchWeightObsVar Catch weight
#'   process model switch and couplings.
#' @slot maturityModel,maturityMean Maturity process model switch and coupling
#'   (SAM: \code{matureModel}, \code{keyMatureMean}).
#' @slot mortalityModel,mortalityMean,mortalityObsVar Natural mortality process
#'   model switch and couplings.
#' @slot XtraSd Integer matrix with 4 columns (fleet, year, age, coupling)
#'   adding extra observation uncertainty (SAM: \code{keyXtraSd}).
#' @slot logNMeanAssumption What the population model corresponds to for
#'   recruitment and other ages: 0 median, 1 mean, 2 mode.
#' @slot initState Add initial-state parameters for the latent processes (0/1).
#' @slot hockeyStickCurve Smoothness of srr 63 (bent hyperbola). \code{NA}
#'   (default) estimates it, a positive value fixes it.
#' @slot simulate Simulate conditional on estimated random effects (SAM:
#'   \code{sim.condRE}).
#' @slot residuals Calculate one-observation-ahead residuals in \code{FLSAM}.
#' @slot sumFleets Names of the catch fleets that are summed in the sum fleet.
#' @slot recapture Tag-recapture data (\code{data.frame} with columns \code{r}
#'   (number recaptured), \code{ReleaseY}, \code{RecaptureY}, \code{Yearclass},
#'   \code{Nscan}, \code{R} (number released) and \code{Type}), modelled as
#'   fleet \code{"Recaptures"} (SAM fleet type 5). Empty when no tag data are used.
#'
#' @seealso \code{\link{FLSAM.control}} (constructor), \code{\link{drop.from.control}},
#'   \code{\link{ctrl2conf}}
#' @rdname FLSAM.control-class
#' @exportClass FLSAM.control
setClass("FLSAM.control",
 representation(
    name            ="character",
    desc            ="character",
    range	    ="numeric",   ## min and max age represented internally in the assessment
    fleets          ="numeric",  ## fleet meta data
    plus.group      ="vector",   ## we model the maximum age as a plus group or not?
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
    fracMixF        = "numeric", ##The fraction of t(3) distribution used in logF increment distribution
    fracMixN        = "vector", ##The fraction of t(3) distribution used in logN increment distribution (for each age group)"rep(0,noAges)
    fracMixObs      = "vector", ##rep(0,length(data$noFleets))
    constRecBreaks  = "vector", ##Vector of break years between which recruitment is at constant level. The break year is included in the left interval. (This option is only used in combination with stock-recruitment code 3)
    predVarObsLink  =  "matrix",##Coupling of parameters used in a prediction-variance link for observations. 
    stockWeightModel  = "logical", ##Integer code describing the treatment of stock weights in the model (0 use as known, 1 use as observations to inform stock weight process (GMRF with cohort and within year correlations))
    stockWeightMean   = "vector",## # Coupling of stock-weight process mean parameters (not used if stockWeightModel==0)
    stockWeightObsVar = "vector", ##Coupling of stock-weight observation variance parameters (not used if stockWeightModel==0)
    catchWeightModel  = "logical", ##Integer code describing the treatment of catch weights in the model (0 use as known, 1 use as observations to inform catch weight process (GMRF with cohort and within year correlations))                                  
    catchWeightMean   = "vector",
    catchWeightObsVar = "vector", ##
    maturityModel     = "logical",##
    maturityMean      = "vector",##
    mortalityModel    = "logical",##
    mortalityMean     = "vector",##
    mortalityObsVar   = "vector",##
    XtraSd            = "matrix",##
    logNMeanAssumption= "vector",##
    initState         = "numeric",##
    hockeyStickCurve  = "numeric", ## Smoothness of srr 63, NA = estimate
    simulate          ="logical",  ## Simulate data
    residuals         ="logical",
    sumFleets         ="vector",     ## Define which catch fleets should be summed over a certain period
    recapture         ="data.frame"), ## Tag-recapture data (SAM fleet type 5)
  prototype=prototype(
    range           =as.numeric(1),   ## minimum age represented internally in the assessment
    plus.group      =as.vector(1),   ## model the maximum age as a plus group?
    states          =as.matrix(0),   ## matrix describing the coupling of the states
    logN.vars       =as.vector(0),   ## vector of coupling of the logN variables
    logP.vars       =as.numeric(0),   ## vector of coupling of the logP variables
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
    fracMixF        = as.numeric(0),
    fracMixN        = as.vector(0),
    fracMixObs      = as.vector(0),
    constRecBreaks  = as.vector(0),
    predVarObsLink  = as.matrix(NA),
    stockWeightModel  = as.logical(FALSE),
    stockWeightMean   = as.vector(NA),
    stockWeightObsVar = as.vector(NA),
    catchWeightModel  = as.logical(FALSE),
    catchWeightMean   = as.vector(NA),
    catchWeightObsVar = as.vector(NA),
    maturityModel     = as.logical(FALSE),
    maturityMean      = as.vector(NA),
    mortalityModel    = as.logical(FALSE),
    mortalityMean     = as.vector(NA),
    mortalityObsVar   = as.vector(NA),
    XtraSd            = as.matrix(0),
    logNMeanAssumption= as.vector(0),
    initState         = as.numeric(0),
    hockeyStickCurve  = as.numeric(NA),
    simulate          =as.logical(FALSE),
    residuals         =as.logical(TRUE),
    sumFleets         =as.vector(0),     ## Define which catch fleets should be summed over a certain period
    recapture         =data.frame()),
  validity=function(object){
               	# Everything is fine
               	return(TRUE)}
)



#' Create a default FLSAM.control object
#'
#' Creates an \code{\link{FLSAM.control-class}} object with a syntactically
#' correct starting configuration, derived from
#' \code{\link[stockassessment]{defcon}}. The default configuration is a first
#' guess only: it is the responsibility of the user to find an appropriate
#' configuration for the stock.
#'
#' @param stcks An \code{FLStock}, or an \code{FLStocks} with elements
#'   \code{residual} and (optionally) \code{sum} for sum-of-fleets catches.
#'   Multiple catch fleets are given as areas of the \code{residual} stock.
#' @param tun An \code{FLIndices} object. Index types \code{"number"} (age
#'   structured), \code{"biomass"} (biomass/SSB index), \code{"con"} and
#'   \code{"partial"} (partial surveys modelled as components) are allowed.
#' @param sumFleets Names of the catch fleets (areas) that are summed into the
#'   \code{sum} stock; required when \code{stcks} has more than one element.
#' @param catch.vars Optional catch observation weights (see \code{\link{FLSAM}}).
#' @param scaleYears Optional vector of years for which a catch multiplier is
#'   estimated.
#' @param recapture Optional tag-recapture data (\code{data.frame} with columns
#'   \code{r}, \code{ReleaseY}, \code{RecaptureY}, \code{Yearclass}, \code{Nscan},
#'   \code{R}, \code{Type}), added as the last fleet, named \code{"Recaptures"}.
#'   See \code{\link{FLSAM.control-class}}.
#' @return An \code{FLSAM.control} object.
#' @seealso \code{\link{FLSAM.control-class}}, \code{\link{FLSAM}},
#'   \code{\link{drop.from.control}}
#' @examples
#' data(NSH)
#' ctrl <- FLSAM.control(NSH, NSH.tun)
#' # bind HERAS catchabilities for ages 5-8 and renumber
#' ctrl@catchabilities["HERAS", as.character(5:8)] <- 101
#' ctrl <- update(ctrl)
#' # treat the weights of the HERAS survey as fixed variances
#' ctrl@fixVarToWeight["HERAS"] <- TRUE
#' @export
FLSAM.control <- function(stcks,tun,sumFleets=vector(),catch.vars=NULL,scaleYears=NULL,recapture=NULL) {
  #Default constructor
  #Create object
  if(length(stcks)>1 & length(sumFleets)==0) stop("sumFleets must be specified when stcks are provided")
  if(is(stcks, "FLStock")) stcks <- FLStocks(residual=stcks)
  ctrl <- new("FLSAM.control")

  if(!is.null(recapture)) recapture <- .checkRecapture(recapture)
  ctrlSAM <- defcon(FLSAM2SAM(stcks,tun,sumFleets,catch.vars,recapture))
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
  ctrl@fleets <-factor(unlist(lapply(tun,function(x){return(x@type)})),levels=c("con","number","biomass","partial"))
  ctrl@fleets <- c(na.omit(c(rep(0,dims(stcks[["residual"]])$area),as.numeric(ctrl@fleets),ifelse(length(sumFleets)==0,numeric(),7))))
  nTag <- ifelse(is.null(recapture),0,1)
  if(length(ctrl@fleets==4)>0)
    ctrl@fleets[ctrl@fleets==4] <- 6
  fleet.names <- c(na.omit(c(paste("catch",dimnames(stcks[["residual"]]@catch.n)$area),names(tun),ifelse(length(sumFleets)==0,character(),"sumFleet"))))
  if(nTag==1){
    ctrl@fleets <- c(ctrl@fleets,5)
    fleet.names <- c(fleet.names,"Recaptures")
  }
  names(ctrl@fleets) <- fleet.names

  ctrl@plus.group <- rep(0,length(ctrl@fleets))
  for(iFleet in names(ctrl@fleets)){
    if(length(grep("catch",iFleet))>0)
      ctrl@plus.group[which(names(ctrl@fleets)==iFleet)] <- ifelse(is.na(stcks[["residual"]]@range["plusgroup"]==stcks[["residual"]]@range["max"])==0,1,0)
    if(length(grep("catch",iFleet))==0 & !iFleet %in% c("sumFleet","Recaptures"))
      ctrl@plus.group[which(names(ctrl@fleets)==iFleet)] <- ifelse(is.na(tun[[iFleet]]@range["plusgroup"]==tun[[iFleet]]@range["max"])==0,1,0)
  }

  #Setup coupling structures - default is for each parameter to be fitted independently
  default.coupling <- matrix(as.integer(NA),nrow=dims(stcks[["residual"]])$area+length(tun)+ifelse(length(sumFleets)==0,0,1)+nTag,ncol=dims(stcks[["residual"]])$age,
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
  ctrl@biomassTreat[which(ctrl@fleets %in% c(3,4,5))] <- rep(0,length(which(ctrl@fleets %in% c(3,4,5))))
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
  ctrl@fixVarToWeight     <- ifelse(ctrlSAM$fixVarToWeight==0,FALSE,TRUE); names(ctrl@fixVarToWeight) <- fleet.names

  ctrl@fracMixF           <- ctrlSAM$fracMixF
  ctrl@fracMixN           <- ctrlSAM$fracMixN
  ctrl@fracMixObs         <- ctrlSAM$fracMixObs; names(ctrl@fracMixObs) <- fleet.names 
  ctrl@constRecBreaks     <- ctrlSAM$constRecBreaks
  ctrl@predVarObsLink     <- default.coupling; ctrl@predVarObsLink[] <- -1
  ctrl@stockWeightModel   <- FALSE
  ctrl@stockWeightMean    <- default.coupling[1,]
  ctrl@stockWeightObsVar  <- default.coupling[1,]
  ctrl@catchWeightModel   <- FALSE
  ctrl@catchWeightMean    <- default.coupling[1,]
  ctrl@catchWeightObsVar  <- default.coupling[1,]
  ctrl@maturityModel        <- FALSE
  ctrl@maturityMean         <- default.coupling[1,]
  ctrl@mortalityModel     <- FALSE
  ctrl@mortalityMean      <- default.coupling[1,]
  ctrl@mortalityObsVar    <- default.coupling[1,]
  ctrl@XtraSd             <- ctrlSAM$keyXtraSd
  ctrl@logNMeanAssumption <- c(0,0)
  ctrl@initState          <- 0
  ctrl@hockeyStickCurve   <- as.numeric(NA)
  ctrl@sumFleets          <- sumFleets
  if(nTag==1) ctrl@recapture <- recapture

  ctrl <- update(ctrl)

  #Finished!
  return(ctrl)
}

#' Drop fleets or ages from an FLSAM.control object
#'
#' Removes one or more fleets and/or age groups from an
#' \code{\link{FLSAM.control-class}} object and renumbers the couplings with
#' \code{update()}. When a partial-survey fleet (type 6) is dropped, all
#' partial-survey fleets are dropped.
#'
#' @param object An \code{FLSAM.control} object.
#' @param fleets Names of the fleets to drop.
#' @param ages Names of the ages to drop.
#' @param ... Not used.
#' @return An \code{FLSAM.control} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' drop.from.control(NSH.ctrl, fleets="IBTS0")
#' drop.from.control(NSH.ctrl, ages=8:9)
#' @rdname drop.from.control
#' @export
    setGeneric("drop.from.control",function(object,...) standardGeneric("drop.from.control"))
#' @rdname drop.from.control
setMethod("drop.from.control",signature(object="FLSAM.control"),
  function(object,fleets="missing",ages="missing") {
    object <- .upgradeControl(object)
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
        if(is(slt, "matrix")){
          remaining.fleets <- setdiff(rownames(slt),fleets)
          slot(object,slt.name) <- slt[remaining.fleets,,drop=FALSE]
        }
      }
      whichFleet      <- which(names(object@fleets) %in% fleets)
      object@fleets   <- object@fleets[setdiff(names(object@fleets),fleets)]
      object@cor.obs.Flag <- object@cor.obs.Flag[-whichFleet]
      object@biomassTreat <- object@biomassTreat[-whichFleet]
      object@likFlag      <- object@likFlag[-whichFleet]
      object@fracMixObs   <- object@fracMixObs[-whichFleet]
      object@plus.group    <- object@plus.group[-whichFleet]
      if(length(object@fixVarToWeight)>1)
        object@fixVarToWeight <- object@fixVarToWeight[-whichFleet]
      if("Recaptures" %in% fleets) object@recapture <- data.frame()
    }
    #Then drop the ages 
    if(!missing("ages")) {
      for(slt.name in slotNames(object)) {
        slt <- slot(object,slt.name)
        if(is(slt, "matrix")) {
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
      object@stockWeightMean   <- object@stockWeightMean[setdiff(names(object@stockWeightMean),ages)]
      object@stockWeightObsVar   <- object@stockWeightObsVar[setdiff(names(object@stockWeightObsVar),ages)]      
      object@catchWeightMean   <- object@catchWeightMean[setdiff(names(object@catchWeightMean),ages)]
      object@catchWeightObsVar   <- object@catchWeightObsVar[setdiff(names(object@catchWeightObsVar),ages)]      
      object@maturityMean   <- object@maturityMean[setdiff(names(object@maturityMean),ages)]      
      object@mortalityMean   <- object@mortalityMean[setdiff(names(object@mortalityMean),ages)]            
      object@mortalityObsVar   <- object@mortalityObsVar[setdiff(names(object@mortalityObsVar),ages)]                  
    }
    #Finally, do an update
    object <- update(object)
    return(object)
})

#' Renumber the couplings of an FLSAM.control object
#'
#' Renumbers all coupling matrices and vectors of an
#' \code{\link{FLSAM.control-class}} object so that the parameter indices are
#' consecutive and start at zero, keeping \code{-1} for unused entries. Run
#' this after modifying couplings by hand.
#'
#' @param object An \code{FLSAM.control} object.
#' @return The renumbered \code{FLSAM.control} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' ctrl <- NSH.ctrl
#' ctrl@catchabilities["HERAS", as.character(5:8)] <- 101
#' update(ctrl)@catchabilities
#' @aliases update,FLSAM.control-method
#' @rdname update-FLSAM.control
#' @export
setMethod("update", signature(object="FLSAM.control"),
  function(object){
  object <- .upgradeControl(object)
  for(iSlt in slotNames(object)){
    if(is(slot(object,iSlt), "matrix")){
      isNA <- which(slot(object,iSlt)==-1)
      slot(object,iSlt)[isNA] <- NA
      slot(object,iSlt)[]  <-  as.numeric(factor(slot(object,iSlt)))-1
      slot(object,iSlt)[isNA] <- -1
    }
    if(iSlt %in% c("logN.vars","scalePars","logP.vars","stockWeightMean","stockWeightObsVar","catchWeightMean","catchWeightObsVar","maturityMean","mortalityMean","mortalityObsVar")){
      isNA <- which(slot(object,iSlt)==-1)
      slot(object,iSlt)[]  <- as.numeric(factor(slot(object,iSlt)))-1
      slot(object,iSlt)[isNA] <- -1
    }
  }
  return(object)}
)

setValidity("FLSAM.control",
  function(object){
    msg <- character()
    #-1 All slots populated
    if(any(object@plus.group==1) && is.na(object@range["plusgroup"]))
      msg <- c(msg,"Specify plusgroup in object@range")
    if(any(is.na(object@range[-which(names(object@range)=="plusgroup")])))
      msg <- c(msg,"Specify values in object@range")
    if(any(is.na(object@logN.vars)))
      msg <- c(msg,"Not all ages in logN.vars specified")

    #-2 Settings within range
    if(!object@srr %in% c(0:3,60:69,90:93,201,202,260,261,263,264,266,267,290,293,401,402,490,493))
      msg <- c(msg,"SRR type outside possible setting range")
    if(any(!object@fleets %in% c(0:3,5:7)))
      msg <- c(msg,paste("Fleet type of",paste(names(object@fleets)[!object@fleets %in% c(0:3,5:7)],collapse=", "),"outside possible setting range"))
    if(any(!object@biomassTreat %in% c(-1,0:6,10)))
      msg <- c(msg,"biomassTreat outside possible setting range (-1, 0-6, 10)")
    if(!length(object@fixVarToWeight) %in% c(1,length(object@fleets)))
      msg <- c(msg,"fixVarToWeight must have length 1 or one value per fleet")

    #-3 Dimension and dimnames check of the fleet x age coupling matrices
    fltAge <- c("states","catchabilities","power.law.exps","f.vars","obs.vars")
    dmns   <- lapply(fltAge,function(x) dimnames(slot(object,x)))
    if(!all(sapply(dmns,identical,dmns[[1]])))
      msg <- c(msg,paste("Dimnames are not the same for",paste(fltAge,collapse=", ")))
    if(!identical(rownames(object@cor.obs),rownames(object@states)) | ncol(object@cor.obs)!=(ncol(object@states)-1))
      msg <- c(msg,"Dimensions of cor.obs do not match the fleet x age-1 structure")

    if(length(msg)==0) return(TRUE) else return(msg)
  }
)

