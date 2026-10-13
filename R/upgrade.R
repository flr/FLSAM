#' Upgrade FLSAM objects created with an older FLSAM version
#'
#' Objects saved with an earlier version of FLSAM lack slots that were added
#' later (e.g. \code{hockeyStickCurve} in \code{\link{FLSAM.control-class}} and
#' \code{fit} in \code{\link{FLSAM-class}}). \code{upgradeFLSAM} adds these
#' slots with their default values. \code{\link{FLSAM}} and
#' \code{\link{ctrl2conf}} do this automatically for control objects.
#'
#' @param object An \code{FLSAM.control}, \code{FLSAM} or \code{FLSAMs} object.
#' @return The upgraded object.
#' @examples
#' data(NSH); data(NSH.sam)
#' NSH.ctrl <- upgradeFLSAM(NSH.ctrl)
#' NSH.sam  <- upgradeFLSAM(NSH.sam)
#' @export
upgradeFLSAM <- function(object){
  if(is(object,"FLSAMs")){
    object@.Data <- lapply(object@.Data,upgradeFLSAM)
    return(object)
  }
  if(is(object,"FLSAM.control"))
    return(.upgradeControl(object))
  if(is(object,"FLSAM")){
    object <- .upgradeSlots(object,"FLSAM")
    object@control <- .upgradeControl(object@control)
    return(object)
  }
  stop("object must be of class FLSAM.control, FLSAM or FLSAMs")
}

#- Copy all available slots into a fresh object of the current class definition
.upgradeSlots <- function(object,cls){
  slts <- slotNames(cls)
  if(all(sapply(slts,function(x) .hasSlot(object,x)))) return(object)
  res  <- new(cls)
  for(iSlt in slts)
    if(.hasSlot(object,iSlt)) slot(res,iSlt,check=FALSE) <- slot(object,iSlt)
  return(res)
}

.upgradeControl <- function(ctrl){
  ctrl <- .upgradeSlots(ctrl,"FLSAM.control")
  #- fixVarToWeight is fleet specific since SAM 0.12.0
  if(length(ctrl@fixVarToWeight)==1 & length(ctrl@fleets)>1)
    ctrl@fixVarToWeight <- rep(ctrl@fixVarToWeight,length(ctrl@fleets))
  if(length(ctrl@fixVarToWeight)==length(ctrl@fleets))
    names(ctrl@fixVarToWeight) <- names(ctrl@fleets)
  return(ctrl)
}

#' Access the underlying stockassessment fit
#'
#' An \code{\link{FLSAM-class}} object stores the \code{sam} fit of the
#' stockassessment package in its \code{fit} slot (unless
#' \code{FLSAM(..., keep.fit=FALSE)} was used). \code{samfit} returns it, so
#' that all functions of the stockassessment package can be used, e.g.
#' \code{stockassessment::modeltable(samfit(sam))}.
#'
#' @param object An \code{FLSAM} or \code{FLSAMs} object.
#' @return A \code{sam} fit, or for \code{FLSAMs} a \code{samset} (list of fits).
#'   An error is raised when no fit is stored.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)
#' fit <- samfit(sam)
#' stockassessment::ssbplot(fit)
#' }
#' @rdname samfit
#' @export
setGeneric("samfit", function(object) standardGeneric("samfit"))
#' @rdname samfit
setMethod("samfit", signature(object="FLSAM"),
  function(object){
    fit <- if(.hasSlot(object,"fit")) object@fit else NULL
    if(is.null(fit))
      stop(paste("No stockassessment fit stored in this FLSAM object.",
                 "Rerun FLSAM() with keep.fit=TRUE (the default)."))
    return(fit)
})
#' @rdname samfit
setMethod("samfit", signature(object="FLSAMs"),
  function(object){
    res <- lapply(object,samfit)
    class(res) <- "samset"
    return(res)
})

#- Printing objects saved with older FLSAM versions (missing slots) upgrades them first
#' @rdname upgradeFLSAM
#' @aliases show,FLSAM.control-method
setMethod("show",signature(object="FLSAM.control"),
  function(object) methods::showDefault(.upgradeControl(object)))

#' @rdname upgradeFLSAM
#' @aliases show,FLSAM-method
setMethod("show",signature(object="FLSAM"),
  function(object){
    methods::showDefault(upgradeFLSAM(object))
})

#' @rdname upgradeFLSAM
#' @param x An \code{FLSAM} object.
#' @param ... Not used.
#' @aliases print,FLSAM-method
setMethod("print",signature(x="FLSAM"),
  function(x,...) show(x))
