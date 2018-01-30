setClass("FLSAM",
         representation(
           "FLComp",
           control    = "FLSAM.control",
           nopar      = "integer",
           n.states   = "integer",
           states     = "matrix",
           components = "matrix",
           nlogl      = "numeric",
           vcov       = "matrix",
           rescov     = "matrix",
           obscov     = "list",
           params     = "data.frame",
           stock.n    = "FLQuant",
           harvest    = "FLQuant",
           residuals  = "data.frame",
           info       = "matrix"),
         prototype=prototype(),
         validity=function(object){
           # Everything is fine
           return(TRUE)}
)

FLSAM <-function(stcks,tun,ctrl,catch.vars=NULL,return.fit=F,starting.values=NULL,...){
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(class(stcks)=="FLStock") stcks <- FLStocks(residual=stcks)
  data  <- FLSAM2SAM(stcks,tun,ctrl@sumFleets,catch.vars)
  conf  <- ctrl2conf(ctrl,data)
  par   <- defpar(data,conf)
  if(!is.null(starting.values))
    par <- updateStart(par,starting.values)

  fit   <- sam.fit(data,conf,par,sim.condRE=ctrl@simulate,...)

  #- Check if you want variance and covariance info
  #- Get FLSAM object back
  if(return.fit){
    res <- fit
  } else {
    res <- try(SAM2FLR(fit,ctrl))
  }

  return(res)
}

# class
setClass("FLSAMs", contains="FLlst")

# constructor
setGeneric("FLSAMs", function(object, ...){
	standardGeneric("FLSAMs")
	}
)

setMethod("FLSAMs", signature(object="ANY"), function(object, ...){
	lst1 <- list(...)
	nlst <- length(lst1)
	lst <- list()
	length(lst) <- nlst + 1
	lst[[1]] <- object
	lst[-1] <- lst1
	FLSAMs(lst)
})

setMethod("FLSAMs", signature(object="FLSAM"), function(object, ...) {
    lst <- c(object, list(...))
    FLSAMs(lst)
})

setMethod("FLSAMs", "missing", function(...){
	if(missing(...)){
		new("FLSAMs")
	} else {
		lst <- list(...)
		FLSAMs(lst)
	}
})

setMethod("FLSAMs", "list", function(object){
        #Check list contents are FLSAM objects
        for(i in seq(object)) {
          if(!is(object[[i]],"FLSAM")) {
            stop(sprintf("Object %i in list is not an FLSAM",i)) }} 
        #Setup unique names
        list.names <- names(object)
        obj.names <- sapply(object,name)
        outnames <- if(is.null(list.names)) {obj.names} else {list.names}
        if(any(table(outnames) !=1) | any(outnames=="") |is.null(outnames)) outnames <- as.character(seq(object))
	new("FLSAMs", .Data=object,names=outnames)
})

setMethod("FLSAMs", "FLSAMs", function(object){
	return(object)
}) # }}}
