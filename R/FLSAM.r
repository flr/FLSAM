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

#FLSAM2parameters <- function(sam


FLSAM.MSE <-function(stcks,tun,ctrl,catch.vars=NULL,starting.values=NULL,...){

  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(class(stcks)=="FLStocks") stop("Not implemented for multi-fleet yet")
  if(class(stcks)=="FLStock") stcks <- FLStocks(residual=stcks)
  if(any(unlist(lapply(stcks,function(x)dims(x)$iter))<=1) & any(unlist(lapply(tun,function(x)dims(x)$iter))<=1))
    stop("Running in MSE mode means supplying FLStock and FLIndices with more iters than 1. Use FLSAM() instead")

  #Count iters
  iters <- unlist(lapply(stcks,function(x)dims(x)$iter))

  #get datasets ready
  data  <- conf <- par <- list()
  for(i in 1:iters){
    iTun <- tun
    for(j in 1:length(tun))
      iTun[[j]]<- iter(iTun[[j]],i)
    data[[i]]  <- FLSAM2SAM(FLStocks("residual"=iter(stcks[["residual"]],i)),iTun,ctrl@sumFleets,catch.vars)
    conf[[i]]  <- ctrl2conf(ctrl,data[[i]])
    par[[i]]   <- stockassessment::defpar(data[[i]],conf[[i]])
  }

  require(doParallel)
  cl <- makeCluster(detectCores()-1) #set up nodes
  clusterEvalQ(cl,library(FLSAM))
  clusterEvalQ(cl,library(stockassessment))
  registerDoParallel(cl)

  if(is.null(starting.values)){
    fit   <- stockassessment::sam.fit(data[[1]],conf[[1]],par[[1]],sim.condRE=ctrl@simulate,...)  #38.76
    starting.values <- fit$pl[-grep("missing",names(fit$pl))]
  }
  # use fit as starting point for sam.fitfast
  for(i in 2:iters)
    par[[i]] <- starting.values

  res <- foreach(i = 2:iters) %dopar% try(sam.fitfast(data[[i]],conf[[i]],par[[i]],silent=T)) #14.14
  stopCluster(cl)
  detach("package:doParallel",unload=TRUE)
  detach("package:foreach",unload=TRUE)
  for(i in 1:iters){
    if(i == 1){
      stcks[["residual"]]@stock.n[,,,,,i] <- exp(fit$rep$logN[,1:dims(stcks[["residual"]]@stock.n)$year])
      stcks[["residual"]]@harvest[,,,,,i] <- exp(fit$rep$totF[,1:dims(stcks[["residual"]]@harvest)$year])
    } else {
      stcks[["residual"]]@stock.n[,,,,,i] <- exp(res[[i-1]]$logN[,1:dims(stcks[["residual"]]@stock.n)$year])
      stcks[["residual"]]@harvest[,,,,,i] <- exp(res[[i-1]]$totF[,1:dims(stcks[["residual"]]@harvest)$year])
    }
  }
  return(stcks)
}

getLowerBounds<-function(parameters){
    list(sigmaObsParUS=rep(-10,length(parameters$sigmaObsParUS)))
}

getUpperBounds<-function(parameters){
    list(sigmaObsParUS=rep(10,length(parameters$sigmaObsParUS)))
}

clean.void.catches<-function(dat, conf){
  cfidx <- which(dat$fleetTypes==0)
  aidx <- unique(dat$aux[dat$aux[,2]%in%cfidx,3]-conf$minAge+1)
  faidx <- as.matrix(expand.grid(cfidx, aidx))
  faidx <- faidx[which(conf$keyLogFsta[faidx]== -1),,drop=FALSE]
  rmidx <- paste0(dat$aux[,2],"x",dat$aux[,3]-conf$minAge+1) %in%  paste0(faidx[,1],"x",faidx[,2])
  dat$aux <- dat$aux[!rmidx,]
  dat$logobs <- dat$logobs[!rmidx]
  dat$weight <- dat$weight[!rmidx]
  dat$nobs<-sum(!rmidx)
  dat$minAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=min))
  dat$maxAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=max))
  newyear<-min(as.numeric(dat$aux[,"year"])):max(as.numeric(dat$aux[,"year"]))
  newfleet<-min(as.numeric(dat$aux[,"fleet"])):max(as.numeric(dat$aux[,"fleet"]))
  mmfun<-function(f,y, ff){idx<-which(dat$aux[,"year"]==y & dat$aux[,"fleet"]==f); ifelse(length(idx)==0, NA, ff(idx)-1)}
  dat$idx1<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=min)
  dat$idx2<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=max)
  dat
}


sam.fitfast <- function(data, conf, parameters, lower=getLowerBounds(parameters), upper=getUpperBounds(parameters), sim.condRE=TRUE, ...){

  data <- clean.void.catches(data,conf)
  tmball <- c(data, conf, simFlag=as.numeric(sim.condRE))
  if(is.null(tmball$resFlag)){tmball$resFlag <- 0}
  nmissing <- sum(is.na(data$logobs))
  parameters$missing <- numeric(nmissing)
  if(length(conf$keyVarLogP)>1){
    ran <- c("logN", "logF","logPS", "missing")
    obj <- TMB::MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",...)
  } else {
    ran <- c("logN", "logF", "missing")
    obj <- TMB::MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",...)
  }

  lower2<-rep(-Inf,length(obj$par))
  upper2<-rep(Inf,length(obj$par))
  for(nn in names(lower)) lower2[names(obj$par)==nn]=lower[[nn]]
  for(nn in names(upper)) upper2[names(obj$par)==nn]=upper[[nn]]

  opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower2,upper=upper2)
  rep <- obj$report()

return(c(opt,rep))}


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
