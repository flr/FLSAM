setClass("FLSAM",
         representation(
           "FLComp",
           control    = "FLSAM.control",
           nopar      = "integer",
           n.states   = "integer",
           states     = "matrix",
           components = "matrix",
           stock.wt   = "FLQuant",
           catch.wt   = "FLQuant",
           mat        = "FLQuant",
           m          = "FLQuant",
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

FLSAM <-function(stcks,tun,ctrl,catch.vars=NULL,return.fit=F,starting.values=NULL,set.pars=NULL,...){
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(is(stcks, "FLStock")) stcks <- FLStocks(residual=stcks)
  dataS  <- FLSAM2SAM(stcks,tun,ctrl@sumFleets,catch.vars)
  confS  <- ctrl2conf(ctrl,dataS)
  parS   <- defpar(dataS,confS)
  if(!is.null(starting.values))
    parS <- updateStart(parS,FLSAM2par(starting.values))

   if(!is.null(set.pars)){
    
    matchNames <- matrix(c("logN.vars","logSdLogN",
                           "logP.vars","logSdLogP",
                           "catchabilities","logFpar",
                           "power.law.exps","logQpow",
                           "f.vars","logSdLogFsta",
                           "obs.vars","logSdLogObs"),ncol=2,byrow=T,dimnames=list(1:6,c("FLSAM","SAM")))
    if(any(!set.pars$par %in% matchNames[,"FLSAM"]))
        warning(paste(set.pars$par[which(!set.pars$par %in% matchNames[,"FLSAM"])],"cannot be set"))
    set.pars <- merge(set.pars,matchNames,by.x="par",by.y="FLSAM")
    

    mapDef        <- as.list(matchNames[,"SAM"]); names(mapDef) <- matchNames[,"SAM"]
    for(i in names(mapDef))
      mapDef[[i]] <- parS[[i]]
    map           <- mapDef

    for(i in 1:nrow(set.pars)){
      parS[[set.pars$SAM[i]]][set.pars$number[i]+1] <- set.pars$value[i]
      map[[set.pars$SAM[[i]]]][set.pars$number[i]+1] <- NA
    }
    map=list(logSdLogN=as.factor(map$logSdLogN),
             logSdLogP=as.factor(map$logSdLogP),
             logFpar=as.factor(map$logFpar),
             logQpow=as.factor(map$logQpow),
             logSdLogFsta=as.factor(map$logSdLogFsta),
             logSdLogObs=as.factor(map$logSdLogObs))
    map <- map[which(names(map) %in% set.pars$SAM)]
    fit   <- sam.fit(dataS,confS,parS,map=map,sim.condRE=ctrl@simulate,...)
  } else {
    fit   <- sam.fit(dataS,confS,parS,sim.condRE=ctrl@simulate,...)
  }

  #- Check if you want variance and covariance info
  #- Get FLSAM object back
  if(return.fit){
    res <- fit
  } else {
    res <- try(SAM2FLR(fit,ctrl))
  }

  return(res)
}


FLSAM.MSE <-function(stcks,tun,ctrl,catch.vars=NULL,starting.sam=NULL,return.sam=F,...){

  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(is(stcks, "FLStocks")) stop("Not implemented for multi-fleet yet")
  if(is(stcks, "FLStock")) stcks <- FLStocks(residual=stcks)
  if(any(unlist(lapply(stcks,function(x)dims(x)$iter))<=1) & any(unlist(lapply(tun,function(x)dims(x)$iter))<=1))
    stop("Running in MSE mode means supplying FLStock and FLIndices with more iters than 1. Use FLSAM() instead")

  #Count iters
  iters <- unlist(lapply(stcks,function(x)dims(x)$iter))

  #Turn residuals off
  ctrl@residuals <- FALSE

  #get datasets ready
  data  <- conf <- par <- list()

  require(doParallel)
  ncores <- detectCores()-1
  ncores <- ifelse(iters<ncores,iters,ncores)
  cl <- makeCluster(ncores) #set up nodes
  clusterEvalQ(cl,library(FLSAM))
  clusterEvalQ(cl,library(stockassessment))
  registerDoParallel(cl)


  data <- foreach(i = 1:iters) %dopar% FLSAM2SAM(FLStocks("residual"=iter(stcks[["residual"]],i)),FLIndices(lapply(tun, function(x) iter(x,i))),ctrl@sumFleets,catch.vars)
  conf <- foreach(i = 1:iters) %dopar% ctrl2conf(ctrl,data[[i]])
  par  <- foreach(i = 1:iters) %dopar% defpar(data[[i]],conf[[i]])
  checkUpdate <- function(i,iSam,iPar){
                 if(is(iSam, "logical")){
                   ret <- updateStart(iPar,FLSAM2par(iSam)) } else {
                   ret <- iPar }
                 return(ret)}
  if(!is.null(starting.sam))
    par <- foreach(i = 1:iters) %dopar% checkUpdate(i,starting.sam,par[[i]])

#  for(i in 1:iters){
#    iTun <- tun
#    for(j in 1:length(tun))
#      iTun[[j]]<- iter(iTun[[j]],i)
#    data[[i]]  <- FLSAM2SAM(FLStocks("residual"=iter(stcks[["residual"]],i)),iTun,ctrl@sumFleets,catch.vars)
#    conf[[i]]  <- ctrl2conf(ctrl,data[[i]])
#    par[[i]]   <- stockassessment::defpar(data[[i]],conf[[i]])
#    if(!is.null(starting.sam)){
#      if(class(starting.sam[[i]])!="logical")
#        par[[i]] <- updateStart(par[[i]],FLSAM2par(starting.sam[[i]]))
#    }
#  }


  #- First run without staring values simply because it's quicker
#  if(!force.starting.sam){
#    system.time(res <- foreach(i = 1:iters) %dopar% try(sam.fitfast(data[[i]],conf[[i]],defpar(data[[i]],conf[[i]]),silent=T))),...))
#  } else {
  res <- foreach(i = 1:iters) %dopar% try(sam.fitfast(data[[i]],conf[[i]],par[[i]],silent=T,...))
#  }

  #- Check for failed runs and do those with starting conditions
#  failed <- which(is.na(unlist(lapply(res,function(x){return(unlist(x$sdrep)[1])}))))
#  if(length(failed)>0 & !is.null(starting.sam)){
#    resFailed <- foreach(i = failed) %dopar% try(sam.fitfast(data[[i]],conf[[i]],par[[i]],silent=T,...))
#    res[failed] <- resFailed
#  }

  #- Return sam objects
  if(return.sam){
    resSAM <- list()
    for(i in 1:iters){
      if(!is.na(unlist(res[[i]]$sdrep)[1])){
        resSAM[[i]] <- SAM2FLR(res[[i]],ctrl)
      } else {
        resSAM[[i]] <- NA
      }
    }
    resSAM <- as(resSAM,"FLSAMs")
  }
  stopCluster(cl)
  
  if(!return.sam){
    for(i in 1:iters){
      if(!is.na(unlist(res[[i]]$sdrep)[1])){
        stcks[["residual"]]@stock.n[,,,,,i] <- exp(res[[i]]$rep$logN[,1:dims(stcks[["residual"]]@stock.n)$year])
        stcks[["residual"]]@harvest[,,,,,i] <- res[[i]]$rep$totF[,1:dims(stcks[["residual"]]@harvest)$year]
      } else {
        stcks[["residual"]]@stock.n[,,,,,i] <- NA
        stcks[["residual"]]@harvest[,,,,,i] <- NA
      }
    }
  }
  if(return.sam)
    ret <- resSAM
  if(!return.sam)
    ret <- stcks
  return(ret)
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
    obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",...)
  } else {
    ran <- c("logN", "logF", "missing")
    obj <- MakeADFun(tmball, parameters, random=ran, DLL="stockassessment",...)
  }

  lower2<-rep(-Inf,length(obj$par))
  upper2<-rep(Inf,length(obj$par))
  for(nn in names(lower)) lower2[names(obj$par)==nn]=lower[[nn]]
  for(nn in names(upper)) upper2[names(obj$par)==nn]=upper[[nn]]

  opt <- try(nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=0, eval.max=1000, iter.max=500),lower=lower2,upper=upper2))
  rep <- obj$report()
#  sdrep <- TMB::sdreport(obj,opt$par)
#
#  # Last two states
#  idx <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
#  sdrep$estY <- sdrep$value[idx]
#  sdrep$covY <- sdrep$cov[idx,idx]
#
#  idx <- c(which(names(sdrep$value)=="beforeLastLogN"),which(names(sdrep$value)=="beforeLastLogF"))
#  sdrep$estYm1 <- sdrep$value[idx]
#  sdrep$covYm1 <- sdrep$cov[idx,idx]
#
#  pl <- as.list(sdrep,"Est")
#  plsd <- as.list(sdrep,"Std")
#
#  #sdrep$cov<-NULL # save memory
#
#  fit <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=conf, opt=opt, obj=obj, rep=rep, low=lower2, hig=upper2)

  if(!is(opt, "try-error")){
    sdrep <- list(par.fixed=length(opt$par),cov.fixed=matrix(NA),cov=matrix(NA))
    pl    <- c(lapply(split(opt$par,names(opt$par)),unname),list(logPS=rep$logP),list(logN=rep$logN),list(logF=log(rep$totF[-which(duplicated(conf$keyLogFsta[1,])),])))
  } else {
    sdrep <- NA
    pl    <- NA
  }
  fit <- list(sdrep=sdrep, pl=pl, plsd=NA, data=data, conf=conf, opt=opt, obj=obj, rep=rep, low=lower2, hig=upper2)
return(fit)}


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
