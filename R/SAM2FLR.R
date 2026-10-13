#' Convert a SAM fit into an FLSAM object
#'
#' \code{SAM2FLR} converts a fitted \code{sam} object (from
#' \code{\link[stockassessment]{sam.fit}} or \code{FLSAM(..., return.fit=TRUE)})
#' into an \code{\link{FLSAM-class}} object. \code{SIM2FLR} converts a matrix
#' of simulated parameter vectors (see \code{\link{simulate}}) into an
#' \code{\link{FLSAMs}} object.
#'
#' @param fit A \code{sam} fit.
#' @param ctrl The \code{\link{FLSAM.control-class}} object used for the fit.
#' @param keep.fit Store \code{fit} in the \code{fit} slot of the result.
#' @param sim Matrix of simulated parameters (rows are simulations).
#' @param set.pars The \code{set.pars} used in the fit, if any.
#' @return An \code{FLSAM} (\code{SAM2FLR}) or \code{FLSAMs} (\code{SIM2FLR}) object.
#' @rdname SAM2FLR
#' @export
SAM2FLR <- function(fit,ctrl,keep.fit=TRUE){
  ctrl        <- .upgradeControl(ctrl)

  #- Create new FLSAM object and fill with ctrl elements
  res         <- new("FLSAM")
  res@control <- ctrl
  res@nopar   <- length(fit$sdrep$par.fixed)
  res@n.states<- as.integer(length(which(unique(c(ctrl@states))>=0))+ c(ctrl@range["max"] - ctrl@range["min"] + 1))
  res@states  <- ctrl@states
  res@nlogl   <- fit$opt$objective
  #- Get overall vcov and the cov of the observations
  res@vcov  <- fit$sdrep$cov.fixed
  if(length(which(fit$conf$keyCorObs>=0))>0){
    res@obscov  <- fit$rep$obsCov[which(apply(ctrl@cor.obs,1,function(x){any(x>=0,na.rm=T)}))]
    names(res@obscov) <- rownames(ctrl@cor.obs)[apply(ctrl@cor.obs,1,function(x){any(x>=0,na.rm=T)})]
  }
  #- get parameters + sd value
  parssummary <- rbind(data.frame(name=ac(rownames(as.data.frame(unlist(fit$pl)))),
                                  value=           as.vector(unlist(fit$pl)),
                                  std.dev=         as.vector(unlist(fit$plsd))),
                       data.frame(name=ac(names(fit$sdrep$value)),value=fit$sdrep$value,std.dev=fit$sdrep$sd))
  parssummary$name <- ac(parssummary$name)
  colnames(parssummary) <- c("name","value","std.dev")
  parssummary[,1]       <- gsub('[0-9]+', '',parssummary[,1])
  res@params  <- parssummary

  rescov    <- sdreport(fit$obj, fit$opt$par)
  res@rescov<- rescov$cov
  paramname<- names(rescov$value)
  rownames(res@rescov) <- paramname
  colnames(res@rescov) <- paramname

  #- Get component information
  if(length(ctrl@logP.vars) > 0){
    res@components <- exp(fit$rep$comps)
    compnames <- names(ctrl@fleets)[which(ctrl@fleets==6)]
    compyears <- range(fit$data$aux[which(fit$data$aux[,2] %in% c(which(ctrl@fleets==6))),1])
    dimnames(res@components) <- list(component=compnames,years=compyears[1]:compyears[2])
  }

  #- Get stock weight information
  if(ctrl@stockWeightModel){
    res@stock.wt <- FLQuant(t(exp(fit$pl$logSW[1:length(ctrl@range["minyear"]:ctrl@range["maxyear"]),])),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                        year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                        unit="unique",
                                                        season="all",
                                                        area="unique",
                                                        iter=1))
  }

  #- Get catch weight information
  if(ctrl@catchWeightModel){
    res@catch.wt <- FLQuant(aperm(exp(fit$pl$logCW[1:length(ctrl@range["minyear"]:ctrl@range["maxyear"]),,]),c(2,1,3)),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                        year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                        unit="unique",
                                                        season="all",
                                                        area=names(ctrl@fleets)[grep("catch",names(ctrl@fleets))],
                                                        iter=1))
  }

  #- Get maturity information
  if(ctrl@maturityModel){
    res@mat <- FLQuant(t(plogis(fit$pl$logitMO[1:length(ctrl@range["minyear"]:ctrl@range["maxyear"]),])),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                        year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                        unit="unique",
                                                        season="all",
                                                        area="unique",
                                                        iter=1))
  }  
  #- Get mortality information
  if(ctrl@mortalityModel){
    res@m <- FLQuant(t(exp(fit$pl$logNM[1:length(ctrl@range["minyear"]:ctrl@range["maxyear"]),])),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                        year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                        unit="unique",
                                                        season="all",
                                                        area="unique",
                                                        iter=1))
  }    


  #- Fill stock and harvest
  res@stock.n <- FLQuant(t(ntable(fit)),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                        year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                        unit="unique",
                                                        season="all",
                                                        area="unique",
                                                        iter=1))
                                                          
  if(length(which(ctrl@fleets==0))==1)
    res@harvest <- FLQuant(t(faytable(fit)),
                                                        dimnames = list(age = ctrl@range["min"]:ctrl@range["max"],
                                                        year = ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                        unit = "unique", season = "all", area = "unique",
                                                        iter = 1))
  if(length(which(ctrl@fleets==0))>1){
    defaultArray<- array(0,dim=c(length(ctrl@range["min"]:ctrl@range["max"]),
                                  length(ctrl@range["minyear"]:ctrl@range["maxyear"]),
                                  1,
                                  1,
                                  length(which(ctrl@fleets==0)),
                                  1),
                            dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                          year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                          unit="unique",
                                          season="all",
                                          area=gsub("catch ","",rownames(ctrl@states[which(ctrl@fleets==0),])),
                                          iter=1))
    for(iArea in 1:dim(defaultArray)[5]){
      ages  <- ctrl@states[which(ctrl@fleets==0)[iArea],]
      ages  <- ages[which(ages>=0)]
      defaultArray[names(ages),,,,iArea,] <- exp(fit$pl$logF[ages+1,])
    }
    res@harvest <- FLQuant(defaultArray)
    if(!all(round(c(areaSums(res@harvest)[,drop=T]/fit$rep$totF),4)==1)) warning("Partial Fs not equal to total F")
  }


  units(res@harvest) <- "f"
  #- Get all the residuals
  if(ctrl@residuals){
    res@residuals <- data.frame(year=fit$data$aux[,"year"],
                                fleet=names(ctrl@fleets)[fit$data$aux[,"fleet"]],
                                age=fit$data$aux[,"age"],
                                log.obs=fit$data$logobs,
                                log.mdl=fit$rep$predObs,
                                std.res=.samResiduals(fit))
    res@residuals$log.obs[which(is.na(res@residuals$log.obs))] <- fit$pl$missing
    res@residuals$log.mdl[which(!is.finite(res@residuals$log.mdl))] <- NA
    res@residuals$std.res[which(!is.finite(res@residuals$std.res))] <- NA
  } else {
    res@residuals <- data.frame(year=fit$data$aux[,"year"],
                                fleet=names(ctrl@fleets)[fit$data$aux[,"fleet"]],
                                age=fit$data$aux[,"age"],
                                log.obs=fit$data$logobs,
                                log.mdl=fit$rep$predObs,
                                std.res=NA)
    if(!is.null(fit$pl$missing))
      res@residuals$log.obs[which(is.na(res@residuals$log.obs))] <- fit$pl$missing
    res@residuals$log.mdl[which(!is.finite(res@residuals$log.mdl))] <- NA
  }

  #- Tag-recapture observations: add the recapture year
  if(any(fit$data$fleetTypes==5))
    res@residuals$RecaptureY <- fit$data$aux[,"RecaptureY"]

  #- Finally some info and save
  info        <- data.frame(FLSAM.version=packageDescription("FLSAM")$Version,
                            FLCore.version=packageDescription("FLCore")$Version,
                            R.version=R.version$version.string,
                            platform=R.version$platform,
                            run.date=Sys.time())
  res@info    <- t(info)
  res@name    <- ctrl@name
  res@desc    <- ctrl@desc
  res@range   <- ctrl@range
  if(keep.fit) res@fit <- fit
return(res)}

#' @rdname SAM2FLR
#' @export
SIM2FLR <- function(sim,fit,ctrl,set.pars=NULL){
  resList       <-list()
  for(i in 1:nrow(sim)){
    #- Create new FLSAM object and fill with ctrl elements
    res         <- new("FLSAM")
    res@control <- ctrl
    res@nopar   <- length(fit$sdrep$par.fixed)
    res@n.states<- as.integer(length(which(unique(c(ctrl@states))>=0))+ c(ctrl@range["max"] - ctrl@range["min"] + 1))
    res@states  <- ctrl@states
    res@nlogl   <- fit$opt$objective
    #- Get overall vcov and the cov of the observations
    res@vcov    <- matrix(0)

    #- get parameter names
    pars2match <- unique(dimnames(sim)[[2]])
    pl          <- fit$pl   #- copy per simulation: the fixed parameters are removed below
    keepPars    <- names(pl)[which(names(pl) %in% pars2match)]

    #- Fix for set.pars: remove the fixed parameters (all at once per parameter)
    if(!is.null(set.pars)){
      set.parsSAM <- .setParsSAM(pl,set.pars)
      for(iPar in unique(set.parsSAM$SAM))
        pl[[iPar]] <- pl[[iPar]][-(set.parsSAM$number[set.parsSAM$SAM==iPar]+1)]
    }

    pars       <- ac(rownames(as.data.frame(unlist(pl))))
    idx       <- character()
    for(iK in keepPars){
      if(iK != "logAlphaSCB"){
        grp <- pars[grepl(iK,pars)]
        if(length(grp)>1){
          idx <- c(idx,grp[which(!is.na(an(substr(grp,nchar(iK)+1,nchar(grp)))))])
        } else {
          idx <- c(idx,grp)
        }
      }
      if(iK == "logAlphaSCB"){
        grp <- pars[grepl(iK,pars)]
        idx <- c(idx,grp)
      }
    }
    #- add names + values
    parssummary <- data.frame(name=idx,
                                    value=           sim[i,],
                                    std.dev=         NA)
    parssummary$name <- ac(parssummary$name)
    colnames(parssummary) <- c("name","value","std.dev")
    parssummary[,1]       <- gsub('[0-9]+', '',parssummary[,1])
    res@params  <- parssummary

    res@rescov<- matrix(0)

    #- Fill stock and harvest
    res@stock.n <- FLQuant(exp(subset(res@params,name=="logN")$value),
                                                          dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                          year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                          unit="unique",
                                                          season="all",
                                                          area="unique",
                                                          iter=1))

    if(length(which(ctrl@fleets==0))==1)
      res@harvest <- FLQuant(exp(matrix(subset(res@params,name=="logF")$value,ncol=length(ctrl@range["minyear"]:ctrl@range["maxyear"])))[ctrl@states[which(ctrl@fleets==0),]+1,],
                                                          dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                          year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                          unit="unique",
                                                          season="all",
                                                          area="unique",
                                                          iter=1))
    if(length(which(ctrl@fleets==0))>1){
      stop("Multifleet simulate not implemented yet")
      defaultArray<- array(0,dim=c(length(ctrl@range["min"]:ctrl@range["max"]),
                                    length(ctrl@range["minyear"]:ctrl@range["maxyear"]),
                                    1,
                                    1,
                                    length(which(ctrl@fleets==0)),
                                    1),
                              dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                            year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                            unit="unique",
                                            season="all",
                                            area=gsub("catch ","",rownames(ctrl@states[which(ctrl@fleets==0),])),
                                            iter=1))
      for(iArea in 1:dim(defaultArray)[5]){
        ages  <- ctrl@states[which(ctrl@fleets==0)[iArea],]
        ages  <- ages[which(ages>=0)]
        defaultArray[names(ages),,,,iArea,] <- exp(fit$pl$logF[ages+1,])
      }
      res@harvest <- FLQuant(defaultArray)
      if(!all(round(c(areaSums(res@harvest)[,drop=T]/fit$rep$totF),4)==1)) warning("Partial Fs not equal to total F")
    }


    units(res@harvest) <- "f"
    #- Finally some info and save
    info        <- data.frame(FLSAM.version=packageDescription("FLSAM")$Version,
                              FLCore.version=packageDescription("FLCore")$Version,
                              R.version=R.version$version.string,
                              platform=R.version$platform,
                              run.date=Sys.time())
    res@info    <- t(info)
    res@name    <- ctrl@name
    res@desc    <- ctrl@desc
    res@range   <- ctrl@range
    resList[[i]]<- res
  }
  res <- as(resList,"FLSAMs")
return(res)}

#- One-observation-ahead residuals; tag-recapture observations (fleet type 5)
#  are discrete and get conditional residuals with oneStepGeneric
.samResiduals <- function(fit){
  dis <- fit$data$fleetTypes[fit$data$aux[,"fleet"]]==5
  if(!any(dis)) return(residuals(fit)$residual)
  res <- rep(NA_real_,length(dis))
  res[!dis] <- residuals(fit,subset=which(!dis))$residual
  res[dis]  <- residuals(fit,discrete=TRUE,conditional=which(!dis),subset=which(dis),
                         method="oneStepGeneric",range=c(0,Inf))$residual
  return(res)
}
