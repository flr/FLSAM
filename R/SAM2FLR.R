SAM2FLR <- function(fit,ctrl){

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
  if(length(ctrl@logP.vars)>0){
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
    require(boot)
    res@mat <- FLQuant(t(inv.logit(fit$pl$logitMO[1:length(ctrl@range["minyear"]:ctrl@range["maxyear"]),])),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
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
                                                          
  if(length(grep("catch",rownames(ctrl@states)))==1)
    res@harvest <- FLQuant(t(faytable(fit)),
                                                        dimnames = list(age = ctrl@range["min"]:ctrl@range["max"],
                                                        year = ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                        unit = "unique", season = "all", area = "unique",
                                                        iter = 1))
  if(length(grep("catch",rownames(ctrl@states)))>1){
    defaultArray<- array(0,dim=c(length(ctrl@range["min"]:ctrl@range["max"]),
                                  length(ctrl@range["minyear"]:ctrl@range["maxyear"]),
                                  1,
                                  1,
                                  length(grep("catch",rownames(ctrl@states))),
                                  1),
                            dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                          year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                          unit="unique",
                                          season="all",
                                          area=gsub("catch ","",rownames(ctrl@states[grep("catch",rownames(ctrl@states)),])),
                                          iter=1))
    for(iArea in 1:dim(defaultArray)[5]){
      ages  <- ctrl@states[grep("catch",rownames(ctrl@states))[iArea],]
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
                                std.res=residuals(fit)$residual)
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
return(res)}

SIM2FLR <- function(sim,fit,ctrl){
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

    #- get parameters + sd value
    parssummary <- data.frame(name=ac(rownames(as.data.frame(unlist(fit$pl)))),
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

    if(length(grep("catch",rownames(ctrl@states)))==1)
      res@harvest <- FLQuant(exp(matrix(subset(res@params,name=="logF")$value,ncol=length(ctrl@range["minyear"]:ctrl@range["maxyear"])))[ctrl@states[grep("catch",rownames(ctrl@states)),]+1,],
                                                          dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                          year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                          unit="unique",
                                                          season="all",
                                                          area="unique",
                                                          iter=1))
    if(length(grep("catch",rownames(ctrl@states)))>1){
      stop("Multifleet simulate not implemented yet")
      defaultArray<- array(0,dim=c(length(ctrl@range["min"]:ctrl@range["max"]),
                                    length(ctrl@range["minyear"]:ctrl@range["maxyear"]),
                                    1,
                                    1,
                                    length(grep("catch",rownames(ctrl@states))),
                                    1),
                              dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                            year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                            unit="unique",
                                            season="all",
                                            area=gsub("catch ","",rownames(ctrl@states[grep("catch",rownames(ctrl@states)),])),
                                            iter=1))
      for(iArea in 1:dim(defaultArray)[5]){
        ages  <- ctrl@states[grep("catch",rownames(ctrl@states))[iArea],]
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







