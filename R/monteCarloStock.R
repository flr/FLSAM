monteCarloStock <- function(stck,tun,sam,realisations,return.sam=FALSE,saveParsDir=NULL,saveSimdataDir=NULL,set.pars=NULL,ncores=detectCores()-1, ...){
  
  require(doParallel)
  ctrl              <- sam@control
  
  #-Create new stock object with nr of realisations in iter slot
  if(class(stck)=="FLStocks"){
    idxStck           <- which.max(unlist(lapply(stck,function(x){x@range["maxyear"]})))
    mcstck            <- propagate(stck[[idxStck]],iter=realisations)
    mcstck            <- window(mcstck,start=range(sam)["minyear"],end=range(sam)["maxyear"])
    mcstck@stock.n[]  <- NA
    mcstck@harvest[]  <- NA
  } else {
    mcstck            <- propagate(stck,iter=realisations)
    mcstck            <- window(mcstck,start=range(sam)["minyear"],end=range(sam)["maxyear"])
    mcstck@stock.n[]  <- NA
    mcstck@harvest[]  <- NA
  }
  
  #- Run first an assessment and generate new data to fit new assessment on
  sam@control@simulate <- TRUE
  sam@control@residuals<- FALSE
  if(!is.null(set.pars)){
    if(class(stck)=="FLStock") stcks <- FLStocks(residual=stck)
    dataS  <- FLSAM2SAM(stcks,tun,ctrl@sumFleets,catch.vars)
    confS  <- ctrl2conf(ctrl,dataS)
    parS   <- defpar(dataS,confS)
    matchNames <- matrix(c("logN.vars","logSdLogN",
                           "logP.vars","logSdLogP",
                           "catchabilities","logFpar",
                           "power.law.exps","logQpow",
                           "f.vars","logSdLogFsta",
                           "obs.vars","logSdLogObs"),ncol=2,byrow=T,dimnames=list(1:6,c("FLSAM","SAM")))
    if(any(!set.pars$par %in% matchNames[,"FLSAM"]))
      warning(paste(set.pars$par[which(!set.pars$par %in% matchNames[,"FLSAM"])],"cannot be set"))
    set.pars.internal <- merge(set.pars,matchNames,by.x="par",by.y="FLSAM")
    
    
    mapDef        <- as.list(matchNames[,"SAM"]); names(mapDef) <- matchNames[,"SAM"]
    for(i in names(mapDef))
      mapDef[[i]] <- jitter(parS[[i]],factor=0.1)
      if(any(duplicated(mapDef[[i]])))
        stop("mapping goes wrong, some parameters get the same starting value, retry")
    map           <- mapDef
    
    for(i in 1:nrow(set.pars.internal)){
      parS[[set.pars.internal$SAM[i]]][set.pars.internal$number[i]+1] <- set.pars.internal$value[i]
      map[[set.pars.internal$SAM[[i]]]][set.pars.internal$number[i]+1] <- NA
    }
    map=list(logSdLogN=as.factor(map$logSdLogN),
             logSdLogP=as.factor(map$logSdLogP),
             logFpar=as.factor(map$logFpar),
             logQpow=as.factor(map$logQpow),
             logSdLogFsta=as.factor(map$logSdLogFsta),
             logSdLogObs=as.factor(map$logSdLogObs))
    map <- map[which(names(map) %in% set.pars$SAM)]
    
    object            <- FLSAM(stcks,tun,sam@control,set.pars=set.pars,return.fit=T)
  } else {
    object            <- FLSAM(stck,tun,sam@control,return.fit=T)
  }
  simdat            <- replicate(realisations,
                                 c(object$data[names(object$data)!="logobs"],#all the old data
                                   object$obj$simulate()["logobs"]),#simulated observations
                                 simplify=FALSE)
  #- Make cluster to speed up process
  ncores <- ifelse(realisations<ncores,realisations,ncores)
  cl <- makeCluster(ncores)
  clusterEvalQ(cl,library(FLSAM))
  clusterEvalQ(cl,library(stockassessment))
  registerDoParallel(cl)
  
  #- Function to run a new assessment on each of the simulated datasets
  for(i in 1:realisations){
    if(length(which(is.na(object$data$logobs)))>0)
      simdat[[i]]$logobs[which(is.na(object$data$logobs))] <- NA
  }
  if(!is.null(saveSimdataDir))
    save(simdat,file=saveSimdataDir)

  #clusterExport(cl,varlist=c("simdat","object"),envir=environment())
  if(!is.null(set.pars)){
    runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdat[[i]],object$conf,object$pl,map=map,silent=T,...))
  } else {
    runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdat[[i]],object$conf,object$pl,silent=T,...))
  }
  
  stopCluster(cl) #shut it down
  
  #- Fill the results of the simulations
  samRuns <- list()
  for(i in 1:realisations){
    if(!is.na(unlist(runs[[i]]$sdrep)[1])){
      samRuns[[i]]           <- SAM2FLR(runs[[i]],sam@control)
      mcstck@stock.n[,,,,,i] <- samRuns[[i]]@stock.n
      mcstck@harvest[,,,,,i] <- samRuns[[i]]@harvest
      #mcstck@catch[,,,,,i]   <- subset(params(samRuns[[i]]),name=="logCatch")$value
    }
    else {
      samRuns[[i]] <- NA
    }
  }
  samRuns <- as(samRuns, "FLSAMs")
  
  #- Save randomly drawn parameters
  pars          <- lapply( samRuns , function(x) params(x)$value)
  random.param  <- matrix(unlist(pars)    , byrow= T, nrow = realisations ,dimnames =list(NULL, params(samRuns[[1]])$name ))
  if(!is.null(saveParsDir))
    save(random.param,file=saveParsDir)
  
  if(return.sam)
    ret <- list(sam=samRuns,stk=mcstck)
  if(!return.sam)
    ret <- mcstck
  return(ret)}

monteCarloStockMscan <- function(stck,tun,sam,realisations,return.sam=FALSE,
                                 saveParsDir=NULL,saveNLLs=NULL,saveSimdataDir=NULL,
                                 set.pars=NULL,mscan=NULL,ncores=detectCores()-1, ...){
  
  require(doParallel)
  ctrl              <- sam@control
  
  #-Create new stock object with nr of realisations in iter slot
  if(class(stck)=="FLStocks"){
    idxStck           <- which.max(unlist(lapply(stck,function(x){x@range["maxyear"]})))
    mcstck            <- propagate(stck[[idxStck]],iter=realisations)
    mcstck            <- window(mcstck,start=range(sam)["minyear"],end=range(sam)["maxyear"])
    mcstck@stock.n[]  <- NA
    mcstck@harvest[]  <- NA
  } else {
    mcstck            <- propagate(stck,iter=realisations)
    mcstck            <- window(mcstck,start=range(sam)["minyear"],end=range(sam)["maxyear"])
    mcstck@stock.n[]  <- NA
    mcstck@harvest[]  <- NA
  }
  
  #- Run first an assessment and generate new data to fit new assessment on
  sam@control@simulate <- TRUE
  sam@control@residuals<- FALSE
  if(!is.null(set.pars)){
    if(class(stck)=="FLStock") stcks <- FLStocks(residual=stck)
    dataS  <- FLSAM2SAM(stcks,tun,ctrl@sumFleets,catch.vars)
    confS  <- ctrl2conf(ctrl,dataS)
    parS   <- defpar(dataS,confS)
    matchNames <- matrix(c("logN.vars","logSdLogN",
                           "logP.vars","logSdLogP",
                           "catchabilities","logFpar",
                           "power.law.exps","logQpow",
                           "f.vars","logSdLogFsta",
                           "obs.vars","logSdLogObs"),ncol=2,byrow=T,dimnames=list(1:6,c("FLSAM","SAM")))
    if(any(!set.pars$par %in% matchNames[,"FLSAM"]))
      warning(paste(set.pars$par[which(!set.pars$par %in% matchNames[,"FLSAM"])],"cannot be set"))
    set.pars.internal <- merge(set.pars,matchNames,by.x="par",by.y="FLSAM")
    
    
    mapDef        <- as.list(matchNames[,"SAM"]); names(mapDef) <- matchNames[,"SAM"]
    for(i in names(mapDef))
      mapDef[[i]] <- jitter(parS[[i]],factor=0.1)
      if(any(duplicated(mapDef[[i]])))
        stop("mapping goes wrong, some parameters get the same starting value, retry")
    map           <- mapDef
    
    for(i in 1:nrow(set.pars.internal)){
      parS[[set.pars.internal$SAM[i]]][set.pars.internal$number[i]+1] <- set.pars.internal$value[i]
      map[[set.pars.internal$SAM[[i]]]][set.pars.internal$number[i]+1] <- NA
    }
    map=list(logSdLogN=as.factor(map$logSdLogN),
             logSdLogP=as.factor(map$logSdLogP),
             logFpar=as.factor(map$logFpar),
             logQpow=as.factor(map$logQpow),
             logSdLogFsta=as.factor(map$logSdLogFsta),
             logSdLogObs=as.factor(map$logSdLogObs))
    map <- map[which(names(map) %in% set.pars$SAM)]
    
    object            <- FLSAM(stcks,tun,sam@control,set.pars=set.pars,return.fit=T)
  } else {
    object            <- FLSAM(stck,tun,sam@control,return.fit=T)
  }
  simdat            <- replicate(realisations,
                                 c(object$data[names(object$data)!="logobs"],#all the old data
                                   object$obj$simulate()["logobs"]),#simulated observations
                                 simplify=FALSE)
  #- Make cluster to speed up process
  ncores <- ifelse(realisations<ncores,realisations,ncores)
  #cl <- makeCluster(ncores)
  #clusterEvalQ(cl,library(FLSAM))
  #clusterEvalQ(cl,library(stockassessment))
  #registerDoParallel(cl)
  
  #- Function to run a new assessment on each of the simulated datasets
  for(i in 1:realisations){
    if(length(which(is.na(object$data$logobs)))>0)
      simdat[[i]]$logobs[which(is.na(object$data$logobs))] <- NA
  }
  if(!is.null(saveSimdataDir))
    save(simdat,file=saveSimdataDir)

  #Store nlls
  simdatM <- simdat
  storeNLL<- matrix(NA,nrow=realisations,ncol=length(mscan),dimnames=list(iter=1:realisations,addM=mscan))
  # Change input data
  for(iM in mscan){
    print(iM)
    for(i in 1:realisations){
      simdatM[[i]]$natMor <- simdat[[i]]$natMor + iM
    }

    cl <- makeCluster(ncores)
    clusterEvalQ(cl,library(stockassessment))
    registerDoParallel(cl)
    
    #run assessments
    if(!is.null(set.pars)){
      capture.output(runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdatM[[i]],object$conf,object$pl[-lengt(object$pl)],map=map,silent=T,newtonsteps=0)))
    } else {
      capture.output(runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdatM[[i]],object$conf,object$pl[-length(object$pl)],silent=T,newtonsteps=0)))
    }

    storeNLL[,ac(iM)] <- unlist(lapply(runs,function(x){try(x$opt$objective)}))
    save(storeNLL,file=saveNLLs)    
    stopCluster(cl)
  }

  #- Get the addm values, and set to default of 0 when non-convergence
  addms <- mscan[apply(storeNLL,1,which.min)]
  if(any(is.na(addms) | !is.finite(addms))){
    idx <- which(is.na(addms) | !is.finite(addms))
    addms[idx] <- 0
  }
  
  #- Final run, adjust the m vectors
  for(i in 1:realisations){
      simdatM[[i]]$natMor <- simdat[[i]]$natMor + addms[i]
  }
  #run assessments
  cl <- makeCluster(ncores)
  clusterEvalQ(cl,library(stockassessment))
  registerDoParallel(cl)
  if(!is.null(set.pars)){
      runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdatM[[i]],object$conf,object$pl[-length(object$pl)],map=map,silent=T,newtonsteps=0,...))
  } else {
      runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdatM[[i]],object$conf,object$pl[-length(object$pl)],silent=T,newtonsteps=0))
  }
  stopCluster(cl) #shut it down
  

  #- Fill the results of the simulations
  samRuns <- list()
  for(i in 1:realisations){
    if(!is.na(unlist(runs[[i]]$sdrep)[1])){
      samRuns[[i]]           <- SAM2FLR(runs[[i]],sam@control)
      mcstck@stock.n[,,,,,i] <- samRuns[[i]]@stock.n
      mcstck@harvest[,,,,,i] <- samRuns[[i]]@harvest
      #mcstck@catch[,,,,,i]   <- subset(params(samRuns[[i]]),name=="logCatch")$value
    }
    else {
      samRuns[[i]] <- NA
    }
  }
  samRuns <- as(samRuns, "FLSAMs")
  
  #- Save randomly drawn parameters
  pars          <- lapply( samRuns , function(x) params(x)$value)
  random.param  <- matrix(unlist(pars)    , byrow= T, nrow = realisations ,dimnames =list(NULL, params(samRuns[[1]])$name ))
  if(!is.null(saveParsDir))
    save(random.param,file=saveParsDir)
  
  if(return.sam)
    ret <- list(sam=samRuns,stk=mcstck,addm=addms)
  if(!return.sam)
    ret <- mcstck
  return(ret)}
