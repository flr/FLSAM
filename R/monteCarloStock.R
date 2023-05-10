monteCarloStock <- function(stck,tun,sam,realisations,return.sam=FALSE,saveParsDir=NULL,set.pars=NULL,ncores=detectCores()-1, ...){
  
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

    object            <- FLSAM(stck,tun,sam@control,set.pars=set.pars,return.fit=T)
  } else {
    object            <- FLSAM(stck,tun,sam@control,return.fit=T)
  }
  simdat            <- replicate(realisations,
                               	c(object$data[names(object$data)!="logobs"],#all the old data
                               	object$obj$simulate(unlist(object$pl))["logobs"]),#simulated observations
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
  #clusterExport(cl,varlist=c("simdat","object"),envir=environment())
  if(!is.null(set.pars)){
    runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdat[[i]],object$conf,object$pl,map=map,silent=T,...))
  } else {
    runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdat[[i]],object$conf,object$pl,silent=T,...))
  }

  stopCluster(cl) #shut it down
  
  if(return.sam){
    resSAM <- list()
    for(i in 1:iters){
      if(!is.na(unlist(res[[i]]$sdrep)[1])){
        resSAM[[i]] <- SAM2FLR(res[[i]],sam@control)
      } else {
        resSAM[[i]] <- NA
      }
    }
    resSAM <- as(resSAM,"FLSAMs")
  }

  #- Fill the results of the simulations
  if(!return.sam){
    samRuns <- list()
    for(i in 1:realisations){
      if(!is.na(unlist(runs[[i]]$sdrep)[1])){
        samRuns[[i]]           <- SAM2FLR(runs[[i]],sam@control)
        mcstck@stock.n[,,,,,i] <- samRuns[[i]]@stock.n
        mcstck@harvest[,,,,,i] <- samRuns[[i]]@harvest
        #mcstck@catch[,,,,,i]   <- subset(params(samRuns[[i]]),name=="logCatch")$value
      }
    }
  }
  
  #- Save randomly drawn parameters
  pars          <- lapply( samRuns , function(x) params(x)$value)
  random.param  <- matrix(unlist(pars)    , byrow= T, nrow = realisations ,dimnames =list(NULL, params(samRuns[[1]])$name ))
  if(!is.null(saveParsDir))
    save(random.param,file=saveParsDir)

  if(return.sam)
    ret <- resSAM
  if(!return.sam)
    ret <- mcstck
return(ret)}
