#' Monte Carlo realisations of a stock assessment
#'
#' Simulates new observations from a fitted assessment (conditional on the
#' estimated random effects), refits the assessment to each simulated data
#' set and returns an \code{FLStock} with one iteration per realisation.
#' \code{monteCarloStockMscan} also profiles an additive change in natural
#' mortality (\code{mscan}) for each realisation and uses the best one.
#' Runs in parallel with \pkg{doParallel}.
#'
#' @param stck An \code{FLStock} (or \code{FLStocks}).
#' @param tun An \code{FLIndices}.
#' @param sam The \code{\link{FLSAM-class}} assessment of \code{stck}.
#' @param realisations Number of realisations.
#' @param return.sam Also return the \code{FLSAMs} of the realisations.
#' @param saveParsDir,saveSimdataDir,saveNLLs Optional file names to save the
#'   parameters, the simulated data and the negative log-likelihoods.
#' @param set.pars Optional \code{set.pars} (see \code{\link{FLSAM}}).
#' @param mscan Vector of additive changes in natural mortality.
#' @param ncores Number of cores.
#' @param catch.vars Optional catch weights (see \code{\link{FLSAM}}).
#' @param ... Further arguments passed to \code{sam.fit}.
#' @return An \code{FLStock} with \code{realisations} iterations, or a list
#'   with elements \code{sam} (\code{FLSAMs}) and \code{stk} (and \code{addm}).
#'   Realisations that do not converge are \code{NA}.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' newStock <- monteCarloStock(NSH, NSH.tun, NSH.sam, 10, ncores=2)
#' }
#' @rdname monteCarloStock
#' @export
monteCarloStock <- function(stck,tun,sam,realisations,return.sam=FALSE,saveParsDir=NULL,saveSimdataDir=NULL,set.pars=NULL,ncores=detectCores()-1,catch.vars=NULL, ...){
  
  require(doParallel)
  ctrl              <- sam@control
  
  #-Create new stock object with nr of realisations in iter slot
  if(is(stck,"FLStocks")){
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
    stcks  <- if(is(stck,"FLStock")) FLStocks(residual=stck) else stck
    dataS  <- FLSAM2SAM(stcks,tun,ctrl@sumFleets,catch.vars)
    confS  <- ctrl2conf(ctrl,dataS)
    parS   <- defpar(dataS,confS)
    map    <- .setPars(parS,set.pars)$map
    
    object            <- FLSAM(stcks,tun,sam@control,catch.vars=catch.vars,set.pars=set.pars,return.fit=T)
  } else {
    object            <- FLSAM(stck,tun,sam@control,catch.vars=catch.vars,return.fit=T)
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
    if(!is(runs[[i]],"try-error")){
      samRuns[[ac(i)]]       <- SAM2FLR(runs[[i]],sam@control,keep.fit=FALSE)
      mcstck@stock.n[,,,,,i] <- samRuns[[ac(i)]]@stock.n
      mcstck@harvest[,,,,,i] <- samRuns[[ac(i)]]@harvest
      #mcstck@catch[,,,,,i]   <- subset(params(samRuns[[i]]),name=="logCatch")$value
    }
  }
  if(length(samRuns)<realisations)
    warning(paste(realisations-length(samRuns),"realisation(s) failed to converge and are set to NA"))
  samRuns <- as(samRuns, "FLSAMs")
  
  #- Save randomly drawn parameters
  pars          <- lapply( samRuns , function(x) params(x)$value)
  random.param  <- matrix(unlist(pars)    , byrow= T, nrow = length(samRuns) ,dimnames =list(NULL, params(samRuns[[1]])$name ))
  if(!is.null(saveParsDir))
    save(random.param,file=saveParsDir)
  
  if(return.sam)
    ret <- list(sam=samRuns,stk=mcstck)
  if(!return.sam)
    ret <- mcstck
  return(ret)}

#' @rdname monteCarloStock
#' @export
monteCarloStockMscan <- function(stck,tun,sam,realisations,return.sam=FALSE,
                                 saveParsDir=NULL,saveNLLs=NULL,saveSimdataDir=NULL,
                                 set.pars=NULL,mscan=NULL,ncores=detectCores()-1,catch.vars=NULL, ...){
  
  require(doParallel)
  ctrl              <- sam@control
  
  #-Create new stock object with nr of realisations in iter slot
  if(is(stck,"FLStocks")){
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
    stcks  <- if(is(stck,"FLStock")) FLStocks(residual=stck) else stck
    dataS  <- FLSAM2SAM(stcks,tun,ctrl@sumFleets,catch.vars)
    confS  <- ctrl2conf(ctrl,dataS)
    parS   <- defpar(dataS,confS)
    map    <- .setPars(parS,set.pars)$map
    
    object            <- FLSAM(stcks,tun,sam@control,catch.vars=catch.vars,set.pars=set.pars,return.fit=T)
  } else {
    object            <- FLSAM(stck,tun,sam@control,catch.vars=catch.vars,return.fit=T)
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
      capture.output(runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdatM[[i]],object$conf,object$pl[-length(object$pl)],map=map,silent=T,newtonsteps=0)))
    } else {
      capture.output(runs <- foreach(i = 1:realisations) %dopar% try(sam.fit(simdatM[[i]],object$conf,object$pl[-length(object$pl)],silent=T,newtonsteps=0)))
    }

    storeNLL[,ac(iM)] <- unlist(lapply(runs,function(x){if(is(x,"try-error")) NA else x$opt$objective}))
    if(!is.null(saveNLLs)) save(storeNLL,file=saveNLLs)
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
    if(!is(runs[[i]],"try-error")){
      samRuns[[ac(i)]]       <- SAM2FLR(runs[[i]],sam@control,keep.fit=FALSE)
      mcstck@stock.n[,,,,,i] <- samRuns[[ac(i)]]@stock.n
      mcstck@harvest[,,,,,i] <- samRuns[[ac(i)]]@harvest
      #mcstck@catch[,,,,,i]   <- subset(params(samRuns[[i]]),name=="logCatch")$value
    }
  }
  if(length(samRuns)<realisations)
    warning(paste(realisations-length(samRuns),"realisation(s) failed to converge and are set to NA"))
  samRuns <- as(samRuns, "FLSAMs")
  
  #- Save randomly drawn parameters
  pars          <- lapply( samRuns , function(x) params(x)$value)
  random.param  <- matrix(unlist(pars)    , byrow= T, nrow = length(samRuns) ,dimnames =list(NULL, params(samRuns[[1]])$name ))
  if(!is.null(saveParsDir))
    save(random.param,file=saveParsDir)
  
  if(return.sam)
    ret <- list(sam=samRuns,stk=mcstck,addm=addms)
  if(!return.sam)
    ret <- mcstck
  return(ret)}
