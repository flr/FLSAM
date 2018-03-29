monteCarloStock <- function(stck,tun,sam,realisations,return.sam=FALSE,...){
  require(doParallel)
  ctrl              <- sam@control
  
  #-Create new stock object with nr of realisations in iter slot
  if(class(stck)=="FLStocks"){
    idxStck           <- which.max(unlist(lapply(NSHs,function(x){x@range["maxyear"]})))
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
  object            <- FLSAM(stck,tun,sam@control,return.fit=T)
  simdat            <- replicate(realisations,
                               	c(object$data[names(object$data)!="logobs"],#all the old data
                               	object$obj$simulate(unlist(object$pl))["logobs"]),#simulated observations
                               	simplify=FALSE)
  #- Make cluster to speed up process
  ncores <- detectCores()-1
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
  runs <- foreach(i = 1:realisations) %dopar% try(sam.fitfast(simdat[[i]],object$conf,object$pl,silent=T,...))
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

  if("doParallel" %in% (.packages()))
    detach("package:doParallel",unload=TRUE)
  if("foreach" %in% (.packages()))
    detach("package:foreach",unload=TRUE)
  if("iterators" %in% (.packages()))
    detach("package:iterators",unload=TRUE)

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
  if(return.sam)
    ret <- resSAM
  if(!return.sam)
    ret <- stcks
return(ret)}
