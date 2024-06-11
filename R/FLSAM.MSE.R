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

  res <- foreach(i = 1:iters) %dopar% try(sam.fit(data[[i]], 
                                                  conf[[i]], par[[i]],
                                                  sim.condRE = ctrl@simulate,
                                                  newtonsteps = 0,
                                                  silent = T))

  stopCluster(cl)
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
  
  if(!return.sam){
    for(i in 1:iters){
      if(!is.na(unlist(res[[i]]$sdrep)[1])){
        stcks[["residual"]]@stock.n[, , , , , i] <- FLQuant(t(ntable(res[[i]])), dimnames = list(age = ctrl@range["min"]:ctrl@range["max"], 
                                                                                                 year = ctrl@range["minyear"]:ctrl@range["maxyear"], unit = "unique", 
                                                                                                 season = "all", area = "unique", iter = 1))
        stcks[["residual"]]@harvest[, , , , , i] <- FLQuant(t(faytable(res[[i]])), dimnames = list(age = ctrl@range["min"]:ctrl@range["max"], 
                                                                                              year = ctrl@range["minyear"]:ctrl@range["maxyear"], 
                                                                                              unit = "unique", season = "all", area = "unique", 
                                                                                              iter = 1))
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
