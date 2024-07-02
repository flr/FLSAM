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
  
  cores <- 1
  if (is.numeric(parallel)) {
    cores <- parallel
    parallel <- TRUE
  }
  else if (getDoParRegistered()) {
    cores <- nbrOfWorkers()
  }
  
  # require(doParallel)
  # ncores <- detectCores()-1
  # ncores <- ifelse(iters<ncores,iters,ncores)
  # cl <- makeCluster(ncores) #set up nodes
  # clusterEvalQ(cl,library(FLSAM))
  # clusterEvalQ(cl,library(stockassessment))
  # registerDoParallel(cl)
  
  message("SAM Running on ", nbrOfWorkers(), " nodes.")
  
  its <- split(seq(it), sort(seq(it)%%cores))
  
  data <- foreach(i = its) %dofuture% FLSAM2SAM(FLStocks("residual"=iter(stcks[["residual"]],i)),FLIndices(lapply(tun, function(x) iter(x,i))),ctrl@sumFleets,catch.vars)
  
  data <- foreach(i = its) %dofuture% FLSAM2SAM(FLStocks("residual"=iter(stcks[["residual"]],i)),FLIndices(lapply(tun, function(x) iter(x,i))),ctrl@sumFleets,catch.vars)
  conf <- foreach(i = its) %dofuture% ctrl2conf(ctrl,data[[i]])
  par  <- foreach(i = its) %dofuture% defpar(data[[i]],conf[[i]])
  
  res <- foreach(i = its) %dofuture% try(sam.fit(data[[i]], 
                                                  conf[[i]], par[[i]],
                                                  sim.condRE = ctrl@simulate,
                                                  newtonsteps = 0,
                                                  silent = T))
  
  #stopCluster(cl)
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

FLSAM.fast <-function(stcks,tun,ctrl,catch.vars=NULL,starting.sam=NULL,return.sam=F,...){
  
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(is(stcks, "FLStocks")) stop("Not implemented for multi-fleet yet")
  if(is(stcks, "FLStock")) stcks <- FLStocks(residual=stcks)
  
  #Turn residuals off
  ctrl@residuals <- FALSE
  
  data  <- FLSAM2SAM(FLStocks("residual"=stcks[["residual"]]),FLIndices(tun),ctrl@sumFleets,catch.vars)
  conf  <- ctrl2conf(ctrl,data)
  par   <- defpar(data,conf)
  
  res <- sam.fit(data, 
                 conf,
                 par,
                 sim.condRE = ctrl@simulate,
                 newtonsteps = 0,
                 silent = T)
  
  if(return.sam){
    resSAM <- SAM2FLR(res,ctrl)
  }
  
  if(!return.sam){
    if(!is.na(unlist(res$sdrep)[1])){
      stcks[["residual"]]@stock.n <- FLQuant(t(ntable(res)), dimnames = list(age = ctrl@range["min"]:ctrl@range["max"], 
                                                                                          year = ctrl@range["minyear"]:ctrl@range["maxyear"], unit = "unique", 
                                                                                          season = "all", area = "unique", iter = 1))
      stcks[["residual"]]@harvest <- FLQuant(t(faytable(res)), dimnames = list(age = ctrl@range["min"]:ctrl@range["max"], 
                                                                                            year = ctrl@range["minyear"]:ctrl@range["maxyear"], 
                                                                                            unit = "unique", season = "all", area = "unique", 
                                                                                            iter = 1))
    } else {
      stcks[["residual"]]@stock.n <- NA
      stcks[["residual"]]@harvest <- NA
    }
  }
  
  if(return.sam)
    ret <- resSAM
  if(!return.sam)
    ret <- stcks
  return(ret)
}