monteCarloStock <- function(stck,tun,sam,realisations){
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
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl,library(FLSAM))
  clusterEvalQ(cl,library(stockassessment))
  registerDoParallel(cl)

  #- Function to run a new assessment on each of the simulated datasets
  runIter <- function(i,object,simdat,sam){
    if(length(which(is.na(object$data$logobs)))>0)
      simdat[[i]]$logobs[which(is.na(object$data$logobs))] <- NA
    res <- try(sam.fit(simdat[[i]],object$conf,defpar(simdat[[i]],object$conf)))
    stockharvest  <- try(SAM2FLR(res,sam@control))
  return(stockharvest)}


  #- Load and run the model
  runs <- foreach(i=1:realisations) %dopar% runIter(i,object=object,simdat=simdat,sam=sam)
  stopCluster(cl) #shut it down
  
  #- Fill the results of the simulations
  for(i in 1:realisations){
    if(class(runs[[i]])!="try-error"){
      mcstck@stock.n[,,,,,i] <- runs[[i]]@stock.n
      mcstck@harvest[,,,,,i] <- runs[[i]]@harvest
      mcstck@catch[,,,,,i]   <- catch(runs[[i]])$value
    }
  }

  return(mcstck)}
