monteCarloStock <- function(stck,tun,sam,realisations){
  require(parallel)
  ctrl              <- sam@control
  
  #-Create new stock object with nr of realisations in iter slot
  mcstck            <- propagate(stck,iter=realisations)
  mcstck            <- window(mcstck,start=range(sam)["minyear"],end=range(sam)["maxyear"])
  mcstck@stock.n[]  <- NA
  mcstck@harvest[]  <- NA

  #- Run first an assessment and generate new data to fit new assessment on
  sam@control@simulate <- TRUE
  object            <- FLSAM(stck,tun,sam@control,return.fit=T)
  simdat            <- replicate(realisations,
                               	c(object$data[names(object$data)!="logobs"],#all the old data
                               	object$obj$simulate(unlist(object$pl))["logobs"]),#simulated observations
                               	simplify=FALSE)
  .get.stem <-function(ctrl) {
    if(length(ctrl@sam.binary)==0) {
	 .stem <- "sam" } else {
	 .stem <- gsub("\\.exe$","",basename(ctrl@sam.binary)) }
  }
  #- Make cluster to speed up process
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl,library(FLSAM))
#  tmb.stem    <- .get.stem(sam@control)
#  #- Find the location of the dll
#  if(length(sam@control@sam.binary)!=0) {
#    tmb.exec <- sam@control@sam.binary
#    if(!file.exists(tmb.exec)) {stop(sprintf("Cannot find specified sam executable (binary) file: %s",tmb.exec))}
#  } else if (.Platform$OS.type=="unix") {
#    tmb.exec <- file.path(system.file("bin", "linux", package="FLSAM",
#                 mustWork=TRUE), tmb.stem)
#  } else if (.Platform$OS.type == "windows") {
#    tmb.exec <- file.path(system.file("bin", "windows", package="FLSAM", mustWork=TRUE),
#    sprintf(tmb.stem))
#  } else {
#    stop(sprintf("Platform type, %s, is not currently supported.",R.version$os))
#  }

  #- Function to run a new assessment on each of the simulated datasets
  runIter <- function(i,object,simdat,sam){
    simdat[[i]]$simFlag <- 0
    simdat[[i]]$resFlag <- 0
    if(length(which(is.na(object$data$logobs)))>0)
      simdat[[i]]$logobs[which(is.na(object$data$logobs))] <- NA
    dyn.load(dynlib("C:/Program Files/R/R-3.3.3/library/FLSAM/bin/windows/sam"))
    obj     <- MakeADFun(simdat[[i]], object$conf, random=c("logN","logF","missing"), DLL="sam")
    opt     <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000))
    rep           <- obj$report()
    sdrep         <- sdreport(obj,opt$par)

    # Last two states
    idx           <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
    sdrep$estY    <- sdrep$value[idx]
    sdrep$covY    <- sdrep$cov[idx,idx]

    idx           <- c(which(names(sdrep$value)=="beforeLastLogN"),which(names(sdrep$value)=="beforeLastLogF"))
    sdrep$estYm1  <- sdrep$value[idx]
    sdrep$covYm1  <- sdrep$cov[idx,idx]

    pl            <- as.list(sdrep,"Est")
    plsd          <- as.list(sdrep,"Std")
    fit           <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=object$conf, opt=opt, obj=obj, rep=rep)
    sam@control@residuals <- FALSE
    stockharvest  <- SAM2FLR(fit,sam@control)
  return(stockharvest)}


  #- Load and run the model
  runs <- parLapply(cl, as.list(1:realisations),runIter,object=object,simdat=simdat,sam=sam)
  stopCluster(cl) #shut it down
  
  #- Fill the results of the simulations
  for(i in 1:realisations){
    mcstck@stock.n[,,,,,i] <- runs[[i]]@stock.n
    mcstck@harvest[,,,,,i] <- runs[[i]]@harvest
    mcstck@catch[,,,,,i]   <- catch(runs[[i]])$value
  }

  return(mcstck)}
