#-------------------------------------------------------------------------------
# Leave one in, leave one out test
#-------------------------------------------------------------------------------

#- Create generic function for 'looi'
if (!isGeneric("looi")) {
  setGeneric('looi', function(stck,tun,ctrl,type) standardGeneric('looi'))
}

setMethod("looi",signature(stck="FLStock",tun="FLIndices",ctrl="FLSAM.control",type="missing"),
  function(stck,tun,ctrl,type){
     looi(stck,tun,ctrl,type="full")
}) 

setMethod("looi",signature(stck="FLStock",tun="FLIndices",ctrl="FLSAM.control",type="character"),
  function(stck,tun,ctrl,type="full"){
    #Check type argument
    if(is.na(pmatch(toupper(type),c("LOI","LOO","FULL")))) {
      stop(sprintf("Invalid 'type' argument: %s",type)) }
    type <- toupper(type)
   
    #- Create run scheme with all possible combinations
    overview          <- matrix(NA,nrow=2^length(tun),ncol=length(tun),dimnames=list(run=1:2^length(tun),fleet=names(tun)))
    for(iFlt in 1:length(tun))
      overview[,iFlt] <- rep(c(rep(1,(2^length(tun))/(2^(iFlt))),rep(0,(2^length(tun))/(2^(iFlt)))),2^(iFlt-1))
    overview          <- overview[-nrow(overview),] #Remove 'no-fleet-included' option
    rownames(overview)<- unlist(lapply(apply(overview,1,function(x){names(which(x==1))}),paste,collapse=" + "))
    LOO.runs <- which(rowSums(overview==0)==1) 
    rownames(overview)[LOO.runs] <-  paste("Drop",colnames(overview)[apply(overview[LOO.runs,]==0,1,which)])
    LOI.runs <- which(rowSums(overview)==1) 
    rownames(overview)[LOI.runs] <-  paste(colnames(overview)[apply(overview[LOI.runs,]==1,1,which)],"Only")
    
    if(type=="LOO"){ overview <- overview[c(1,which(rowSums(overview) == (length(tun)-1))),];  cat("Running LOO analyses\n")}
    if(type=="LOI"){ overview <- overview[c(1,which(rowSums(overview) == 1)),];               cat("Running LOI analyses\n")}
    
    #Perform the initial base assessment i.e. "everybody-in"
    rundir <- tempdir()
    base.run <- FLSAM(stck,tun,ctrl,run.dir=rundir,batch.mode=FALSE)  #Has to work
    result <- FLSAMs("All fleets"=base.run)

    #- Now use the update functionality to do the LOO, LOI, LOO trickery really quickly
    #  We achieve this using the reduced.cfg file in ADMB together with the pin functionality - 
    #  in this way we don't need to stress about adjusting the parameters of the pin file
    #  for changes in the number of fleets etc
    for(iRun in rownames(overview)[-1]){
      #

      #Write the outputs
      FLR2SAM(stck,tun,ctrl,run.dir=rundir)

      #Write the pin file
      params2pin(base.run,save.dir=rundir)

      #Mark which surveys to use via the reduced.cfg file
      flts <- ctrl@fleets
      flts[] <- 0     #Everybody in
      flts[colnames(overview)] <- overview[iRun,]-1   #-1 is out for reduced.cfg
      cat("# This allows to exclude data in ways customized for retrospective runs\n",
         "# The following specifies one integer for each fleet (catches and surveys)\n",
         "# if:\n",
         "#     0: all data for that fleet is used\n",
         "#    -1: all data for that fleet is excluded\n",
         "#     n: n is a positive integer the corresponding fleet's data is reduced\n",
         "#        by n years starting from the most recent year.\n",
         flts,"\n", 
         "\n\n# Checksums to ensure correct reading of input data \n",3142, 3142,"\n",
         file=file.path(rundir,"reduced.cfg"),append=FALSE) 

      #- Run the assessment
      cat(sprintf('\nRunning "%s" assessment...\n',iRun))
      rtn <- runSAM(ctrl, run.dir=rundir,use.pin=TRUE)
      if(rtn!=0) {
        warning(sprintf("%s run failed\n",iRun))
        next  #Loop 
      }

      #- Collect the results
      res <- SAM2FLR(ctrl,run.dir=rundir)    
      result[[iRun]] <- res
      }
  return(result)
})

#- Create generic function for 'loo' and 'loi'
if (!isGeneric("loo")) {
  setGeneric('loo', function(stck,tun,ctrl) standardGeneric('loo'))
}

if (!isGeneric("loi")) {
  setGeneric('loi', function(stck,tun,ctrl) standardGeneric('loi'))
}
setMethod("loo",signature(stck="FLStock",tun="FLIndices",ctrl="FLSAM.control"),
  function(stck,tun,ctrl){
     looi(stck,tun,ctrl,type="loo")
}) 
setMethod("loi",signature(stck="FLStock",tun="FLIndices",ctrl="FLSAM.control"),
  function(stck,tun,ctrl){
     looi(stck,tun,ctrl,type="loi")
}) 

#-------------------------------------------------------------------------------
# Retro function
#-------------------------------------------------------------------------------

if (!isGeneric("retro"))
	setGeneric("retro", function(stock, indices, control, retro, year.range,base.assess)
    	standardGeneric("retro"))

setMethod('retro', signature(stock='FLStock', indices='FLIndices', control='FLSAM.control',retro='numeric'),
function(stock, indices, control, retro, year.range="missing",base.assess="missing"){
  # ---------- Checks ----------------------
    if (!inherits(stock, "FLStock"))
      stop("stock must be an 'FLStock' object!")
    if (!validObject(stock))
      stop("stock is not a valid object!")
    if (!inherits(indices, "FLIndices"))
      stop("indices must be an 'FLIndices' object!")
    if (!validObject(indices))
      stop("indices is not a valid object!")
    # Check we have at least a usable retro or years.range.
    if ((missing(year.range) || !is.numeric(year.range)) && (!is.numeric(retro) || retro < 0 || is.na(retro) || length(retro) > 1))
      stop("Either 'retro' argument must be a non-negative integer or 'year.range' must be a numeric vector.")
    # ------------ Sort out retrospective years over which to perform assessment ------------
    stck.min.yr <- stock@range["minyear"]
    stck.max.yr <- stock@range["maxyear"]
    if(missing(year.range))
      year.range <- (stck.max.yr-retro):stck.max.yr
    # Calculate hangover years - but only consider positive hangovers
    hangover.yrs <- sapply(indices,function(x) dims(x)$maxyear) - dims(stock)$maxyear
    hangover.yrs <- ifelse(hangover.yrs < 0,0,hangover.yrs)
    # Check year.range is sensible
    if(min(year.range) < stck.min.yr || max(year.range) > stck.max.yr)
      stop("Year range outside stock object range")
    # ---------- Run that retrospective -------------
    cat("Running retrospective...\n")
    res <- new("FLSAMs")
    for (yr in rev(year.range)){  #yr is the year in which the assessment is being simulated
      Stock <- trim(stock, year=stck.min.yr:yr)
      Indices.temp<-indices
      for (j in 1:length(indices)) {
        idx.min.yr <- dims(indices[[j]])$minyear
        if (yr < idx.min.yr) {stop(paste("Year requested (",yr,
              ") is earlier than the first year (",
              idx.min.yr,") of the ",indices[[j]]@name," index.",sep="")) }
        idx.max.yr <- min(dims(indices[[j]])$maxyear, yr+hangover.yrs[j])
        Indices.temp[[j]] <- trim(indices[[j]],year=idx.min.yr:idx.max.yr)
      }
      control@name <- as.character(yr)
      control@range["maxyear"] <- max(Stock@range["maxyear"],
                max(sapply(Indices.temp,function(x) max(x@range[c("maxyear")]))))
      if(yr == rev(year.range)[1]) {  #First run
        if(missing(base.assess)){
            assess  <- FLSAM(Stock, Indices.temp,control,run.dir=tempdir(),batch.mode=TRUE)
            base.assess <- assess
        } else { assess <- update(base.assess,Stock,Indices.temp,run.dir=tempdir()); base.assess <- assess}
      } else {
          assess  <- update(base.assess,Stock,Indices.temp,run.dir=tempdir()) }
      if(is.null(assess)) {
        warning(sprintf("Retrospective for year %i failed\n",yr))
      } else {
        assess@desc    <-  paste(as.character(yr), "Retrospective")
        res[[as.character(yr)]] <- assess
      }
    }
    res@desc   <- paste("Retrospective analysis from object", stock@desc)
    return(res) 
  })  #End setMethod

