#-------------------------------------------------------------------------------
# Leave one in, leave one out test
#-------------------------------------------------------------------------------

#- Create generic function for 'looi'
if (!isGeneric("looi")) {
  setGeneric('looi', function(e1,e2,e3,type="full") standardGeneric('looi'))
}

setMethod("looi",signature(e1="FLStock",e2="FLIndices",e3="FLSAM.control",type="character"),
  function(e1,e2,e3,type="full"){
  
    #- Create run scheme with all possible combinations
    overview          <- matrix(NA,nrow=2^length(e2),ncol=length(e2),dimnames=list(run=1:2^length(e2),fleet=names(e2)))
    for(iFlt in 1:length(e2))
      overview[,iFlt] <- rep(c(rep(1,(2^length(e2))/(2^(iFlt))),rep(0,(2^length(e2))/(2^(iFlt)))),2^(iFlt-1))
    overview          <- overview[-nrow(overview),] #Remove 'no-fleet-included' option
    rownames(overview)<- unlist(lapply(apply(overview,1,function(x){names(which(x==1))}),paste,collapse=" + "))
    
    if(type=="LOO"){ overview <- overview[c(1,which(rowSums(overview) == (length(e2)-1))),];  cat("Running LOO analyses\n")}
    if(type=="LOI"){ overview <- overview[c(1,which(rowSums(overview) == 1)),];               cat("Running LOI analyses\n")}
    
    #- Create object to save outcomes
    result <- new("FLSAMs")
    for(iRun in rownames(overview)){
      stck                  <- e1
      tun                   <- e2[which(overview[iRun,]==1)]
      ctrl                  <- FLSAM.control(stck,tun)
      #- Update the ctrl and stck object given the changed tuning object
      for(iSlot in slotNames(e3))
        if(class(slot(ctrl,iSlot))=="matrix") {
            slot(ctrl,iSlot) <- slot(e3,iSlot)[c(names(which(e3@fleets==0)),names(which(overview[iRun,]==1))),]}
      ctrl@logN.vars        <- e3@logN.vars
      ctrl@srr              <- e3@srr
      ctrl@name             <- iRun
      ctrl                  <- update(ctrl)

      #- Run the assessment
      cat(sprintf('\nRunning "%s" assessment...\n',iRun))
      if(iRun == rownames(overview)[1]) {  #If first run
           res <- FLSAM(stck,tun,ctrl,run.dir=tempdir(),batch.mode=TRUE)
           first.res <- res  #Store first result separately, and use to initialise other runs 
      } else {
           res <- update(first.res,stck,tun,run.dir=tempdir())}
      if(is.null(res)) {
        warning(sprintf('"%s" run failed',iRun))
      } else {
        result[[iRun]] <- res
      }}
  return(result)
})

#-------------------------------------------------------------------------------
# Retro function
#-------------------------------------------------------------------------------

if (!isGeneric("retro"))
	setGeneric("retro", function(stock, indices, control, retro, ...)
    	standardGeneric("retro"))

setMethod('retro', signature(stock='FLStock', indices='FLIndices', control='FLSAM.control',retro='numeric'),
function(stock, indices, control, retro, year.range="missing"){
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
          assess  <- FLSAM(Stock, Indices.temp,control,run.dir=tempdir(),batch.mode=TRUE)
          base.assess <- assess
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

