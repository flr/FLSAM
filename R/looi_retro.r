#-------------------------------------------------------------------------------
# Leave one in, leave one out test
#-------------------------------------------------------------------------------

#- Create generic function for 'looi'
if (!isGeneric("looi")) {
  setGeneric('looi', function(stck,tun,ctrl,type,base.run) standardGeneric('looi'))
}

setMethod("looi",signature(stck="FLStock",tun="FLIndices",ctrl="FLSAM.control",type="missing",base.run="FLSAM"),
  function(stck,tun,ctrl,type,base.run){
     looi(stck,tun,ctrl,type="full",base.run="missing")
}) 

setMethod("looi",signature(stck="FLStock",tun="FLIndices",ctrl="FLSAM.control",type="character",base.run="FLSAM"),
  function(stck,tun,ctrl,type="full",base.run="missing"){
    #Check type argument
    if(is.na(pmatch(toupper(type),c("LOI","LOO","FULL")))) {
      stop(sprintf("Invalid 'type' argument: %s",type)) }
    type <- toupper(type)

    #- Create run scheme with all possible combinations
    overview          <- matrix(NA,nrow=2^length(tun),ncol=length(tun),dimnames=list(run=1:2^length(tun),fleet=names(tun)))
    if(any(ctrl@fleets == 6)){
      idxPart         <- which(ctrl@fleets == 6)
      overview        <- matrix(NA,nrow=2^(length(tun)-(length(idxPart)-1)),ncol=(length(tun)-(length(idxPart)-1)),dimnames=list(run=1:2^(length(tun)-(length(idxPart)-1)),fleet=names(ctrl@fleets[c(which(!ctrl@fleets %in% c(0,6,7)),idxPart[1])])))
    }
    for(iFlt in 1:ncol(overview))
      overview[,iFlt] <- rep(c(rep(1,(2^ncol(overview))/(2^(iFlt))),rep(0,(2^ncol(overview))/(2^(iFlt)))),2^(iFlt-1))
    overview          <- overview[-nrow(overview),] #Remove 'no-fleet-included' option
    rownames(overview)<- unlist(lapply(apply(overview,1,function(x){names(which(x==1))}),paste,collapse=" + "))
    LOO.runs <- which(rowSums(overview==0)==1)
    rownames(overview)[LOO.runs] <-  paste("Drop",colnames(overview)[apply(overview[LOO.runs,]==0,1,which)])
    LOI.runs <- which(rowSums(overview)==1)
    rownames(overview)[LOI.runs] <-  paste(colnames(overview)[apply(overview[LOI.runs,]==1,1,which)],"Only")

    if(type=="LOO"){ overview <- overview[c(1,which(rowSums(overview) == (length(tun)-1))),];  cat("Running LOO analyses\n")}
    if(type=="LOI"){ overview <- overview[c(1,which(rowSums(overview) == 1)),];               cat("Running LOI analyses\n")}

    #Perform the initial base assessment i.e. "everybody-in"
    if(missing(base.run))
      base.run <- FLSAM(stck,tun,ctrl)  #Has to work
    result <- FLSAMs("All fleets"=base.run)

    #- Now use the update functionality to do the LOO, LOI, LOO trickery really quickly
    #  We achieve this using the reduced.cfg file in ADMB together with the pin functionality -
    #  in this way we don't need to stress about adjusting the parameters of the pin file
    #  for changes in the number of fleets etc
    for(iRun in rownames(overview)[-1]){
      keepTun       <- names(which(overview[iRun,]==1))
      keepPart      <- which(names(ctrl@fleets) %in% keepTun & ctrl@fleets == 6)
      if(length(keepPart)>0)
        keepTun     <- unique(c(keepTun,names(which(ctrl@fleets==6))))
      Indices.temp  <- tun[keepTun]
      dropTun       <- names(tun)[which(!names(tun) %in% keepTun)]
      Control.temp  <- drop.from.control(ctrl,fleets=dropTun)
      Control.temp@range["maxyear"] <- max(stck@range["maxyear"],
                                           max(sapply(Indices.temp,function(x) max(x@range[c("maxyear")]))))
      Control.temp@range["minyear"] <- min(stck@range["minyear"],
                                           min(sapply(Indices.temp,function(x) min(x@range[c("minyear")]))))
      res           <- try(FLSAM(stck,Indices.temp,Control.temp))
      result[[iRun]]<- res
    }
  return(result)
})

setMethod("looi",signature(stck="FLStocks",tun="FLIndices",ctrl="FLSAM.control",type="character",base.run="FLSAM"),
  function(stck,tun,ctrl,type="full",base.run="missing"){
    #Check type argument
    if(is.na(pmatch(toupper(type),c("LOI","LOO","FULL")))) {
      stop(sprintf("Invalid 'type' argument: %s",type)) }
    type <- toupper(type)

    #- Create run scheme with all possible combinations
    overview          <- matrix(NA,nrow=2^length(tun),ncol=length(tun),dimnames=list(run=1:2^length(tun),fleet=names(tun)))
    if(any(ctrl@fleets == 6)){
      idxPart         <- which(ctrl@fleets == 6)
      overview        <- matrix(NA,nrow=2^(length(tun)-(length(idxPart)-1)),ncol=(length(tun)-(length(idxPart)-1)),dimnames=list(run=1:2^(length(tun)-(length(idxPart)-1)),fleet=names(ctrl@fleets[c(which(!ctrl@fleets %in% c(0,6,7)),idxPart[1])])))
    }
    for(iFlt in 1:ncol(overview))
      overview[,iFlt] <- rep(c(rep(1,(2^ncol(overview))/(2^(iFlt))),rep(0,(2^ncol(overview))/(2^(iFlt)))),2^(iFlt-1))
    overview          <- overview[-nrow(overview),] #Remove 'no-fleet-included' option
    rownames(overview)<- unlist(lapply(apply(overview,1,function(x){names(which(x==1))}),paste,collapse=" + "))
    LOO.runs <- which(rowSums(overview==0)==1)
    rownames(overview)[LOO.runs] <-  paste("Drop",colnames(overview)[apply(overview[LOO.runs,]==0,1,which)])
    LOI.runs <- which(rowSums(overview)==1)
    rownames(overview)[LOI.runs] <-  paste(colnames(overview)[apply(overview[LOI.runs,]==1,1,which)],"Only")

    if(type=="LOO"){ overview <- overview[c(1,which(rowSums(overview) == (length(tun)-1))),];  cat("Running LOO analyses\n")}
    if(type=="LOI"){ overview <- overview[c(1,which(rowSums(overview) == 1)),];               cat("Running LOI analyses\n")}

    #Perform the initial base assessment i.e. "everybody-in"
    if(missing(base.run))
      base.run <- FLSAM(stck,tun,ctrl)  #Has to work
    result <- FLSAMs("All fleets"=base.run)

    #- Now use the update functionality to do the LOO, LOI, LOO trickery really quickly
    #  We achieve this using the reduced.cfg file in ADMB together with the pin functionality -
    #  in this way we don't need to stress about adjusting the parameters of the pin file
    #  for changes in the number of fleets etc
    for(iRun in rownames(overview)[-1]){
      keepTun       <- names(which(overview[iRun,]==1))
      keepPart      <- which(names(ctrl@fleets) %in% keepTun & ctrl@fleets == 6)
      if(length(keepPart)>0)
        keepTun     <- unique(c(keepTun,names(which(ctrl@fleets==6))))
      Indices.temp  <- tun[keepTun]
      dropTun       <- names(tun)[which(!names(tun) %in% keepTun)]
      Control.temp  <- drop.from.control(ctrl,fleets=dropTun)
      Control.temp@range["maxyear"] <- max(sapply(stck,function(x) max(x@range["maxyear"])),
                                           max(sapply(Indices.temp,function(x) max(x@range[c("maxyear")]))))
      Control.temp@range["minyear"] <- min(sapply(stck,function(x) min(x@range["minyear"])),
                                           min(sapply(Indices.temp,function(x) min(x@range[c("minyear")]))))
      res           <- try(FLSAM(stck,Indices.temp,Control.temp))
      result[[iRun]]<- res
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

setMethod("loo",signature(stck="FLStocks",tun="FLIndices",ctrl="FLSAM.control"),
  function(stck,tun,ctrl){
     looi(stck,tun,ctrl,type="loo")
})
setMethod("loi",signature(stck="FLStocks",tun="FLIndices",ctrl="FLSAM.control"),
  function(stck,tun,ctrl){
     looi(stck,tun,ctrl,type="loi")
})

#-------------------------------------------------------------------------------
# Retro function
#-------------------------------------------------------------------------------

if (!isGeneric("retro"))
	setGeneric("retro", function(stock, indices, control, retro, year.range)
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
    res <- list()
    starting.vals <- NULL
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
      if(is.null(starting.vals)){
        assess  <- try(FLSAM(Stock, Indices.temp,control,return.fit=T))
        starting.vals <- assess$pl
      } else {
        assess  <- try(FLSAM(Stock, Indices.temp,control,return.fit=T,starting.values=starting.vals))
        if(class(assess)=="try-error"){
          starting.vals <- NULL
        } else {
          starting.vals <- assess$pl
        }
      }
      if(class(assess)=="try-error") {
        warning(sprintf("Retrospective for year %i failed\n",yr))
      } else {
        res[[as.character(yr)]] <- SAM2FLR(assess,control)
        #res[[as.character(yr)]]@desc    <-  paste(as.character(yr), "Retrospective")
      }
    }
    res        <- as(res,"FLSAMs")
    for(yr in rev(year.range))
      res[[ac(yr)]]@desc <- paste(as.character(yr), "Retrospective")
    res@desc   <- paste("Retrospective analysis from object", stock@desc)
    return(res) 
  })  #End setMethod

setMethod('retro', signature(stock='FLStocks', indices='FLIndices', control='FLSAM.control',retro='numeric'),
function(stock, indices, control, retro, year.range="missing"){

  # ---------- Checks ----------------------
    if (!inherits(stock, "FLStocks"))
      stop("stock must be an 'FLStocks' object!")
    if (!validObject(stock))
      stop("stocks is not a valid object!")
    if (!inherits(indices, "FLIndices"))
      stop("indices must be an 'FLIndices' object!")
    if (!validObject(indices))
      stop("indices is not a valid object!")
    # Check we have at least a usable retro or years.range.
    if ((missing(year.range) || !is.numeric(year.range)) && (!is.numeric(retro) || retro < 0 || is.na(retro) || length(retro) > 1))
      stop("Either 'retro' argument must be a non-negative integer or 'year.range' must be a numeric vector.")
    # ------------ Sort out retrospective years over which to perform assessment ------------
    
    stck.min.yr <- min(unlist(lapply(stock,function(x){return(x@range["minyear"])})))
    stck.max.yr <- max(unlist(lapply(stock,function(x){return(x@range["maxyear"])})))
    stcks.range <- lapply(stock,function(x){return(x@range[c("minyear","maxyear")])})
    trim.stock  <- which.max(do.call(rbind,stcks.range)[,"maxyear"])
    
    if(missing(year.range))
      year.range <- (stck.max.yr-retro):stck.max.yr
    if(min(year.range)<= stcks.range[[trim.stock]]["minyear"])
      stop("Running retrospective beyond combination of sum and residual fleets")
    # Calculate hangover years - but only consider positive hangovers
    hangover.yrs <- sapply(indices,function(x) dims(x)$maxyear) - stck.max.yr
    hangover.yrs <- ifelse(hangover.yrs < 0,0,hangover.yrs)
    # Check year.range is sensible
    if(min(year.range) < stck.min.yr || max(year.range) > stck.max.yr)
      stop("Year range outside stock object range")
    # ---------- Run that retrospective -------------
    cat("Running retrospective...\n")
    res <- list()
    starting.vals <- NULL
    for (yr in rev(year.range)){  #yr is the year in which the assessment is being simulated
      Stock <- new("FLStocks")
      for(iStk in 1:length(stock)){
        if(iStk == trim.stock)
          Stock[[iStk]] <- trim(stock[[trim.stock]], year=stcks.range[[trim.stock]]["minyear"]:yr)
        if(iStk != trim.stock)
          Stock[[iStk]] <- stock[[iStk]]
      }
      names(Stock) <- names(stock)
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
      control@range["maxyear"] <- max(max(unlist(lapply(Stock,function(x){return(x@range["maxyear"])}))),
                max(sapply(Indices.temp,function(x) max(x@range[c("maxyear")]))))

      if(is.null(starting.vals)){
        assess  <- try(FLSAM(Stock, Indices.temp,control,return.fit=T))
        starting.vals <- assess$pl
      } else {
        assess  <- try(FLSAM(Stock, Indices.temp,control,return.fit=T,starting.values=starting.vals))
        if(class(assess)=="try-error"){
          starting.vals <- NULL
        } else {
          starting.vals <- assess$pl
        }
      }
      if(class(assess)=="try-error") {
        warning(sprintf("Retrospective for year %i failed\n",yr))
      } else {
        res[[as.character(yr)]] <- SAM2FLR(assess,control)
        #res[[as.character(yr)]]@desc    <-  paste(as.character(yr), "Retrospective")
      }
    }
    res        <- as(res,"FLSAMs")
    for(yr in rev(year.range))
      res[[ac(yr)]]@desc <- paste(as.character(yr), "Retrospective")
    res@desc   <- paste("Retrospective analysis from object", stock@desc)
    return(res)
  })  #End setMethod