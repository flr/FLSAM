ctrl2conf <- function(ctrl,data){
  conf              <- defcon(data)
  conf$minAge[]       <- ctrl@range["min"];
  conf$maxAge[]       <- ctrl@range["max"];
  conf$maxAgePlusGroup[] <- ifelse(ctrl@plus.group,1,0);
  conf$keyLogFsta[]   <- ctrl@states
  conf$corFlag[]      <- ctrl@cor.F
  conf$keyLogFpar[]   <- ctrl@catchabilities
  conf$keyQpow[]      <- ctrl@power.law.exps
  conf$keyVarF[]      <- ctrl@f.vars
  conf$keyVarLogN[]   <- ctrl@logN.vars;
  conf$keyVarLogP[]   <- if(is.na(ctrl@logP.vars[1])){numeric(0)}else{ctrl@logP.vars}
  conf$keyVarObs[]    <- ctrl@obs.vars
  conf$obsCorStruct[] <- factor(ctrl@cor.obs.Flag,levels=c("ID","AR","US"))
  conf$keyCorObs[]    <- ctrl@cor.obs
  conf$stockRecruitmentModelCode[] <- ctrl@srr
  conf$noScaledYears[] <- ctrl@scaleNoYears
  if(conf$noScaledYears>0){
    conf$keyScaledYears[]  <- ctrl@scaleYears
    conf$keyParScaledYA[] <- ctrl@scalePars
  } else {
    conf$keyScaledYears <- numeric()
    conf$keyParScaledYA <- matrix(numeric(),nrow=0,ncol=0)
  }
  conf$fbarRange[]    <- c(ctrl@range[c("minfbar","maxfbar")])
  conf$keyBiomassTreat[] <- rep(-1,data$noFleets)
  conf$keyBiomassTreat[which(data$fleetTypes %in% c(3,4))] <- ctrl@biomassTreat[which(ctrl@fleets %in% c(3,4))]
  conf$obsLikelihoodFlag[] <- factor(ctrl@likFlag,levels=c("LN","ALN"))
  conf$fixVarToWeight[] <- 0
return(conf)}


FLSAM2SAM <- function(stcks,tun,sumFleets=NULL,catch.vars=NULL){

  allMax <- max(sapply(stcks,function(x) max(x@range["maxyear"])),
            max(sapply(tun,function(x) max(x@range[c("maxyear")]))))
  stcksMax <- max(sapply(stcks,function(x) max(x@range["maxyear"])))

  if (stcksMax<allMax)  {
    stck <- stcks[[which.max(sapply(stcks,function(x) max(x@range["maxyear"])))]]
    stck<- window(stck,end=allMax,
              frequency=1,extend=TRUE)  # extend by one year
    stck@mat[,as.character((stcksMax+1):allMax),,,,]         <- stck@mat[,as.character(stcksMax),,,,]      # MV (there must be an easier way!!
    stck@stock.wt[,as.character((stcksMax+1):allMax),,,,]    <- stck@stock.wt[,as.character(stcksMax),,,,]
    stck@catch.wt[,as.character((stcksMax+1):allMax),,,,]    <- stck@catch.wt[,as.character(stcksMax),,,,]
    stck@discards.wt[,as.character((stcksMax+1):allMax),,,,] <- stck@discards.wt[,as.character(stcksMax),,,,]
    stck@landings.wt[,as.character((stcksMax+1):allMax),,,,] <- stck@landings.wt[,as.character(stcksMax),,,,]
    stck@m[,as.character((stcksMax+1):allMax),,,,]           <- stck@m[,as.character(stcksMax),,,,]
    stck@landings.n[,as.character((stcksMax+1):allMax),,,,]  <- stck@landings.n[,as.character(stcksMax),,,,]
    stck@harvest.spwn[,as.character((stcksMax+1):allMax),,,,]<- stck@harvest.spwn[,as.character(stcksMax),,,,]
    stck@m.spwn[,as.character((stcksMax+1):allMax),,,,]      <- stck@m.spwn[,as.character(stcksMax),,,,]
    stcks[[which.max(sapply(stcks,function(x) max(x@range["maxyear"])))]] <- stck
  }

  #- Prepare survey data
  surveysFLR  <- lapply(tun,function(x){ret <- t(x@index[,drop=T]);if(dim(ret)[1]==1) ret <- t(ret); return(ret)})
  dmnsFLR     <- lapply(tun,dims)
  typeFLR     <- lapply(tun,function(x){return(x@type)})
  for(iSurv in names(surveysFLR)){
    colnames(surveysFLR[[iSurv]]) <- seq(dmnsFLR[[iSurv]]$min,dmnsFLR[[iSurv]]$max,1)
    attr(surveysFLR[[iSurv]],"time") <- c(dmnsFLR[[iSurv]]$startf,dmnsFLR[[iSurv]]$endf)
    if(typeFLR[[iSurv]]=="partial")
      attr(surveysFLR[[iSurv]],"part") <- which(names(which(unlist(typeFLR)=="partial"))==iSurv)
  }

  #- Prepare catch data
  if(dims(stcks[["residual"]])$area>1){
    residual.fleets <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(t(stcks[["residual"]]@catch.n[,,,,x,drop=T]))})
  } else {
    residual.fleets <- list(t(stcks[["residual"]]@catch.n[,drop=T]))
  }
  if("sum" %in% names(stcks)){
    sum.fleets      <- list(t(stcks[["sum"]]@catch.n[,drop=T]))
    attr(sum.fleets[[1]], "sumof") <- which(dimnames(stcks[["residual"]]@catch.n)$area %in% sumFleets)
  } else {
    sum.fleets <- NULL
  }

  #- Prepare props, natmort and stock weigth
  if("sum" %in% names(stcks))
    propMat   <- rbind(t(stcks[["residual"]]@mat[,,,,1,drop=T]),t(stcks[["sum"]]@mat[,drop=T]))
  if(!"sum" %in% names(stcks))
    propMat   <- t(stcks[["residual"]]@mat[,,,,1,drop=T])
  propMat     <- propMat[ac(sort(rownames(propMat))),]

  if("sum" %in% names(stcks))
    propF     <- rbind(t(stcks[["residual"]]@harvest.spwn[,,,,1,drop=T]),t(stcks[["sum"]]@harvest.spwn[,drop=T]))
  if(!"sum" %in% names(stcks))
    propF     <- t(stcks[["residual"]]@harvest.spwn[,,,,1,drop=T])
  propF       <- propF[ac(sort(rownames(propF))),]

  if("sum" %in% names(stcks))
    propM     <- rbind(t(stcks[["residual"]]@m.spwn[,,,,1,drop=T]),t(stcks[["sum"]]@m.spwn[,drop=T]))
  if(!"sum" %in% names(stcks))
    propM     <- t(stcks[["residual"]]@m.spwn[,,,,1,drop=T])
  propM       <- propM[ac(sort(rownames(propM))),]

  if("sum" %in% names(stcks))
    stockWeight <- rbind(t(stcks[["residual"]]@stock.wt[,,,,1,drop=T]),t(stcks[["sum"]]@stock.wt[,drop=T]))
  if(!"sum" %in% names(stcks))
    stockWeight <- t(stcks[["residual"]]@stock.wt[,,,,1,drop=T])
  stockWeight <- stockWeight[ac(sort(rownames(stockWeight))),]

  if("sum" %in% names(stcks))
    natMort   <- rbind(t(stcks[["residual"]]@m[,,,,1,drop=T]),t(stcks[["sum"]]@m[,drop=T]))
  if(!"sum" %in% names(stcks))
    natMort   <- t(stcks[["residual"]]@m[,,,,1,drop=T])
  natMort     <- natMort[ac(sort(rownames(natMort))),]

  #- Prepare catch, discards, landings and landing fraction
  if("sum" %in% names(stcks)){
    catchWeight   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(rbind(t(stcks[["residual"]]@catch.wt[,,,,x,drop=T]),t(stcks[["sum"]]@catch.wt[,drop=T])))})
    catchWeight   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(catchWeight[[x]][ac(sort(rownames(catchWeight[[x]]))),])})
  }
  if(!"sum" %in% names(stcks)){
    catchWeight   <- t(stcks[["residual"]]@catch.wt[,,,,1,drop=T])
    catchWeight   <- catchWeight[ac(sort(rownames(catchWeight))),]
  }

  if("sum" %in% names(stcks)){
    discardWeight   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(rbind(t(stcks[["residual"]]@discards.wt[,,,,x,drop=T]),t(stcks[["sum"]]@discards.wt[,drop=T])))})
    discardWeight   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(discardWeight[[x]][ac(sort(rownames(discardWeight[[x]]))),])})
  }
  if(!"sum" %in% names(stcks)){
    discardWeight   <- t(stcks[["residual"]]@discards.wt[,,,,1,drop=T])
    discardWeight   <- discardWeight[ac(sort(rownames(discardWeight))),]
  }

  if("sum" %in% names(stcks)){
    landingWeight   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(rbind(t(stcks[["residual"]]@landings.wt[,,,,x,drop=T]),t(stcks[["sum"]]@landings.wt[,drop=T])))})
    landingWeight   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(landingWeight[[x]][ac(sort(rownames(landingWeight[[x]]))),])})
  }
  if(!"sum" %in% names(stcks)){
    landingWeight   <- t(stcks[["residual"]]@landings.wt[,,,,1,drop=T])
    landingWeight   <- landingWeight[ac(sort(rownames(landingWeight))),]
  }

  if("sum" %in% names(stcks)){
    landFrac   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(rbind(t(stcks[["residual"]]@landings.n[,,,,x,drop=T])/t(stcks[["residual"]]@catch.n[,,,,x,drop=T]),t(stcks[["sum"]]@landings.n[,drop=T])/t(stcks[["sum"]]@catch.n[,drop=T])))})
    landFrac   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){return(landFrac[[x]][ac(sort(rownames(landFrac[[x]]))),])})
    landFrac   <- lapply(as.list(1:dims(stcks[["residual"]])$area),function(x){landFrac[[x]][is.na(landFrac[[x]])] <- 1; return(landFrac[[x]])})
  }
  if(!"sum" %in% names(stcks)){
    landFrac   <- t(stcks[["residual"]]@landings.n[,,,,1,drop=T]/stcks[["residual"]]@catch.n[,,,,1,drop=T])
    landFrac   <- landFrac[ac(sort(rownames(landFrac))),]
    landFrac[is.na(landFrac)] <- 1
  }


    sam.dat <-setup.sam.data(surveys=surveysFLR,
                              residual.fleets=residual.fleets, # Notice list
                              sum.residual.fleets=sum.fleets,
                              prop.mature=propMat,
                              stock.mean.weight=stockWeight,
                              catch.mean.weight=catchWeight,
                              dis.mean.weight=discardWeight,
                              land.mean.weight=landingWeight,
                              prop.f=propF,
                              prop.m=propM,
                              natural.mortality=natMort,
                              land.frac=landFrac)

  # Get data weighting
  if(is.null(catch.vars)==F){
    if("sum" %in% names(stcks)){
      idxSum      <- which(sam.dat$aux[,"fleet"] == which(sam.dat$fleetTypes == 7))
      sam.dat$weight[idxSum] <- c(t(catch.vars[["sum"]][,drop=T]))
    }
    counter     <- 1
    for(iFleet in which(sam.dat$fleetTypes==0)){
      idxRes    <- which(sam.dat$aux[,"fleet"] == iFleet)
      sam.dat$weight[idxRes] <- c(t(catch.vars[["residual"]][,,,,counter,drop=T]))
      counter   <- counter + 1
    }
  }

  tun.var       <- lapply(tun,index.var)
  for(iFleet in which(!sam.dat$fleetTypes %in% c(0,7))){
    iTun        <- which(names(tun)==attr(sam.dat,"fleetNames")[iFleet])
    idxTun      <- which(sam.dat$aux[,"fleet"] == iFleet)
    tunWeight   <- as.data.frame(tun.var[[iTun]])$data
    if(sam.dat$fleetTypes[iFleet]==6)
      tunWeight <- tunWeight[which(as.data.frame(tun[[iTun]]@index)$data>0)]
    sam.dat$weight[idxTun] <- tunWeight
  }
return(sam.dat)}



FLR2SAM <-function(stck,tun,ctrl,map=NULL,catch.var=NULL) {
  #---------------------------------------------------
  # Setup for output
  #---------------------------------------------------
  #General Setup
  sam.stem <- .get.stem(ctrl)
  run.time <- Sys.time()
  miss.val <- 0

  #Setup meta data
  samp.times <- c(rep(miss.val,dims(stck)$area),sapply(tun,function(x) mean(x@range[c("startf","endf")])))
  samp.times <- ifelse(is.na(samp.times),miss.val,samp.times)
  names(samp.times) <- names(ctrl@fleets)
  yrs <- ctrl@range["minyear"]:ctrl@range["maxyear"]
  lastYear <- ctrl@range["maxyear"]
  nyrs <- length(yrs)
  
  #Add a year to the stock object when tuning data are newer than catch data
  #Other combinations are not taken into account. Credit to Morten Vinther for this patch
  extraYr  <- F
  if (stck@range["maxyear"]<lastYear)  {
   extraYr <- T
   stck<- window(stck,start=stck@range["minyear"],end=ctrl@range["maxyear"],
              frequency=1,extend=TRUE)  # extend by one year
   stck@mat[,as.character(lastYear),,,,]          <- stck@mat[,as.character(lastYear-1),,,,]      # MV (there must be an easier way!!
   stck@stock.wt[,as.character(lastYear),,,,]     <- stck@stock.wt[,as.character(lastYear-1),,,,]
   stck@catch.wt[,as.character(lastYear),,,,]     <- stck@catch.wt[,as.character(lastYear-1),,,,]
   stck@discards.wt[,as.character(lastYear),,,,]  <- stck@discards.wt[,as.character(lastYear-1),,,,]
   stck@landings.wt[,as.character(lastYear),,,,]  <- stck@landings.wt[,as.character(lastYear-1),,,,]
   stck@m[,as.character(lastYear),,,,]            <- stck@m[,as.character(lastYear-1),,,,]
   stck@landings.n[,as.character(lastYear),,,,]   <- stck@landings.n[,as.character(lastYear-1),,,,]
   stck@harvest.spwn[,as.character(lastYear),,,,] <- stck@harvest.spwn[,as.character(lastYear-1),,,,]
   stck@m.spwn[,as.character(lastYear),,,,]       <- stck@m.spwn[,as.character(lastYear-1),,,,]
  }

  #Generate observation matrix
  # Separate sumFleets and separate fleets
  if(dim(ctrl@sumFleets)[1] > 1){
    sumof <- unlist(strsplit(ctrl@sumFleets["sumFleets","sumof"],","))
    catch.dat.sum <- as.data.frame(areaSums(stck@catch.n[,ac(ctrl@sumFleets["sumFleets","startyrsum"]:ctrl@sumFleets["sumFleets","endyrsum"]),,,sumof],na.rm=T))
    catch.dat.fleet <- numeric()
    for(iArea in dimnames(stck@catch.n)$area)
      catch.dat.fleet <- rbind(catch.dat.fleet,as.data.frame(stck@catch.n[,ac(ctrl@sumFleets[paste("catch",iArea),"startyrfleet"]:ctrl@sumFleets[paste("catch",iArea),"endyrfleet"])]))
    catch.dat   <- rbind(catch.dat.sum,catch.dat.fleet)
  } else {
    catch.dat     <- as.data.frame(stck@catch.n)
  }
  catch.dat$qname <- paste("catch",catch.dat$area)
  tun.dat       <- as.data.frame(lapply(tun,index))
  obs.dat       <- rbind(catch.dat,tun.dat)
  obs.dat       <- subset(obs.dat,obs.dat$year %in% yrs)
  if(dim(ctrl@sumFleets)[1] >1)
    obs.dat$fleet <- as.numeric(factor(obs.dat$qname,levels=c(names(ctrl@fleets),"catch unique")))
  if(dim(ctrl@sumFleets)[1] == 1)
      obs.dat$fleet <- as.numeric(factor(obs.dat$qname,levels=names(ctrl@fleets)))
  obs.dat$age[which(ctrl@fleets[obs.dat$fleet]%in%c(3,4))] <- -1
  obs.dat       <- obs.dat[,c("year","fleet","age","data")]
  obs.dat       <- obs.dat[order(obs.dat$year,obs.dat$fleet,obs.dat$age),]
  if(dim(ctrl@sumFleets)[1] >1){
    obs.dat$fleet.name <- paste("#",c("catch unique",names(ctrl@fleets))[obs.dat$fleet],sep="")
  }
  if(dim(ctrl@sumFleets)[1] == 1){
   obs.dat$fleet.name <- paste("#",names(ctrl@fleets)[obs.dat$fleet],sep="")
  }
  #obs.dat       <- subset(obs.dat,!(obs.dat$data<=0 | is.na(obs.dat$data)))
  idx.start     <-which(!duplicated(obs.dat$year))
  idx.end       <-c(idx.start[-1]-1,nrow(obs.dat))
  nobs          <- nrow(obs.dat)

  if(dim(ctrl@sumFleets)[1] >1){
    idx1          <- matrix(NA,nrow=length(unique(obs.dat$fleet)),ncol=nyrs,dimnames=list(c("catch unique",names(ctrl@fleets)),yrs))
    idx2          <- matrix(NA,nrow=length(unique(obs.dat$fleet)),ncol=nyrs,dimnames=list(c("catch unique",names(ctrl@fleets)),yrs))
  }
  if(dim(ctrl@sumFleets)[1] == 1){
    idx1          <- matrix(NA,nrow=length(unique(obs.dat$fleet)),ncol=nyrs,dimnames=list(names(ctrl@fleets),yrs))
    idx2          <- matrix(NA,nrow=length(unique(obs.dat$fleet)),ncol=nyrs,dimnames=list(names(ctrl@fleets),yrs))
  }

  for(i in 1:nrow(obs.dat)){
    iFlt        <- obs.dat[i,"fleet"]
    iYr         <- obs.dat[i,"year"]
    idx1[iFlt,ac(iYr)] <- min(idx1[iFlt,ac(iYr)],i-1,na.rm=T)
    idx2[iFlt,ac(iYr)] <- max(idx2[iFlt,ac(iYr)],i-1,na.rm=T)
  }
  
  #Generate observation.var matrix
  if(is.null(catch.var)==F){
    if(!identical(dims(stck@catch.n),dims(catch.var))) stop("Dimensions of catch.var not identical to catch.n(stock)")
    if(dim(ctrl@sumFleets)[1] > 1){
      catch.var.sum <- as.data.frame(areaSums(catch.var[,ac(ctrl@sumFleets["sumFleets","startyrsum"]:ctrl@sumFleets["sumFleets","endyrsum"])],na.rm=T))
      catch.var.fleet <- numeric()
      for(iArea in dimnames(catch.var)$area)
        catch.var.fleet <- rbind(catch.var.fleet,as.data.frame(catch.var[,ac(ctrl@sumFleets[paste("catch",iArea),"startyrfleet"]:ctrl@sumFleets[paste("catch",iArea),"endyrfleet"])]))
      catch.var   <- rbind(catch.var.sum,catch.var.fleet)
    } else {
      catch.var     <- as.data.frame(catch.var)
    }
    catch.var$qname <- paste("catch",catch.var$area)
  } else {
    catch.var       <- catch.dat
    catch.var$data  <- NA
  }
  tun.var       <- as.data.frame(lapply(tun,index.var))
  obs.var       <- rbind(catch.var,tun.var)
  obs.var       <- subset(obs.var,obs.var$year %in% yrs)

  #---------------------------------------------------
  # Create data file
  #---------------------------------------------------

  #- data from the stock object
  dat             <- list()
  if(dim(ctrl@sumFleets)[1] > 1){
    dat$noFleets    <- length(ctrl@fleets)+1
    dat$fleetTypes  <- c(ctrl@fleets,7);      attr(dat$fleetTypes,"names")    <- NULL
    dat$sampleTimes <- c(samp.times,0);          attr(dat$sampleTimes,"names")   <- NULL
    dat$minAgePerFleet <- c(rep(range(stck)["min"],dims(stck)$area),do.call(rbind,lapply(tun,range))[,"min"],range(stck)["min"])
      dat$minAgePerFleet[is.na(dat$minAgePerFleet)] <- -1
    names(dat$minAgePerFleet)[length(dat$minAgePerFleet)] <- "catch unique"
    dat$maxAgePerFleet <- c(rep(range(stck)["max"],dims(stck)$area),do.call(rbind,lapply(tun,range))[,"max"],range(stck)["max"])
      dat$maxAgePerFleet[is.na(dat$maxAgePerFleet)] <- -1
    names(dat$maxAgePerFleet)[length(dat$maxAgePerFleet)] <- "catch unique"
  }
  if(dim(ctrl@sumFleets)[1] == 1){
    dat$noFleets    <- length(ctrl@fleets)
    dat$fleetTypes  <- c(ctrl@fleets);      attr(dat$fleetTypes,"names")    <- NULL
    dat$sampleTimes <- samp.times;          attr(dat$sampleTimes,"names")   <- NULL
    dat$minAgePerFleet <- c(rep(range(stck)["min"],dims(stck)$area),do.call(rbind,lapply(tun,range))[,"min"])
      dat$minAgePerFleet[is.na(dat$minAgePerFleet)] <- -1
    names(dat$minAgePerFleet)[1:(dims(stck)$area)] <- names(ctrl@fleets)[which(ctrl@fleets==0)]
    dat$maxAgePerFleet <- c(rep(range(stck)["max"],dims(stck)$area),do.call(rbind,lapply(tun,range))[,"max"])
      dat$maxAgePerFleet[is.na(dat$maxAgePerFleet)] <- -1
     names(dat$maxAgePerFleet)[1:(dims(stck)$area)] <- names(ctrl@fleets)[which(ctrl@fleets==0)]
  }
  dat$noYears     <- nyrs
  dat$years       <- yrs
  dat$nobs        <- nobs
  dat$idx1        <- idx1;                attr(dat$idx1,"dimnames")       <- NULL
  dat$idx2        <- idx2;                attr(dat$idx2,"dimnames")       <- NULL
  dat$minWeek     <- -1
  dat$maxWeek     <- -1
  if(length(ctrl@fleets==6)>0){
    idxSurvP <- which(names(tun) %in% names(ctrl@fleets[which(ctrl@fleets==6)]))
    supP     <- cumsum(unlist(lapply(tun[idxSurvP],function(x){return(max(as.integer(dimnames(x@index)$age))+1)})))
    minWeek <- c(1,supP[-length(supP)]+1)-1
    maxWeek <- supP-1
    names(minWeek) <- names(maxWeek)
    dat$minWeek <- minWeek
    dat$maxWeek <- maxWeek
  }
  dat$aux         <- data.matrix(obs.dat[,c("year","fleet","age")])
  dat$logobs      <- log(obs.dat$data)
  dat$logobs[is.infinite(dat$logobs)] <- NA
  dat$logobs[is.nan(dat$logobs)]      <- NA
  dat$weight      <- obs.var$data
  dat$propMat     <- t(stck@mat[,,,,1,drop=T]);
  dat$stockMeanWeight <- t(stck@stock.wt[,,,,1,drop=T])
  dms             <- dims(stck)
  dat$catchMeanWeight <- aperm(array(stck@catch.wt[,drop=T],dim=c(dms$age,dms$year,dms$area)),perm=c(2,1,3))
  dat$natMor      <- t(stck@m[,,,,1,drop=T])
  dat$landFrac    <- aperm(array((stck@landings.n/stck@catch.n)[,drop=T],dim=c(dms$age,dms$year,dms$area)),perm=c(2,1,3))
  dat$landFrac[is.na(dat$landFrac)] <- 1; dat$landFrac[is.infinite(dat$landFrac)] <- 0
  if(extraYr)
    dat$landFrac[nrow(dat$landFrac),,] <- dat$landFrac[nrow(dat$landFrac)-1,,]
  dat$disMeanWeight   <- aperm(array(stck@discards.wt[,drop=T],dim=c(dms$age,dms$year,dms$area)),perm=c(2,1,3))
  dat$landMeanWeight  <- aperm(array(stck@landings.wt[,drop=T],dim=c(dms$age,dms$year,dms$area)),perm=c(2,1,3))
  dat$propF       <- aperm(array(stck@harvest.spwn[,drop=T],dim=c(dms$age,dms$year,dms$area)),perm=c(2,1,3))
  dat$propM       <- t(stck@m.spwn[,,,,1,drop=T])

  #- data on the control file
  setNAtominone               <- function(x){for(i in slotNames(x)){idx <- which(is.na(slot(x,i)));if(length(idx)>0){if(class(slot(x,i))=="data.frame"){ y <- slot(x,i); y <- as.matrix(y); y[idx] <- -1; y <- as.data.frame(y,stringsAsFactors=F); slot(x,i)<-y } else {slot(x,i)[idx] <- -1}}};return(x)}
  sumFleets       <- ifelse(dim(ctrl@sumFleets)[1]>1,T,F)
  addRowForSumFleet           <- function(x,TF){x <- rbind(x,x[1,]); x[nrow(x),] <- -1; rownames(x)[nrow(x)] <- "catch unique"; ret <- if(TF){x} else {x[-nrow(x),]}; return(ret)}
  ctrlOrig        <- ctrl
  ctrl            <- setNAtominone(ctrl)
  dat$minAge      <- ctrl@range["min"];
  dat$maxAge      <- ctrl@range["max"];
  dat$maxAgePlusGroup <- ifelse(ctrl@plus.group,1,0);
  dat$keyLogFsta  <- addRowForSumFleet(ctrl@states,sumFleets);
  dat$corFlag     <- ctrl@cor.F
  dat$keyLogFpar  <- addRowForSumFleet(ctrl@catchabilities,sumFleets);
  dat$keyQpow     <- addRowForSumFleet(ctrl@power.law.exps,sumFleets);
  dat$keyVarF     <- addRowForSumFleet(ctrl@f.vars,sumFleets);
  dat$keyVarLogN  <- ctrl@logN.vars;
  dat$keyVarLogP  <- if(is.na(ctrl@logP.vars[1])){numeric(0)}else{ctrl@logP.vars}
  dat$keyVarObs   <- addRowForSumFleet(ctrl@obs.vars,sumFleets)
  dat$obsCorStruct<- factor(ctrl@cor.obs.Flag,levels=c("ID","AR","US"))
  dat$keyCorObs   <- addRowForSumFleet(ctrl@cor.obs,sumFleets)
  dat$stockRecruitmentModelCode <- ctrl@srr
  dat$noScaledYears <- ctrl@scaleNoYears
  if(dat$noScaledYears>0){
    dat$keyScaledYears  <- ctrl@scaleYears
    dat$keyParScaledYA <- ctrl@scalePars
  } else {
    dat$keyScaledYears <- numeric()
    dat$keyParScaledYA <- matrix(numeric(),nrow=0,ncol=0)
  }
  dat$fbarRange   <- c(ctrl@range[c("minfbar","maxfbar")])
  dat$keyBiomassTreat <- rep(-1,dat$noFleets)
  dat$keyBiomassTreat[which(dat$fleetTypes %in% c(3,4))] <- ctrl@biomassTreat[which(ctrl@fleets %in% c(3,4))]
  dat$obsLikelihoodFlag <- factor(ctrl@likFlag,levels=c("LN","ALN"))
  dat$fixVarToWeight <- 0
  dat$simFlag     <- ifelse(ctrl@simulate,1,0)
  dat$resFlag     <- ifelse(ctrl@residuals,1,0)
  dat$sumKey      <- matrix(0,nrow=nrow(dat$keyLogFpar),ncol=ncol(dat$keyLogFpar))
  if(sumFleets)
    dat$sumKey[nrow(dat$sumKey),which(names(ctrl@fleets) %in% paste("catch",unlist(strsplit(ctrl@sumFleets["sumFleets","sumof"],","))))] <- 1
  ctrl            <- ctrlOrig

  #---------------------------------------------------
  # Create initialisation file
  #---------------------------------------------------
  parameters <- list(logFpar = numeric(), logQpow = numeric(), logSdLogFsta = numeric(), logSdLogN = numeric(), logSdLogP = numeric(),
                     logSdLogObs = numeric(), logSdLogTotalObs = numeric(), transfIRARdist = numeric(), sigmaObsParUS = numeric(),
                     rec_loga = numeric(), rec_logb = numeric(), itrans_rho = numeric(), rhop = numeric(), logScale = numeric(), logitReleaseSurvival = numeric(),
                     logitRecapturePhi = numeric(), logAlphaSCB = numeric(), logF = numeric(), logN = numeric(),logPS=numeric())
  if(length(unique(na.omit(c(ctrl@catchabilities)))) > 0)
    parameters$logFpar        <- rep(-5,length(unique(na.omit(c(ctrl@catchabilities)))))
  if(length(unique(na.omit(c(ctrl@power.law.exps)))) > 0)
    parameters$logQpow        <- rep(1,unique(na.omit(c(ctrl@power.law.exps))))
  if(length(unique(na.omit(c(ctrl@f.vars)))) > 0)
    parameters$logSdLogFsta   <- rep(-0.7,length(unique(na.omit(c(ctrl@f.vars)))))
  if(length(unique(na.omit(c(ctrl@logN.vars)))) > 0)
    parameters$logSdLogN      <- rep(-0.35,length(unique(na.omit(c(ctrl@logN.vars)))))
  if(length(unique(na.omit(c(ctrl@logP.vars)))) > 0)
    parameters$logSdLogP      <- rep(-0.7,length(unique(na.omit(c(ctrl@logP.vars)))))
  if(length(unique(na.omit(c(ctrl@obs.vars)))) > 0)
    parameters$logSdLogObs    <- rep(-0.35,length(unique(na.omit(c(ctrl@obs.vars)))))
  if(any(ctrl@likFlag %in% "ALN"))
    parameters$logSdLogTotalObs <- 0
  if(length(unique(na.omit(c(ctrl@cor.obs))))>0){
    if(!any(ctrl@cor.obs.Flag=="AR")) warning("Trying to estimate cor.obs without cor.obs.Flag set")
    parameters$transfIRARdist <- rep(0.05,length(unique(na.omit(c(ctrl@cor.obs)))))
  }
  if(any(ctrl@cor.obs.Flag %in% "US")){
    idx                       <- which(ctrl@cor.obs.Flag == "US")
    len                       <- length(na.omit(c((ctrl@cor.obs.Flag[idx]))))
    parameters$sigmaObsParUS  <- rep(0,len*(len-1)/2)
  }
  if(ctrl@srr > 0){
    parameters$rec_loga       <- 1
    parameters$rec_logb       <- 1
  }
  parameters$itrans_rho     <- unlist(lapply(as.list(ctrl@cor.F),function(x){if(x==0){ ret <- numeric()} else { ret <- numeric(1)+.5}; return(ret)}))
  if(length(unique(na.omit(c(ctrl@logP.vars)))) > 0)
    parameters$rhop           <- 0.5
  if(length(unique(na.omit(c(ctrl@scalePars))))>0)
    parameters$logScale       <- rep(0,length(unique(na.omit(c(ctrl@scalePars)))))
  if(length(unique(na.omit(c(ctrl@logP.vars)))) > 0)
    parameters$logAlphaSCB    <- unlist(lapply(mapply(seq,dat$minWeek,dat$maxWeek,SIMPLIFY=F),function(x){log(rep(1/length(x),length(x)-1))}))
  parameters$logF             <- matrix(0,nrow=length(unique(na.omit(c(ctrl@states)))),ncol=nyrs)
  parameters$logN             <- matrix(0,nrow=length(dat$minAge:dat$maxAge),ncol=nyrs)
  if(length(unique(na.omit(c(ctrl@logP.vars)))) > 0){
    idxPart <- which(dat$fleetTypes==6)
    colP    <- length(unique(dat$aux[which(dat$aux[,"fleet"] %in% idxPart),"year"]))
    parameters$logPS=matrix(0, nrow=length(idxPart)-1, ncol=colP)
  }

  return(list(data=dat,parameters=parameters,map=map))
}

.get.stem <-function(ctrl) {
  if(length(ctrl@sam.binary)==0) {
	.stem <- "stockassessment" } else {
	.stem <- gsub("\\.exe$","",basename(ctrl@sam.binary)) }
}
