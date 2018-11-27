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
  oldTunNames <- names(tun)
  newTunNames <-lapply(as.list(names(tun)),function(x){gsub(" ","",x)})
  names(tun) <- newTunNames

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
    if(is.na(dmnsFLR[[iSurv]]$min) & dimnames(tun[[iSurv]]@index)$age[1]=="all"){
      colnames(surveysFLR[[iSurv]]) <- -1
    } else {
      colnames(surveysFLR[[iSurv]]) <- seq(dmnsFLR[[iSurv]]$min,dmnsFLR[[iSurv]]$max,1)
    }
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

