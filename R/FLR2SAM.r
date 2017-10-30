FLR2SAM <-function(stck,tun,ctrl,pin.sam=NULL,map=NULL) {
  #---------------------------------------------------
  # Setup for output
  #---------------------------------------------------
  #General Setup
  sam.stem <- .get.stem(ctrl)
  run.time <- Sys.time()
  miss.val <- 0

  #Setup meta data
  samp.times <- c(miss.val,sapply(tun,function(x) mean(x@range[c("startf","endf")])))
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
  catch.dat     <- as.data.frame(FLQuants(catch=stck@catch.n))
  tun.dat       <- as.data.frame(lapply(tun,index))
  obs.dat       <- rbind(catch.dat,tun.dat)
  obs.dat       <- subset(obs.dat,obs.dat$year %in% yrs)
  obs.dat$fleet <- as.numeric(factor(obs.dat$qname,levels=names(ctrl@fleets)))
  obs.dat$age[which(ctrl@fleets[obs.dat$fleet]%in%c(3,4))] <- -1
  obs.dat       <- obs.dat[,c("year","fleet","age","data")]
  obs.dat       <- obs.dat[order(obs.dat$year,obs.dat$fleet,obs.dat$age),]
  obs.dat$fleet.name <- paste("#",names(ctrl@fleets)[obs.dat$fleet],sep="")
  #obs.dat       <- subset(obs.dat,!(obs.dat$data<=0 | is.na(obs.dat$data)))
  idx.start     <-which(!duplicated(obs.dat$year))
  idx.end       <-c(idx.start[-1]-1,nrow(obs.dat))
  nobs          <- nrow(obs.dat)

  idx1          <- matrix(NA,nrow=length(unique(obs.dat$fleet)),ncol=nyrs,dimnames=list(names(ctrl@fleets),yrs))
  idx2          <- matrix(NA,nrow=length(unique(obs.dat$fleet)),ncol=nyrs,dimnames=list(names(ctrl@fleets),yrs))
  for(i in 1:nrow(obs.dat)){
    iFlt        <- obs.dat[i,"fleet"]
    iYr         <- obs.dat[i,"year"]
    idx1[iFlt,ac(iYr)] <- min(idx1[iFlt,ac(iYr)],i-1,na.rm=T)
    idx2[iFlt,ac(iYr)] <- max(idx2[iFlt,ac(iYr)],i-1,na.rm=T)
  }
  

  #---------------------------------------------------
  # Create data file
  #---------------------------------------------------

  #- data from the stock object
  dat             <- list()
  dat$noFleets    <- length(ctrl@fleets)
  dat$fleetTypes  <- c(ctrl@fleets);      attr(dat$fleetTypes,"names")    <- NULL
  dat$sampleTimes <- samp.times;          attr(dat$sampleTimes,"names")   <- NULL
  dat$noYears     <- nyrs
  dat$years       <- yrs
  dat$minAgePerFleet <- c(catch=c(range(stck)["min"]),do.call(rbind,lapply(tun,range))[,"min"])
    dat$minAgePerFleet[is.na(dat$minAgePerFleet)] <- -1
  dat$maxAgePerFleet <- c(catch=c(range(stck)["max"]),do.call(rbind,lapply(tun,range))[,"max"])
    dat$maxAgePerFleet[is.na(dat$maxAgePerFleet)] <- -1
  dat$nobs        <- nobs
  dat$idx1        <- idx1;                attr(dat$idx1,"dimnames")       <- NULL
  dat$idx2        <- idx2;                attr(dat$idx2,"dimnames")       <- NULL
  dat$aux         <- data.matrix(obs.dat[,c("year","fleet","age")])
  dat$logobs      <- log(obs.dat$data)
  dat$logobs[is.infinite(dat$logobs)] <- NA
  dat$logobs[is.nan(dat$logobs)]      <- NA
  weight          <- numeric(nobs)
  for(i in 1:nobs)
    weight[i]     <- ctrl@obs.weight[dat$aux[i,2],ac(abs(dat$aux[i,3]))]
  dat$weight      <- weight
  dat$propMat     <- t(stck@mat[,drop=T]);
  dat$stockMeanWeight <- t(stck@stock.wt[,drop=T])
  dat$catchMeanWeight <- t(stck@catch.wt[,drop=T])
  dat$natMor      <- t(stck@m[,drop=T])
  dat$landFrac    <- t((stck@landings.n/stck@catch.n)[,drop=T])
  if(extraYr)
    dat$landFrac[nrow(dat$landFrac),] <- dat$landFrac[nrow(dat$landFrac)-1,]
  dat$disMeanWeight   <- t(stck@discards.wt[,drop=T])
  dat$landMeanWeight  <- t(stck@landings.wt[,drop=T])
  dat$propF       <- t(stck@harvest.spwn[,drop=T])
  dat$propM       <- t(stck@m.spwn[,drop=T])

  #- data on the control file
  setNAtominone               <- function(x){for(i in slotNames(x)){idx <- which(is.na(slot(x,i)));if(length(idx)>0){slot(x,i)[idx] <- -1}};return(x)}
  ctrlOrig        <- ctrl
  ctrl            <- setNAtominone(ctrl)
  dat$minAge      <- ctrl@range["min"];   attr(dat$minAge,"names")            <- NULL
  dat$maxAge      <- ctrl@range["max"];   attr(dat$maxAge,"names")            <- NULL
  dat$maxAgePlusGroup <- ifelse(ctrl@plus.group,1,0);   attr(dat$maxAgePlusGroup,"names")        <- NULL
  dat$keyLogFsta  <- ctrl@states;         attr(dat$keyLogFsta,"dimnames")     <- NULL
  dat$corFlag     <- ifelse(is.logical(ctrl@cor.F),ifelse(ctrl@cor.F,1,0),ctrl@cor.F)
  dat$keyLogFpar  <- ctrl@catchabilities; attr(dat$keyLogFpar,"dimnames")     <- NULL
  dat$keyQpow     <- ctrl@power.law.exps; attr(dat$keyQpow,"dimnames")        <- NULL
  dat$keyVarF     <- ctrl@f.vars;         attr(dat$keyVarF,"dimnames")        <- NULL
  dat$keyVarLogN  <- ctrl@logN.vars;      attr(dat$keyVarLogN,"dimnames")     <- NULL
  dat$keyVarObs   <- ctrl@obs.vars;       attr(dat$keyVarObs,"dimnames")      <- NULL
  dat$obsCorStruct<- factor(ctrl@cor.obs.Flag,levels=c("ID","AR","US"))
  dat$keyCorObs   <- ctrl@cor.obs
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
  dat$keyBiomassTreat[which(dat$fleetTypes %in% c(3,4))] <- 0
  dat$obsLikelihoodFlag <- factor(ctrl@likFlag,levels=c("LN","ALN"))
  dat$fixVarToWeight <- 0
  dat$simFlag     <- ifelse(ctrl@simulate,1,0)
  dat$resFlag     <- 0
  ctrl            <- ctrlOrig

  #---------------------------------------------------
  # Create initialisation file
  #---------------------------------------------------
  parameters <- list(logFpar = numeric(), logQpow = numeric(), logSdLogFsta = numeric(), logSdLogN = numeric(),
                     logSdLogObs = numeric(), logSdLogTotalObs = numeric(), transfIRARdist = numeric(), sigmaObsParUS = numeric(),
                     rec_loga = numeric(), rec_logb = numeric(), itrans_rho = numeric(), logScale = numeric(), logitReleaseSurvival = numeric(),
                     logitRecapturePhi = numeric(), logF = numeric(), logN = numeric())
  if(length(unique(na.omit(c(ctrl@catchabilities)))) > 0)
    parameters$logFpar        <- rep(-5,length(unique(na.omit(c(ctrl@catchabilities)))))
  if(length(unique(na.omit(c(ctrl@power.law.exps)))) > 0)
    parameters$logQpow        <- rep(1,unique(na.omit(c(ctrl@power.law.exps))))
  if(length(unique(na.omit(c(ctrl@f.vars)))) > 0)
    parameters$logSdLogFsta   <- rep(-0.7,length(unique(na.omit(c(ctrl@f.vars)))))
  if(length(unique(na.omit(c(ctrl@logN.vars)))) > 0)
    parameters$logSdLogN      <- rep(-0.35,length(unique(na.omit(c(ctrl@logN.vars)))))
  if(length(unique(na.omit(c(ctrl@obs.vars)))) > 0)
    parameters$logSdLogObs    <- rep(-0.35,length(unique(na.omit(c(ctrl@obs.vars)))))
  if(any(ctrl@likFlag %in% "ALN"))
    parameters$logSdLogTotalObs <- 0
  if(length(unique(na.omit(c(ctrl@cor.obs))))>0)
    parameters$transfIRARdist <- rep(0.05,length(unique(na.omit(c(ctrl@cor.obs)))))
  if(any(ctrl@cor.obs.Flag %in% "US")){
    idx                       <- which(ctrl@cor.obs.Flag == "US")
    len                       <- length(na.omit(c((ctrl@cor.obs.Flag[idx]))))
    parameters$sigmaObsParUS  <- rep(0,len*(len-1)/2)
  }
  if(ctrl@srr > 0){
    parameters$rec_loga       <- 1
    parameters$rec_logb       <- 1
  }
  if(ctrl@cor.F)
    parameters$itrans_rho     <- 0.5
  if(length(unique(na.omit(c(ctrl@scalePars))))>0)
    parameters$logScale       <- rep(0,length(unique(na.omit(c(ctrl@scalePars)))))
  parameters$logF             <- matrix(0,nrow=length(unique(na.omit(c(ctrl@states)))),ncol=nyrs)
  parameters$logN             <- matrix(0,nrow=length(dat$minAge:dat$maxAge),ncol=nyrs)

  #---------------------------------------------------
  # If pin file is given, write files to disk
  #---------------------------------------------------
  if(is.null(pin.sam)==FALSE){
    parameters                <- createPin(pin.sam,ctrl)
  }

  return(list(data=dat,parameters=parameters,map=map))
}

.get.stem <-function(ctrl) {
  if(length(ctrl@sam.binary)==0) {
	.stem <- "sam" } else {
	.stem <- gsub("\\.exe$","",basename(ctrl@sam.binary)) }
}
