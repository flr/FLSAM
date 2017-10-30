SAM2FLR <- function(fit,ctrl){

  #- Create new FLSAM object and fill with ctrl elements
  res         <- new("FLSAM")
  res@control <- ctrl
  res@nohess  <- ctrl@nohess
  if(!ctrl@nohess)
    res@nopar <- length(fit$sdrep$par.fixed)
  res@n.states<- length(unique(na.omit(c(ctrl@states))))+length(seq(ctrl@range["min"],ctrl@range["max"]))
  res@states  <- ctrl@states
  res@nlogl   <- fit$opt$objective
  #- Get overall vcov and the cov of the observations
  if(!ctrl@nohess){
    res@vcov  <- fit$sdrep$cov.fixed
    res@rescov<- fit$sdrep$cov
  }
  if(length(fit$conf$transfIRARdist)>0){
    res@obscov  <- fit$rep$obsCov[apply(ctrl@cor.obs,1,function(x){any(x>=0,na.rm=T)})]
    names(res@obscov) <- rownames(ctrl@cor.obs)[apply(ctrl@cor.obs,1,function(x){any(x>=0,na.rm=T)})]
  }
  #- get parameters + sd value
  if(!ctrl@nohess){
    parssummary <- rbind(data.frame(name=ac(rownames(as.data.frame(unlist(fit$pl)))),
                                    value=           as.vector(unlist(fit$pl)),
                                    std.dev=         as.vector(unlist(fit$plsd))),
                         data.frame(name=ac(names(fit$sdrep$value)),value=fit$sdrep$value,std.dev=fit$sdrep$sd))
    parssummary$name <- ac(parssummary$name)
    colnames(parssummary) <- c("name","value","std.dev")
    parssummary[,1]       <- gsub('[0-9]+', '',parssummary[,1])
    res@params  <- parssummary
  }
  if(ctrl@nohess)
    res@params  <- cbind(name=names(fit$opt$par),value=as.data.frame(fit$opt$par),std.dev=NA)

  #- Fill stock and harvest
  if(!ctrl@nohess){
    res@stock.n <- FLQuant(exp(fit$pl$logN),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                          year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                          unit="unique",
                                                          season="all",
                                                          area="unique",
                                                          iter=1))
    res@harvest <- FLQuant(exp(fit$pl$logF[ctrl@states["catch",]+1,]),dimnames=list(age=ctrl@range["min"]:ctrl@range["max"],
                                                          year=ctrl@range["minyear"]:ctrl@range["maxyear"],
                                                          unit="unique",
                                                          season="all",
                                                          area="unique",
                                                          iter=1))
    units(res@harvest) <- "f"
  }
  #- Get all the residuals
  if(ctrl@residuals){
    res@residuals <- data.frame(year=fit$data$aux[,"year"],
                                fleet=names(ctrl@fleets)[fit$data$aux[,"fleet"]],
                                age=fit$data$aux[,"age"],
                                log.obs=fit$data$logobs,
                                log.mdl=fit$rep$predObs,
                                std.res=oneStepPredict(fit$obj,observation.name="logobs",data.term.indicator="keep",discrete=FALSE)[,"residual"])
    res@residuals$log.obs[which(is.na(res@residuals$log.obs))] <- fit$pl$missing
  }
  #- Finally some info and save
  info        <- data.frame(FLSAM.version=packageDescription("FLSAM")$Version,
                            FLCore.version=packageDescription("FLCore")$Version,
                            R.version=R.version$version.string,
                            platform=R.version$platform,
                            run.date=Sys.time())
  res@info    <- t(info)
  res@name    <- ctrl@name
  res@desc    <- ctrl@desc
  res@range   <- ctrl@range
return(res)}




