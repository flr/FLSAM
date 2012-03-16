FLR2SAM <-function(stck,tun,ctrl,run.dir="missing") {
  #---------------------------------------------------
  # Setup for output
  #---------------------------------------------------
  #General Setup
  if(missing(run.dir)) {run.dir <- tempdir() }
  admb.stem <- "sam" 
  run.time <- Sys.time()
  miss.val <- -99999

  #Internal Helper functions
  .format.matrix.ADMB <- function(mat,na.replace="missing") {
     if(na.replace!="missing") {mat[is.na(mat)] <- na.replace}
     colnames(mat)[1] <- paste("#",colnames(mat)[1])
     tbl <- capture.output(write.table(mat, row.names=FALSE, col.names=TRUE, quote=FALSE))
     tbl <- paste(tbl,"#",c(" ",rownames(mat)),"\n")
     return(tbl)
  }
  .file.header <- function(fname,ftime) {
     cat("# Auto generated file\n", file=fname)
     cat(sprintf("# Datetime : %s\n\n",ftime),file=fname,append=TRUE)
   }
  .flqout <- function(desc,flq,na.replace="missing") { #Local export function
                  cat("#",desc,"\n",file=dat.file, append=TRUE)
                  cat(.format.matrix.ADMB(t(flq[,,drop=TRUE]@.Data),na.replace=na.replace), 
                     file=dat.file, append=TRUE)
                  return(invisible(NULL))}

  #Setup meta data
  samp.times <- c(miss.val,sapply(tun,function(x) mean(x@range[c("startf","endf")])))
  samp.times <- ifelse(is.na(samp.times),miss.val,samp.times)
  names(samp.times) <- names(ctrl@fleets)
  yrs <- ctrl@range["minyear"]:ctrl@range["maxyear"]
  lastYear <- ctrl@range["maxyear"]
  nyrs <- length(yrs)
  
  #Add a year to the stock object when tuning data are newer than catch data
  #Other combinations are not taken into account. Credit to Morten Vinther for this patch
  if (stck@range["maxyear"]<lastYear)  {  
   stck<- window(stck,start=stck@range["minyear"],end=ctrl@range["maxyear"],
              frequency=1,extend=TRUE)  # extend by one year
   stck@mat[,as.character(lastYear),,,,]<-stck@mat[,as.character(lastYear-1),,,,]      # MV (there must be an easier way!!
   stck@stock.wt[,as.character(lastYear),,,,]<-stck@stock.wt[,as.character(lastYear-1),,,,]
   stck@catch.wt[,as.character(lastYear),,,,]<-stck@catch.wt[,as.character(lastYear-1),,,,]
   stck@discards.wt[,as.character(lastYear),,,,]<-stck@discards.wt[,as.character(lastYear-1),,,,]
   stck@landings.wt[,as.character(lastYear),,,,]<-stck@landings.wt[,as.character(lastYear-1),,,,]
   stck@m[,as.character(lastYear),,,,]<-stck@m[,as.character(lastYear-1),,,,]
   stck@landings.n[,as.character(lastYear),,,,]<-stck@landings.n[,as.character(lastYear-1),,,,]
   stck@harvest.spwn[,as.character(lastYear),,,,]<-stck@harvest.spwn[,as.character(lastYear-1),,,,]
   stck@m.spwn[,as.character(lastYear),,,,]<-stck@m.spwn[,as.character(lastYear-1),,,,]
  }

  #Generate observation matrix
  catch.dat <- as.data.frame(FLQuants(catch=stck@catch.n))
  tun.dat <- as.data.frame(lapply(tun,index))
  obs.dat <- rbind(catch.dat,tun.dat)
  obs.dat <- subset(obs.dat,obs.dat$year %in% yrs)
  obs.dat$fleet <- as.numeric(factor(obs.dat$qname,levels=names(ctrl@fleets)))
  obs.dat$age[which(ctrl@fleets[obs.dat$fleet]%in%c(3,4))] <- ctrl@range["max"]
  obs.dat <- obs.dat[,c("year","fleet","age","data")]
  obs.dat <- obs.dat[order(obs.dat$year,obs.dat$fleet,obs.dat$age),]
  obs.dat$fleet.name <- paste("#",names(ctrl@fleets)[obs.dat$fleet],sep="")
  obs.dat <- subset(obs.dat,!(obs.dat$data<=0 | is.na(obs.dat$data)))
  idx.start <-which(!duplicated(obs.dat$year))
  idx.end   <-c(idx.start[-1]-1,nrow(obs.dat))
  nobs  <- nrow(obs.dat)

  #---------------------------------------------------
  # Create data file
  #---------------------------------------------------
  # Now write the data file!
  dat.file <- file.path(run.dir,sprintf("%s.dat",admb.stem))
  .file.header(dat.file,run.time)

  #Write meta data
  cat("# Number of fleets (res+con+sur)\n",length(ctrl@fleets),"\n", file=dat.file, append=TRUE)
  cat("# Fleet types (res=0, con=1, sur=2, ssb=3)\n",.format.matrix.ADMB(t(ctrl@fleets)),file=dat.file, append=TRUE)
  cat("# Sample times (only relevent for sur)\n",.format.matrix.ADMB(t(samp.times)), file=dat.file, append=TRUE)
  cat("# Number of years\n",nyrs,"\n", file=dat.file, append=TRUE)
  cat("# Years\n",yrs,"\n", file=dat.file, append=TRUE)
  cat("# Number of observations \n",nobs,"\n", file=dat.file, append=TRUE)
  cat("# Index1 (index of first obs in each year) \n",idx.start,"\n", file=dat.file, append=TRUE)
  cat("# Index2 (index of last obs in each year) \n",idx.end,"\n", file=dat.file, append=TRUE)

  #Observations
  cat("# The observation matrix \n#",colnames(obs.dat),"\n", file=dat.file, append=TRUE)
  write.table(obs.dat, row.names=FALSE, col.names=FALSE, quote=FALSE, file=dat.file, append=TRUE)

  #Now write the rest of the stock information
  .flqout("Proportion mature",stck@mat)
  .flqout("Stock mean weights",stck@stock.wt)
  .flqout("Catch mean weights",stck@catch.wt)
  .flqout("Natural Mortality", stck@m)
  .flqout("Landing Fraction L/(L+D)",
       stck@landings.n/(stck@landings.n + stck@discards.n),na.replace=1)
  .flqout("Catch mean weights DISCARD",stck@discards.wt)
  .flqout("Catch mean weights LANDINGS",stck@landings.wt)
  .flqout("Fprop",stck@harvest.spwn)
  .flqout("Mprop",stck@m.spwn)
  
  #Finally, write the checksums
  cat("# Checksums to ensure correct reading of input data \n",42,42,"\n", file=dat.file, append=TRUE)
  
  #---------------------------------------------------
  # Create configuration file
  #---------------------------------------------------
  #Write configuration file
  cfg.file <- file.path(run.dir,"model.cfg")
  .file.header(cfg.file,run.time)

  #Write the headers
  cat("# Min, max age represented internally in model \n",ctrl@range[c("min","max")],"\n", file=cfg.file, append=TRUE)
  cat("# Max age considered a plus group? (0 = No, 1= Yes)\n",as.numeric(ctrl@plus.group[1]),"\n", file=cfg.file, append=TRUE)
  
  #Coupling Matrices
  cat("\n# Coupling of fishing mortality STATES (ctrl@states)\n",.format.matrix.ADMB(ctrl@states,na.replace=0),file=cfg.file,append=TRUE)
  cat("\n# Use correlated random walks for the fishing mortalities\n# ( 0 = independent, 1 = correlation estimated)\n",as.integer(ctrl@cor.F),file=cfg.file,append=TRUE)
  cat("\n\n# Coupling of catchability PARAMETERS (ctrl@catchabilities)\n",.format.matrix.ADMB(ctrl@catchabilities,na.replace=0),file=cfg.file,append=TRUE)
  cat("\n# Coupling of power law model EXPONENTS (ctrl@power.law.exps)\n",.format.matrix.ADMB(ctrl@power.law.exps,na.replace=0),file=cfg.file,append=TRUE)
  cat("\n# Coupling of fishing mortality RW VARIANCES (ctrl@f.vars)\n",.format.matrix.ADMB(ctrl@f.vars,na.replace=0),file=cfg.file,append=TRUE)
  cat("\n# Coupling of log N RW VARIANCES (ctrl@logN.vars)\n",ctrl@logN.vars,file=cfg.file,append=TRUE)
  cat("\n\n# Coupling of OBSERVATION VARIANCES (ctrl@obs.vars)\n",.format.matrix.ADMB(ctrl@obs.vars,na.replace=0),file=cfg.file,append=TRUE)
  
  #Final values
  cat("\n# Stock recruitment model code (0=RW, 1=Ricker, 2=BH, ... more in time\n",ctrl@srr,"\n",file=cfg.file,append=TRUE)
  cat("# Years in which catch data are to be scaled by an estimated parameter \n",0,"\n",file=cfg.file,append=TRUE)
  cat("# Fbar range \n",ctrl@range[c("minfbar","maxfbar")],"\n",file=cfg.file,append=TRUE)
  cat("# Model timeout \n",ctrl@timeout,"\n",file=cfg.file,append=TRUE)
  
  #Finally, write the checksums
  cat("\n\n# Checksums to ensure correct reading of input data \n",123456,123456,"\n", file=cfg.file, append=TRUE)

  #---------------------------------------------------
  # Create initialisation file
  #---------------------------------------------------
  init.file <- file.path(run.dir,"model.init")
  .file.header(init.file,run.time)
  cat("#varLogFsta\n 0.5\n",file=init.file,append=TRUE)
  cat("#varLogN\n0.5 \n",file=init.file,append=TRUE)
  cat("#varLogObs \n 0.5\n",file=init.file,append=TRUE)
  cat("#logFpar \n 0.3 \n",file=init.file,append=TRUE)
  cat("#rec_loga \n 1 \n",file=init.file,append=TRUE)
  cat("#rec_logb \n -12 \n",file=init.file,append=TRUE)
  cat("\n\n# Checksums to ensure correct reading of input data \n",90210, 90210,"\n", file=init.file, append=TRUE)


  #---------------------------------------------------
  # Create reduced run file
  #---------------------------------------------------
  red.file <- file.path(run.dir,"reduced.cfg")
  .file.header(red.file,run.time)
  cat("# This allows to exclude data in ways customized for retrospective runs\n",
      "# The following specifies one integer for each fleet (catches and surveys)\n",
      "# if:\n",
      "#     0: all data for that fleet is used\n",
      "#    -1: all data for that fleet is excluded\n",
      "#     n: n is a positive integer the corresponding fleet's data is reduced\n",
      "#        by n years starting from the most recent year.\n",
      rep(0,length(ctrl@fleets)),"\n", 
      "\n\n# Checksums to ensure correct reading of input data \n",3142, 3142,"\n",
      file=red.file,append=TRUE) 
  return(invisible(NULL))
}

