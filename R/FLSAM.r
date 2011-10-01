FLSAM <-function(stck,tun,ctrl,run.dir="missing") {
  #---------------------------------------------------
  # Setup for output
  #---------------------------------------------------
  #General Setup
  if(missing(run.dir)) {run.dir <- tempdir() }
  admb.stem <- "ssass" 
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

  #Setup meta data
  fleet.names <- c(catch="catch",sapply(tun,name))
  fleet.types <- factor(sapply(tun,type),levels=c("con","number","biomass"))
  fleet.types <- c(0,as.numeric(fleet.types))
  names(fleet.types) <- fleet.names
  samp.times <- c(miss.val,sapply(tun,function(x) mean(x@range[c("startf","endf")])))
  samp.times <- ifelse(is.na(samp.times),miss.val,samp.times)
  names(samp.times) <- fleet.names
  yrs <- seq(stck@range["minyear"],stck@range["maxyear"])
  nyrs <- length(yrs)

  #Generate observation matrix
  obs.dat <- as.data.frame(lapply(tun,index))
  obs.dat <- rbind(obs.dat,as.data.frame(FLQuants(catch=stck@catch.n)))
  obs.dat <- subset(obs.dat,obs.dat$year %in% yrs)
  obs.dat$fleet <- as.numeric(factor(obs.dat$qname,levels=fleet.names))
  obs.dat$age[which(fleet.types[obs.dat$fleet]==3)] <- median(as.numeric(obs.dat$age),na.rm=TRUE)   #Set ssb indices equal to median age
  obs.dat <- obs.dat[,c("year","fleet","age","data")]
  obs.dat <- obs.dat[order(obs.dat$year,obs.dat$fleet,obs.dat$age),]
  obs.dat$fleet.name <- paste("#",fleet.names[obs.dat$fleet],sep="")
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
  cat("# Number of fleets (res+con+sur)\n",length(fleet.names),"\n", file=dat.file, append=TRUE)
  cat("# Fleet types (res=0, con=1, sur=2, ssb=3)\n",.format.matrix.ADMB(t(fleet.types)),file=dat.file, append=TRUE)
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
  flqout <- function(desc,flq) { #Local export function
                  cat("#",desc,"\n",file=dat.file, append=TRUE)
                  cat(.format.matrix.ADMB(t(flq[,,drop=TRUE]@.Data)), file=dat.file, append=TRUE)
                  return(invisible(NULL))}
  flqout("Proportion mature",stck@mat)
  flqout("Stock mean weights",stck@stock.wt)
  flqout("Catch mean weights",stck@catch.wt)
  flqout("Natural Mortality", stck@m)
  flqout("Landing Fraction L/(L+D)",stck@landings.n*0+1)
  flqout("Fprop",stck@harvest.spwn)
  flqout("Mprop",stck@m.spwn)
  
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
  cat("\n# Coupling of catchability PARAMETERS (ctrl@catchabilities)\n",.format.matrix.ADMB(ctrl@catchabilities,na.replace=0),file=cfg.file,append=TRUE)
  cat("\n# Coupling of power law model EXPONENTS (ctrl@power.law.exps)\n",.format.matrix.ADMB(ctrl@power.law.exps,na.replace=0),file=cfg.file,append=TRUE)
  cat("\n# Coupling of fishing mortality RW VARIANCES (ctrl@f.vars)\n",.format.matrix.ADMB(ctrl@f.vars,na.replace=0),file=cfg.file,append=TRUE)
  cat("\n# Coupling of log N RW VARIANCES (ctrl@logN.vars)\n",ctrl@logN.vars,file=cfg.file,append=TRUE)
  cat("\n\n# Coupling of OBSERVATION VARIANCES (ctrl@obs.vars)\n",.format.matrix.ADMB(ctrl@obs.vars,na.replace=0),file=cfg.file,append=TRUE)
  
  #Final values
  cat("\n# Stock recruitment model code (0=RW, 1=Ricker, 2=BH, ... more in time\n",ctrl@srr,"\n",file=cfg.file,append=TRUE)
  cat("# Years in which catch data are to be scaled by an estimated parameter (mainly cod related)\n",0,"\n",file=cfg.file,append=TRUE)
  cat("# Fbar range \n",ctrl@range[c("minfbar","maxfbar")],"\n",file=cfg.file,append=TRUE)
  
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
  cat("#logFpar \n, 0.3 \n",file=init.file,append=TRUE)
  cat("#rec_loga \n, 1 \n",file=init.file,append=TRUE)
  cat("#rec_logb \n, -12 \n",file=init.file,append=TRUE)

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
      "0 0 0 0 0\n", file=red.file,append=TRUE)

  #---------------------------------------------------
  # We're ready! Run the executable
  #---------------------------------------------------
  if(!ctrl@simulate) {
    admb.args <-  "-nr 2 -noinit -iprint 1"
    #Platform specific issues
    if (.Platform$OS.type=="unix") {
      admb.exec <- file.path(system.file("bin", "linux", package="FLSAM", 
                     mustWork=TRUE), admb.stem)
      file.copy(admb.exec, run.dir)
    } else if (.Platform$OS.type == "windows") {
      admb.exec <- file.path(system.file("bin", "windows", package="FLSAM", mustWork=TRUE),
        sprintf("%s.exe",admb.stem))
      file.copy(admb.exec, run.dir)
    } else {
      stop(sprintf("Platform type, %s, is not currently supported.",R.version$os))
    }

    #Run!
    cmd <- sprintf("./%s -nr 2 -noinit -iprint 1" ,basename(admb.exec))
    olddir <- setwd(run.dir)
    rtn <- system(cmd)
    setwd(olddir)
    if(rtn!=0) {
      stop(sprintf("An error occurred while running ADMB. Return code %s.",rtn))
    }
   } else {
    cat("Simulated run.\n")
   }

  #---------------------------------------------------
  # Now read the results from the assessment
  #---------------------------------------------------
  #Create return object
  res       <- new("FLSAM")
  res@name  <- ctrl@name
  res@desc  <- ctrl@desc
  
  #Read parameter file
  par.fname    <- file.path(run.dir,sprintf("%s.par",admb.stem))
  parfile     <-  as.numeric(scan(par.fname,what="", n=16, quiet=TRUE)[c(6,11,16)])
  res@nopar   <-as.integer(parfile[1])
  res@nlogl   <-parfile[2]
  res@maxgrad <-parfile[3]
  
  #Read report file
  rept.fname    <- file.path(run.dir,sprintf("%s.rep",admb.stem))
  rept          <-  scan(rept.fname,what=integer(), quiet=TRUE)
  stateDim  <-  rept[1]
  yrs     <-  rept[-1]

  #Read residual files
  res.fname    <- file.path(run.dir,sprintf("%s.res",admb.stem))
  res@residuals  <-  read.table(res.fname,header=FALSE,
                       col.names=c("year","fleet","age","log.obs","log.mdl","std.res"))
  res@residuals$fleet <- factor(res@residuals$fleet)
  levels(res@residuals$fleet) <- fleet.names 

  #Read (parts of) the correlation matrix - first the determinant of the hessian
  cor.fname    <- file.path(run.dir,sprintf("%s.cor",admb.stem))
  lin             <-  readLines(cor.fname,n=1)
  res@logDetHess  <-  as.numeric(gsub("^.*=(.+)$","\\1",lin[1]))
  
  #Then the correlation matrix
  npar           <-  as.integer(length(lin)-2)
  cor.dat        <- read.table(cor.fname,skip=1,header=TRUE,fill=TRUE)
  cor.df         <- cor.dat[,c(1:4)]
  cor.mat        <- as.matrix(cor.dat[,-c(1:4)])

  #Am not sure what is going on here. Will leave this out until we get clarification
  #Problem is is that cor.mat is not square
#  ret@cor<-matrix(NA, ret@npar, ret@npar)
#  for(i in 1:ret@npar){
#    ret@cor[1:i,i]<-as.numeric(unlist(lapply(sublin[i],
#      function(x)x[5:(4+i)])))
#    ret@cor[i,1:i]<-ret@cor[1:i,i]
#  }
#  ret@cov<-ret@cor*(ret@std%o%ret@std)

  #Now extract parameter estimates together with uncertainties
  #These should really be put into FLQuants or similar
  mslh<-function(variable){
    sub.df <- subset(cor.df,name==variable,select=c("value","std"))
    sub.df$low.bnd <- sub.df$value - 1.96 * sub.df$std
    sub.df$up.bnd <- sub.df$value + 1.96 * sub.df$std
    if(nrow(sub.df)==length(yrs)) rownames(sub.df) <- yrs
    return(sub.df)
  }
  res@ssb<-mslh('ssb')
  res@fbar<-mslh('fbar')
  res@tsb<-mslh('tsb')
  res@logssb<-mslh('logssb')
  res@logfbar<-mslh('logfbar')
  res@logtsb<-mslh('logtsb')
  res@logscale<-mslh('logScale')
  res@logFpar<-mslh('logFpar')
  res@logCatch<-mslh('logCatch')

  #Extract the state variables
  x<-mslh('U')
  stateEst<-matrix(x$value,ncol=stateDim, byrow=TRUE,dimnames=list(year=yrs,state=NULL))
  stateStd<-matrix(x$std,ncol=stateDim, byrow=TRUE,dimnames=list(year=yrs,state=NULL))
  stateLow<-matrix(x$low.bnd,ncol=stateDim, byrow=TRUE,dimnames=list(year=yrs,state=NULL))
  stateHigh<-matrix(x$up.bnd,ncol=stateDim, byrow=TRUE,dimnames=list(year=yrs,state=NULL))

  #And copy into the appropriate slots as data.frames or FLQuants
  n.stateEst <- stateEst[,1:dims(stck)$age]
  f.stateEst <- exp(t(stateEst[,-c(1:dims(stck)$age)]))
  res@stock.n <- stck@stock.n
  res@stock.n@.Data[,,,,,] <- exp(t(stateEst[,1:dims(stck)$age]))
  res@harvest <- stck@harvest
  for(a in  dimnames(ctrl@states)$age){
    res@harvest@.Data[a,,,,,] <- f.stateEst[ctrl@states["catch",a],]
  }
  res@recruitment<-data.frame(value=exp(stateEst[,1]), std=NA, low.bnd=exp(stateLow[,1]),up.bnd= exp(stateHigh[,1]))
  rownames(res@recruitment) <- yrs

  #Finished! Remove the temporary directory
  if(missing(run.dir)) {unlink(run.dir)}
  return(res)
}

