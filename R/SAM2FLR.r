SAM2FLR <-function(ctrl,run.dir="missing",admb.stem="ssass") {
  #---------------------------------------------------
  # Read the results from the assessment
  #---------------------------------------------------
  #Create return object
  res       <- new("FLSAM")
  res@name  <- ctrl@name
  res@desc  <- ctrl@desc
  res@range <- ctrl@range
  res@control <- ctrl
  
  #Read parameter file
  par.fname    <- file.path(run.dir,sprintf("%s.par",admb.stem))
  parfile     <-  as.numeric(scan(par.fname,what="", n=16, quiet=TRUE)[c(6,11,16)])
  res@nopar   <-as.integer(parfile[1])
  res@nlogl   <-parfile[2]
  
  #Read report file
  rept.fname    <- file.path(run.dir,sprintf("%s.rep",admb.stem))
  rept          <-  scan(rept.fname,what=integer(), quiet=TRUE)
  n.states  <-  rept[1]
  yrs       <-  rept[-1]

  #Read residual files
  res.fname    <- file.path(run.dir,sprintf("%s.res",admb.stem))
  res@residuals  <-  read.table(res.fname,header=FALSE,
                       col.names=c("year","fleet","age","log.obs","log.mdl","std.res"))
  res@residuals$fleet <- factor(res@residuals$fleet)
  levels(res@residuals$fleet) <- names(ctrl@fleets)

  #Read standard deviation report
  std.fname <- file.path(run.dir,sprintf("%s.std",admb.stem))
  res@params <- read.table(std.fname,header=FALSE,skip=1,
                  col.names=c("index","name","value","std.dev"))

  #TODO: Following code reads the correlation matrix. I have left this
  #      out for the meantime, but it can be added back in later if required
  #Read (parts of) the correlation matrix - first the determinant of the hessian
  #cor.fname    <- file.path(run.dir,sprintf("%s.cor",admb.stem))
  #lin             <-  readLines(cor.fname,n=1)
  #res@logDetHess  <-  as.numeric(gsub("^.*=(.+)$","\\1",lin[1]))
  
  #Then the correlation matrix
  #npar           <-  as.integer(length(lin)-2)
  #cor.dat        <- read.table(cor.fname,skip=1,header=TRUE,fill=TRUE)
  #cor.df         <- cor.dat[,c(1:4)]
  #cor.mat        <- as.matrix(cor.dat[,-c(1:4)])

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
#  mslh<-function(variable){
#    sub.df <- subset(cor.df,name==variable,select=c("value","std"))
#    sub.df$low.bnd <- sub.df$value - 1.96 * sub.df$std
#    sub.df$up.bnd <- sub.df$value + 1.96 * sub.df$std
#    if(nrow(sub.df)==length(yrs)) rownames(sub.df) <- yrs
#    return(sub.df)
#  }

  #Extract the state variables
  u<-subset(res@params,name=="U")
  stateEst<-matrix(u$value,nrow=n.states, byrow=FALSE,dimnames=list(state=NULL,year=yrs))
  stateSd <-matrix(u$std.dev,nrow=n.states, byrow=FALSE,dimnames=list(state=NULL,year=yrs))

  #And copy into the appropriate slots as data.frames or FLQuants
  flq <- FLQuant(NA, dimnames=list(age=ctrl@range["min"]:ctrl@range["max"], 
                                   year=ctrl@range["minyear"]:ctrl@range["maxyear"]))
  n.ages <- dims(flq)$age
  res@stock.n <- flq
  res@stock.n <- exp(stateEst[1:n.ages,])
  res@harvest <- flq
  res@harvest@units <- "f"
  f.stateEst <- stateEst[-c(1:n.ages),]
  for(a in  dimnames(ctrl@states)$age){
    states.key <- ctrl@states["catch",a] 
    res@harvest[a,] <- exp(f.stateEst[states.key,]
  }

  #Finished! 
  return(res)
}

