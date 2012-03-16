SAM2FLR <-function(ctrl,run.dir=tempdir()) {
  #---------------------------------------------------
  # Read the results from the assessment
  #---------------------------------------------------
  admb.stem <- "sam"

  #Create return object
  res       <- new("FLSAM")
  res@name  <- ctrl@name
  res@desc  <- ctrl@desc
  res@range <- ctrl@range
  res@control <- ctrl
  
  #Read parameter file
  par.fname    <- file.path(run.dir,sprintf("%s.par",admb.stem))
  par.hdr     <-  as.numeric(scan(par.fname,what="", n=16, quiet=TRUE)[c(6,11,16)])
  res@nopar   <-  as.integer(par.hdr[1])
  res@nlogl   <-  par.hdr[2]
  
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


  #---------------------------
  #Handle missing std deviation if nohess=TRUE
  #---------------------------
  if(ctrl@nohess) {
     #Load and parse par file
     par.txt     <- readLines(par.fname)[-1]
     n.lines     <- length(par.txt)
     par.chunks  <- paste(par.txt[seq(1,n.lines,by=2)],par.txt[seq(2,n.lines,by=2)])
     par.lst     <- lapply(par.chunks,function(x) {
                       var.name <- gsub("# (.*?):.*","\\1",x) 
                       var.vals <- gsub("# .*?:(.*)","\\1",x) 
                       var.dat  <- scan(con <- textConnection(var.vals),quiet=TRUE) 
                       close(con)
                       data.frame(name=var.name,value=var.dat,std.dev=NA)
                    })
     par.df <- do.call(rbind,par.lst)
     res@params <- par.df

     #Create synthetic (empty) vcov matrix
     res@vcov <- matrix(NA)

  } else { #Else read the results

     #Read standard deviation report (if it exists)
     std.fname <- file.path(run.dir,sprintf("%s.std",admb.stem))
     res@params <- read.table(std.fname,header=FALSE,skip=1,
                     col.names=c("index","name","value","std.dev"))
     res@params$index <- NULL
   
     #Now read the variance-covariance matrix. However, the problem is
     #is that the result return by ADMB is not a square file, but only
     #the lower triangle
     cor.fname    <- file.path(run.dir,sprintf("%s.cor",admb.stem))
     cor.lines.in <- readLines(cor.fname)[-c(1:2)] #First lines are headers
     cor.lines    <- strsplit(cor.lines.in," +") #Split into elements 
     vcov.dim     <- length(cor.lines) #Size of (squqare) vcov matrix 
     #Write individual lines into a matrix for further processing. I can't
     #see a cleaner way to do this apart from a for-loop :-(
     cor.file.mat      <- matrix(as.character(NA), vcov.dim, vcov.dim+5)
     for(i in 1:vcov.dim){ 
       cor.file.mat[i,seq(cor.lines[[i]])]<- cor.lines[[i]]
     }
     #Extract correlation matrix
     cor.mat.lt  <- apply(cor.file.mat[,-c(1:5)],2,as.numeric)
     cor.mat.lt[upper.tri(cor.mat.lt,diag=FALSE)] <- 0
     cor.mat <- cor.mat.lt + t(cor.mat.lt)
     diag(cor.mat) <- diag(cor.mat)/2   #Diagonals get included twice, so halve them
   
     #Convert correlation matrix to vcov matrix
     stds <- as.numeric(cor.file.mat[,5]) 
     vcov.mat <- cor.mat*(stds%o%stds)
     dimnames(vcov.mat) <- list(cor.file.mat[,3],cor.file.mat[,3])
     res@vcov <- vcov.mat
   
  }

  #Extract the state variables
  u<-subset(res@params,name=="U")
  stateEst<-matrix(u$value,nrow=n.states, byrow=FALSE,dimnames=list(state=NULL,year=yrs))
  stateSd <-matrix(u$std.dev,nrow=n.states, byrow=FALSE,dimnames=list(state=NULL,year=yrs))

  #And copy into the appropriate slots as data.frames or FLQuants
  flq <- FLQuant(NA, dimnames=list(age=ctrl@range["min"]:ctrl@range["max"], 
                                   year=ctrl@range["minyear"]:ctrl@range["maxyear"]))
  n.ages <- dims(flq)$age
  res@stock.n <- flq
  res@stock.n[,] <- exp(stateEst[1:n.ages,])
  res@harvest <- flq
  res@harvest@units <- "f"
  f.stateEst <- stateEst[-c(1:n.ages),,drop=FALSE]
  for(a in  dimnames(ctrl@states)$age){
    states.key <- ctrl@states["catch",a] 
    res@harvest[a,] <- exp(f.stateEst[states.key,])
  }

  #Populate the info slot as last step
  info <- data.frame(FLSAM.version=packageDescription("FLSAM")$Version,
                     FLCore.version=packageDescription("FLCore")$Version,
                     R.version=R.version$version.string,
                     platform=R.version$platform,
                     run.date=Sys.time())
  res@info <- t(info)
  colnames(res@info) <- "_"

  #Finished! 
  return(res)
}

