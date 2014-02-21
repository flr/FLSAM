SAM2FLR <-function(ctrl="missing",admb.stem="sam",run.dir=tempdir()) {
  #It may not always be possible to create a ctrl object e.g.
  #e.g. when trying to read the output of a model from
  #stockassessment.org. We start by adding some code to 
  #handle this situation and recreate an object from
  #the contents of the run.dir
  #Create object
  res       <- new("FLSAM")
  
  #Read report file
  rept.fname    <- file.path(run.dir,sprintf("%s.rep",admb.stem))
  rept          <-  scan(rept.fname,what=integer(), quiet=TRUE)
  res@n.states  <-  rept[1]
  yrs           <-  rept[-1]
  
  if(missing(ctrl)) {
    res@name <- sprintf("Read from %s",run.dir)
    res@desc <- date()
    #Figure out whether the hessian has been estimated
    res@nohess <- !file.exists(sprintf("%s/%s.cor",run.dir,admb.stem))
    #Work out the contents of the range object from the
    #contents of the directory. First the ages
    mdl.cfg <- readLines(sprintf("%s/model.cfg",run.dir))
    res@range["min"] <- scan(text=mdl.cfg[grep("Min Age",mdl.cfg,
                                               ignore.case=TRUE)+1],quiet=TRUE)
    max.ages <- scan(text=mdl.cfg[grep("Max Age",mdl.cfg,
                                       ignore.case=TRUE)+1],quiet=TRUE)
    res@range["max"] <- max.ages[1]
    res@range["plusgroup"] <- ifelse(max.ages[2]==1,max.ages[1],NA)
    #Extract the bindings
    comment.line <- grepl("^ *#",mdl.cfg)
    num.block <- ifelse(!comment.line,cumsum(comment.line),NA)
    state.line    <- grep("STATES",mdl.cfg,ignore.case=FALSE)
    state.block <- min(num.block[-(1:state.line)],na.rm=TRUE)
    res@states <- as.matrix(read.table(text=mdl.cfg[which(num.block==state.block)]))
    colnames(res@states) <- res@range["min"]:res@range["max"]
    rownames(res@states) <- c("catch",seq(nrow(res@states))[-1])
    #And now the years (already read in from sam.rep)
    res@range["minyear"] <- min(yrs)
    res@range["maxyear"] <- max(yrs)
    #And now the fbar
    fbar.vals <- scan(text=mdl.cfg[grep("Fbar",mdl.cfg)+1],quiet=TRUE)
    res@range["minfbar"] <- min(fbar.vals)
    res@range["maxfbar"] <- max(fbar.vals)
    #Work out the fishing mortality couplings
    
  } else {
    #Extract stem
    admb.stem <- .get.admb.stem(ctrl)
    #Create return object
    res@name  <- ctrl@name
    res@desc  <- ctrl@desc
    res@range <- ctrl@range
    res@nohess  <- ctrl@nohess
    res@control <- ctrl
    res@states <- ctrl@states
  }
  
  #---------------------------------------------------
  # Read the results from the assessment
  #---------------------------------------------------
  #Read parameter file
  par.fname    <- file.path(run.dir,sprintf("%s.par",admb.stem))
  par.hdr     <-  as.numeric(scan(par.fname,what="", n=16, quiet=TRUE)[c(6,11,16)])
  res@nopar   <-  as.integer(par.hdr[1])
  res@nlogl   <-  par.hdr[2]
  res@maxgrad <-  par.hdr[3]

  #Read residual files
  res.fname    <- file.path(run.dir,"sam.res") #Not ideal, but easier than changing it in tpl
  res@residuals  <-  read.table(res.fname,header=FALSE,
                       col.names=c("year","fleet","age","log.obs","log.mdl","std.res"))
  res@residuals$fleet <- factor(sprintf("Fleet %s",res@residuals$fleet))

  #---------------------------
  #The model can occasionally be run in "nohess" mode, where 
  #the hessian is not estimated. This can create problems as 
  #neither the vcov matrix nor the std.errs are estimated
  #---------------------------
  if(res@nohess) {
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
  stateEst<-matrix(u$value,nrow=res@n.states, byrow=FALSE,
                   dimnames=list(state=NULL,
                                 year=res@range["minyear"]:res@range["maxyear"]))
  stateSd <-matrix(u$std.dev,nrow=res@n.states, byrow=FALSE,
                   dimnames=list(state=NULL,
                                 year=res@range["minyear"]:res@range["maxyear"]))

  #And copy into the appropriate slots as data.frames or FLQuants
  flq <- FLQuant(NA, dimnames=list(age=res@range["min"]:res@range["max"], 
                                   year=res@range["minyear"]:res@range["maxyear"]))
  n.ages <- dims(flq)$age
  res@stock.n <- flq
  res@stock.n[,] <- exp(stateEst[1:n.ages,])
  
  #Populate the harvest slot using the state bindings
  n.ages <- dims(flq)$age
  colnames(res@states) <- dimnames(flq)$age
  res@harvest <- flq
  res@harvest@units <- "f"
  f.stateEst <- stateEst[-c(1:n.ages),,drop=FALSE]
  for(a in  colnames(res@states)){
    states.key <- res@states["catch",a] 
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

