  ### ======================================================================================================
  ### write.FLSAM.out - Used to generate a "pseudo-ica.out" file from FLSAM
  ### ======================================================================================================
  FLSAM.out   <-  function(FLStock,FLIndices,FLSAM,format="TABLE %i.") {
    #Check  validity of objects
    if(!validObject(FLStock) | !is.FLStock(FLStock)) {stop("\"FLStock\" argument is not a valid FLStock object")}
    if(!validObject(FLIndices)) {stop("\"FLIndices\" argument is not a valid FLIndices object")}
    if(!validObject(FLSAM)) {stop("\"FLSAM\" argument is not a valid FLSAM object")}
    #Also check that the FLSAM object corresponds to the FLStock and FLIndices object (somehow)
    
    #Initialise key values
    counter <- 1 
    opt     <- NULL
    
    #First go through the input data
    output.structure <- rbind(c("catch.n","CATCH IN NUMBER"),
                          c("catch.wt","WEIGHTS AT AGE IN THE CATCH"),
                          c("stock.wt","WEIGHTS AT AGE IN THE STOCK"),
                          c("m","NATURAL MORTALITY"),
                          c("mat","PROPORTION MATURE"),
                          c("harvest.spwn","FRACTION OF HARVEST BEFORE SPAWNING"),
                          c("m.spwn","FRACTION OF NATURAL MORTALITY BEFORE SPAWNING"))
    for(i in 1:nrow(output.structure)) {
      #First the title line, then the data
      opt <- c(opt, paste(sprintf(format,counter),output.structure[i,2]),"")
      opt <- c(opt, format.quant(slot(FLStock,output.structure[i,1])),"","")
      counter <- counter +1
    }  
  
    #Then the indices 
    opt <- c(opt, paste(sprintf(format,counter),"SURVEY INDICES"),"")
    for(idx in FLIndices) {
      #First the survey configuration, then the index values, then the variances
      opt <- c(opt, paste(idx@name,"- Configuration"),"")
      opt <- c(opt, paste(idx@desc,collapse=" "))
      opt <- c(opt, capture.output(print(idx@range)))
      opt <- c(opt, paste("Index type :",idx@type),"")
      opt <- c(opt, paste(idx@name,"- Index Values"),"")
      opt <- c(opt, format.quant(idx@index),"")
    }  
    opt <- c(opt,"")
    counter <- counter +1
    
    #Now information about how the stock object is configured
    opt <- c(opt, paste(sprintf(format,counter),"STOCK OBJECT CONFIGURATION"),"")
    opt <- c(opt, capture.output(FLStock@range),"","")
    counter <- counter + 1
  
    #Now the configuration of FLSAM
    opt <- c(opt, paste(sprintf(format,counter),"FLSAM CONFIGURATION SETTINGS"),"")
    longest.name  <-  max(sapply(slotNames(FLSAM@control),nchar))
    matrix.type   <-  c("range","fleets","states","catchabilities","power.law.exps","f.vars","obs.vars")
    for(sltName in slotNames(FLSAM@control)) {
      padded.name <-  paste(c(sltName,rep(" ",longest.name-nchar(sltName))),collapse="")
      if(sltName %in% matrix.type){
        opt <- c(opt,  paste(padded.name,":",paste(format.matrix(slot(FLSAM@control,sltName)))))
      } else {
          opt <- c(opt,paste(padded.name,":",paste(slot(FLSAM@control,sltName),collapse=" ")))
      }
    }
    opt <- c(opt,"","")
    counter <- counter + 1
    
    #FLR and R Package Information
    descrips <- list(packageDescription("FLSAM"),packageDescription("FLAssess"),packageDescription("FLCore"))
    descrips <- lapply(descrips, function(b) {
      rbind(paste("Package  :",b$Package),
            paste("Version  :",b$Version),
            paste("Packaged :",b$Packaged),
            paste("Built    :",b$Built),
            "")
    })
    descrip.str <- do.call(rbind,c(R.Version()$version.string,"",descrips))
    opt <- c(opt, paste(sprintf(format,counter),"FLR, R SOFTWARE VERSIONS"),"",
              capture.output(write.table(descrip.str,row.names=FALSE,quote=FALSE,col.names=FALSE)),"")
    counter <- counter + 1
  
    #Stock summary - the tricky one!
    opt <- c(opt, paste(sprintf(format,counter),"STOCK SUMMARY"),"")
    
    if(length(ssb(FLSAM)$year)!=length(as.data.frame(ssb(FLStock))$year)){
      landings=c(format(round(as.data.frame(FLStock@landings)$data,0)),"")
      sop=     c(format(round(as.data.frame(sop(FLStock,"landings"))$data,4)),"")
    } else {
      landings=format(round(as.data.frame(FLStock@landings)$data,0))
      sop=     format(round(as.data.frame(sop(FLStock,"landings"))$data,4))
    }

    ss.df <-  cbind(year =format(as.data.frame(ssb(FLSAM))$year),
                    recsm=format(round(as.data.frame(rec(FLSAM))$value,0)),
                    recsl=format(round(as.data.frame(rec(FLSAM))$lbnd,0)),
                    recsu=format(round(as.data.frame(rec(FLSAM))$ubnd,0)),
                    tsbm=format(round(as.data.frame(tsb(FLSAM))$value,0)),
                    tsbl=format(round(as.data.frame(tsb(FLSAM))$lbnd,0)),
                    tsbu=format(round(as.data.frame(tsb(FLSAM))$ubnd,0)),
                    ssbm=format(round(as.data.frame(ssb(FLSAM))$value,0)),
                    ssbl=format(round(as.data.frame(ssb(FLSAM))$lbnd,0)),
                    ssbu=format(round(as.data.frame(ssb(FLSAM))$ubnd,0)),
                    fbarm=format(round(as.data.frame(fbar(FLSAM))$value,4)),
                    fbarl=format(round(as.data.frame(fbar(FLSAM))$lbnd,4)),
                    fbaru=format(round(as.data.frame(fbar(FLSAM))$ubnd,4)),
                    landings=landings,
                    sop=sop)

    ss.units  <-  c("",rep(units(FLStock@stock.n),3),rep("",3),rep("",3),rep(units(fbar(FLStock)),3),units(FLStock@landings),"")
    ss.units  <-  ifelse(ss.units=="NA","",ss.units)
    ss.df <-  rbind(c("Year","Recruitment","Low","High","TSB","Low","High","SSB","Low","High","Fbar","Low","High","Landings","Landings"),
                    c("",paste("Age",dims(FLStock)$min),rep("",2),rep("",3),rep("",3),paste("(Ages ",FLStock@range["minfbar"],"-",FLStock@range["maxfbar"],")",sep=""),rep("",2),"","SOP"),
                    ss.units,ss.df)
    ss.df <-  apply(ss.df,2,"format",justify="right")
    opt  <- c(opt,capture.output(write.table(ss.df,row.names=FALSE,quote=FALSE,col.names=FALSE)),"","")
    counter <- counter + 1

    #Now, the estimated f and N (taken from the stock object, rather than the ica object
    #so that short term forcasts can be included 
    output.structure <- rbind(c("harvest","ESTIMATED FISHING MORTALITY"),
                              c("stock.n","ESTIMATED POPULATION ABUNDANCE"))
    for(i in 1:nrow(output.structure)) {
      #First the title line, then the data
      opt <- c(opt, paste(sprintf(format,counter),output.structure[i,2]),"")
      opt <- c(opt, format.quant(slot(FLSAM,output.structure[i,1])),"","")
      counter <- counter +1
    }  
    
  #And now, the rest of the outputs from the assessment. Here just needs pred catch and catch resid

  #Subtracting pred catch and catch resid from NSH.sam  
  x           <- subset(FLSAM@residuals,fleet=="catch")
  year        <- x$year
  age         <- x$age
  pred.catch  <- FLQuant(exp(x$log.mdl),dimnames=list(age=sort(unique(age)),year=sort(unique(year)),unit="unique",season="all",area="unique",iter="1"))
  catch.resid <- FLQuant(x$std.res,dimnames=list(age=sort(unique(age)),year=sort(unique(year)),unit="unique",season="all",area="unique",iter="1"))
  output.structure <- rbind(c("pred.catch","PREDICTED CATCH NUMBERS AT AGE"),
                              c("catch.resid","CATCH AT AGE RESIDUALS"))
    for(i in 1:nrow(output.structure)) {
      #First the title line, then the data
      opt <- c(opt, paste(sprintf(format,counter),output.structure[i,2]),"")
      opt <- c(opt, format.quant(get(output.structure[i,1])),"","")
      counter <- counter +1
    }

    opt <- c(opt,"")
    counter <- counter +1

  
  #Predicted index values
  flts  <- ac(unique(FLSAM@residuals$fleet))
  for(j in flts[which(!flts %in% "catch")]){
    x           <- subset(FLSAM@residuals,fleet==j)
    year        <- x$year
    if(j %in% names(FLSAM@control@fleets[which(FLSAM@control@fleets == 3)])){
      age         <- "all"
    } else { age <- x$age }
    pred.catch  <- FLQuant(exp(x$log.mdl),dimnames=list(age=sort(unique(age)),year=sort(unique(year)),unit="unique",season="all",area="unique",iter="1"))
    catch.resid <- FLQuant(x$std.res,dimnames=list(age=sort(unique(age)),year=sort(unique(year)),unit="unique",season="all",area="unique",iter="1"))
    output.structure <- rbind(c("pred.catch","PREDICTED INDEX AT AGE"),
                              c("catch.resid","INDEX AT AGE RESIDUALS"))

    for(i in 1:nrow(output.structure)) {
      #First the title line, then the data
      opt <- c(opt, paste(sprintf(format,counter),output.structure[i,2],j),"")
      opt <- c(opt, format.quant(get(output.structure[i,1])),"","")
      counter <- counter +1
    }
  }
    opt <- c(opt,"")
    counter <- counter +1
  

    #Parameter estimates
    opt <- c(opt, paste(sprintf(format,counter),"FIT PARAMETERS"),"")
    opt <-  c(opt,capture.output(coef(FLSAM)),"")
    opt <- c(opt,"")
    counter <- counter + 1
    
    #Negative log-likelihood
    opt <- c(opt, paste(sprintf(format,counter),"NEGATIVE LOG-LIKELIHOOD"),"")
    opt <- c(opt, FLSAM@nlogl,"","")
    counter <- counter + 1
    
    #And finish  
    return(opt)
  
  }
  
  sam.out   <-  function(...) {
    FLSAM.out(...)
  }
  
  
  format.quant <- function(x) {
    #Drop extraneous dimensions
    mat.obj <-  drop(x@.Data)   
    #Check that we haven't "overdropped", and are left with at least quant vs year
    if(is.null(dim(mat.obj))) {
      mat.obj <-  array(as.vector(mat.obj),dim=dim(x@.Data)[1:2],dimnames=dimnames(x@.Data)[1:2])
    }
    
    #Write output
    opt <-  capture.output({
      cat(paste("Units  : ",units(x),"\n"))
      print.default(mat.obj,quote=FALSE,right=TRUE)
      })  #End  capture output
    
    return(opt)
      
  }   #End function, end setMethod

  format.matrix <- function(x){
    #Write output
    opt <-  capture.output({
      print.default(x,quote=FALSE,right=TRUE)
      })  #End  capture output
    return(opt)
  }