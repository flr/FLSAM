# merge {{{
# merge results from FLSAM into FLStock
if (!isGeneric("merge")) {
    setGeneric("merge", useAsDefault = merge)
}
setMethod("merge", signature(x="FLStock", y="FLSAM"),
  function(x, y, ...)
  {
    quant <- quant(stock.n(x))
    dnx <- dimnames(stock.n(x))
    dny <- dimnames(y@stock.n)

    # check dimensions match
    if(!all.equal(dnx[c(-2,-6)], dny[c(-2,-6)]))
      stop("Mismatch in dimensions: only year can differ between stock and assess")

    # same plusgroup
    if(x@range['plusgroup'] != x@range['plusgroup'])
      stop("Mismatch in plusgroup: x and y differ")

    # year ranges match?
    if(!all(dny[['year']] %in% dnx[['year']]))
    {
      x@stock.n <- window(y@stock.n,start=x@range["minyear"],end=x@range["maxyear"])
      x@harvest <- window(y@harvest,start=x@range["minyear"],end=x@range["maxyear"])
    } else {
       x@stock.n <- y@stock.n
       x@harvest <- y@harvest
    } 
    x@desc <- paste(x@desc, "+ FLSAM:", y@name)
    x@harvest@units <- y@harvest@units 
    x@range=c(unlist(dims(x)[c('min', 'max', 'plusgroup','minyear', 'maxyear')]),
      x@range[c('minfbar', 'maxfbar')])
        
    return(x)
  }
)   # }}}

# "+"      {{{
setMethod("+", signature(e1="FLStock", e2="FLSAM"),
  function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(merge(e1, e2))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)
setMethod("+", signature(e1="FLSAM", e2="FLStock"),
	function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(merge(e2, e1))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)   # }}}

setMethod("+", signature(e1="FLStock", e2="FLSAMs"),
  function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(FLStocks(lapply(e2,merge,x=e1)))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)
setMethod("+", signature(e1="FLSAMs", e2="FLStock"),
	function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(FLStocks(lapply(e1,merge,x=e2)))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)   # }}}


#General helper function to extract a given parameter from an FLSAM object
#and return it as a data.frame
.params2flq <- function(object,param) {
   dat <- subset(object@params,name==param)
   flq <- FLQuant(dat$value,dimnames=list(age="all",
              year=object@range["minyear"]:object@range["maxyear"]))
   return(flq)
}

.extract.params <- function(object,params) {
  #Extract numbers at age parameters
  ps <- subset(object@params,name%in%params)[,c("name","value","std.dev")]
  colnames(ps)[which(colnames(ps)=="name")] <- "param.name"
  ps$ubnd <- exp(ps$value + 1.96*ps$std.dev)
  ps$lbnd <- exp(ps$value - 1.96*ps$std.dev)
  ps$value <- exp(ps$value)
  ps$std.dev <- NULL
  return(ps)
}
.extract.states <- function(object) {
          Us <- .extract.params(object,"U")
          yrs <- seq(object@control@range["minyear"],
                     object@control@range["maxyear"]) 
          ages <- seq(object@control@range["min"],
                      object@control@range["max"])
          n.states <- nrow(Us)/length(yrs)
          states <- c(ages,seq(-1,by=-1,length.out=n.states-length(ages)))
          Us <- cbind(expand.grid(age=states,year=yrs),Us)
          return(Us)
}

#ssb          {{{
setMethod("ssb", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logssb")   
          res <- cbind(year=seq(object@control@range["minyear"],
                            object@control@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

setMethod("ssb", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object[[i]]@name,ssb(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

setMethod("fbar", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logfbar")   
          res <- cbind(year=seq(object@control@range["minyear"],
                            object@control@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

setMethod("fbar", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object[[i]]@name,fbar(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

setMethod("tsb", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logtsb")   
          res <- cbind(year=seq(object@control@range["minyear"],
                            object@control@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

setMethod("tsb", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object[[i]]@name,tsb(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

setMethod("n", signature(object="FLSAM"),
        function(object, ...) {
	  return(subset(.extract.states(object),age>=0))
        }
)       # }}}

setMethod("n", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object[[i]]@name,n(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

if (!isGeneric("f")) {
    setGeneric("f", function(object) standardGeneric("f"))
}
setMethod("f", signature(object="FLSAM"),
        function(object) {
          f.states <- subset(.extract.states(object),age<0)
          f.states$param <- -f.states$age 
          f.bindings <- object@control@states["catch",]
          f.bindings <- as.data.frame(as.table(f.bindings),responseName="param")
          res <- merge(f.states,f.bindings,by="param") 
          res$age <- res$Var1
          res$Var1 <- NULL
          res$param <- NULL
          res <- res[order(res$year,res$age),]
          return(res)
        }
)       # }}}

setMethod("f", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object[[i]]@name,f(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

setMethod("rec", signature(object="FLSAM"),
        function(object, ...) {
          return(subset(n(object),age==object@control@range["min"]))
        }
)       # }}}

setMethod("rec", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object[[i]]@name,rec(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

setMethod("AIC",signature(object="FLSAM"),
        function(object, ...) {
        return(2*object@nopar +2*object@nlogl)
        }
)

setMethod("AIC",signature(object="FLSAMs"),
        function(object, ...) {
        return(sapply(object,AIC,...))
        }
)

#- Create generic function for 'looi'
if (!isGeneric("looi")) {
  setGeneric('looi', function(e1,e2,e3,type="full") standardGeneric('looi'))
}

setMethod("looi",signature(e1="FLStock",e2="FLIndices",e3="FLSAM.control",type="character"),
  function(e1,e2,e3,type="full"){
  
    #- Create run scheme with all possible combinations
    overview          <- matrix(NA,nrow=2^length(e2),ncol=length(e2),dimnames=list(run=1:2^length(e2),fleet=names(e2)))
    for(iFlt in 1:length(e2))
      overview[,iFlt] <- rep(c(rep(1,(2^length(e2))/(2^(iFlt))),rep(0,(2^length(e2))/(2^(iFlt)))),2^(iFlt-1))
    overview          <- overview[-nrow(overview),] #Remove 'no-fleet-included' option
    rownames(overview)<- unlist(lapply(apply(overview,1,function(x){names(which(x==1))}),paste,collapse="+"))
    
    if(type=="LOO"){ overview <- overview[c(1,which(rowSums(overview) == (length(e2)-1))),];  print("Running LOO analyses")}
    if(type=="LOI"){ overview <- overview[c(1,which(rowSums(overview) == 1)),];               print("Running LOI analyses")}
    
    #- Create object to save outcomes
    result <- new("FLSAMs")
    for(iRun in rownames(overview)){
      stck                  <- e1
      tun                   <- e2[which(overview[iRun,]==1)]
      ctrl                  <- FLSAM.control(stck,tun)
      #- Update the ctrl and stck object given the changed tuning object
      for(iSlot in slotNames(e3))
        if(class(slot(ctrl,iSlot))=="matrix") slot(ctrl,iSlot) <- slot(e3,iSlot)[c(names(which(e3@fleets==0)),names(which(overview[iRun,]==1))),]
      ctrl@logN.vars        <- e3@logN.vars
      ctrl@srr              <- e3@srr

      #- Run the assessment
      result[[iRun]]        <- FLSAM(stck,tun,ctrl)
    }
  return(result)}
)


if (!isGeneric("catchabilities")) {
  setGeneric('catchabilities', function(object) standardGeneric('catchabilities'))
}
setMethod("catchabilities",signature(object="FLSAM"),
  function(object) {
       #Extract numbers at age parameters
       params <- .extract.params(object,"logFpar")
       if(nrow(params)==0) stop("No catchabiltities fitted in model.")
       params$no <- 1:nrow(params)
       bindings <-  as.data.frame(as.table(object@control@catchabilities),responseName="no")
       #Merge
       bindings <- subset(bindings,!is.na(bindings$no))
       res <- merge(bindings,params)
       #Tidy up
       res <- res[order(res$fleet,res$age),]
       res$no <- NULL
       return(res)
})

obs.var <- function(object) {
       #Extract data
       params <- .extract.params(object,"logSdLogObs")
       if(nrow(params)==0) stop("No observation variance fitted in model.")
       params$no <- 1:nrow(params)
       bindings <-  as.data.frame(as.table(object@control@obs.vars),responseName="no")
       #Merge
       bindings <- subset(bindings,!is.na(bindings$no))
       res <- merge(bindings,params)
       #Tidy up
       res <- res[order(res$fleet,res$age),]
       res$no <- NULL
       return(res)
}

power.law.exps <- function(object) {
       #Extract data
       params <- .extract.params(object,"logQpow")
       if(nrow(params)==0) { stop("No power law exponents fitted in model.")}
       params$no <- 1:nrow(params)
       bindings <-  as.data.frame(as.table(object@control@power.law.exps),responseName="no")
       #Merge
       bindings <- subset(bindings,!is.na(bindings$no))
       res <- merge(bindings,params)
       #Tidy up
       res <- res[order(res$fleet,res$age),]
       res$no <- NULL
       return(res)
}


#- Create generic function for 'lr.test'
if (!isGeneric("lr.test")) {
  setGeneric('lr.test', function(object,...,type="sequential") standardGeneric('lr.test'))
}

setMethod("lr.test",signature("FLSAMs"),
  function(object,...,type="sequential"){
    sams    <- object
    
    #- Get negative log likelihood and number of parameter values
    dat.l <- lapply(sams,function(mdl) {data.frame(nll=mdl@nlogl,npar=mdl@nopar)})

    #- Check if models are ordered by increasing or decreasing complexity
    if(!all(abs(diff(do.call(rbind,lapply(sams,function(mdl) {data.frame(nll=mdl@nlogl,npar=mdl@nopar)}))$npar))>0)) stop("Models not ordered by increasing or decreasing complexity")

    #- Give each model a name
    modNames <- unlist(lapply(sams,name))
    if(any(nchar(modNames)==0) | any(length(nchar(modNames))==0))    modNames[which(nchar(modNames)==0 | length(nchar(modNames))==0)] <- paste("FLSAM",which(nchar(modNames)==0 | length(nchar(modNames))==0))
    if(any(duplicated(modNames)==T)) modNames[which(duplicated(modNames)==T)] <- paste(modNames[which(duplicated(modNames)==T)],which(duplicated(modNames)==T))
    modNames <- paste(1:length(modNames),modNames)

    tbl <- matrix(NA,nrow=length(sams),ncol=6,dimnames=list(models=modNames,statistics=c("Comparison","Neg. log likel","# Parameters","Likel difference","Degrees of freedom","P value")))
    #- Perform Log-ratio test in sequential mode
    if(type == "sequential"){
      tbl[2:length(modNames),"Comparison"]      <- paste(2:length(modNames),"vs.",1:(length(modNames)-1))
      tbl[,c("Neg. log likel","# Parameters")]  <- as.matrix(do.call(rbind,dat.l))
      for(i in 2:length(modNames)){
        if(as.numeric(tbl[i,"# Parameters"]) >  as.numeric(tbl[i-1,"# Parameters"])){
          tbl[i,  "Likel difference"]           <- round(as.numeric(tbl[i-1,"Neg. log likel"])            - as.numeric(tbl[i,  "Neg. log likel"]),2)
          tbl[i,  "Degrees of freedom"]         <- as.numeric(      tbl[i,"# Parameters"])                - as.numeric(tbl[i-1,"# Parameters"])
          tbl[i,  "P value"]                    <- round(1-pchisq(2*(as.numeric(tbl[i-1,"Neg. log likel"])- as.numeric(tbl[i,  "Neg. log likel"])),as.numeric(tbl[i,"Degrees of freedom"])),4)
        }
        if(as.numeric(tbl[i,"# Parameters"]) <= as.numeric(tbl[i-1,"# Parameters"])){
          tbl[i-1,"Comparison"]                 <- paste(i-1,"vs.",i); tbl[i,"Comparison"] <- NA
          tbl[i-1,"Likel difference"]           <- round(as.numeric(tbl[i,"Neg. log likel"])              - as.numeric(tbl[i-1,"Neg. log likel"]),2)
          tbl[i-1,"Degrees of freedom"]         <- as.numeric(      tbl[i-1,"# Parameters"])              - as.numeric(tbl[i,  "# Parameters"])
          tbl[i-1,"P value"]                    <- round(1-pchisq(2*(as.numeric(tbl[i,"Neg. log likel"])  - as.numeric(tbl[i-1,"Neg. log likel"])),as.numeric(tbl[i-1,"Degrees of freedom"])),4)
        }
      }
    }

    #- Perform Log-ratio test between each models and the baseline (first) model
    if(type == "first"){
      tbl[2:length(modNames),"Comparison"]      <- paste(2:length(modNames),"vs.",1)
      tbl[,c("Neg. log likel","# Parameters")]  <- as.matrix(do.call(rbind,dat.l))
      for(i in 2:length(modNames)){
        if(as.numeric(tbl[i,"# Parameters"]) >  as.numeric(tbl[1,"# Parameters"])){
          tbl[i,  "Likel difference"]           <- round(as.numeric(tbl[1,"Neg. log likel"])              - as.numeric(tbl[i,"Neg. log likel"]),2)
          tbl[i,  "Degrees of freedom"]         <- as.numeric(      tbl[i,"# Parameters"])                - as.numeric(tbl[1,"# Parameters"])
          tbl[i,  "P value"]                    <- round(1-pchisq(2*(as.numeric(tbl[1,"Neg. log likel"])  - as.numeric(tbl[i,"Neg. log likel"])),as.numeric(tbl[i,"Degrees of freedom"])),4)
        }
        if(as.numeric(tbl[i,"# Parameters"]) <= as.numeric(tbl[1,"# Parameters"])){
          tbl[i,"Likel difference"]             <- round(as.numeric(tbl[i,"Neg. log likel"])              - as.numeric(tbl[1,"Neg. log likel"]),2)
          tbl[i,"Degrees of freedom"]           <- as.numeric(      tbl[1,"# Parameters"])                - as.numeric(tbl[i,"# Parameters"])
          tbl[i,"P value"]                      <- round(1-pchisq(2*(as.numeric(tbl[i,"Neg. log likel"])  - as.numeric(tbl[1,"Neg. log likel"])),as.numeric(tbl[i,"Degrees of freedom"])),4)
        }
      }
    }
    return(as.table(tbl))
  }
)
    
    
    
setMethod("lr.test",signature("FLSAM"),
  function(object,...,type="sequential"){
    mdls    <- list(object,...)
    sams    <- new("FLSAMs")
    for(i in 1:length(mdls)) sams[[i]] <- mdls[[i]]

    #- Get negative log likelihood and number of parameter values
    dat.l <- lapply(sams,function(mdl) {data.frame(nll=mdl@nlogl,npar=mdl@nopar)})

    #- Check if models are ordered by increasing or decreasing complexity
    if(!all(abs(diff(do.call(rbind,lapply(sams,function(mdl) {data.frame(nll=mdl@nlogl,npar=mdl@nopar)}))$npar))>0)) stop("Models not ordered by increasing or decreasing complexity")

    #- Give each model a name
    modNames <- unlist(lapply(sams,name))
    if(any(nchar(modNames)==0) | any(length(nchar(modNames))==0))    modNames[which(nchar(modNames)==0 | length(nchar(modNames))==0)] <- paste("FLSAM",which(nchar(modNames)==0 | length(nchar(modNames))==0))
    if(any(duplicated(modNames)==T)) modNames[which(duplicated(modNames)==T)] <- paste(modNames[which(duplicated(modNames)==T)],which(duplicated(modNames)==T))
    modNames <- paste(1:length(modNames),modNames)

    tbl <- matrix(NA,nrow=length(sams),ncol=6,dimnames=list(models=modNames,statistics=c("Comparison","Neg. log likel","# Parameters","Likel difference","Degrees of freedom","P value")))
    #- Perform Log-ratio test in sequential mode
    if(type == "sequential"){
      tbl[2:length(modNames),"Comparison"]      <- paste(2:length(modNames),"vs.",1:(length(modNames)-1))
      tbl[,c("Neg. log likel","# Parameters")]  <- as.matrix(do.call(rbind,dat.l))
      for(i in 2:length(modNames)){
        if(as.numeric(tbl[i,"# Parameters"]) >  as.numeric(tbl[i-1,"# Parameters"])){
          tbl[i,  "Likel difference"]           <- round(as.numeric(tbl[i-1,"Neg. log likel"])            - as.numeric(tbl[i,  "Neg. log likel"]),2)
          tbl[i,  "Degrees of freedom"]         <- as.numeric(      tbl[i,"# Parameters"])                - as.numeric(tbl[i-1,"# Parameters"])
          tbl[i,  "P value"]                    <- round(1-pchisq(2*(as.numeric(tbl[i-1,"Neg. log likel"])- as.numeric(tbl[i,  "Neg. log likel"])),as.numeric(tbl[i,"Degrees of freedom"])),4)
        }
        if(as.numeric(tbl[i,"# Parameters"]) <= as.numeric(tbl[i-1,"# Parameters"])){
          tbl[i-1,"Comparison"]                 <- paste(i-1,"vs.",i); tbl[i,"Comparison"] <- NA
          tbl[i-1,"Likel difference"]           <- round(as.numeric(tbl[i,"Neg. log likel"])              - as.numeric(tbl[i-1,"Neg. log likel"]),2)
          tbl[i-1,"Degrees of freedom"]         <- as.numeric(      tbl[i-1,"# Parameters"])              - as.numeric(tbl[i,  "# Parameters"])
          tbl[i-1,"P value"]                    <- round(1-pchisq(2*(as.numeric(tbl[i,"Neg. log likel"])  - as.numeric(tbl[i-1,"Neg. log likel"])),as.numeric(tbl[i-1,"Degrees of freedom"])),4)
        }
      }
    }

    #- Perform Log-ratio test between each models and the baseline (first) model
    if(type == "first"){
      tbl[2:length(modNames),"Comparison"]      <- paste(2:length(modNames),"vs.",1)
      tbl[,c("Neg. log likel","# Parameters")]  <- as.matrix(do.call(rbind,dat.l))
      for(i in 2:length(modNames)){
        if(as.numeric(tbl[i,"# Parameters"]) >  as.numeric(tbl[1,"# Parameters"])){
          tbl[i,  "Likel difference"]           <- round(as.numeric(tbl[1,"Neg. log likel"])              - as.numeric(tbl[i,  "Neg. log likel"]),2)
          tbl[i,  "Degrees of freedom"]         <- as.numeric(tbl[i,"# Parameters"])                      - as.numeric(tbl[1,  "# Parameters"])
          tbl[i,  "P value"]                    <- round(1-pchisq(2*(as.numeric(tbl[1,"Neg. log likel"])  - as.numeric(tbl[i,"Neg. log likel"])),as.numeric(tbl[i,"Degrees of freedom"])),4)
        }
        if(as.numeric(tbl[i,"# Parameters"]) <= as.numeric(tbl[1,"# Parameters"])){
          tbl[i,"Likel difference"]             <- round(as.numeric(tbl[i,"Neg. log likel"])              - as.numeric(tbl[1,"Neg. log likel"]),2)
          tbl[i,"Degrees of freedom"]           <- as.numeric(tbl[1,"# Parameters"])                      - as.numeric(tbl[i,  "# Parameters"])
          tbl[i,"P value"]                      <- round(1-pchisq(2*(as.numeric(tbl[i,"Neg. log likel"])  - as.numeric(tbl[1,"Neg. log likel"])),as.numeric(tbl[i,"Degrees of freedom"])),4)
        }
      }
    }
    return(as.table(tbl))
  }
)

