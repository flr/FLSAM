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

    # year ranges match? If not, trim the longest to equal the shortest
    if(!identical(dny[['year']],dnx[['year']]))
    {
      #Get common range
      common.rng <- range(as.numeric(intersect(dnx$year,dny$year)))
      x <- window(x,start=common.rng[1],end=common.rng[2])
      x@stock.n <- window(y@stock.n,start=common.rng[1],end=common.rng[2])
      x@harvest <- window(y@harvest,start=common.rng[1],end=common.rng[2])
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
.extract.params <- function(object,params) {
  #Check that this is sensible first
  if(object@control@nohess & params %in% c("logssb","logtsb","logfbar","logCatch")) {
    stop(paste("Cannot extract SSB, Fbar, TSB or modelled catch from an FLSAM object that",
               "has been run with the nohess=TRUE option. To calculate these variables",
               "please update the corresponding FLStock object e.g. <FLStock> <-",
               "<FLStock> + <FLSAM> and the use the corresponding function"))}

  #Extract numbers at age parameters
  ps <- subset(object@params,name%in%params)[,c("name","value","std.dev")]
  ps$CV <- ps$std.dev
  ps$lbnd <- exp(ps$value - 1.96*ps$std.dev)
  ps$ubnd <- exp(ps$value + 1.96*ps$std.dev)
  ps$value <- exp(ps$value)
  ps$std.dev <- NULL
  ps$name <- NULL
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
            res[[i]] <- cbind(name=names(object)[i],ssb(object[[i]]))
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
            res[[i]] <- cbind(name=names(object)[i],fbar(object[[i]]))
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
            res[[i]] <- cbind(name=names(object)[i],tsb(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

setMethod("catch", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logCatch")   
          res <- cbind(year=seq(object@control@range["minyear"],
                            object@control@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

setMethod("catch", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=names(object)[i],catch(object[[i]]))
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
            res[[i]] <- cbind(name=names(object)[i],n(object[[i]]))
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
          res$age <- as.numeric(as.character(res$Var1))
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
            res[[i]] <- cbind(name=names(object)[i],f(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

setMethod("rec", signature(object="FLSAM"),
        function(object, ...) {
          rec.rtn <- subset(n(object),age==object@control@range["min"])
          rec.rtn$age <- NULL
          return(rec.rtn)
        }
)       # }}}

setMethod("rec", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=names(object)[i],rec(object[[i]]))
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

setGeneric("coef",useAsDefault=coef)
setGeneric("coefficients",useAsDefault=coefficients)
setMethod("coef",signature(object="FLSAM"),
        function(object,...) {
          return(object@params[1:object@nopar,])
        }
)
setMethod("coef", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(object.name=names(object)[i],coef(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}
setMethod("coefficients",signature(object="FLSAM"),
        function(object, ...) {
          return(coef(object))
        })
setMethod("coefficients",signature(object="FLSAMs"),
        function(object, ...) {
          return(coef(object))
        }
)

#-------------------------------------------------------------------------------
# Extract fleet parameters
# There are three types of fleet parameters - observation variances,
# catchabilities, and power law exponents. Due to the way SAM is 
# constructed, there is also a distinction between the parameters for
# age structured fleets and SSB fleets. Here we use a common function
# to try and deal with all these issues.
#-------------------------------------------------------------------------------
.extract.fleet.parameters <- function(object,type) {
       #Setup parameter type to extract
       ext <- switch(type,
         "observation variance"=list(age="logSdLogObs",ssb="logSdSSB",
                          slt="obs.vars",ssb.fleet=c(3,4)),
         "catchability"=list(age="logFpar",ssb="logScaleSSB",
                          slt="catchabilities",ssb.fleet=c(3,4)),
         "power law"=list(age="logQpow",ssb="logPowSSB",
                          slt="power.law.exps",ssb.fleet=4))
  
       #Extract data
       age.params <- .extract.params(object,ext$age)
       if(nrow(age.params)!=0) {
         age.params$no <- 1:nrow(age.params)
         bindings.mat <- slot(object@control,ext$slt) 
         bindings <-  as.data.frame(as.table(bindings.mat),responseName="no")
         #Merge
         bindings <- subset(bindings,!is.na(bindings$no))
         bindings$age <- as.numeric(as.character(bindings$age))
         age.res <- merge(bindings,age.params)
         age.res$no <- NULL
       } else {
         age.res <- .extract.params(object,"asfasfsdf") #space filler
       }

       #Add in SSB indices (if any)
       ssb.params <- .extract.params(object,ext$ssb)
       if(nrow(ssb.params)!=0) {
         ssb.fleets <- names(object@control@fleets)[object@control@fleets%in%ext$ssb.fleet]
         if(nrow(ssb.params)!=length(ssb.fleets)) {
           stop(paste("Number of fleets does not match number of",type,"parameters.",
                      "This is a bug. Please report it to the FLSAM users mailing",
                      "list, FLSAM@googlegroups.com"))}
         ssb.res <- cbind(fleet=ssb.fleets,age=NA,ssb.params)
       } else {
         ssb.res <- .extract.params(object,"asfasfsdf") #space filler
       }
       #Tidy up
       res <- rbind(age.res,ssb.res)
       if(nrow(res)==0) {
         stop(paste("FLSAM object does not contain",type,"parameters"))
       } else {
         return(res[order(res$fleet,res$age),])
       }
}


if (!isGeneric("catchabilities")) {
  setGeneric('catchabilities', function(object) standardGeneric('catchabilities'))
}
setMethod("catchabilities",signature(object="FLSAM"),
  function(object) {
    .extract.fleet.parameters(object,"catchability")
})

setMethod("catchabilities", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],catchabilities(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

if (!isGeneric("obs.var")) {
  setGeneric('obs.var', function(object) standardGeneric('obs.var'))
}
setMethod("obs.var",signature(object="FLSAM"),
    function(object) {
       .extract.fleet.parameters(object,"observation variance")
})

setMethod("obs.var", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=objects@name[i],obs.var(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

if (!isGeneric("power.law.exps")) {
  setGeneric('power.law.exps', function(object) standardGeneric('power.law.exps'))
}
setMethod("power.law.exps",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"power law")
})

setMethod("power.law.exps", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@name[i],power.law.exps(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

#-------------------------------------------------------------------------------
# Log-ratio test
#-------------------------------------------------------------------------------

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
    lr.test(FLSAMs(object,...),type=type)
  }
)

