# merge {{{
# merge results from FLSAM into FLStock
#' @rdname FLSAM-merge
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

#' @rdname FLSAM-merge
setMethod("+", signature(e1="FLStock", e2="FLSAM"),
  function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(merge(e1, e2))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)
#' @rdname FLSAM-merge
setMethod("+", signature(e1="FLSAM", e2="FLStock"),
	function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(merge(e2, e1))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)   # }}}

#' @rdname FLSAM-merge
setMethod("+", signature(e1="FLStock", e2="FLSAMs"),
  function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(FLStocks(lapply(e2,merge,x=e1)))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)
#' @rdname FLSAM-merge
setMethod("+", signature(e1="FLSAMs", e2="FLStock"),
	function(e1, e2) {
    if(validObject(e1) & validObject(e2))
      return(FLStocks(lapply(e1,merge,x=e2)))
    else
      stop("Input objects are not valid: validObject == FALSE")
    }
)   # }}}

#' @rdname FLSAM-merge
setMethod("+", signature(e1="FLSAM", e2="FLStocks"),
	function(e1, e2) {
    if(validObject(e1) & validObject(e2)){
    yrranges <- lapply(e2,function(x){return(x@range[c("minyear","maxyear")])})
    for(iStk in 1:length(e2)){
      x <- e2[[iStk]]
      y <- window(e1,start=yrranges[[iStk]]["minyear"],end=yrranges[[iStk]]["maxyear"])
      if(dims(x)$area > dims(y)$area){
        x@stock.n[] <- y@stock.n
        x@harvest[] <- y@harvest
      }
      if(dims(x)$area <= dims(y)$area){
        x@stock.n[] <- y@stock.n
        x@harvest[] <- areaSums(y@harvest,na.rm=T)
      }
      e2[[iStk]] <- x
    }
    return(e2)
    } else {
      stop("Input objects are not valid: validObject == FALSE")}
    }
)   # }}}

#' @rdname FLSAM-merge
setMethod("+", signature(e1="FLStocks", e2="FLSAM"),
	function(e1, e2) {
    eone <- e2;    etwo <- e1
    e1 <- eone;    e2 <- etwo
    if(validObject(e1) & validObject(e2)){
    yrranges <- lapply(e2,function(x){return(x@range[c("minyear","maxyear")])})
    for(iStk in 1:length(e2)){
      x <- e2[[iStk]]
      y <- window(e1,start=yrranges[[iStk]]["minyear"],end=yrranges[[iStk]]["maxyear"])
      if(dims(x)$area > dims(y)$area){
        x@stock.n[] <- y@stock.n
        x@harvest[] <- y@harvest
      }
      if(dims(x)$area <= dims(y)$area){
        x@stock.n[] <- y@stock.n
        x@harvest[] <- areaSums(y@harvest,na.rm=T)
      }
      e2[[iStk]] <- x
    }
    return(e2)
    } else {
      stop("Input objects are not valid: validObject == FALSE")}
    }
)   # }}}





#General helper function to extract a given parameter from an FLSAM object
#and return it as a data.frame
.extract.params <- function(object,params) {
  #Check that this is sensible first

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
  Ns <- .extract.params(object,c("logN"))
  Fs <- .extract.params(object,c("logF"))
  yrs <- seq(object@range["minyear"],
             object@range["maxyear"]) 
  ages <- seq(object@range["min"],
              object@range["max"])
  states <- c(ages,seq(-1,by=-1,length.out=object@n.states-length(ages)))
  storeNs <- expand.grid(age=states[states>=0],year=yrs)
  storeFs <- expand.grid(age=states[states<0],year=yrs)
  Us <- rbind(cbind(storeNs,Ns),cbind(storeFs,Fs))
  return(Us)
}

#ssb          {{{
#' @rdname accessors
setMethod("ssb", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logssb")   
          res <- cbind(year=seq(object@range["minyear"],
                            object@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

#' @rdname accessors
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

#' @rdname accessors
    setGeneric("components", function(object) standardGeneric("components"))
#' @rdname accessors
setMethod("components", signature(object="FLSAM"),
        function(object) {
        return(object@components)}
)       # }}}

#' @rdname accessors
setMethod("components", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=names(object)[i],components(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

#' @rdname accessors
setMethod("fbar", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logfbar")   
          res <- cbind(year=seq(object@range["minyear"],
                            object@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

#' @rdname accessors
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

#' @rdname accessors
setMethod("tsb", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logtsb")   
          res <- cbind(year=seq(object@range["minyear"],
                            object@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

#' @rdname accessors
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

#' @rdname accessors
setMethod("catch", signature(object="FLSAM"),
        function(object, ...) {
          res <- .extract.params(object,"logCatch")   
          res <- cbind(year=seq(object@range["minyear"],
                            object@range["maxyear"]),
                       res)
          return(res)
        }
)       # }}}

#' @rdname accessors
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

#' @rdname accessors
setMethod("n", signature(object="FLSAM"),
        function(object, ...) {
	  return(subset(.extract.states(object),age>=0))
        }
)       # }}}

#' @rdname accessors
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

#' @rdname accessors
setMethod("catch.n", signature(object="FLSAM"),
        function(object, ...) {
	  return(subset(residuals(object),age>=0 & fleet %in% names(which(object@control@fleets == 0))))
        }
)       # }}}

#' @rdname accessors
setMethod("catch.n", signature(object="FLSAMs"),
        function(object, ...) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=names(object)[i],catch.n(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}


#' @rdname accessors
    setGeneric("f", function(object) standardGeneric("f"))
#' @rdname accessors
setMethod("f", signature(object="FLSAM"),
        function(object) {
          f.states <- subset(.extract.states(object),age<0)
          f.states$param <- -f.states$age-1
          if(length(grep("catch",rownames(object@control@states)))>1){
            f.bindings <- object@control@states[grep("catch",rownames(object@control@states)),]
            f.bindings <- as.data.frame(as.table(f.bindings),responseName="param")
            f.bindings  <- subset(f.bindings,param>=0)
            res <- merge(f.states,f.bindings,by="param")
            res$age <- as.numeric(as.character(res$age.y))
          } else {
            f.bindings <- object@control@states[grep("catch",rownames(object@control@states)),]
            f.bindings <- as.data.frame(as.table(f.bindings),responseName="param")
            f.bindings  <- subset(f.bindings,param>=0)
            res <- merge(f.states,f.bindings,by="param")
            res$age <- as.numeric(as.character(res$Var1))
            res$fleet <- "catch unique"
          }
          res$param <- NULL
          res <- res[order(res$year,res$age),c("fleet","age","year","value","CV","lbnd","ubnd")]
          return(res)
        }
)       # }}}

#' @rdname accessors
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

#' @rdname accessors
setMethod("rec", signature(object="FLSAM"),
        function(object, ...) {
          rec.rtn <- subset(n(object),age==object@range["min"])
          rec.rtn$age <- NULL
          return(rec.rtn)
        }
)       # }}}

#' @rdname accessors
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

#' @rdname AIC-FLSAM
setMethod("AIC",signature(object="FLSAM"),
        function(object, ...) {
        return(2*object@nopar +2*object@nlogl)
        }
)

#' @rdname AIC-FLSAM
setMethod("AIC",signature(object="FLSAMs"),
        function(object, ...) {
        return(sapply(object,AIC,...))
        }
)

#' @rdname parameter.accessors
setGeneric("coef",useAsDefault=coef)
#' @rdname parameter.accessors
setGeneric("coefficients",useAsDefault=coefficients)
#' @rdname parameter.accessors
setMethod("coef",signature(object="FLSAM"),
        function(object,...) {
          return(object@params[1:object@nopar,])
        }
)
#' @rdname parameter.accessors
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
#' @rdname parameter.accessors
setMethod("coefficients",signature(object="FLSAM"),
        function(object, ...) {
          return(coef(object))
        })
#' @rdname parameter.accessors
setMethod("coefficients",signature(object="FLSAMs"),
        function(object, ...) {
          return(coef(object))
        }
)

#' @rdname simulate-FLSAM
setGeneric("simulate", function(x,y,z,n,set.pars,set.seed) standardGeneric("simulate"))

#' @rdname simulate-FLSAM
setMethod("simulate",signature(x="FLStock",y="FLIndices",z="FLSAM.control"),
          function(x,y,z,n=100,set.pars=NULL,set.seed=NULL){
            if(missing(n))        n        <- 100
            if(missing(set.pars)) set.pars <- NULL
            if(missing(set.seed)) set.seed <- NULL
            if(!is.null(set.pars))
              fit   <- FLSAM(x,y,z,return.fit=T,set.pars=set.pars)
            if(is.null(set.pars))
              fit   <- FLSAM(x,y,z,return.fit=T)
            if(!is.null(set.seed))
              set.seed(set.seed)
            sdrep <- sdreport(fit$obj,getJointPrecision=T)
            sigma <- as.matrix(solve(sdrep$jointPrecision))
            mu    <- c(sdrep$par.fixed,sdrep$par.random)
            sim   <- mvrnorm(n,mu=mu,Sigma=sigma)
            colnames(sim) <- rownames(sigma)
            rownames(sim) <- 1:n
            
            sim   <- SIM2FLR(sim,fit,z,set.pars)

          return(sim)})


#-------------------------------------------------------------------------------
# Extract fleet parameters
# There are three types of fleet parameters - observation variances,
# catchabilities, and power law exponents. Due to the way SAM is 
# constructed, there is also a distinction between the parameters for
# age structured fleets and SSB fleets. Here we use a common function
# to try and deal with all these issues.
#-------------------------------------------------------------------------------
.extract.fleet.parameters <- function(object,type) {
       #
  
       #Setup parameter type to extract
       ext <- switch(type,
         "observation variance"=list(age="logSdLogObs",
                          slt="obs.vars"),
         "catchability"=list(age="logFpar",
                          slt="catchabilities"),
         "power law"=list(age="logQpow",
                          slt="power.law.exps"),
         "catch multiplier"=list(age="logScale",
                          slt="scalePars"),
         "number variance"=list(age="logSdLogN",
                          slt="logN.vars"),
         "harvest variance"=list(age="logSdLogFsta",
                          slt="f.vars"),
         "component variance"=list(age="logSdLogP",
                          slt="logP.vars"),
         "observation correlation"=list(age="transfIRARdist",
                          slt="cor.obs"),
         "recruitment"=list(val=c("rec_pars"),
                          slt="srr"),
         "harvest correlation"=list(val="itrans_rho",
                          slt="cor.F"))

       #Extract data
       age.params <- .extract.params(object,ext$age)
       if(nrow(age.params)!=0) {
         age.params$no <- 0:(nrow(age.params)-1)
         bindings.mat <- slot(object@control,ext$slt) 
         bindings <-  as.data.frame(as.table(bindings.mat),responseName="no")
         if(length(grep("age",colnames(bindings)))==0)
          colnames(bindings)[1] <- "age"
         if(type %in% c("number variance","harvest variance")){
          if(length(grep("fleet",colnames(bindings)))>0){
            fleets   <- sort(unique(bindings$fleet))
            #bindings <- cbind(fleet=strsplit(type," ")[[1]][1],bindings[,-grep("fleet",colnames(bindings))])
            bindings <- cbind(fleet=paste(strsplit(type," ")[[1]][1],fleets),bindings[,-grep("fleet",colnames(bindings))])
          } else {
            bindings <- cbind(fleet=strsplit(type," ")[[1]][1],bindings)
          }
         }
         #Merge
         bindings <- subset(bindings,!is.na(bindings$no))
         if(type!="observation correlation")
           bindings$age <- as.numeric(as.character(bindings$age))
         age.res <- merge(bindings,age.params)
         age.res$no <- NULL
       } else {
         age.res <- .extract.params(object,"asfasfsdf") #space filler
       }
       
       #Extract data
       val.params <- .extract.params(object,ext$val)
       if(nrow(val.params)!=0){
        rownames(val.params) <- NULL
         if(type=="recruitment"){
           #rec_pars are reported both as parameter and as derived quantity: keep them once
           nRec     <- length(object@params$name[object@params$name=="rec_pars"])
           if(.hasSlot(object,"fit") && !is.null(object@fit)){
             nRec <- length(object@fit$pl$rec_pars)
           } else if(nRec %% 2 == 0 && isTRUE(all.equal(val.params$value[seq_len(nRec/2)],val.params$value[nRec/2+seq_len(nRec/2)]))){
             nRec <- nRec/2
           }
           val.params <- val.params[seq_len(min(nRec,nrow(val.params))),,drop=FALSE]
           val.res <- cbind(expand.grid(fleet=paste("srr",letters[seq_len(nrow(val.params))]),age=c(object@control@range["min"])),val.params)
         }
         if(type=="harvest correlation")
           val.res <- cbind(expand.grid(fleet=names(object@control@fleets[which(object@control@fleets==0)]),age="all"),val.params)
       } else {
         val.res <- .extract.params(object,"asfasfsdf") #space filler
       }


       #Tidy up
       res <- rbind(age.res,val.res)
       if(nrow(res)==0) {
         stop(paste("FLSAM object does not contain",type,"parameters"))
       } else {
         return(res[order(res$fleet,res$age),])
       }
}


#' @rdname accessors
  setGeneric('catchabilities', function(object) standardGeneric('catchabilities'))
#' @rdname accessors
setMethod("catchabilities",signature(object="FLSAM"),
  function(object) {
    .extract.fleet.parameters(object,"catchability")
})

#' @rdname accessors
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

#' @rdname accessors
  setGeneric('obs.var', function(object) standardGeneric('obs.var'))
#' @rdname accessors
setMethod("obs.var",signature(object="FLSAM"),
    function(object) {
       .extract.fleet.parameters(object,"observation variance")
})

#' @rdname accessors
setMethod("obs.var", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],obs.var(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

#' @rdname accessors
  setGeneric('power.law.exps', function(object) standardGeneric('power.law.exps'))
#' @rdname accessors
setMethod("power.law.exps",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"power law")
})

#' @rdname accessors
setMethod("power.law.exps", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],power.law.exps(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

#' @rdname accessors
  setGeneric('catch.scale', function(object) standardGeneric('catch.scale'))
#' @rdname accessors
setMethod("catch.scale",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"catch multiplier")
})

#' @rdname accessors
setMethod("catch.scale", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],catch.scale(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}


#' @rdname accessors
  setGeneric('n.var', function(object) standardGeneric('n.var'))
#' @rdname accessors
setMethod("n.var",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"number variance")
})

#' @rdname accessors
setMethod("n.var", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],n.var(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

#' @rdname accessors
  setGeneric('f.var', function(object) standardGeneric('f.var'))
#' @rdname accessors
setMethod("f.var",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"harvest variance")
})

#' @rdname accessors
setMethod("f.var", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],f.var(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}


#' @rdname accessors
  setGeneric('cor.F', function(object) standardGeneric('cor.F'))
#' @rdname accessors
setMethod("cor.F",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"harvest correlation")
})

#' @rdname accessors
setMethod("cor.F", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],cor.F(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}


#' @rdname accessors
  setGeneric('cor.obs', function(object) standardGeneric('cor.obs'))
#' @rdname accessors
setMethod("cor.obs",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"observation correlation")
})

#' @rdname accessors
setMethod("cor.obs", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],cor.obs(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

#' @rdname accessors
  setGeneric('rec.pars', function(object) standardGeneric('rec.pars'))
#' @rdname accessors
setMethod("rec.pars",signature(object="FLSAM"),
   function(object) {
    .extract.fleet.parameters(object,"recruitment")
})

#' @rdname accessors
setMethod("rec.pars", signature(object="FLSAMs"),
        function(object) {
          res <- list()
          length(res) <- length(object)
          for(i in seq(object)) {
            res[[i]] <- cbind(name=object@names[i],rec.pars(object[[i]]))
          }
          return(do.call(rbind,res))
        }
)       # }}}

#-------------------------------------------------------------------------------
# Log-ratio test
#-------------------------------------------------------------------------------

#- Create generic function for 'lr.test'
#' @rdname lr.test
  setGeneric('lr.test', function(object,...,type="sequential") standardGeneric('lr.test'))

#' @rdname lr.test
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
    
    
    
#' @rdname lr.test
setMethod("lr.test",signature("FLSAM"),
  function(object,...,type="sequential"){
    lr.test(FLSAMs(object,...),type=type)
  }
)

#' Mohn's rho of a retrospective analysis
#'
#' Calculates the relative difference between the terminal-year estimate of
#' each retrospective peel and the estimate for the same year in the
#' reference (most recent) assessment, following Mohn (1999) and
#' \code{\link[stockassessment]{mohn}}: \code{100 * (peel - ref) / ref}.
#' A positive rho means that the peels overestimated the quantity.
#'
#' @param retro An \code{\link{FLSAMs}} object, as returned by \code{\link{retro}},
#'   named by assessment year.
#' @param ref.year The reference assessment year (default: the most recent).
#' @param span Number of peels to include (default 5).
#' @param type Quantity: \code{"fbar"}, \code{"ssb"} or \code{"rec"}.
#' @param ... Not used.
#' @return A \code{data.frame} with columns \code{rho} (\%) and \code{year},
#'   one row per peel (the reference year has rho 0). The mean rho over the
#'   peels (Mohn's rho) is stored in the attribute \code{"mohn"}.
#' @note Versions of FLSAM before 3.1 returned \code{100 * (1 - peel/ref)},
#'   i.e. the opposite sign, and \code{ref.year} had no working default.
#' @references Mohn, R. (1999). The retrospective problem in sequential
#'   population analysis. ICES Journal of Marine Science 56: 473-488.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' r <- retro(NSH, NSH.tun, NSH.ctrl, retro=3)
#' mohns.rho(r, type="ssb")
#' attr(mohns.rho(r, type="fbar"), "mohn")
#' }
#' @rdname mohns.rho
#' @export
  setGeneric('mohns.rho', function(retro,ref.year,span,type,...) standardGeneric('mohns.rho'))
#' @rdname mohns.rho
setMethod("mohns.rho",signature("FLSAMs"),
function(retro,ref.year,span,type,...){
  if(missing(ref.year)) ref.year <- max(an(names(retro)),na.rm=T)
  if(missing(span))     span     <- 5
  if(missing(type))     type     <- "fbar"
  if(type == "fbar")
    retro.fbar  <- lapply(retro,fbar)
  if(type == "ssb")
    retro.fbar  <- lapply(retro,ssb)
  if(type == "rec")
    retro.fbar  <- lapply(retro,rec)
  if(ref.year>max(an(names(retro)),na.rm=T)) stop("Reference year set higher than the latest assessment year")
  if(span >= length(retro)){
    span <- length(retro)-1
    warning("Span longer than number of retro years available")
  }
  retro.fbar  <- retro.fbar[which(names(retro.fbar) %in% ref.year:(ref.year-span))]
  rho         <- data.frame(rho=NA, year=c((ref.year-span):ref.year))
  termFs      <- do.call(rbind,lapply(as.list(names(retro.fbar)),
                         function(x){return(retro.fbar[[ac(x)]][which(retro.fbar[[ac(x)]]$year==an(x)),])}))
  refFs       <- retro.fbar[[ac(ref.year)]]
  refFs       <- refFs[which(refFs$year %in% termFs$year),]
  colnames(refFs) <- c("year","refvalue","refCV","reflbnd","refubnd")
  combFs      <- merge(refFs,termFs,by="year")
  combFs      <- combFs[order(combFs$year),]
  rhos        <- (combFs$value/combFs$refvalue-1)*100
  rho$rho     <- rhos
  attr(rho,"mohn") <- mean(rho$rho[rho$year!=ref.year])
return(rho)})
