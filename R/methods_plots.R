#' Plot FLSAM assessment results
#'
#' Stock summary plots (SSB, Fbar and recruitment with 95\% confidence
#' intervals) of one or more assessments, optionally together with an
#' \code{FLStock}. The plots are made with \pkg{ggplot2} and return a
#' \code{ggplot} object that can be modified further.
#'
#' @param x An \code{\link{FLSAM-class}}, \code{\link{FLSAMs}} or \code{FLStock} object.
#' @param y Missing, or an \code{FLSAM}, \code{FLSAMs} or \code{FLStock} object
#'   to plot together with \code{x}.
#' @param futureYrs Include the last year (in which only survey information
#'   is available)? Default \code{TRUE}.
#' @param what Not used (kept for backward compatibility).
#' @param main Plot title.
#' @param ... Not used.
#' @return A \code{ggplot} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' plot(NSH.sam)
#' data(HERAS.sams)
#' plot(HERAS.sams)
#' @name plot-FLSAM
#' @aliases plot,FLSAM,missing-method plot,FLSAM,FLSAM-method plot,FLSAMs,missing-method plot,FLSAMs,FLStock-method plot,FLSAM,FLStock-method plot,FLStock,FLSAM-method plot,FLStock,FLSAMs-method
NULL

.plotLevels <- c("Spawning stock biomass","Fishing mortality","Recruitment")

#' @rdname plot-FLSAM
#' @export
setMethod("plot",signature(x="FLSAM",y="missing"),
 function(x,y,futureYrs=T,what=c("ssb","fbar","rec"),...) {
    #Extract data
    ssb.dat <- ssb(x)
    ssb.dat$name <- "Spawning stock biomass"
    rec.dat <- rec(x)
    rec.dat$name <- "Recruitment"
    rec.dat$age <- NULL
    fbar.dat <- fbar(x)
    fbar.dat$name <- "Fishing mortality"
    plot.dat    <- rbind(ssb.dat, fbar.dat,rec.dat)
    plot.dat$name <- factor(plot.dat$name,levels=.plotLevels)
    if(!futureYrs) plot.dat    <- subset(plot.dat,year %in% range(x)["minyear"]:(range(x)["maxyear"]-1))
    #Plot it
    ggplot(plot.dat,aes(x=year)) +
      geom_ribbon(aes(ymin=lbnd,ymax=ubnd),fill="grey",colour=NA) +
      geom_line(aes(y=value),colour="black",linewidth=1) +
      facet_wrap(~name,ncol=1,scales="free_y") +
      expand_limits(y=0) +
      labs(x="Year",y=NULL,title=x@name) +
      theme_bw()
})

#' @rdname plot-FLSAM
#' @export
setMethod("plot",signature(x="FLSAM",y="FLSAM"),
 function(x,y,main="",...) {
    plot(FLSAMs(x,y),...)
})

#' @rdname plot-FLSAM
#' @export
setMethod("plot",signature(x="FLSAMs",y="missing"),
 function(x,y,main="",futureYrs=T,...) {
    #Extract data - ssb and fbar are easy but recs harder
    ssb.dat <- ssb(x)
    if(!futureYrs) ssb.dat <- ssb(x)[-c(which(diff(ssb(x)$year)<0),nrow(ssb(x))),]
    ssb.dat$var <- "Spawning stock biomass"
    rec.dat <- rec(x)
    if(!futureYrs) rec.dat <- rec(x)[-c(which(diff(rec(x)$year)<0),nrow(rec(x))),]
    rec.dat$var <- "Recruitment"
    fbar.dat <- fbar(x)
    if(!futureYrs) fbar.dat <- fbar(x)[-c(which(diff(fbar(x)$year)<0),nrow(fbar(x))),]
    fbar.dat$var <- "Fishing mortality"
    plot.dat    <- rbind(ssb.dat, fbar.dat,rec.dat)
    plot.dat$var <- factor(plot.dat$var,levels=.plotLevels)
    plot.dat$name <- factor(plot.dat$name,levels=unique(plot.dat$name))

    #Plot it
    ggplot(plot.dat,aes(x=year,y=value,colour=name)) +
      geom_line() +
      facet_wrap(~var,ncol=1,scales="free_y") +
      expand_limits(y=0) +
      labs(x="Year",y=NULL,colour=NULL,title=main) +
      theme_bw()
})

#' @rdname plot-FLSAM
#' @export
setMethod("plot",signature(x="FLSAMs",y="FLStock"),
 function(x,y,main="",futureYrs=TRUE,...) {
    #extract data from FLSAMs first
    ssb.dat <- ssb(x)
    if(!futureYrs) ssb.dat <- ssb(x)[-c(which(diff(ssb(x)$year)<0),nrow(ssb(x))),]
    ssb.dat$var <- "Spawning stock biomass"
    rec.dat <- rec(x)
    if(!futureYrs) rec.dat <- rec(x)[-c(which(diff(rec(x)$year)<0),nrow(rec(x))),]
    rec.dat$var <- "Recruitment"
    fbar.dat <- fbar(x)
    if(!futureYrs) fbar.dat <- fbar(x)[-c(which(diff(fbar(x)$year)<0),nrow(fbar(x))),]
    fbar.dat$var <- "Fishing mortality"
    sams.dat    <- rbind(ssb.dat, fbar.dat,rec.dat)

    #Then add the data from the stock
    stck.dat <- rbind(data.frame(as.data.frame(rec(y)),var="Recruitment"),
                  data.frame(as.data.frame(ssb(y)),var="Spawning stock biomass"),
                  data.frame(as.data.frame(fbar(y)),var="Fishing mortality"))
    stck.dat$name <- y@name
    stck.dat$value <- stck.dat$data

    #Now combine
    common.cols <- intersect(names(stck.dat),names(sams.dat))
    plot.dat <- rbind(stck.dat[,common.cols],sams.dat[,common.cols])
    plot.dat$var <- factor(plot.dat$var,levels=.plotLevels)
    plot.dat$name <- factor(plot.dat$name,levels=unique(plot.dat$name))

    if(!futureYrs) plot.dat    <- subset(plot.dat,year %in% min(unlist(lapply(x,function(y){range(y)["minyear"]}))) :
                                                            max(unlist(lapply(x,function(y){range(y)["maxyear"]-1}))))
    #plot it
    ggplot(plot.dat,aes(x=year,y=value,colour=name)) +
      geom_line() +
      geom_point(size=0.8) +
      facet_wrap(~var,ncol=1,scales="free_y") +
      expand_limits(y=0) +
      labs(x="Year",y=NULL,colour=NULL,title=main) +
      theme_bw()
})

#' @rdname plot-FLSAM
#' @export
setMethod("plot",signature(x="FLSAM",y="FLStock"),
 function(x,y,main="",...) {
    plot(FLSAMs(x),y,main=main,...)
})

#' @rdname plot-FLSAM
#' @export
setMethod("plot",signature(x="FLStock",y="FLSAM"),
 function(x,y,main="",...) {
    plot(y,x,main=main,...)
})

#' @rdname plot-FLSAM
#' @export
setMethod("plot",signature(x="FLStock",y="FLSAMs"),
 function(x,y,main="",...) {
    plot(y,x,main=main,...)
})

#' Components plot
#'
#' Stacked area plot of the estimated proportions of the components of a
#' partial survey (fleet type 6).
#'
#' @param sam An \code{\link{FLSAM-class}} object with components.
#' @return Called for its side effect (base graphics).
#' @export
comp.plot <- function(sam){
  require(RColorBrewer)
  if(length(sam@components)==0) stop("No component information in assessment output")
  yrs <- as.numeric(range(dimnames(sam@components)$years)[1]):as.numeric(range(dimnames(sam@components)$years)[2])
  plot(0,0,type="n",xlim=range(yrs),ylim=c(0,1),yaxs="i",xaxs="i",xlab="Year",ylab="Fraction")
  area.plot.dat <- apply(sam@components,2,function(x) 1-cumsum(c(0,x))/sum(x))
  cols<- rev(brewer.pal(nrow(sam@components),"PuOr"))
  for(i in 1:nrow(sam@components)) {
    x.to.plot <- c(yrs,rev(yrs))
    y.to.plot <- c(area.plot.dat[i,],rev(area.plot.dat[i+1,]))
    polygon(x.to.plot,y.to.plot,col=cols[i])
  }
  legend("topleft",bg="white",legend=rownames(sam@components),pch=22,pt.bg=cols,pt.cex=2)
}

#' Partial fishing mortality plot
#'
#' Plots the mean fishing mortality over the Fbar age range for each catch
#' fleet of a multi-fleet assessment.
#'
#' @param sam An \code{\link{FLSAM-class}} object.
#' @return A \code{ggplot} object.
#' @export
fpart.plot <- function(sam){
  f.dat <- f(sam)
  fbars <- sam@control@range["minfbar"]:sam@control@range["maxfbar"]
  f.dat <- subset(f.dat,age %in% fbars)
  f.dat <- aggregate(f.dat[,"value"],by=as.list(f.dat[,c("fleet","year")]),FUN=mean,na.rm=T)
  ggplot(f.dat,aes(x=year,y=x,colour=fleet)) +
    geom_line() +
    geom_point() +
    labs(x="Year",y="Fbar-partial",colour=NULL,title=sam@name) +
    theme_bw()
}

#' Parameter correlation plot
#'
#' Heat map of the correlation matrix of the estimated parameters (derived
#' from the \code{vcov} slot). This is the FLSAM equivalent of
#' \code{\link[stockassessment]{parplot}} / \code{corplot}.
#'
#' @param sam An \code{\link{FLSAM-class}} object.
#' @param cols Three colours for correlations of -1, 0 and 1.
#' @param full Show the full matrix (\code{TRUE}) or only the fitted
#'   parameters (\code{FALSE}, default).
#' @return A \code{ggplot} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' cor.plot(NSH.sam)
#' @export
cor.plot <- function(sam,cols=c("red","white","blue"),full=FALSE) {
   n.coefs <- nrow(coef(sam))
   cor.mat <- cov2cor(sam@vcov)
   if(!full) {cor.mat <- cor.mat[1:n.coefs,1:n.coefs]} #Don't show the full matrix
   var.names <- colnames(cor.mat)
   var.names <- factor(var.names,unique(var.names))
   var.names <- paste(var.names,do.call(c,lapply(table(var.names),seq,from=1)),sep=".")
   dimnames(cor.mat) <- list(Var2=var.names,Var1=var.names)
   plt.dat <- as.data.frame(as.table(cor.mat),responseName="cor")
   ggplot(plt.dat,aes(x=Var2,y=Var1,fill=cor)) +
     geom_tile() +
     scale_fill_gradient2(low=cols[1],mid=cols[2],high=cols[3],midpoint=0,limits=c(-1,1)) +
     labs(x=NULL,y=NULL,fill="Correlation",title=sam@name) +
     theme_bw() +
     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
}

#Otolith plot (resample from the variance co-variance matrix and create a plot of all the resampled value)
#' Joint uncertainty of SSB and Fbar ("otolith" plot)
#'
#' Draws SSB and Fbar in a given year from their estimated joint
#' (multivariate normal) distribution and shows the two-dimensional
#' confidence region with contours, optionally with marginal densities.
#'
#' @param object An \code{\link{FLSAM-class}} object.
#' @param year Year of the plot (default the last year).
#' @param plot Draw the plot.
#' @param show.points Show the individual draws.
#' @param do.contours Draw contours of the cumulative density.
#' @param margin.plots Add marginal density plots.
#' @param show.estimate Show the point estimate.
#' @param n Number of draws.
#' @param pch Plot symbol of the draws.
#' @param alpha Confidence level of the marginal intervals (default 0.05).
#' @param show.grid Draw a grid.
#' @param n.grid Grid size of the kernel density estimate.
#' @param contour.args List of arguments passed to \code{contour}.
#' @param ... Not used.
#' @return The contour object (or the draws) invisibly.
#' @examples
#' data(NSH); data(NSH.sam)
#' otolith(NSH.sam, year=2010, n=1000)
#' @rdname otolith
  setGeneric("otolith", function(object, ...)
    standardGeneric("otolith"))
#' @rdname otolith
setMethod("otolith", signature(object="FLSAM"),
function(object, year=object@range["maxyear"], plot=TRUE, show.points=FALSE, do.contours=TRUE,
                      margin.plots=TRUE, show.estimate=TRUE,
                      n=100, pch=".", alpha=0.05, show.grid=TRUE,
                      n.grid=50, contour.args=list(),...){

  debug <- FALSE
  filled.contours <- FALSE

  #Get the range of ages for fbar
  f.ages        <-  object@range["minfbar"]:object@range["maxfbar"]

  #set up year ranges and identify year which plot is being done for
  years <- object@range["minyear"]:object@range["maxyear"]
  year.index <- which(years==year)
  
  # Generate n iterations of all the parameters
  if(n > 200){
    Fbar.all  <- numeric()
    SSB.all   <- numeric()
    paramvalue<- subset(object@params,name%in%c("logssb","logfbar","logCatch","logCatchByFleet","logLand",
                                                "logtsb","logR","logLagR","lastLogN","lastLogF",
                                                "beforeLastLogN","beforeLastLogF"))   
    idx       <- is.finite(paramvalue$value)
    for(i in 1:ceiling(n/200)){
      d <- mvrnorm(n=200,paramvalue$value[idx], object@rescov[idx,idx])
      colnames(d) <- paramvalue$name[idx]
      Fbar.all  <- rbind(Fbar.all,d[,colnames(d)=="logfbar"])
      SSB.all   <- rbind(SSB.all,d[,colnames(d)=="logssb"])
    }
    Fbar.all  <- Fbar.all[1:n,]
    SSB.all   <- SSB.all[ 1:n,]
  } else {
      pars <- subset(object@params,name%in%c("logssb","logfbar","logCatch","logCatchByFleet","logLand",
                                                "logtsb","logR","logLagR","lastLogN","lastLogF",
                                                "beforeLastLogN","beforeLastLogF")) 
      idx <- is.finite(pars$value)
      d <- mvrnorm(n=n,pars$value[idx], object@rescov[idx,idx])
      colnames(d) <- pars$name[idx]
      Fbar.all <- d[,colnames(d)=="logfbar"]
      SSB.all <-  d[,colnames(d)=="logssb"]
    }

  #Set these points up to be returned and return them to normal space
  Fbar <- exp(Fbar.all[,year.index])
  SSB <- exp(SSB.all[,year.index])
  return.obj  <- data.frame(Fbar,SSB)

  #Extract SSB and Fbar estimates from the assessment outputs for the year requested
  SSB.est   <-  ssb(object)$value[ssb(object)$year==year]
  Fbar.est  <-  fbar(object)$value[fbar(object)$year==year]


  #Now calculate the contour lines, if requested
  if(do.contours) {
    #Calculate kernel density estimate
    kern  <-  kde2d(Fbar,SSB,n=n.grid)
    #Calculate cumulative distribution function
    kz    <-  as.vector(kern$z)
    ord   <-  order(kz)
    cumfrac <-  cumsum(kz[ord])/sum(kz)
    cumfrac.matrix  <-  matrix(cumfrac[rank(kz)],nrow=nrow(kern$z),ncol=ncol(kern$z))
    if(is.null(contour.args$levels)) contour.args$levels <-   c(0.01,0.05,0.25,0.5,0.75)
#      if(is.null(contour.args$lty))   contour.args$lty <-   c(1,1,2,3,4)
#      if(is.null(contour.args$lwd))   contour.args$lwd <-   c(1,3,1,1,1)
    if(is.null(contour.args$method))contour.args$method <-    "edge"
#      if(filled.contours) {
#       do.call(filled.contour,c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,nlevels=100,color.palette=heat.colors))
#    }
    otolith.obj  <-  c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,contour.args)
    return.obj    <- otolith.obj
  }

  #Now do the plot
  if(plot) {
    if(!show.points) pch <- NA
    if(margin.plots) {
      par(mar=c(0.5,0.5,0.5,0.5),oma=c(4,4,1.5,0))
      layout(matrix(c(1,4,3,2),2,2,byrow=TRUE), c(3,lcm(4)), c(lcm(4),3), respect=FALSE)
    } else { #Else force to a 1x1 plot
      par(mfrow=c(1,1))
    }

    x.lab <-  paste("Fbar (",ifelse(is.null(f.ages),"All Ages",paste(range(f.ages),collapse="-")),")",sep="")
    xlim  <- range(pretty(Fbar))
    ylim  <- range(pretty(SSB))

    #First the horizontal plot
    if(margin.plots) {
      densF   <-  density(Fbar)
      plot(densF,ann=FALSE,xaxt="n",yaxt="n",type="l",xlim=xlim)
      if(show.grid) grid()
      title(ylab="Probability\nDensity",xpd=NA,mgp=c(1,1,0))

      #Calculate 95% confidence intervals
      cumsumF.fun <-  approxfun(cumsum(densF$y)/sum(densF$y),densF$x)
      densF.fun   <-  approxfun(densF$x,densF$y)
      ul.F    <-  cumsumF.fun(1-alpha/2)
      ul.dens <-  densF.fun(ul.F)
      ll.F    <-  cumsumF.fun(alpha/2)
      ll.dens <-  densF.fun(ll.F)
      points(c(ll.F,ul.F),c(ll.dens,ul.dens),pch="|",cex=1.5)
      text(c(ll.F,ul.F),c(ll.dens,ul.dens),label=sprintf("%.3f",round(c(ll.F,ul.F),3)),pos=4,xpd=NA)
      if(show.estimate) {
        points(Fbar.est,densF.fun(Fbar.est),pch=19,cex=1.5)
        text(Fbar.est,densF.fun(Fbar.est),label=sprintf("%.3f",round(Fbar.est,3)),pos=4,xpd=NA)
      }
    }

    #Now the vertical plot
    if(margin.plots) {
      densSSB <-  density(SSB)
      plot(densSSB$y,densSSB$x,xaxt="n",yaxt="n",type="l",ylim=ylim)
      abline(v=0,col="grey")
      if(show.grid) grid()
      title(xlab="Probability\nDensity",xpd=NA,mgp=c(2,1,0))
      #Calculate 95% confidence intervals
      cumsumSSB.fun <-  approxfun(cumsum(densSSB$y)/sum(densSSB$y),densSSB$x)
      densSSB.fun   <-  approxfun(densSSB$x,densSSB$y)
      ul.SSB    <-  cumsumSSB.fun(1-alpha/2)
      ul.dens <-  densSSB.fun(ul.SSB)
      ll.SSB    <-  cumsumSSB.fun(alpha/2)
      ll.dens <-  densSSB.fun(ll.SSB)
      points(c(ll.dens,ul.dens),c(ll.SSB,ul.SSB),pch="-",cex=2)
      text(c(ll.dens,ul.dens),c(ll.SSB,ul.SSB),label=round(c(ll.SSB,ul.SSB),0),pos=4,xpd=NA)
      if(show.estimate) {
        points(densSSB.fun(SSB.est),SSB.est,pch=19,cex=1.5)
        text(densSSB.fun(SSB.est),SSB.est,label=round(SSB.est,0),pos=2,xpd=NA)
      }
    }

    #Now the main plot
    plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="")
    if(show.points) points(Fbar,SSB,pch=pch)
    title(xlab=x.lab,ylab="SSB",xpd=NA)
    if(show.estimate) points(Fbar.est,SSB.est,pch=19,cex=1.5)
    if(show.grid) grid()
    if(do.contours) {
      do.call(contour,otolith.obj)
    }

  }

  return(invisible(return.obj))
})

#-Observation variance plot
#' Observation variance plots
#'
#' \code{obsvar.plot} shows the estimated observation standard deviations by
#' data source; \code{obscv.plot} plots them against their CV;
#' \code{obscor.plot} shows the estimated observation correlation matrices of
#' fleets with an \code{"AR"} or \code{"US"} correlation structure. See
#' \code{\link{obssd.plot}} and \code{\link{obscor.empirical.plot}} for
#' ggplot2 alternatives.
#'
#' @param sam An \code{\link{FLSAM-class}} object.
#' @return Called for its side effect (base graphics).
#' @examples
#' data(NSH); data(NSH.sam)
#' obsvar.plot(NSH.sam)
#' obscv.plot(NSH.sam)
#' @rdname obsvar.plot
#' @export
obsvar.plot <- function(sam){
                  obv <- obs.var(sam)
                  obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
                  obv <- obv[order(obv$value),]
                  bp <- barplot(obv$value,ylab="Observation Variance",
                         main="Observation variances by data source",col=factor(obv$fleet))
                  axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
                  legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)
              }

#-CV versus observation variance plot
#' @rdname obsvar.plot
#' @export
obscv.plot  <- function(sam){

  obv <- obs.var(sam)
  obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
  obv <- obv[order(obv$value),]
  plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
    pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
  text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
}

#-Observation correlation plot
#' @rdname obsvar.plot
#' @export
obscor.plot <- function(sam){
  oc <- sam@obscov
  idx <- apply(sam@control@cor.obs,1,function(x){any(x>=0,na.rm=T)})
  for(i in 1:length(oc)){
    ridx <- which(!is.na(sam@control@cor.obs[idx,]))
    dimnames(oc[[i]])[[1]] <- ridx[1]:(nrow(oc[[i]])+ridx[1]-1)
    dimnames(oc[[i]])[[2]] <- dimnames(oc[[i]])[[1]]
  }
  cols <- ifelse(length(oc)==1,1,2)
  par(mfrow=c(ceiling(length(oc)/2),cols))
  for(i in 1:length(oc))
    plotcorr(cov2cor(oc[[i]]),main=names(oc)[i])
}

#-Processes error plot: expressed as deviation in type = "mort","n","tsb"
#' Process error plots
#'
#' Shows the process error of an assessment, i.e. the difference between the
#' estimated numbers at age and the numbers predicted from last year's numbers
#' with the survivor equation, expressed in numbers (\code{type="n"}),
#' mortality (\code{"mort"}) or total biomass (\code{"tsb"}). Use an
#' \code{FLStock} updated with the assessment results (\code{stock + sam}).
#'
#' @param stck An \code{FLStock} with estimated \code{stock.n} and \code{harvest}.
#' @param weight Weights used for \code{type="tsb"}: \code{"stock.wt"} or \code{"catch.wt"}.
#' @param type \code{"n"}, \code{"mort"} or \code{"tsb"}.
#' @param rel Express the deviations relative (\%) instead of absolute.
#' @return The \code{ggplot} object (printed, returned invisibly).
#' @examples
#' data(NSH); data(NSH.sam)
#' NSH <- NSH + NSH.sam
#' procerr.plot(NSH, type="n")
#' procerr.plot(NSH, type="mort", rel=TRUE)
#' procerr.plot(NSH, type="tsb", weight="stock.wt")
#' @export
procerr.plot <- function(stck,weight="stock.wt",type="n",rel=F){

      #- convert N-at-age, F-at-age and m-at-age
      ns      <- as.data.frame(stck@stock.n) #N-at-age predicted by SAM
      colnames(ns)[7]       <- "n"
      fs      <- as.data.frame(stck@harvest) #F-at-age predicted by SAM
      colnames(fs)[7]       <- "f"
      ms      <- as.data.frame(stck@m)
      colnames(ms)[7]       <- "m"
      zs      <- merge(fs,ms,by=c("age","year","unit","season","area","iter"))
      zs$z    <- zs$f+zs$m
      if(weight=="stock.wt") ws      <- as.data.frame(stck@stock.wt)
      if(weight=="catch.wt") ws      <- as.data.frame(stck@catch.wt)
      colnames(ws)[7]       <- "ws"

      #- Calculate difference in Ns
      predns  <- ns; predns$year <- predns$year-1; predns$age <- predns$age-1   #N-at-age predicted by SAM, shift one year to align
      colnames(predns)[7]   <- "predn"
      ncomb   <- merge(ns,predns,by=c("age","year","unit","season","area","iter"))
      sens    <- as.data.frame(stck@stock.n * exp(-stck@harvest - stck@m))      #N-at-age predicted by survivor equation
      colnames(sens)[7]     <- "sens"
      ncomb   <- merge(ncomb,sens,by=c("age","year","unit","season","area","iter"))
      ncomb$ndiff <- ncomb$predn - ncomb$sens

      #- Calculate difference in mortality
      nmcomb  <- ncomb
      nmcomb$predmort  <- log(nmcomb$n/nmcomb$predn)
      nmcomb  <- merge(nmcomb,zs,by=c("age","year","unit","season","area","iter"))
      nmcomb$mdiff <- nmcomb$predmort - nmcomb$z
      
      #- Calculate difference in biomass
      nmwcomb <- nmcomb
      nmwcomb <- merge(nmwcomb,ws,by=c("age","year","unit","season","area","iter"))
      nmwcomb$bdiff <- nmwcomb$ndiff * nmwcomb$ws
      nmwcomb$bnorm <- nmwcomb$n * nmwcomb$ws
      bdiff   <- aggregate(nmwcomb$bdiff,by=as.list(nmwcomb[,c("year","unit","season","area","iter")]),FUN=sum,na.rm=T)
      bnorm   <- aggregate(nmwcomb$bnorm,by=as.list(nmwcomb[,c("year","unit","season","area","iter")]),FUN=sum,na.rm=T)
      colnames(bdiff)[6] <- "bdiff"
      colnames(bnorm)[6] <- "bnorm"
      
      #- Make the plots
      if(type=="n"){
        #- Deviance from n
        if(rel)
          nmwcomb[,"ndiff"] <- nmwcomb[,"ndiff"] / nmwcomb[,"n"] * 100
        inp   <- nmwcomb[,c(1:6,grep("ndiff",colnames(nmwcomb)))]
        colnames(inp)[7] <- "data"
        if(rel)
          yTxt <- "% difference numbers at age"
        if(rel==F)
          yTxt <- "Numbers at age"
        pic   <- ggplot(inp,aes(x=year,y=data)) +
                   geom_col(fill="black") +
                   facet_wrap(~age,scales="free_y") +
                   labs(x="Year",y=yTxt,title="Process error deviation in N") +
                   theme_bw()
      }
      if(type=="mort"){
        #- Deviance from mortality
        if(rel){
          inp   <- nmwcomb[,c(1:6,grep("mdiff",colnames(nmwcomb)))]
          inp$mdiff <- nmwcomb[,"mdiff"] / nmwcomb[,"z"] * 100
        }
        if(rel==F)
          inp   <- nmwcomb[,c(1:6,grep("mdiff",colnames(nmwcomb)))]
        colnames(inp)[7] <- "data"
        inp$sign <- factor(ifelse(inp$data<0,"negative","positive"),levels=c("negative","positive"))
        if(rel)
          yTxt <- "% difference mortality (year)"
        if(rel==F)
          yTxt <- "Mortality (year)"
        pic   <- ggplot(inp,aes(x=year,y=age)) +
                   geom_point(aes(size=abs(data),colour=sign),shape=1) +
                   scale_colour_manual(values=c(negative="black",positive="red"),drop=FALSE) +
                   scale_size_area(max_size=8) +
                   labs(x="Years",y="Age",size=yTxt,colour=NULL,title="Process error deviation in M") +
                   theme_bw()
      }
      if(type=="tsb"){
        if(rel)
          bdiff$bdiff <- bdiff$bdiff/bnorm$bnorm *100
        if(rel)
          yTxt <- "% difference in biomass (t)"
        if(rel==F)
          yTxt <- "Biomass (t)"
        pic <- ggplot(bdiff,aes(x=year,y=bdiff)) +
                 geom_col(fill="black") +
                 labs(x="Year",y=yTxt,title="Process error deviation in biomass") +
                 theme_bw()
      }
      print(pic)
      return(invisible(pic))
}

#' Retrospective diagnostic plots
#'
#' \code{retroResiduals} plots the residual pattern of one fleet across the
#' peels of a retrospective analysis, \code{retroSelectivity} the fishing
#' mortality at age, and \code{retroParams} the estimated parameters by
#' assessment year. The plots are printed and returned invisibly as
#' \code{ggplot} objects.
#'
#' @param x An \code{\link{FLSAMs}} object from \code{\link{retro}}.
#' @param fleet Name of the fleet.
#' @param yrs Years to include.
#' @param cohort Plot residuals by cohort instead of by year.
#' @return A \code{ggplot} object (invisibly).
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' r <- retro(NSH, NSH.tun, NSH.ctrl, retro=3)
#' retroResiduals(r, "HERAS", 2005:2011)
#' retroSelectivity(r, 2005:2010)
#' retroParams(r)
#' }
#' @rdname retroPlots
#' @export
retroResiduals <- function(x,fleet,yrs,cohort=F){

res <- lapply(x,residuals)
res <- lapply(res,function(y){y[which(y$fleet == fleet & y$year %in% yrs),]})
res <- lapply(res,function(y){cbind(y,retro=max(y$year))})
res <- do.call(rbind,res)
res$cohort <- res$year - res$age
res$retro  <- factor(res$retro)

if(length(unique(res$age))==1){
      pic <- ggplot(res,aes(x=year,y=std.res,colour=retro)) +
               labs(x="Years",title=paste("Residual pattern in",fleet,"at age"))
    } else {
      if(cohort){
        pic <- ggplot(res,aes(x=age,y=std.res,colour=retro)) +
                 facet_wrap(~cohort) +
                 labs(x="Ages",title=paste("Residual pattern in",fleet,"at cohort"))
      } else  {
        pic <- ggplot(res,aes(x=age,y=std.res,colour=retro)) +
                 facet_wrap(~year) +
                 labs(x="Ages",title=paste("Residual pattern in",fleet,"at age"))
      }
   }
pic <- pic +
         geom_hline(yintercept=0,linetype=3) +
         geom_line() +
         labs(y="Standardized residuals",colour="Assessment") +
         theme_bw()
print(pic)
return(invisible(pic))
}

#' @rdname retroPlots
#' @export
retroSelectivity <- function(x,yrs){

      res <- lapply(x,f)
      res <- lapply(res,function(y){y[which(y$year %in% (max(y$year)-20):(max(y$year)-1)),]})
      res <- lapply(res,function(y){cbind(y,retro=max(y$year))})
      res <- do.call(rbind,res)

      res <- subset(res,year %in% yrs)
      res$retro <- factor(res$retro)

      pic <- ggplot(res,aes(x=an(age),y=value,colour=retro)) +
               geom_line() +
               facet_wrap(~year,scales="free_y") +
               labs(x="Ages",y="F",colour="Assessment",title="Retrospective pattern in F at age") +
               theme_bw()
      print(pic)
      return(invisible(pic))
}

#' @rdname retroPlots
#' @export
retroParams <- function(x){
  retroPars     <- lapply(x,params)
  subretroPars  <- lapply(retroPars,function(y){return(subset(y,name %in% c("logFpar","logQpow","logSdLogFsta","logSdLogN","logSdLogObs",
                                                                            "logSdLogTotalObs","transfIRARdist","rec_pars","itrans_rho","logScale")))})
  subretroPars  <- lapply(as.list(1:length(subretroPars)),function(y){return(cbind(year=names(subretroPars)[y],subretroPars[[y]]))})
  subretroPars  <- lapply(subretroPars,function(y){
                          lapply(as.list(names(table(ac(y$name)))),
                                 function(z){tmpy <- subset(y,name==z)
                                             tmpy$nameOrig <- tmpy$name
                                             if(nrow(tmpy)>1)
                                              tmpy$name <- paste(tmpy$name,1:nrow(tmpy))
                                             return(tmpy)})})
  subretroPars  <- do.call(rbind,lapply(subretroPars,function(y){do.call(rbind,y)}))

  pic <- ggplot(subretroPars,aes(x=an(year),y=exp(value),colour=name,group=name)) +
           geom_line() +
           geom_point() +
           facet_wrap(~nameOrig,scales="free_y") +
           labs(x="Assessment year",y="Parameter value",colour=NULL) +
           theme_bw() +
           theme(legend.position="none")
  print(pic)
  return(invisible(pic))
  }

#' Kobe plot of an FLSAM assessment
#'
#' Plots the SSB/Bmsy versus F/Fmsy trajectory using the \pkg{kobe} package.
#'
#' @param object An \code{\link{FLSAM-class}} object.
#' @param fmsy,bmsy Reference points used to scale F and SSB.
#' @param nits Number of random draws.
#' @param prob Probability level.
#' @param ref.year Year to highlight (default the last assessment year).
#' @return A \code{ggplot} object.
#' @export
plotkobe <- function(object,fmsy=1,bmsy=1,nits=1000,prob=c(0.95),ref.year=NULL){
  if(is.null(ref.year)) ref.year <- object@range["maxyear"]
  require(kobe)

  yft <- kobeFLSAM(object,bmsy=bmsy,fmsy=fmsy,what="trks",prob=prob,nits=nits)
  pts <- kobeFLSAM(object,bmsy=bmsy,fmsy=fmsy,what="pts",prob=prob,nits=nits)
  kobePhase(subset(yft,quantity=="value"))+
    #geom_point(aes(stock,harvest),data=pts,lty=2,pch=19,col="blue",cex=0.5) +
    geom_path(aes(stock,harvest)) +
    geom_point(aes(stock,harvest),data=subset(yft,year==ref.year & quantity=="value"),pch=19,cex=2)
}
