setMethod("plot",signature(x="FLSAM",y="missing"),
 function(x,y,futureYrs=T,...) {
    #Extract data - ssb and fbar are easy but recs harder
    ssb.dat <- ssb(x)
    ssb.dat$name <- "Spawning stock biomass"
    rec.dat <- rec(x)
    rec.dat$name <- "Recruitment"
    rec.dat$age <- NULL
    fbar.dat <- fbar(x)
    fbar.dat$name <- "Fishing mortality"
    plot.dat    <- rbind(ssb.dat, fbar.dat,rec.dat) 
    plot.dat$name <- factor(plot.dat$name,
                       levels=c("Spawning stock biomass","Fishing mortality","Recruitment"))
    if(!futureYrs) plot.dat    <- subset(plot.dat,year %in% range(x)["minyear"]:(range(x)["maxyear"]-1))
    #Plot it
    xyplot(value+ubnd+lbnd~year|name,data=plot.dat,...,
              prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
              main=list(x@name,cex=0.9),
              ylab=rev(c("SSB","Fbar","Rec")),xlab="Year",
              layout=c(1,3),as.table=TRUE,
              panel=function(...) {
                   panel.grid(h=-1,v=-1)
                   arg <- list(...)
                   d<- split(data.frame(x=arg$x,y=arg$y),arg$groups[arg$subscripts])
                   panel.polygon(c(d$ubnd$x,rev(d$lbnd$x)),c(d$ubnd$y,rev(d$lbnd$y)),
                      col="grey",border=NA)
                   panel.xyplot(d$value$x,d$value$y,col="black",lwd=2,type="l")
                     },
              scales=list(alternating=1,y=list(relation="free",rot=0)))
})

setMethod("plot",signature(x="FLSAM",y="FLSAM"),
 function(x,y,main="",...) {
    plot(FLSAMs(x,y),...)
})

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
    plot.dat$var <- factor(plot.dat$var,
                       levels=c("Spawning stock biomass","Fishing mortality","Recruitment"))

    #Plot it
    xyplot(value~year|var,data=plot.dat,...,
              prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
              main=main,groups=name,type="l",
              ylab=rev(c("SSB","Fbar","Rec")),xlab="Year",
              layout=c(1,3),as.table=TRUE,
              auto.key=list(space="right",lines=TRUE,points=FALSE),
              scales=list(alternating=1,y=list(relation="free",rot=0)),
              panel=function(...) {
                   panel.grid(h=-1,v=-1)
                   panel.superpose.2(...)})
})

setMethod("plot",signature(x="FLSAMs",y="FLStock"),
 function(x,y,main="",futureYrs,...) {
    #extract data from FLSAMs first
    ssb.dat <- ssb(x)
    if(!futureYrs) ssb.dat <- ssb(x)[-c(which(diff(ssb(x)$year)<0),nrow(ssb(x))),]
    ssb.dat$var <- "spawning stock biomass"
    rec.dat <- rec(x)
    if(!futureYrs) rec.dat <- rec(x)[-c(which(diff(rec(x)$year)<0),nrow(rec(x))),]
    rec.dat$var <- "recruitment"
    fbar.dat <- fbar(x)
    if(!futureYrs) fbar.dat <- fbar(x)[-c(which(diff(fbar(x)$year)<0),nrow(fbar(x))),]
    fbar.dat$var <- "fishing mortality"
    sams.dat    <- rbind(ssb.dat, fbar.dat,rec.dat) 

    #Then add the data from the stock
    stck.dat <- rbind(data.frame(as.data.frame(rec(y)),var="recruitment"),
                  data.frame(as.data.frame(ssb(y)),var="spawning stock biomass"),
                  data.frame(as.data.frame(fbar(y)),var="fishing mortality"))
    stck.dat$name <- y@name
    stck.dat$value <- stck.dat$data
    
    #Now combine
    common.cols <- intersect(names(stck.dat),names(sams.dat))
    plot.dat <- rbind(stck.dat[,common.cols],sams.dat[,common.cols])
    plot.dat$var <- factor(plot.dat$var,
                       levels=c("spawning stock biomass","fishing mortality","recruitment"))

    if(!futureYrs) plot.dat    <- subset(plot.dat,year %in% min(unlist(lapply(x,function(y){range(y)["minyear"]}))) :
                                                            max(unlist(lapply(x,function(y){range(y)["maxyear"]-1}))))
    #plot it
    xyplot(value~year|var,data=plot.dat,...,
              prepanel=function(...) {list(ylim=range(pretty(c(0,list(...)$y))))},
              main=main,groups=name,type="l",
              ylab=rev(c("ssb","fbar","rec")),xlab="year",
              layout=c(1,3),as.table=TRUE,
              auto.key=list(space="right",lines=TRUE,points=TRUE),
              scales=list(alternating=1,y=list(relation="free",rot=0)),
              panel=function(...) {
                   panel.grid(h=-1,v=-1)
                   panel.superpose.2(...)})
})

setMethod("plot",signature(x="FLSAM",y="FLStock"),
 function(x,y,main="",...) {
    plot(FLSAMs(x),y,...)
})

setMethod("plot",signature(x="FLStock",y="FLSAM"),
 function(x,y,main="",...) {
    plot(y,x,...)
})

setMethod("plot",signature(x="FLStock",y="FLSAMs"),
 function(x,y,main="",...) {
    plot(y,x,...)
})

#Correlation plot
cor.plot <- function(sam,cols=c("red","white","blue"),full=FALSE) {
   if(sam@control@nohess) {
    stop(paste("Cannot generate a correlation plot for an FLSAM object that has been run",
               "with the nohess=TRUE option. Please rerun the model with nohess=FALSE and",
               "try again."))}
   n.coefs <- nrow(coef(sam))   
   cor.mat <- cov2cor(sam@vcov)
   if(!full) {cor.mat <- cor.mat[1:n.coefs,1:n.coefs]} #Don't show the full matrix 
   var.names <- colnames(cor.mat)
   var.names <- factor(var.names,unique(var.names))
   var.names <- paste(var.names,do.call(c,lapply(table(var.names),seq,from=1)),sep=".")
   dimnames(cor.mat) <- list(Var2=var.names,Var1=var.names)
   plt.dat <- as.data.frame(as.table(cor.mat),responseName="cor")
   levelplot(cor ~ Var2 * Var1,plt.dat,
       at=seq(-1,1,by=0.02),zlim=c(-1,1),
       scales=list(x=list(rot=90)),
       xlab="",ylab="",main=sam@name,
       col.regions=colorRampPalette(cols,space="Lab")(102))
}

#Otolith plot (resample from the variance co-variance matrix and create a plot of all the resampled value)
#- Create generic function for 'looi'
if (!isGeneric("otolith")) {
  setGeneric("otolith", function(sam.out,year, ...)
    standardGeneric("otolith"))
}
setMethod("otolith", signature(sam.out="FLSAM", year="numeric"),
function(sam.out, year=2011, plot=TRUE, show.points=FALSE, do.contours=TRUE,
                      margin.plots=TRUE, show.estimate=TRUE,
                      n=100, pch=".", alpha=0.05, show.grid=TRUE,
                      n.grid=50, contour.args=list(),...){

  require(MASS)
  debug <- FALSE
  filled.contours <- FALSE

  #Get the range of ages for fbar
  f.ages        <-  sam.out@control@range["minfbar"]:sam.out@control@range["maxfbar"]

  #set up year ranges and identify year which plot is being done for
  years <- sam.out@control@range["minyear"]:sam.out@control@range["maxyear"]
  year.index <- which(years==year)


  # Generate n iterations of all the parameters
  if(n > 200){
    Fbar.all  <- numeric()
    SSB.all   <- numeric()
    paramvalue<- sam.out@params$value
    samvcov   <- sam.out@vcov
    for(i in 1:ceiling(n/200)){
      d <- mvrnorm(n=200,paramvalue, samvcov)
      Fbar.all  <- rbind(Fbar.all,d[,colnames(d)=="logfbar"])
      SSB.all   <- rbind(SSB.all,d[,colnames(d)=="logssb"])
    }
    Fbar.all  <- Fbar.all[1:n,]
    SSB.all   <- SSB.all[ 1:n,]
  } else {
      d <- mvrnorm(n=n,sam.out@params$value, sam.out@vcov)
      Fbar.all <- d[,colnames(d)=="logfbar"]
      SSB.all <-  d[,colnames(d)=="logssb"]
    }

  #Set these points up to be returned and return them to normal space
  Fbar <- exp(Fbar.all[,year.index])
  SSB <- exp(SSB.all[,year.index])
  return.obj  <- data.frame(Fbar,SSB)

  #Extract SSB and Fbar estimates from the assessment outputs for the year requested
  SSB.est   <-  ssb(sam.out)$value[ssb(sam.out)$year==year]
  Fbar.est  <-  fbar(sam.out)$value[fbar(sam.out)$year==year]


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
obsvar.plot <- function(FLSAM){
                  obv <- obs.var(FLSAM)
                  obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
                  obv <- obv[order(obv$value),]
                  bp <- barplot(obv$value,ylab="Observation Variance",
                         main="Observation variances by data source",col=factor(obv$fleet))
                  axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
                  legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)
              }

#-CV versus observation variance plot
obscv.plot  <- function(FLSAM){

  obv <- obs.var(FLSAM)
  obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
  obv <- obv[order(obv$value),]
  plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
    pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
  text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
}
