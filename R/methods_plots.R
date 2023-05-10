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
    plot.dat$name <- factor(plot.dat$name,
                       levels=c("Spawning stock biomass","Fishing mortality","Recruitment"))
    if(!futureYrs) plot.dat    <- subset(plot.dat,year %in% range(x)["minyear"]:(range(x)["maxyear"]-1))
    #Plot it
    xyplot(value+ubnd+lbnd~year|name,data=plot.dat,#...,
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

#Components plot
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

fpart.plot <- function(sam){
  f.dat <- f(sam)
  fbars <- sam@control@range["minfbar"]:sam@control@range["maxfbar"]
  f.dat <- subset(f.dat,age %in% fbars)
  f.dat <- aggregate(f.dat[,"value"],by=as.list(f.dat[,c("fleet","year")]),FUN=mean,na.rm=T)
  xyplot(x ~ year,data=f.dat,groups=fleet,type="b",pch=19,
         main=list(sam@name,cex=0.9),
         ylab="Fbar-partial",xlab="Year",
         panel=function(...) {
          panel.grid(h=-1,v=-1)
          panel.xyplot(...)
         },
         scales=list(alternating=1,y=list(relation="free",rot=0)))
}

#Correlation plot
cor.plot <- function(sam,cols=c("red","white","blue"),full=FALSE) {
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
  setGeneric("otolith", function(object, ...)
    standardGeneric("otolith"))
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
    paramvalue<- subset(object@params,name%in%c("beforeLastLogF","beforeLastLogN","lastLogF","lastLogN","logCatch",
                                                 "logCatchByFleet","logfbar","logR","logssb","logtsb","comps"))
    for(i in 1:ceiling(n/200)){
      d <- mvrnorm(n=200,paramvalue$value, object@rescov)
      colnames(d) <- paramvalue$name
      Fbar.all  <- rbind(Fbar.all,d[,colnames(d)=="logfbar"])
      SSB.all   <- rbind(SSB.all,d[,colnames(d)=="logssb"])
    }
    Fbar.all  <- Fbar.all[1:n,]
    SSB.all   <- SSB.all[ 1:n,]
  } else {
      pars <- subset(object@params,name%in%c("beforeLastLogF","beforeLastLogN","lastLogF","lastLogN","logCatch",
                                                 "logCatchByFleet","logfbar","logR","logssb","logtsb","comps"))
      browser()
      d <- mvrnorm(n=n,c(pars$value), object@rescov)
      colnames(d) <- pars$name
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
obscv.plot  <- function(sam){

  obv <- obs.var(sam)
  obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
  obv <- obv[order(obv$value),]
  plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
    pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
  text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
}

#-Observation correlation plot
obscor.plot <- function(sam){
  require(ellipse)
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
        inpQ  <- as(inp,"FLQuant")
        x_at <- pretty(1:(max(as.data.frame(inpQ)$year)-min(as.data.frame(inpQ)$year))+1)-1
        x_labels <- formatC(pretty(as.data.frame(inpQ)$year), digits = 0, format = "f")
        if(rel)
          yTxt <- "% difference numbers at age"
        if(rel==F)
          yTxt <- "Numbers at age"
        pic   <- barchart(data ~ year | as.factor(age),data=as.data.frame(inpQ),origin = 0, scales = list(x=list(at=x_at,labels=x_labels),y="free"),col=1,horizontal=F,
                          xlab="Year",ylab=yTxt,main="Process error deviation in N",
                          panel=function(...) {
                            panel.grid(h=-1,v=0);
                            panel.barchart(...)})
        print(pic)
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
        inpQ  <- as(inp,"FLQuant")
        if(rel)
          yTxt <- "% difference mortality (year)"
        if(rel==F)
          yTxt <- "Mortality (year)"
        pic   <- bubbles(age ~ year,data=inpQ,bub.scale=5,col=c("black","red"),xlab="Years",ylab=yTxt,main="Process error deviation in M")
        print(pic)
      }
      if(type=="tsb"){
        x_at <- pretty(1:(max(bdiff$year)-min(bdiff$year)+1))-1
        x_labels <- formatC(pretty(bdiff$year), digits = 0, format = "f")
        
        if(rel)
          bdiff$bdiff <- bdiff$bdiff/bnorm$bnorm *100
        if(rel)
          yTxt <- "% difference in biomass (t)"
        if(rel==F)
          yTxt <- "Biomass (t)"
        pic <- barchart(bdiff ~ year,data=bdiff,origin = 0, col=1,horizontal=F,scales = list(x=list(at=x_at,labels=x_labels)),
                          xlab="Year",ylab=yTxt,main="Process error deviation in biomass",
                          panel=function(...) {
                            panel.grid(h=-1,v=0);
                            panel.barchart(...)})
        print(pic)
      }
}
      
#- Plot the retrospective residual pattern in any fleet
retroResiduals <- function(x,fleet,yrs){

res <- lapply(x,residuals)
res <- lapply(res,function(y){y[which(y$fleet == fleet & y$year %in% yrs),]})
res <- lapply(res,function(y){cbind(y,retro=max(y$year))})
res <- do.call(rbind,res)

print(xyplot(std.res ~ age | as.factor(year),data=res,type="l",groups=retro,
auto.key=list(space="right",points=FALSE,lines=TRUE,type="l"),main=paste("Residual pattern in",fleet,"at age"),
ylab="Standardized residuals",xlab="Ages",
panel = panel.superpose,
 panel.groups = function(...) {
    panel.grid(v=-1,h=-1,lty=3)
    panel.xyplot(...)
},
scales=list(alternating=1,y=list(relation="free",rot=0))))
}
#- Plot the retrospective selectivity pattern
retroSelectivity <- function(x,yrs){

      res <- lapply(x,f)
      res <- lapply(res,function(y){y[which(y$year %in% (max(y$year)-20):(max(y$year)-1)),]})
      res <- lapply(res,function(y){cbind(y,retro=max(y$year))})
      res <- do.call(rbind,res)

      res <- subset(res,year %in% yrs)

      print(xyplot(value ~ an(age) | as.factor(year),data=res,type="l",groups=retro,
      auto.key=list(space="right",points=FALSE,lines=TRUE,type="l"),main=paste("Retrospective pattern in F at age"),
      ylab="F",xlab="Ages",
      panel = panel.superpose,
       panel.groups = function(...) {
          panel.grid(v=-1,h=-1,lty=3)
          panel.xyplot(...)
      },
      scales=list(alternating=1,y=list(relation="free",rot=0))))
}

retroParams <- function(x){
  retroPars     <- lapply(x,params)
  subretroPars  <- lapply(retroPars,function(y){return(subset(y,name %in% c("logFpar","lowQpow","logSdLogFsta","logSdLogN","logSdLogObs","rec_loga",
                                                                            "rec_logb","rho","logScale","logScaleSSB","logPowSSB","logSdSSB")))})
  subretroPars  <- lapply(as.list(1:length(subretroPars)),function(y){return(cbind(year=names(subretroPars)[y],subretroPars[[y]]))})
  subretroPars  <- lapply(subretroPars,function(y){
                          lapply(as.list(names(table(ac(y$name)))),
                                 function(z){tmpy <- subset(y,name==z)
                                             tmpy$nameOrig <- tmpy$name
                                             if(nrow(tmpy)>1)
                                              tmpy$name <- paste(tmpy$name,1:nrow(tmpy))
                                             return(tmpy)})})
  subretroPars  <- do.call(rbind,lapply(subretroPars,function(y){do.call(rbind,y)}))

  print(xyplot(exp(value) ~ an(year) | as.factor(nameOrig),data=subretroPars,groups=name,scale=list(y="free"),type="b",pch=19,xlab="Assessment year",ylab="Parameter value"))
  }

#kobeplot
plotkobe <- function(object,fmsy=1,bmsy=1,nits=1000,prob=c(0.95),ref.year=NULL){
  if(is.null(ref.year)) ref.year <- dims(NSH.sam)$maxyear
  require(kobe)

  yft <- kobeFLSAM(object,bmsy=bmsy,fmsy=fmsy,what="trks",prob=prob,nits=nits)
  pts <- kobeFLSAM(object,bmsy=bmsy,fmsy=fmsy,what="pts",prob=prob,nits=nits)
  kobePhase(subset(yft,quantity=="value"))+
    #geom_point(aes(stock,harvest),data=pts,lty=2,pch=19,col="blue",cex=0.5) +
    geom_path(aes(stock,harvest)) +
    geom_point(aes(stock,harvest),data=subset(yft,year==ref.year & quantity=="value"),pch=19,cex=2)
}

      
        
