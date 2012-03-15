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
    rec.dat$age <- NULL
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
    rec.dat$age <- NULL
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
cor.plot <- function(sam,cols=c("#D7191C","#FDAE61","#FFFFBF","#ABDDA4","#2B83BA"),
                      full=FALSE) {
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
