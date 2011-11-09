setMethod("plot",signature(x="FLSAM",y="missing"),
 function(x,y,...) {
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

    #Plot it
    xyplot(value+ubnd+lbnd~year|name,data=plot.dat,
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
