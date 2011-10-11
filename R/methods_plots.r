setMethod("plot",signature(x="FLSAM",y="missing"),
 function(x,y,...) {
    #Extract data - ssb and fbar are easy but recs harder
    ssb.fbar.dat <- subset(x@params,name %in% c("logssb","logfbar"))
    states <- subset(x@params,name=="U")
    yrs    <- x@range["minyear"]:x@range["maxyear"]
    n.states <- nrow(states)/length(yrs)
    recs.dat    <- states[seq(start=1,by=n.states,length.out=length(yrs)),]
    plot.dat    <- rbind(ssb.fbar.dat,recs.dat) 
    plot.dat$year <- yrs

    #Lets roll
    plot.dat$ubnd <- exp(plot.dat$value + 1.96*plot.dat$std.dev)
    plot.dat$lbnd <- exp(plot.dat$value - 1.96*plot.dat$std.dev)
    plot.dat$est  <- exp(plot.dat$value)
    plot.dat$name <- factor(plot.dat$name,levels=c("logssb","logfbar","U"))
    levels(plot.dat$name) <- c("Spawning stock biomass","Fishing mortality","Recruitment")

    #Plot it
    xyplot(est+ubnd+lbnd~year|name,data=plot.dat,
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
                   panel.xyplot(d$est$x,d$est$y,col="black",lwd=2,type="l")
                     },
              scales=list(alternating=1,y=list(relation="free",rot=0)))
})
