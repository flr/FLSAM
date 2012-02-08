  #######################################################################
  ## Need to do:
  ##
  ## Add smoother to plot d
  ## Find nicer way to do x axis in plot a,c,e 
  ## 
  #######################################################################


residual.diagnostics <- function(x,title=x@name) {
# extracts residuals dataframe from x
  index.res <- x@residuals


# Back transform log transformed observed and modelled to normal space
  index.res$obs <- exp(index.res$log.obs)
  index.res$mdl <- exp(index.res$log.mdl)


# split dataframe by age and survey (fleet)
  index.res.l   <- split(index.res,list(index.res$age,index.res$fleet),drop=TRUE)


#Setup plots
  oldpar <-  par(mfrow=c(3,2),las=0,oma=c(0,0,3,0),mgp=c(1.75,0.5,0),
                mar=c(3,3,2.5,1),cex.main=1,tck=-0.01)

# Run through each combination of survey and age                                                                         

  for (i in 1:length(index.res.l)){

# create working titel to identify age and survey

  ind.age   <- unlist(strsplit(names(index.res.l[i]), "\\."))
  ttl.fmt   <- ifelse(x@control@fleets[ind.age[2]]%in%c(3,4),
                 "Diagnostics - %s",   #SSB fleey
                 "Diagnostics - %s, age %s")  #other fleet
  ttl       <- sprintf(ttl.fmt,ind.age[2],ind.age[1])
  ttl       <- paste(title,ttl)

#Scale index axes
  idx.rng       <-  range(c(index.res.l[[i]]$obs,index.res.l[[i]]$mdl),na.rm=TRUE)
  idx.lim       <-  range(c(0,idx.rng))
  idx.exp       <-  floor(log10(max(pretty(idx.lim)))/3)*3
  idx.div       <-  10^idx.exp
  idx.label     <-  ifelse(idx.exp>1,paste("values ","[","10^",idx.exp,"]",
                    sep=""),paste("values "))

  index.res.l[[i]]$obs <-  index.res.l[[i]]$obs/idx.div
  index.res.l[[i]]$mdl <-  index.res.l[[i]]$mdl/idx.div
  idx.rng       <-  idx.rng/idx.div
  idx.lim       <-  idx.lim/idx.div
  idx.min       <- min(c(index.res.l[[i]]$obs,index.res.l[[i]]$mdl))
  idx.max       <- max(c(index.res.l[[i]]$obs,index.res.l[[i]]$mdl))
  idx.lim       <- c(idx.min/1.15,ceiling(idx.max/0.12)) #scaling keeps data points clear of legend and x-axis


# scale for residuals axes
  res.range     <- abs(range(index.res.l[[i]]$std.res))
  res.lim       <- c(-max(res.range),max(res.range))


#plot 1 obs and mdl time series plotted on a log scale

  plot(obs ~ year, index.res.l[[i]], log="y",ylim=idx.lim,xlab="Year", 
        ylab=idx.label,pch=16)
  points(mdl ~ year, index.res.l[[i]],pch=4)
  points(mdl ~ year, index.res.l[[i]],type="l")
  legend("topleft",c("Observed","Fitted"),pch=c(16,4),lty=c(NA,1),horiz=TRUE)
  title("a) Observed and fitted values time series")


# plot 2 observed against modelled directly, both axes log scale.

  plot(obs~mdl,index.res.l[[i]],log="xy",ylim=idx.lim,xlim=idx.lim,pch=16,
          ylab=paste("Observed ",idx.label,sep=""),
          xlab=paste("Fitted ",idx.label,sep=""))
  abline(0,1,col="black")                    
  legend("topleft","1:1 line",lty=1,horiz=TRUE)
  title("b) Observed vs fitted values")


# Plot 3 Index residuals over time

  plot(std.res ~year, index.res.l[[i]],ylim=res.lim, ylab="Standardised Residuals", xlab="Year")
          points(std.res ~year, index.res.l[[i]], type="h")
          points(std.res ~year, index.res.l[[i]], pch=19,cex=0.75)
  abline(h=0)  
  title("c) Standardised residuals over time")


# Plot 4 Tukey-Anscombe plot, Standardised residuals vs Fitted values, log x-axis

  plot(std.res ~mdl, index.res.l[[i]], ylab="Standardised Residuals", 
        ylim=res.lim, xlim=idx.lim,pch=19,cex=0.75,
        xlab=paste("Fitted ",idx.label,sep=""),log="x")
  plt.dat <- index.res.l[[i]]
  abline(h=0)  
  smoother <- loess.smooth(plt.dat$mdl,plt.dat$std.res)
  lines(smoother,col="red")
  title("d) Tukey-Anscombe plot")

#plot 5 Normal Q-Q plot

  qqnorm(index.res.l[[i]]$std.res,ylim=res.lim,xlab="Quantiles of the Normal Distribution",
        ylab="Standardised Residuals",pch=19,main="")
  qqline(index.res.l[[i]]$std.res,col="red")
  abline(0,1,lty=2)
  legend("topleft",c("qqline","1:1 line"),lty=c(1,2),col=c("red","black"),horiz=TRUE)
  title("e) Normal Q-Q plot")

#Plot 6 Autocorrelation function plot
  
  acf(as.ts(index.res.l[[i]]$std.res),ylab="ACF",xlab="Lag (yrs)",type=c("partial"),
        ci.col="black",main="")
  legend("topright",legend=c("95% Conf. Int."),lty=c(2),pch=c(NA),horiz=TRUE,box.lty=0)
  
  title("f) Autocorrelation of Residuals")  


  # Add a main titel to the plots with Survey and age information
  title(main=ttl,outer=TRUE)


  }

### ============================================================================
### Finishing up
### ============================================================================

  par(oldpar)
}


#Now run the code
#residual.diagnostics(NSH.sam)
