#-------------------------------------------------------------------------------
# Diagnostics following the stockassessment package, drawn with ggplot2
#-------------------------------------------------------------------------------

#- Build a 'samres' object from the residuals stored in an FLSAM object
.samres <- function(sam){
  res <- sam@residuals
  if(nrow(res)==0 || all(is.na(res$std.res))) return(NULL)
  ret <- list(year=res$year,fleet=match(res$fleet,names(sam@control@fleets)),age=res$age,
              observation=res$log.obs,mean=res$log.mdl,residual=res$std.res)
  attr(ret,"fleetNames") <- names(sam@control@fleets)
  class(ret) <- "samres"
  return(ret)
}

#' Statistical tests on the residuals of an FLSAM assessment
#'
#' Tests the one-observation-ahead residuals per fleet and age, using
#' \code{\link[stockassessment]{residplot}}: Student's t-test for bias,
#' chi-square test for a variance of 1, Ljung-Box tests for autocorrelation in
#' the time and age direction, an F-test for equal variance of observations
#' below and above the median prediction (mean-variance relationship), and the
#' Shapiro-Wilk test for normality. The p-values are shown as a heat map.
#'
#' @param sam An \code{\link{FLSAM-class}} object with a stored fit (see \code{\link{samfit}}).
#' @param plot Print the plot (default \code{TRUE}).
#' @param p.breaks Breaks for the p-value classes.
#' @return A \code{data.frame} with columns \code{test}, \code{fleet}, \code{age}
#'   (\code{NA} for tests per fleet) and \code{p}, with the \code{ggplot} in
#'   the attribute \code{"plot"} (invisibly).
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)
#' residual.tests(sam)
#' }
#' @export
residual.tests <- function(sam,plot=TRUE,p.breaks=c(0,0.001,0.05,0.1,1)){
  fit     <- samfit(sam)
  pv      <- stockassessment::residplot(fit,resid=.samres(sam),plot=FALSE,p.breaks=p.breaks,
                                        palette=rep("white",length(p.breaks)-1))
  fleets  <- names(sam@control@fleets)
  ages    <- sam@control@range["min"]:sam@control@range["max"]
  byAge   <- c(bias="Bias",variance="Variance",correlation.time="Correlation time direction")
  byFleet <- c(correlation.age="Correlation age direction",meanvar="Mean-variance relationship",normality="Normality (Shapiro-Wilk)")
  res     <- rbind(do.call(rbind,lapply(names(byAge),function(i)
                     data.frame(test=byAge[i],fleet=rep(fleets,length(ages)),age=rep(ages,each=length(fleets)),p=as.vector(pv[[i]]),row.names=NULL))),
                   do.call(rbind,lapply(names(byFleet),function(i)
                     data.frame(test=byFleet[i],fleet=fleets,age=NA,p=as.vector(pv[[i]]),row.names=NULL))))
  res     <- res[!is.na(res$p),]
  rownames(res) <- NULL
  res$test  <- factor(res$test,levels=c(byAge,byFleet))
  res$fleet <- factor(res$fleet,levels=rev(fleets))
  res$class <- cut(res$p,breaks=p.breaks,include.lowest=TRUE)
  pdat      <- res
  pdat$age  <- ifelse(is.na(pdat$age),"all",as.character(pdat$age))
  pdat$age  <- factor(pdat$age,levels=c(ages,"all"))
  pic <- ggplot(pdat,aes(x=age,y=fleet,fill=class)) +
           geom_tile(colour="white") +
           geom_text(aes(label=sprintf("%.3f",p)),size=2.5) +
           scale_fill_manual(values=grDevices::colorRampPalette(c(.col("main"),.col("light"),"white"))(length(p.breaks)-1),drop=FALSE) +
           facet_wrap(~test,scales="free_x") +
           labs(x="Age",y=NULL,fill="p-value",title=paste("Residual tests",sam@name)) +
           theme_bw()
  if(plot) print(pic)
  res$class <- NULL
  attr(res,"plot") <- pic
  return(invisible(res))
}

#' Jitter test for an FLSAM assessment
#'
#' Refits the assessment from randomly perturbed starting values (as
#' \code{\link[stockassessment]{jit}}) to check that the optimum is unique and
#' the model is identifiable. Parameters fixed through \code{set.pars} stay
#' fixed.
#'
#' @param sam An \code{\link{FLSAM-class}} object with a stored fit (see \code{\link{samfit}}).
#' @param nojit Number of jittered refits.
#' @param sd Standard deviation of the normal noise added to the default
#'   starting values (\code{stockassessment::defpar}).
#' @param ncores Number of cores; 1 runs sequentially.
#' @param seed Optional random seed.
#' @return An \code{\link{FLSAMs}} object with the original fit (\code{"base"})
#'   and the converged jittered fits. The attribute \code{"table"} holds per
#'   run the negative log-likelihood, its difference to the base fit, the
#'   maximum absolute difference in parameter estimates and the convergence flag.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)
#' jt  <- jitter.test(sam, nojit=5, seed=1)
#' attr(jt, "table")
#' plot(jt)
#' }
#' @export
jitter.test <- function(sam,nojit=10,sd=0.25,ncores=1,seed=NULL){
  fit     <- samfit(sam)
  if(!is.null(seed)) set.seed(seed)
  par     <- defpar(fit$data,fit$conf)
  mp      <- fit$obj$env$map
  userMap <- mp[setdiff(names(mp),c("logFScaleMSY","implicitFunctionDelta","splinePenalty"))]
  parv    <- unlist(par)
  pars    <- lapply(seq_len(nojit),function(i){
               p <- relist(parv+rnorm(length(parv),sd=sd),par)
               for(nm in names(userMap)){
                 fx <- is.na(userMap[[nm]])
                 p[[nm]][fx] <- fit$pl[[nm]][fx]
               }
               p})
  runFit  <- function(p) try(sam.fit(fit$data,fit$conf,p,map=userMap,silent=TRUE),silent=TRUE)
  if(ncores>1){
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl,suppressMessages(library(stockassessment)))
    fits <- parallel::parLapply(cl,pars,runFit)
  } else {
    fits <- lapply(pars,runFit)
  }
  ok      <- sapply(fits,function(x) !is(x,"try-error"))
  tab     <- data.frame(run=seq_len(nojit),
                        nlogl=sapply(fits,function(x) if(is(x,"try-error")) NA else x$opt$objective),
                        max.par.diff=sapply(fits,function(x) if(is(x,"try-error")) NA else max(abs(x$opt$par-fit$opt$par))),
                        converged=sapply(fits,function(x) if(is(x,"try-error")) FALSE else x$opt$convergence==0))
  tab$dnlogl <- tab$nlogl - fit$opt$objective
  tab     <- tab[,c("run","nlogl","dnlogl","max.par.diff","converged")]
  ctrl    <- sam@control; ctrl@residuals <- FALSE
  res     <- c(list(base=sam),lapply(fits[ok],SAM2FLR,ctrl=ctrl,keep.fit=FALSE))
  names(res) <- c("base",paste("jit",which(ok)))
  for(i in names(res)) res[[i]]@name <- i
  res     <- FLSAMs(res)
  attr(res,"table") <- tab
  return(res)
}

#' Empirical correlation of residuals between ages
#'
#' For each fleet with more than one age, the correlation between the
#' residuals of the different ages across years (as
#' \code{\link[stockassessment]{empirobscorrplot}}). Strong correlations
#' suggest that an \code{"AR"} or \code{"US"} observation correlation
#' structure (\code{cor.obs.Flag}) should be considered.
#'
#' @param sam An \code{\link{FLSAM-class}} object.
#' @return A \code{ggplot} object; the data are in its \code{data} element.
#' @examples
#' data(NSH); data(NSH.sam)
#' obscor.empirical.plot(NSH.sam)
#' @export
obscor.empirical.plot <- function(sam){
  res   <- subset(sam@residuals,!is.na(std.res))
  dat   <- do.call(rbind,lapply(split(res,res$fleet),function(x){
             if(length(unique(x$age))<2) return(NULL)
             tmp <- xtabs(std.res ~ age + year,data=x)
             cm  <- cor(t(tmp))
             data.frame(fleet=x$fleet[1],age1=as.numeric(rownames(cm))[row(cm)],
                        age2=as.numeric(colnames(cm))[col(cm)],cor=as.vector(cm))}))
  if(is.null(dat)) stop("No fleets with more than one age with residuals")
  ggplot(dat,aes(x=factor(age1),y=factor(age2),fill=cor)) +
    geom_tile() +
    geom_text(aes(label=sprintf("%.2f",cor)),size=2.5) +
    scale_fill_gradient2(low=.col("accent"),mid="white",high=.col("main"),midpoint=0,limits=c(-1,1)) +
    facet_wrap(~fleet,scales="free") +
    labs(x="Age",y="Age",fill="Correlation",title="Empirical correlation of residuals") +
    theme_bw()
}

#' Observation standard deviations
#'
#' Bar plot of the estimated observation standard deviations by fleet and age,
#' ordered from small to large (as \code{\link[stockassessment]{sdplot}}),
#' or of the relative weights (\code{1/sd^2} scaled to a maximum of 1).
#'
#' @param sam An \code{\link{FLSAM-class}} object.
#' @param show.rel.w Plot relative weights instead of standard deviations.
#' @return A \code{ggplot} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' obssd.plot(NSH.sam)
#' obssd.plot(NSH.sam, show.rel.w=TRUE)
#' @export
obssd.plot <- function(sam,show.rel.w=FALSE){
  obv        <- obs.var(sam)
  obv$name   <- paste(obv$fleet,ifelse(is.na(obv$age),"",paste("age",obv$age)))
  obv$rel.w  <- (1/obv$value^2)/max(1/obv$value^2)
  obv        <- obv[order(obv$value),]
  obv$name   <- factor(obv$name,levels=unique(obv$name))
  obv$y      <- if(show.rel.w) obv$rel.w else obv$value
  ggplot(obv,aes(x=name,y=y,fill=fleet)) +
    geom_col() +
    scale_fill_manual(values=FLSAM.colours(length(unique(obv$fleet)))) +
    labs(x=NULL,y=ifelse(show.rel.w,"Relative weight","SD"),fill=NULL,title=sam@name) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
}

#' Prediction-variance relationship
#'
#' Plots the estimated relationship between the predicted observation and its
#' standard deviation for a fleet with a prediction-variance link
#' (\code{predVarObsLink}), with a 95\% confidence band from 20000 draws of the
#' fixed-effect parameters (as \code{\link[stockassessment]{predstdplot}}).
#'
#' @param sam An \code{\link{FLSAM-class}} object with a stored fit (see \code{\link{samfit}}).
#' @param fleet Fleet name or number.
#' @param age Age(s), needed when the fleet has more than one relationship.
#' @param type \code{"log"} (sd of log observations) or \code{"natural"}.
#' @param N Number of parameter draws for the confidence band.
#' @return A \code{ggplot} object.
#' @export
predsd.plot <- function(sam,fleet,age=NULL,type="log",N=20000){
  fit   <- samfit(sam)
  if(is.character(fleet)) fleet <- match(fleet,names(sam@control@fleets))
  link  <- fit$conf$predVarObsLink[fleet,]
  use   <- link>=0 & !is.na(link)
  b     <- unique(link[use])
  a     <- unique(fit$conf$keyVarObs[fleet,use])
  if(length(b)==0) stop("The prediction-variance relation is not used for this fleet.")
  if(length(b)>1 | length(a)>1){
    if(is.null(age)) stop("Multiple relations in the fleet; need age to plot the relation for.")
    ageIndex <- age - fit$conf$minAge+1
    if(length(unique(link[ageIndex]))!=1 | length(unique(fit$conf$keyVarObs[fleet,ageIndex]))!=1)
      stop("Several or no relations in the provided age span; only one relation is needed.")
  } else {
    age <- which(link==b) + fit$conf$minAge-1
  }
  ageIndex <- age[1] - fit$conf$minAge+1
  sel      <- fit$data$aux[,"fleet"]==fleet & fit$data$aux[,"age"] %in% age
  predObs  <- fit$rep$predObs[sel]
  beta     <- exp(fit$pl$predVarObs[link[ageIndex]+1])+1
  alpha    <- exp(fit$pl$logSdLogObs[fit$conf$keyVarObs[fleet,ageIndex]+1])
  pred     <- exp(seq(min(predObs),max(predObs),length.out=30))
  ss       <- stockassessment::rmvnorm(N,mu=fit$opt$par,Sigma=fit$sdrep$cov.fixed)
  alphaSim <- exp(ss[,which(names(fit$opt$par)=="logSdLogObs"),drop=FALSE])[,fit$conf$keyVarObs[fleet,ageIndex]+1]
  betaSim  <- exp(ss[,which(names(fit$opt$par)=="predVarObs"),drop=FALSE])[,link[ageIndex]+1]+1
  sdFun    <- if(type=="log") function(al,be,p) sqrt(log(al*p^(be-2)+1)) else function(al,be,p) sqrt(al*p^be)
  sigSim   <- sapply(pred,function(p) sdFun(alphaSim,betaSim,p))
  qi       <- apply(sigSim,2,quantile,c(0.025,0.975))
  xTrans   <- if(type=="log") log else identity
  curve    <- data.frame(x=xTrans(pred),sd=sdFun(alpha,beta,pred),lbnd=qi[1,],ubnd=qi[2,])
  pts      <- data.frame(x=xTrans(exp(predObs)),sd=sdFun(alpha,beta,exp(predObs)))
  ggplot(curve,aes(x=x)) +
    geom_ribbon(aes(ymin=lbnd,ymax=ubnd),fill=.col("light"),alpha=0.4) +
    geom_line(aes(y=sd),colour=.col("main")) +
    geom_point(aes(y=sd),data=pts,shape=3,colour=.col("accent")) +
    labs(x=ifelse(type=="log","Prediction (log)","Prediction"),y="Standard deviation",
         title=paste0("Pred-std relation: ",names(sam@control@fleets)[fleet]," age ",paste(range(age),collapse="-"))) +
    theme_bw()
}

#' Data overview
#'
#' Shows which fleets have observations in which years, and for how many ages
#' (as \code{\link[stockassessment]{dataplot}}).
#'
#' @param sam An \code{\link{FLSAM-class}} object.
#' @return A \code{ggplot} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' data.plot(NSH.sam)
#' @export
data.plot <- function(sam){
  if(.hasSlot(sam,"fit") && !is.null(sam@fit)){
    dat <- data.frame(year=sam@fit$data$aux[,"year"],fleet=names(sam@control@fleets)[sam@fit$data$aux[,"fleet"]],
                      obs=!is.na(sam@fit$data$logobs))
  } else {
    dat <- data.frame(year=sam@residuals$year,fleet=sam@residuals$fleet,obs=!is.na(sam@residuals$std.res))
  }
  dat   <- aggregate(obs ~ year + fleet,data=dat,FUN=sum)
  dat   <- subset(dat,obs>0)
  dat$fleet <- factor(dat$fleet,levels=rev(names(sam@control@fleets)))
  ggplot(dat,aes(x=year,y=fleet,fill=obs)) +
    geom_tile(colour="white") +
    scale_fill_gradient(low=.col("light"),high=.col("main")) +
    labs(x="Year",y=NULL,fill="No. ages\nobserved",title=sam@name) +
    theme_bw()
}

#' Selectivity plot
#'
#' Fishing mortality at age scaled to the maximum F in each year, by year and
#' catch fleet (comparable to \code{\link[stockassessment]{fselectivityplot}}).
#'
#' @param sam An \code{\link{FLSAM-class}} object.
#' @return A \code{ggplot} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' fselectivity.plot(NSH.sam)
#' @export
fselectivity.plot <- function(sam){
  fs      <- f(sam)
  fs$sel  <- fs$value/ave(fs$value,fs$fleet,fs$year,FUN=max)
  ggplot(fs,aes(x=age,y=sel,group=year,colour=year)) +
    geom_line() +
    facet_wrap(~fleet) +
    scale_colour_gradient(low=.col("light"),high=.col("main")) +
    labs(x="Age",y="F / max(F)",colour="Year",title=sam@name) +
    theme_bw()
}
