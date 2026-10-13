#-------------------------------------------------------------------------------
# Parameters with readable names, parameter comparison plot and the analytical
# SSB-Fbar confidence ellipse
#-------------------------------------------------------------------------------

#- Readable description of the cells with coupling number k in a fleet x age matrix
.couplingLabel <- function(mat,k,noAge=character()){
  if(is.null(dim(mat))){                         # vector by age (e.g. logN.vars)
    ages <- names(mat)[which(mat==k)]
    return(.ageRange(ages))
  }
  cells <- which(mat==k,arr.ind=TRUE)
  if(nrow(cells)==0) return(as.character(k))
  flts  <- unique(rownames(mat)[cells[,1]])
  paste(sapply(flts,function(f){
    if(f %in% noAge) return(f)                   # biomass indices: no ages
    ages <- colnames(mat)[cells[rownames(mat)[cells[,1]]==f,2]]
    paste(f,.ageRange(ages))}),collapse=", ")
}
.ageRange <- function(ages){
  if(length(ages)==0) return("")
  if(any(grepl("-",ages))) return(paste(ages[1],ages[length(ages)],sep=" to "))  # age pairs (cor.obs)
  a <- suppressWarnings(as.numeric(ages))
  if(any(is.na(a))) return(paste(ages,collapse=","))
  if(length(a)==1) return(as.character(a))
  if(all(diff(sort(a))==1)) return(paste(min(a),max(a),sep="-"))
  return(paste(sort(a),collapse=","))
}

#' Parameters with readable names
#'
#' Returns all estimated (fixed-effect) parameters of an assessment with a
#' readable name derived from the couplings in the control object (e.g.
#' \code{"catchability HERAS 5-8"}, \code{"obs sd catch unique 0-1"},
#' \code{"N process sd 1-8"}), their type, the estimate on the scale of the
#' model and the back-transformed value with 95\% confidence interval:
#' \code{exp()} for \code{log*} parameters, the inverse logit for
#' \code{logit*} parameters, the correlation for \code{itrans_rho}, and
#' \code{exp()+1} for the prediction-variance power (\code{predVarObs});
#' other parameters are reported as estimated.
#'
#' @param object An \code{\link{FLSAM-class}} or \code{\link{FLSAMs}} object.
#' @return A \code{data.frame} with columns \code{name} (SAM parameter),
#'   \code{number}, \code{label}, \code{type}, \code{estimate}, \code{std.dev},
#'   \code{value}, \code{lbnd}, \code{ubnd} (and \code{model} for \code{FLSAMs}).
#' @seealso \code{\link{params.plot}}, \code{\link{params}}
#' @examples
#' data(NSH); data(NSH.sam)
#' parameter.table(NSH.sam)
#' @export
parameter.table <- function(object){
  if(is(object,"FLSAMs")){
    res <- do.call(rbind,lapply(seq_along(object),function(i)
             cbind(model=names(object)[i],parameter.table(object[[i]]))))
    res$model <- factor(res$model,levels=names(object))
    return(res)
  }
  ctrl  <- .upgradeControl(object@control)
  p     <- object@params
  fixed <- c("logFpar","logQpow","logSdLogFsta","logSdLogN","logSdLogP","logSdLogObs","logSdLogTotalObs",
             "transfIRARdist","sigmaObsParUS","rec_pars","itrans_rho","rhop","logScale","logitReleaseSurvival",
             "logitRecapturePhi","logAlphaSCB","sepFalpha","sepFlogitRho","sepFlogSd","predVarObs",
             "logPhiSW","logSdProcLogSW","meanLogSW","logSdLogSW","logPhiCW","logSdProcLogCW","meanLogCW",
             "logSdLogCW","logPhiMO","logSdProcLogitMO","meanLogitMO","logSdMO","logPhiNM","logSdProcLogNM",
             "meanLogNM","logSdLogNM","logXtraSd","initF","initN")
  res <- do.call(rbind,lapply(intersect(unique(p$name),fixed),function(nm){
    x <- p[p$name==nm,c("value","std.dev")]
    n <- nrow(x)
    #- parameters reported twice (as parameter and as derived quantity)
    if(n %% 2==0 && n>0 && isTRUE(all.equal(x$value[1:(n/2)],x$value[n/2+1:(n/2)]))) x <- x[1:(n/2),]
    data.frame(name=nm,number=seq_len(nrow(x)),estimate=x$value,std.dev=x$std.dev)
  }))
  if(is.null(res)) return(NULL)
  catchF <- names(ctrl@fleets)[ctrl@fleets==0]
  bioF   <- names(ctrl@fleets)[ctrl@fleets==3]
  lab <- function(nm,k){
    switch(nm,
      logFpar         = paste("catchability",.couplingLabel(ctrl@catchabilities,k-1,bioF)),
      logQpow         = paste("q power",.couplingLabel(ctrl@power.law.exps,k-1,bioF)),
      logSdLogObs     = paste("obs sd",.couplingLabel(ctrl@obs.vars,k-1,bioF)),
      logSdLogFsta    = paste("F process sd",.couplingLabel(ctrl@f.vars,k-1)),
      logSdLogN       = paste("N process sd",.couplingLabel(ctrl@logN.vars,k-1)),
      logSdLogP       = paste("component process sd",k),
      transfIRARdist  = paste("obs correlation",.couplingLabel(ctrl@cor.obs,k-1)),
      predVarObs      = paste("pred-var power",.couplingLabel(ctrl@predVarObsLink,k-1,bioF)),
      itrans_rho      = paste("F correlation",catchF[pmin(k,length(catchF))]),
      rec_pars        = paste("SR parameter",letters[k]),
      logScale        = paste("catch scaling",k),
      logSdLogTotalObs= paste("obs sd total",names(ctrl@fleets)[ctrl@likFlag=="ALN"][k]),
      sigmaObsParUS   = paste("obs covariance (US)",k),
      logitReleaseSurvival = paste("tag release survival",k),
      logitRecapturePhi    = paste("tag recapture phi",k),
      logXtraSd       = paste("extra obs sd",k),
      paste(nm,k))
  }
  res$label <- mapply(lab,res$name,res$number)
  res$type  <- ifelse(res$name %in% c("logFpar","logQpow","logScale"),"catchability / scaling",
               ifelse(res$name %in% c("logSdLogObs","logSdLogTotalObs","logXtraSd","predVarObs"),"observation sd",
               ifelse(res$name %in% c("logSdLogFsta","logSdLogN","logSdLogP"),"process sd",
               ifelse(res$name %in% c("transfIRARdist","itrans_rho","sigmaObsParUS","rhop"),"correlation",
               ifelse(res$name %in% c("rec_pars"),"stock-recruitment",
               ifelse(grepl("Release|Recapture",res$name),"tags","other"))))))
  tf <- function(nm,x){
    if(nm=="itrans_rho") return(2/(1+exp(-2*x))-1)
    if(nm=="predVarObs") return(exp(x)+1)
    if(grepl("^logit",nm)) return(plogis(x))
    if(grepl("^log",nm)) return(exp(x))
    return(x)
  }
  sd0 <- ifelse(is.na(res$std.dev),0,res$std.dev)
  res$value <- mapply(tf,res$name,res$estimate)
  res$lbnd  <- mapply(tf,res$name,res$estimate-1.96*sd0)
  res$ubnd  <- mapply(tf,res$name,res$estimate+1.96*sd0)
  rownames(res) <- NULL
  res[,c("name","number","label","type","estimate","std.dev","value","lbnd","ubnd")]
}

#' Parameter comparison plot
#'
#' Plots the estimated parameters with readable names (see
#' \code{\link{parameter.table}}) and their 95\% confidence intervals, by
#' type of parameter. With an \code{FLSAMs} object (several assessments or
#' retrospective peels) the estimates are shown side by side.
#'
#' @param object An \code{\link{FLSAM-class}} or \code{\link{FLSAMs}} object.
#' @param type Optional selection of parameter types (see \code{parameter.table}).
#' @return A \code{ggplot} object.
#' @examples
#' data(HERAS.sams)
#' params.plot(HERAS.sams)
#' params.plot(HERAS.sams, type="catchability / scaling")
#' @export
params.plot <- function(object,type=NULL){
  tab <- parameter.table(object)
  if(!is.null(type)) tab <- tab[tab$type %in% type,]
  if(!"model" %in% names(tab)) tab$model <- factor(object@name)
  tab$label <- factor(tab$label,levels=unique(tab$label))
  nm  <- nlevels(tab$model)
  ggplot(tab,aes(x=label,y=value,colour=model)) +
    geom_pointrange(aes(ymin=lbnd,ymax=ubnd),position=position_dodge(width=0.6),size=0.3) +
    scale_colour_manual(values=FLSAM.colours(nm),name=NULL) +
    facet_wrap(~type,scales="free") +
    labs(x=NULL,y="Estimate (95% CI)") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          legend.position=if(nm>1) "right" else "none")
}

#' Analytical confidence ellipse of SSB and Fbar
#'
#' Confidence ellipse of SSB and Fbar in one or more years, based on the
#' estimated covariance of log(SSB) and log(Fbar) (assumed bivariate normal on
#' the log scale), for one or more assessments. A fast analytical alternative
#' to \code{\link{otolith}}.
#'
#' @param object An \code{\link{FLSAM-class}} or \code{\link{FLSAMs}} object.
#' @param years Years (default: the last assessment year of each object).
#' @param level Confidence level (default 0.95).
#' @param npoints Number of points on each ellipse.
#' @return A \code{ggplot} object; the ellipses are in its \code{data} element.
#' @examples
#' data(NSH.sam); data(HERAS.sams)
#' ssbfbar.ellipse(NSH.sam, years=2008:2011)
#' ssbfbar.ellipse(HERAS.sams)
#' @export
ssbfbar.ellipse <- function(object,years=NULL,level=0.95,npoints=100){
  sams <- if(is(object,"FLSAMs")) object else FLSAMs(object)
  if(is.null(names(sams)) || any(names(sams)=="")) names(sams) <- sapply(sams,name)
  ell <- list(); est <- list()
  for(i in seq_along(sams)){
    s    <- sams[[i]]
    rc   <- s@rescov
    yrs  <- s@range["minyear"]:s@range["maxyear"]
    iS   <- which(rownames(rc)=="logssb"); iF <- which(rownames(rc)=="logfbar")
    if(length(iS)!=length(yrs) || length(iF)!=length(yrs))
      stop("Covariance of log(SSB) and log(Fbar) not available in ",names(sams)[i])
    yy   <- if(is.null(years)) max(yrs) else intersect(years,yrs)
    lS   <- log(ssb(s)$value); lF <- log(fbar(s)$value)
    for(y in yy){
      j   <- which(yrs==y)
      cv  <- rc[c(iF[j],iS[j]),c(iF[j],iS[j])]
      e   <- ellipse::ellipse(cv,centre=c(lF[j],lS[j]),level=level,npoints=npoints)
      ell[[length(ell)+1]] <- data.frame(model=names(sams)[i],year=y,fbar=exp(e[,1]),ssb=exp(e[,2]))
      est[[length(est)+1]] <- data.frame(model=names(sams)[i],year=y,fbar=exp(lF[j]),ssb=exp(lS[j]))
    }
  }
  ell <- do.call(rbind,ell); est <- do.call(rbind,est)
  ell$model <- factor(ell$model,levels=names(sams)); est$model <- factor(est$model,levels=names(sams))
  byYear <- length(sams)==1 && length(unique(ell$year))>1
  ell$grp <- if(byYear) factor(ell$year) else ell$model
  est$grp <- if(byYear) factor(est$year) else est$model
  ng <- nlevels(ell$grp)
  p <- ggplot(ell,aes(x=fbar,y=ssb,colour=grp)) +
         geom_path(aes(group=interaction(model,year))) +
         geom_point(data=est,size=2) +
         scale_colour_manual(values=FLSAM.colours(ng),name=if(byYear) "Year" else NULL) +
         labs(x=paste0("Fbar (",sams[[1]]@range["minfbar"],"-",sams[[1]]@range["maxfbar"],")"),y="SSB",
              title=paste0(100*level,"% confidence ellipse of SSB and Fbar")) +
         theme_bw()
  if(!byYear && length(unique(ell$year))>1)
    p <- p + geom_text(aes(label=year),data=est,vjust=-0.8,size=3,show.legend=FALSE)
  if(ng==1) p <- p + theme(legend.position="none")
  p
}
