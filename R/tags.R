#-------------------------------------------------------------------------------
# Tag-recapture data (SAM fleet type 5)
#-------------------------------------------------------------------------------

#' Tag-recapture parameters
#'
#' Estimates of the tag-recapture parameters of an assessment with tag data
#' (see the \code{recapture} argument of \code{\link{FLSAM.control}}): the
#' survival of released fish (\code{ReleaseSurvival}) and the recapture
#' (overdispersion) parameter \code{RecapturePhi}. Both are estimated on the
#' logit scale and back-transformed with the inverse logit, as in SAM.
#'
#' @param object An \code{\link{FLSAM-class}} object fitted with tag data.
#' @return A \code{data.frame} with columns \code{name}, \code{number},
#'   \code{value}, \code{lbnd} and \code{ubnd} (95\% confidence limits).
#' @seealso \code{\link{tags.plot}}
#' @export
tag.pars <- function(object){
  p <- subset(object@params,name %in% c("logitReleaseSurvival","logitRecapturePhi"))
  if(nrow(p)==0) stop("FLSAM object does not contain tag-recapture parameters")
  p$number <- ave(seq_len(nrow(p)),p$name,FUN=seq_along)
  data.frame(name=sub("^logit","",p$name),number=p$number,
             value=plogis(p$value),
             lbnd=plogis(p$value-1.96*p$std.dev),
             ubnd=plogis(p$value+1.96*p$std.dev),row.names=NULL)
}

#' Observed and predicted tag recaptures
#'
#' Plots the observed (points) and predicted (lines) number of recaptures by
#' release year, for each age at release, coloured by the number of years at
#' liberty (comparable to the \code{fitplot.tags} of the WGWIDE mackerel
#' assessment).
#'
#' @param sam An \code{\link{FLSAM-class}} object fitted with tag data and with
#'   a stored fit (see \code{\link{samfit}}).
#' @param type Optional tag types (column \code{Type}) to plot.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{tag.pars}}
#' @export
tags.plot <- function(sam,type=NULL){
  fit <- samfit(sam)
  aux <- as.data.frame(fit$data$aux)
  idx <- which(fit$data$fleetTypes[aux$fleet]==5)
  if(length(idx)==0) stop("No tag-recapture data in this assessment")
  dat <- data.frame(year=aux$year[idx],age=aux$age[idx],Type=aux$Type[idx],
                    liberty=aux$RecaptureY[idx]-aux$year[idx],
                    obs=fit$data$logobs[idx],pred=fit$rep$predObs[idx])
  if(!is.null(type)) dat <- dat[dat$Type %in% type,]
  dat$liberty <- factor(dat$liberty,levels=sort(unique(dat$liberty)))
  dat$relAge  <- factor(paste("Release age",dat$age),levels=paste("Release age",sort(unique(dat$age))))
  p <- ggplot(dat,aes(x=year,colour=liberty,linetype=liberty)) +
         geom_line(aes(y=pred)) +
         geom_point(aes(y=obs),size=1) +
         .discreteScales(nlevels(dat$liberty),name="Years at\nliberty") +
         labs(x="Release year",y="Recaptures (line: predicted, points: observed)",title=sam@name) +
         theme_bw()
  if(length(unique(dat$Type))>1){
    p + facet_grid(paste("Type",Type)~relAge,scales="free_y")
  } else {
    p + facet_wrap(~relAge,scales="free_y")
  }
}
