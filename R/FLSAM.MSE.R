#' Fast FLSAM assessments for MSE simulations
#'
#' \code{FLSAM.MSE} runs an assessment for every iteration of an
#' \code{FLStock} (and \code{FLIndices}) and fills the \code{stock.n} and
#' \code{harvest} slots of the stock with the estimates. \code{FLSAM.fast}
#' does the same for a single iteration. Both skip the extra Newton steps and
#' the residual calculation to save time.
#'
#' @param stcks An \code{FLStock} (multi-fleet \code{FLStocks} are not
#'   implemented).
#' @param tun An \code{FLIndices} object. Indices with a single iteration are
#'   used for all stock iterations.
#' @param ctrl An \code{\link{FLSAM.control-class}} object.
#' @param catch.vars Optional catch weights (see \code{\link{FLSAM}}).
#' @param starting.sam Not used (kept for backward compatibility).
#' @param return.sam If \code{TRUE}, also return the \code{FLSAM} results.
#' @param ncores Number of cores used to run the iterations in parallel
#'   (\code{FLSAM.MSE} only); 1 runs sequentially.
#' @param ... Further arguments passed to \code{\link[stockassessment]{sam.fit}}.
#' @return \code{FLSAM.MSE}: an \code{FLStocks} (element \code{residual}) with
#'   estimated numbers and fishing mortality, or (\code{return.sam=TRUE}) an \code{\link{FLSAMs}} of the
#'   converged iterations. \code{FLSAM.fast}: the \code{FLStocks}, or
#'   (\code{return.sam=TRUE}) a list with elements \code{sam.fit} and \code{stk}.
#'   Iterations that fail are set to \code{NA}.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' NSH.iter <- propagate(NSH, 3)
#' res <- FLSAM.MSE(NSH.iter, NSH.tun, NSH.ctrl)
#' }
#' @rdname FLSAM.MSE
#' @export
FLSAM.MSE <-function(stcks,tun,ctrl,catch.vars=NULL,starting.sam=NULL,return.sam=F,ncores=1,...){

  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(is(stcks, "FLStocks")) stop("Not implemented for multi-fleet yet")
  if(is(stcks, "FLStock")) stcks <- FLStocks(residual=stcks)
  if(any(unlist(lapply(stcks,function(x)dims(x)$iter))<=1) & any(unlist(lapply(tun,function(x)dims(x)$iter))<=1))
    stop("Running in MSE mode means supplying FLStock and FLIndices with more iters than 1. Use FLSAM() instead")

  #Count iters
  iters <- dims(stcks[["residual"]])$iter

  #Turn residuals off
  ctrl@residuals <- FALSE
  ctrl           <- .upgradeControl(ctrl)
  dots           <- list(...)

  #Run one iteration
  runIter <- function(i){
    stckI <- FLStocks(residual=iter(stcks[["residual"]],i))
    tunI  <- FLIndices(lapply(tun,function(x) if(dims(x)$iter>1) iter(x,i) else x))
    data  <- FLSAM2SAM(stckI,tunI,ctrl@sumFleets,catch.vars)
    conf  <- ctrl2conf(ctrl,data)
    par   <- defpar(data,conf)
    try(do.call(sam.fit,c(list(data=data,conf=conf,parameters=par,sim.condRE=ctrl@simulate,
                               newtonsteps=0,silent=TRUE),dots)))
  }

  if(ncores>1){
    cl <- parallel::makeCluster(min(ncores,iters))
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl,suppressMessages(library(FLSAM)))
    res <- parallel::parLapply(cl,seq_len(iters),runIter)
  } else {
    res <- lapply(seq_len(iters),runIter)
  }
  ok <- sapply(res,function(x) !is(x,"try-error"))
  if(any(!ok)) warning(paste("Assessment failed for iteration(s)",paste(which(!ok),collapse=", ")))

  #- Return sam objects
  if(return.sam){
    resSAM <- lapply(res[ok],SAM2FLR,ctrl=ctrl,keep.fit=FALSE)
    names(resSAM) <- which(ok)
    return(as(resSAM,"FLSAMs"))
  }

  for(i in seq_len(iters)){
    if(ok[i]){
      stcks[["residual"]] <- .fillStock(stcks[["residual"]],res[[i]],i)
    } else {
      stcks[["residual"]]@stock.n[, , , , , i] <- NA
      stcks[["residual"]]@harvest[, , , , , i] <- NA
    }
  }
  return(stcks)
}

#' @rdname FLSAM.MSE
#' @export
FLSAM.fast <-function(stcks,tun,ctrl,catch.vars=NULL,starting.sam=NULL,return.sam=F,...){

  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  if(is(stcks, "FLStocks")) stop("Not implemented for multi-fleet yet")
  if(is(stcks, "FLStock")) stcks <- FLStocks(residual=stcks)

  #Turn residuals off
  ctrl@residuals <- FALSE
  ctrl  <- .upgradeControl(ctrl)

  data  <- FLSAM2SAM(FLStocks("residual"=stcks[["residual"]]),FLIndices(tun),ctrl@sumFleets,catch.vars)
  conf  <- ctrl2conf(ctrl,data)
  par   <- defpar(data,conf)

  res <- try(sam.fit(data,
                 conf,
                 par,
                 sim.condRE = ctrl@simulate,
                 newtonsteps = 0,
                 silent = T,...))

  # fill in stock objects
  if(!is(res,"try-error")){
    stcks[["residual"]] <- .fillStock(stcks[["residual"]],res,1)
  } else {
    stcks[["residual"]]@stock.n[] <- NA
    stcks[["residual"]]@harvest[] <- NA
  }

  if(return.sam){
    resSAM <- if(is(res,"try-error")) NA else SAM2FLR(res,ctrl)
    ret <- list(sam.fit=resSAM,stk=stcks)
  } else {
    ret <- stcks
  }
  return(ret)
}

#- Fill stock.n and harvest of iteration i of an FLStock with a sam fit, only for the
#  years of the stock (surveys may run beyond the last catch year)
.fillStock <- function(stck,fit,i){
  nt   <- ntable(fit)
  ft   <- faytable(fit)
  yrs  <- intersect(rownames(nt),dimnames(stck@stock.n)$year)
  ages <- intersect(colnames(nt),dimnames(stck@stock.n)$age)
  stck@stock.n[ages,yrs,,,,i] <- t(nt[yrs,ages,drop=FALSE])
  stck@harvest[ages,yrs,,,,i] <- t(ft[yrs,ages,drop=FALSE])
  units(stck@harvest) <- "f"
  return(stck)
}
