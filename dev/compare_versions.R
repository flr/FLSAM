#-------------------------------------------------------------------------------
# Compare FLSAM 3.0.0.9005 (original, tag v3.0.0.9005) with FLSAM 3.1.0
# (branch dev-sam-align) and demonstrate the new functionality.
#
# Part A runs the same assessments with both versions, each in a fresh R
#        process with its own library, and compares all results.
# Part B shows the new functionality and the bug fixes (new version only).
#
# Requirements: callr, remotes (only to install the versions if needed)
# Output: dev/compare_output/ (comparison table, plots, R objects)
#-------------------------------------------------------------------------------

#- Settings --------------------------------------------------------------------
libRoot   <- "C:/Users/NielsHintzen/Repositories/FLSAM_libs"  # one sub-library per version
lib.orig  <- file.path(libRoot,"orig")
lib.new   <- file.path(libRoot,"new")
repo.path <- "C:/Users/NielsHintzen/Repositories/FLSAM"       # local checkout of dev-sam-align (NULL: install from GitHub)
outDir    <- file.path(repo.path,"dev","compare_output")
reinstall.new <- FALSE   # TRUE: reinstall the new version (e.g. after changes)
rerun.A   <- TRUE        # FALSE: reuse the saved results of part A
tol       <- 1e-8        # relative tolerance for "equal"

dir.create(outDir,showWarnings=FALSE,recursive=TRUE)

#- Install both versions in their own library (once) ---------------------------
if(!file.exists(file.path(lib.orig,"FLSAM","DESCRIPTION"))){
  dir.create(lib.orig,showWarnings=FALSE,recursive=TRUE)
  remotes::install_github("flr/FLSAM",ref="v3.0.0.9005",lib=lib.orig,upgrade="never",dependencies=FALSE)
}
if(reinstall.new | !file.exists(file.path(lib.new,"FLSAM","DESCRIPTION"))){
  dir.create(lib.new,showWarnings=FALSE,recursive=TRUE)
  if(is.null(repo.path)){
    remotes::install_github("flr/FLSAM",ref="dev-sam-align",lib=lib.new,upgrade="never",dependencies=FALSE)
  } else {
    install.packages(repo.path,repos=NULL,type="source",lib=lib.new)
  }
}

#===============================================================================
# PART A: identical results?
#===============================================================================

#- Function run in a fresh R session for each version
runCases <- function(lib){
  suppressMessages(library(FLSAM,lib.loc=lib))
  data(NSH); data(NSH.sam); data(SOL)
  fq  <- function(x) as.data.frame(x)[,c("age","year","area","data")]
  ext <- function(s) list(nlogl=s@nlogl,nopar=s@nopar,params=s@params,vcov=unname(s@vcov),
                          rescov=unname(s@rescov),stock.n=fq(s@stock.n),harvest=fq(s@harvest),
                          residuals=s@residuals,ssb=ssb(s),fbar=fbar(s),rec=rec(s),tsb=tsb(s),
                          catch=catch(s),f=f(s),catchabilities=catchabilities(s),
                          obs.var=obs.var(s),n.var=n.var(s),f.var=f.var(s))
  tryNull <- function(expr) tryCatch(expr,error=function(e) paste("ERROR:",conditionMessage(e)))
  out <- list(version=as.character(packageVersion("FLSAM",lib.loc=lib)))

  #- Data and configuration that is passed to SAM
  out$NSH.data <- FLSAM2SAM(FLStocks(residual=NSH),NSH.tun,NSH.ctrl@sumFleets)
  out$NSH.conf <- ctrl2conf(NSH.ctrl,out$NSH.data)
  out$SOL.conf <- ctrl2conf(SOL.ctrl,FLSAM2SAM(FLStocks(residual=SOL),SOL.tun,SOL.ctrl@sumFleets))

  #- Assessments
  out$NSH     <- tryNull(ext(FLSAM(NSH,NSH.tun,NSH.ctrl,silent=TRUE)))
  out$SOL     <- tryNull(ext(FLSAM(SOL,SOL.tun,SOL.ctrl,silent=TRUE)))
  r           <- retro(NSH,NSH.tun,NSH.ctrl,retro=2)
  out$retro   <- lapply(r,ext)
  out$mohn    <- tryNull(mohns.rho(r,ref.year=2011,span=2,type="ssb"))
  l           <- loo(NSH,NSH.tun,NSH.ctrl)
  out$loo     <- lapply(l,function(x) list(nlogl=x@nlogl,ssb=ssb(x),fbar=fbar(x)))
  #- Fix the observation sd of coupling 1 (HERAS youngest ages) at 0.3
  sp          <- data.frame(par="obs.vars",number=1,value=log(0.3))
  out$setpars <- tryNull(ext(FLSAM(NSH,NSH.tun,NSH.ctrl,set.pars=sp,silent=TRUE)))
  return(out)
}

if(rerun.A | !file.exists(file.path(outDir,"partA.rds"))){
  cat("Running original version ...\n")
  resOrig <- callr::r(runCases,args=list(lib=lib.orig),libpath=c(lib.orig,.libPaths()),show=FALSE)
  cat("Running new version ...\n")
  resNew  <- callr::r(runCases,args=list(lib=lib.new), libpath=c(lib.new, .libPaths()),show=FALSE)
  saveRDS(list(orig=resOrig,new=resNew),file.path(outDir,"partA.rds"))
} else {
  tmp <- readRDS(file.path(outDir,"partA.rds")); resOrig <- tmp$orig; resNew <- tmp$new
}

#- Comparison helpers
numOf  <- function(x){
  if(is.data.frame(x)) x <- x[,sapply(x,is.numeric),drop=FALSE]
  if(is.list(x) && !is.data.frame(x)) return(unlist(lapply(x,numOf)))
  as.numeric(unlist(x))
}
compare <- function(a,b,case,item,expected=""){
  if(identical(a,b)){
    st <- "identical"; md <- 0
  } else {
    na <- numOf(a); nb <- numOf(b)
    if(length(na)!=length(nb)){
      st <- "DIFFERENT (structure)"; md <- NA
    } else {
      md <- suppressWarnings(max(abs(na-nb)/pmax(abs(nb),1e-12),na.rm=TRUE))
      st <- ifelse(isTRUE(all.equal(na,nb,tolerance=tol)),"equal within tolerance","DIFFERENT")
    }
  }
  data.frame(case=case,item=item,status=st,max.rel.diff=signif(md,3),expected=expected)
}

items <- c("nlogl","nopar","params","vcov","rescov","stock.n","harvest","residuals","ssb","fbar","rec",
           "tsb","catch","f","catchabilities","obs.var","n.var","f.var")
cmpTab <- rbind(
  compare(resOrig$NSH.conf,resNew$NSH.conf,"NSH","SAM configuration (ctrl2conf)"),
  compare(resOrig$SOL.conf,resNew$SOL.conf,"SOL","SAM configuration (ctrl2conf)"),
  do.call(rbind,lapply(setdiff(names(resOrig$NSH.data),c("landMeanWeight","disMeanWeight")),function(i)
    compare(resOrig$NSH.data[[i]],resNew$NSH.data[[i]],"NSH",paste("SAM data:",i)))),
  compare(resOrig$NSH.data$landMeanWeight,resNew$NSH.data$landMeanWeight,"NSH","SAM data: landMeanWeight",
          "may differ: landings.wt is now used instead of catch.wt (bug fix)"),
  compare(resOrig$NSH.data$disMeanWeight,resNew$NSH.data$disMeanWeight,"NSH","SAM data: disMeanWeight",
          "may differ: discards.wt is now used instead of catch.wt (bug fix)"),
  do.call(rbind,lapply(items,function(i) compare(resOrig$NSH[[i]],resNew$NSH[[i]],"NSH assessment",i))),
  do.call(rbind,lapply(items,function(i) compare(resOrig$SOL[[i]],resNew$SOL[[i]],"SOL assessment",i))),
  do.call(rbind,lapply(names(resOrig$retro),function(y) do.call(rbind,lapply(items,function(i)
    compare(resOrig$retro[[y]][[i]],resNew$retro[[y]][[i]],paste("retro peel",y),i))))),
  do.call(rbind,lapply(names(resOrig$loo),function(y) do.call(rbind,lapply(c("nlogl","ssb","fbar"),function(i)
    compare(resOrig$loo[[y]][[i]],resNew$loo[[y]][[i]],paste("loo:",y),i))))),
  compare(abs(resOrig$mohn$rho),abs(resNew$mohn$rho),"mohns.rho","|rho|"),
  compare(resOrig$mohn$rho,-resNew$mohn$rho,"mohns.rho","rho (original) vs -rho (new)",
          "sign is reversed in the new version (bug fix, now as Mohn 1999)"),
  do.call(rbind,lapply(setdiff(items,c("vcov","residuals")),function(i) compare(resOrig$setpars[[i]],resNew$setpars[[i]],"set.pars (obs.vars 1 fixed)",i,
          "rounding differences: parameter order in the TMB map is now deterministic (was random via jitter)"))),
  compare(sort(resOrig$setpars$vcov),sort(resNew$setpars$vcov),"set.pars (obs.vars 1 fixed)","vcov (sorted values)",
          "same values; the order of the logSdLogObs parameters was random in the original version"),
  compare(subset(resOrig$setpars$residuals,year>min(year)),subset(resNew$setpars$residuals,year>min(year)),
          "set.pars (obs.vars 1 fixed)","residuals excluding first year",
          "one-observation-ahead residuals of the first data year are numerically unstable (depend on parameter order)")
)
rownames(cmpTab) <- NULL
write.csv(cmpTab,file.path(outDir,"comparison_table.csv"),row.names=FALSE)

cat("\n================ PART A: original",resOrig$version,"vs new",resNew$version,"================\n")
print(table(cmpTab$status))
cat("\nItems that are not identical:\n")
print(subset(cmpTab,status!="identical"),row.names=FALSE)
cat("\nOriginal mohns.rho:\n"); print(resOrig$mohn)
cat("New mohns.rho (Mohn's rho =",attr(resNew$mohn,"mohn"),"):\n"); print(resNew$mohn)

#===============================================================================
# PART B: new functionality and bug fixes (new version)
#===============================================================================
suppressMessages(library(FLSAM,lib.loc=lib.new))
library(ggplot2)
cat("\n================ PART B: new functionality (FLSAM",as.character(packageVersion("FLSAM",lib.loc=lib.new)),") ================\n")
data(NSH); data(NSH.sam); data(HERAS.sams)
pdf(file.path(outDir,"new_functionality.pdf"),width=10,height=8)

#- B1. The sam fit is stored in the FLSAM object --------------------------------
sam <- FLSAM(NSH,NSH.tun,NSH.ctrl,silent=TRUE)
fit <- samfit(sam)                      # class 'sam' from stockassessment
class(fit)
stockassessment::modeltable(fit)        # any stockassessment function works now
stockassessment::ssbplot(fit)
model.description(sam)

#- B2. Configuration extensions -----------------------------------------------
# fixVarToWeight per fleet: index.var of HERAS used as fixed variance
ctrlW <- NSH.ctrl
ctrlW@fixVarToWeight["HERAS"] <- TRUE
NSH.tunW <- NSH.tun
index.var(NSH.tunW[["HERAS"]])[] <- 0.3^2
datW  <- FLSAM2SAM(FLStocks(residual=NSH),NSH.tunW,ctrlW@sumFleets)
ctrl2conf(ctrlW,datW)$fixVarToWeight     # was always 0 in the original version
samW  <- FLSAM(NSH,NSH.tunW,ctrlW,silent=TRUE)
c(original.weights=nlogl(sam),HERAS.fixed.variance=nlogl(samW))

# hockeyStickCurve: smooth hockey stick (srr 63) with fixed smoothness
ctrl63 <- NSH.ctrl; ctrl63@srr <- 63L; ctrl63@hockeyStickCurve <- 20
sam63  <- try(FLSAM(NSH,NSH.tun,ctrl63,silent=TRUE))
if(!is(sam63,"try-error")) print(rec.pars(sam63))

# set.pars now accepts all fixed-effect parameters
setParsTable()
sp   <- data.frame(par=c("obs.vars","cor.F"),number=c(1,0),value=c(log(0.3),2))
samF <- FLSAM(NSH,NSH.tun,NSH.ctrl,set.pars=sp,silent=TRUE)
c(free=nopar(sam),fixed=nopar(samF))
lr.test(FLSAMs(fixed=samF,free=sam))

#- B3. ggplot2 plots (were lattice) ---------------------------------------------
print(plot(sam))
print(plot(HERAS.sams))
print(plot(FLSAMs(free=sam,fixed=samF),NSH))
print(cor.plot(sam))
NSHs <- NSH + sam
procerr.plot(NSHs,type="n")
procerr.plot(NSHs,type="mort",rel=TRUE)
procerr.plot(NSHs,type="tsb")

#- B4. New diagnostics (as in stockassessment) ---------------------------------
rt <- residual.tests(sam)                 # bias, variance, autocorrelation, mean-variance, normality
subset(rt,p<0.05)
print(obscor.empirical.plot(sam))
print(obssd.plot(sam))
print(obssd.plot(sam,show.rel.w=TRUE))
print(data.plot(sam))
print(fselectivity.plot(sam))
jt <- jitter.test(sam,nojit=5,seed=1)     # jitter test of the optimum
attr(jt,"table")
print(plot(jt))

#- B5. Forecast, reference points, yield per recruit -----------------------------
fc <- FLSAM.forecast(sam,fscale=c(1,1,1),ave.years=2009:2011,rec.years=2000:2011)
fc
print(ggplot(fc,aes(year,value))+geom_ribbon(aes(ymin=lbnd,ymax=ubnd),fill="grey")+geom_line()+
        facet_wrap(~quantity,scales="free_y")+theme_bw()+ggtitle("FLSAM.forecast"))
mf <- FLSAM.modelforecast(sam,fscale=c(1,1,1))
head(mf)
FLSAM.refpts(sam,c("Max","0.1dYPR","0.35SPR"))
yp <- FLSAM.ypr(sam); attr(yp,"refpts")
head(FLSAM.catchbyfleet(sam,obs.show=TRUE))

#- B6. Retrospective diagnostics -----------------------------------------------
r <- retro(NSH,NSH.tun,NSH.ctrl,retro=3)
mohns.rho(r,type="ssb")                   # ref.year and span have working defaults now
attr(mohns.rho(r,type="fbar"),"mohn")
retroResiduals(r,"HERAS",2005:2011)
retroSelectivity(r,2006:2010)
retroParams(r)

#- B7. Bug fixes ----------------------------------------------------------------
validObject(NSH.ctrl)                                     # crashed before
ctrl <- NSH.ctrl; ctrl@stockWeightMean[] <- c(0,2,2,5,5,7,7,7,7)
update(ctrl)@stockWeightMean                              # renumbered now
class(conf2ctrl(ctrl2conf(NSH.ctrl,datW),datW)@logN.vars)   # vector (was matrix)
names(looi(NSH,NSH.tun[c("HERAS","IBTS-Q1")],
           drop.from.control(NSH.ctrl,fleets=c("SCAI","IBTS0"))))  # looi without type failed
ctrl2 <- NSH.ctrl; ctrl2@srr <- 2L
samBH <- FLSAM(NSH,NSH.tun,ctrl2,silent=TRUE)
rm(list=intersect(ls(),"ctrl"))
rec.pars(samBH)                                           # no global 'ctrl' needed, no duplicate rows
fst <- FLSAM.fast(NSH,NSH.tun,NSH.ctrl)                   # failed when surveys run beyond the catch data
all.equal(c(fst[["residual"]]@stock.n),c(window(stock.n(sam),end=2011)))
NSH.iter <- propagate(NSH,2)
mse <- FLSAM.MSE(NSH.iter,NSH.tun,NSH.ctrl)               # did not run (undefined objects)
iter(mse[["residual"]]@stock.n[,"2011"],1:2)

#- B8. Objects saved with the original version (missing the new slots) ----------
old <- new.env()
data("NSH","NSH.sam","HERAS.sams",package="FLSAM",lib.loc=lib.orig,envir=old)  # old class definitions
.hasSlot(old$NSH.ctrl,"hockeyStickCurve")                 # FALSE: old object
old$NSH.ctrl                                              # printing upgrades automatically
old$NSH.sam
samOld <- FLSAM(old$NSH,old$NSH.tun,old$NSH.ctrl,silent=TRUE)   # fitting upgrades automatically
all.equal(nlogl(samOld),nlogl(sam))
names(retro(old$NSH,old$NSH.tun,old$NSH.ctrl,retro=1))
drop.from.control(old$NSH.ctrl,fleets="IBTS0")@fixVarToWeight
head(FLSAM.out(old$NSH,old$NSH.tun,old$NSH.sam))
ctrlUp <- upgradeFLSAM(old$NSH.ctrl)                      # or upgrade explicitly ...
validObject(ctrlUp)
samsUp <- upgradeFLSAM(old$HERAS.sams)                    # ... also FLSAM and FLSAMs objects
is.null(samsUp[[1]]@fit)
# objects that have been upgraded can be saved again with save()

dev.off()
cat("\nPlots written to",file.path(outDir,"new_functionality.pdf"),"\n")
