#-------------------------------------------------------------------------------
# New functionality of FLSAM 3.1 applied to the 2026 herring assessment.
#   Rscript run_extras.R <FLSAM library (new)> <output dir> <assessment dir>
# Uses the fits saved by run_version.R (label "new").
#-------------------------------------------------------------------------------
args      <- commandArgs(TRUE)
lib       <- args[1]; outRoot <- args[2]; assessDir <- args[3]
outDir    <- file.path(outRoot,"new_extra"); dir.create(outDir,showWarnings=FALSE,recursive=TRUE)
suppressMessages({library(FLSAM,lib.loc=lib); library(ggplot2)})
load(file.path(outRoot,"new","assessment.RData"))
pal   <- FLSAM.colours(4)
pngf  <- function(name,width=1600,height=1200) png(file.path(outDir,paste0(name,".png")),width=width,height=height,res=150)
plotTry <- function(name,expr,...){
  pngf(name,...)
  r <- try({p <- eval.parent(substitute(expr)); if(inherits(p,"ggplot")) print(p)},silent=TRUE)
  if(inherits(r,"try-error")){ plot.new(); title(paste(name,"\nFAILED:",conditionMessage(attr(r,"condition")))) }
  dev.off(); invisible(r)
}
log <- file(file.path(outDir,"extras_log.txt"),open="wt"); sink(log); sink(log,type="message")
lastYr <- range(NSH)["maxyear"]

#- Residual tests (stockassessment::residplot) ---------------------------------
plotTry("sf_residual_tests", {rt <- residual.tests(NSH.sam,plot=FALSE); write.csv(rt,file.path(outDir,"sf_residual_tests.csv"),row.names=FALSE); attr(rt,"plot")},width=2600,height=1800)
plotTry("mf_residual_tests", {rt <- residual.tests(NSH3f.sam,plot=FALSE); write.csv(rt,file.path(outDir,"mf_residual_tests.csv"),row.names=FALSE); attr(rt,"plot")},width=2600,height=1800)

#- Other new diagnostics ----------------------------------------------------------
plotTry("sf_obscor_empirical",  obscor.empirical.plot(NSH.sam),width=2400,height=1800)
plotTry("mf_obscor_empirical",  obscor.empirical.plot(NSH3f.sam),width=2400,height=1800)
plotTry("sf_obs_sd",            obssd.plot(NSH.sam))
plotTry("sf_obs_relative_weight", obssd.plot(NSH.sam,show.rel.w=TRUE))
plotTry("mf_obs_sd",            obssd.plot(NSH3f.sam))
plotTry("sf_data_overview",     data.plot(NSH.sam))
plotTry("mf_data_overview",     data.plot(NSH3f.sam))
plotTry("sf_selectivity",       fselectivity.plot(NSH.sam))
plotTry("mf_selectivity",       fselectivity.plot(NSH3f.sam),width=2400,height=1200)

#- Jitter test -------------------------------------------------------------------
cat("\n== Jitter test (sf, 5 runs)\n")
jt <- try(jitter.test(NSH.sam,nojit=5,seed=1))
if(!inherits(jt,"try-error")){
  print(attr(jt,"table")); write.csv(attr(jt,"table"),file.path(outDir,"sf_jitter_table.csv"),row.names=FALSE)
  plotTry("sf_jitter",plot(jt))
}

#- Mohn's rho (with mean) -----------------------------------------------------------
cat("\n== Mohn's rho\n")
mr <- do.call(rbind,lapply(c("sf","mf"),function(m){
        r <- get(ifelse(m=="sf","NSH.retro","NSH3f.retro"))
        do.call(rbind,lapply(c("ssb","fbar","rec"),function(tp){
          x <- mohns.rho(r,span=length(r)-1,type=tp)
          data.frame(model=m,type=tp,mohn=attr(x,"mohn"))}))}))
print(mr); write.csv(mr,file.path(outDir,"mohns_rho.csv"),row.names=FALSE)

#- Forecast, reference points, per recruit (sf) ------------------------------------
cat("\n== Forecasts (sf, status quo F)\n")
fc <- try(FLSAM.forecast(NSH.sam,fscale=c(1,1,1),ave.years=(lastYr-2):lastYr,rec.years=(lastYr-9):lastYr))
if(!inherits(fc,"try-error")){
  print(fc); write.csv(fc,file.path(outDir,"sf_forecast.csv"),row.names=FALSE)
  plotTry("sf_forecast", ggplot(fc,aes(year,value))+
            geom_ribbon(aes(ymin=lbnd,ymax=ubnd),fill=pal[3],alpha=0.4)+geom_line(colour=pal[1],linewidth=1)+
            facet_wrap(~quantity,scales="free_y")+theme_bw()+labs(x="Year",y=NULL,title="FLSAM.forecast: status quo F"))
}
mf <- try(FLSAM.modelforecast(NSH.sam,fscale=c(1,1,1)))
if(!inherits(mf,"try-error")){ print(mf); write.csv(mf,file.path(outDir,"sf_modelforecast.csv"),row.names=FALSE) }

cat("\n== Reference points (sf, per recruit)\n")
rp <- try(FLSAM.refpts(NSH.sam,c("Max","0.1dYPR","0.35SPR")))
if(!inherits(rp,"try-error")){ print(rp); write.csv(rp,file.path(outDir,"sf_refpts.csv"),row.names=FALSE) }

cat("\n== Yield per recruit (sf)\n")
yp <- try(FLSAM.ypr(NSH.sam))
if(!inherits(yp,"try-error")){
  print(attr(yp,"refpts"))
  refs <- data.frame(name=c("Fmax","F0.1","F35%SPR"),f=attr(yp,"refpts")[c("fmax","f01","fsprProp")])
  plotTry("sf_ypr", ggplot(yp,aes(fbar,yield))+geom_line(colour=pal[1],linewidth=1)+
            geom_vline(aes(xintercept=f,colour=name),data=refs,linetype=2)+
            scale_colour_manual(values=pal[2:4])+theme_bw()+labs(x="Fbar",y="Yield per recruit",colour=NULL,title="FLSAM.ypr"))
}

#- Catch by fleet (mf) -----------------------------------------------------------------
cat("\n== Catch by fleet (mf)\n")
cbf <- try(FLSAM.catchbyfleet(NSH3f.sam,obs.show=TRUE))
if(!inherits(cbf,"try-error")){
  print(head(cbf)); write.csv(cbf,file.path(outDir,"mf_catch_by_fleet.csv"),row.names=FALSE)
  plotTry("mf_catch_by_fleet", ggplot(cbf,aes(year,value))+
            geom_ribbon(aes(ymin=lbnd,ymax=ubnd,fill=fleet),alpha=0.3)+geom_line(aes(colour=fleet))+
            geom_point(aes(y=obs,colour=fleet),size=0.8)+
            scale_colour_manual(values=FLSAM.colours(length(unique(cbf$fleet))))+
            scale_fill_manual(values=FLSAM.colours(length(unique(cbf$fleet))))+
            facet_wrap(~fleet,scales="free_y")+theme_bw()+labs(x="Year",y="Catch (t)",title="FLSAM.catchbyfleet (lines: model, points: observed)"),
          width=2400,height=1400)
}

#- Model description and the stored stockassessment fit ------------------------------
cat("\n== Model description (sf)\n"); model.description(NSH.sam)
cat("\n== stockassessment::modeltable on the stored fits\n")
print(stockassessment::modeltable(c(sf=samfit(NSH.sam),mf=samfit(NSH3f.sam))))

#- Old saved model objects of the assessment repository ----------------------------------
cat("\n== Objects saved with the old FLSAM in",file.path(assessDir,"model"),"\n")
old <- new.env()
load(file.path(assessDir,"model","NSAS_HAWG2026_sf.RData"),envir=old)
cat("saved with FLSAM",old$NSH.sam@info["FLSAM.version",1],"; has slot 'fit':",.hasSlot(old$NSH.sam,"fit"),
    "; control has 'hockeyStickCurve':",.hasSlot(old$NSH.ctrl,"hockeyStickCurve"),"\n")
invisible(capture.output(print(old$NSH.sam))); cat("printing old FLSAM: ok\n")
plotTry("old_saved_sf_stock_trajectory",plot(FLSAMs(saved.HAWG2026=upgradeFLSAM(old$NSH.sam),rerun.FLSAM3.1=NSH.sam)))
refit <- FLSAM(old$NSH,old$NSH.tun,old$NSH.ctrl,silent=TRUE)
cat("refit old objects with FLSAM 3.1: nlogl saved =",old$NSH.sam@nlogl,"; refit =",refit@nlogl,"\n")

sink(type="message"); sink(); close(log)
cat("extras done\n")
