#-------------------------------------------------------------------------------
# Run the 2026 North Sea herring assessment (her.27.3a47d, HAWG 2026) with one
# FLSAM version and make the diagnostic plots of report_plots.R.
# Called by compare_her27.R in a fresh R session per version:
#   Rscript run_version.R <FLSAM library> <label> <output dir> <assessment dir> <n retro> [reuse]
# reuse = TRUE: load the fits saved earlier and only redo the plots
# Reproduces model_sf.R and model_mf.R (data_sms2023, addM = 0.02, IBP 2021 configs).
#-------------------------------------------------------------------------------
args      <- commandArgs(TRUE)
lib       <- args[1]; label <- args[2]; outRoot <- args[3]; assessDir <- args[4]
nRetro    <- as.numeric(args[5])
reuse     <- length(args)>5 && as.logical(args[6])
outDir    <- file.path(outRoot,label); dir.create(outDir,showWarnings=FALSE,recursive=TRUE)

suppressMessages({library(FLSAM,lib.loc=lib); library(ggplot2)})
cat("FLSAM",as.character(packageVersion("FLSAM",lib.loc=lib)),"from",lib,"\n")
source(file.path(assessDir,"utilities_model_config.R"))

addM  <- 0.02
pngf  <- function(name,width=1600,height=1200) png(file.path(outDir,paste0(name,".png")),width=width,height=height,res=150)
show1 <- function(p) if(inherits(p,c("trellis","ggplot"))) print(p)
plotTry <- function(name,expr,...){
  pngf(name,...)
  r <- try(show1(eval.parent(substitute(expr))),silent=TRUE)
  if(inherits(r,"try-error")){ plot.new(); title(paste(name,"\nFAILED:",conditionMessage(attr(r,"condition")))) }
  dev.off()
  invisible(r)
}
timing <- list()

if(reuse){
  load(file.path(outDir,"assessment.RData"))
} else {
#===============================================================================
# Single fleet (model_sf.R)
#===============================================================================
load(file.path(assessDir,"data","data_sms2023_sf.RData"))      # NSH, NSH.tun
NSH@m    <- NSH@m + addM
NSH.ctrl <- config_sf_IBPNSherring2021(NSH,NSH.tun,range(NSH)[2])
NSH.ctrl@residuals <- TRUE
timing$sf <- system.time(NSH.sam <- FLSAM(NSH,NSH.tun,NSH.ctrl,silent=TRUE))["elapsed"]
NSH@stock.n <- NSH.sam@stock.n[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]
NSH@harvest <- NSH.sam@harvest[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]
ctrlR <- NSH.ctrl; ctrlR@residuals <- FALSE
timing$sf.retro <- system.time(NSH.retro <- retro(NSH,NSH.tun,ctrlR,retro=nRetro))["elapsed"]

#===============================================================================
# Multi fleet (model_mf.R)
#===============================================================================
load(file.path(assessDir,"data","data_sms2023_mf.RData"))      # NSH3f, NSH.tun (same indices)
NSH3f@m   <- NSH3f@m + addM
NSH3.ctrl <- config_mf_IBPNSherring2021_alt(NSH3f,NSH.tun,range(NSH3f)[2])
NSH3.ctrl@residuals <- TRUE
timing$mf <- system.time(NSH3f.sam <- FLSAM(NSH3f,NSH.tun,NSH3.ctrl,silent=TRUE))["elapsed"]
ctrl3R <- NSH3.ctrl; ctrl3R@residuals <- FALSE
timing$mf.retro <- system.time(NSH3f.retro <- retro(NSH3f,NSH.tun,ctrl3R,retro=nRetro))["elapsed"]
}

#===============================================================================
# Diagnostic plots as in report_plots.R (same calls for both versions)
#===============================================================================
pdf(file.path(outDir,"sf_fit.pdf")); residual.diagnostics(NSH.sam); dev.off()
#- residual.diagnostics: one png per fleet and age (sf and mf)
dir.create(file.path(outDir,"resdiag"),showWarnings=FALSE)
unlink(list.files(file.path(outDir,"resdiag"),full.names=TRUE))
png(file.path(outDir,"resdiag","sf_%03d.png"),width=1500,height=1300,res=130); residual.diagnostics(NSH.sam); dev.off()
png(file.path(outDir,"resdiag","mf_%03d.png"),width=1500,height=1300,res=130); residual.diagnostics(NSH3f.sam); dev.off()
plotTry("sf_stock_trajectory",          plot(NSH.sam))
plotTry("sf_cor_params",                cor.plot(NSH.sam),width=2000,height=2000)
plotTry("sf_procErr_N",                 procerr.plot(NSH+NSH.sam,weight="stock.wt",type="n",rel=T))
plotTry("sf_procErr_M",                 procerr.plot(NSH+NSH.sam,weight="stock.wt",type="mort",rel=T))
plotTry("sf_observation_var_by_source", obsvar.plot(NSH.sam))
plotTry("sf_variance_vs_uncertainty",   obscv.plot(NSH.sam))
plotTry("sf_otolith",                   {set.seed(1); otolith(NSH.sam,n=1000)})
plotTry("sf_components",                comp.plot(NSH.sam))
plotTry("sf_retrospective_stock",       plot(NSH.retro))
plotTry("sf_retrospective_SAM_parameters", retroParams(NSH.retro),width=2000,height=1600)
plotTry("sf_retrospective_selectivity", retroSelectivity(NSH.retro,(range(NSH)["maxyear"]-10):range(NSH)["maxyear"]),width=2000,height=1600)
plotTry("sf_retrospective_residuals_HERAS", retroResiduals(NSH.retro,"HERAS",(range(NSH)["maxyear"]-10):range(NSH)["maxyear"]),width=2000,height=1600)
plotTry("mf_stock_trajectory",          plot(NSH3f.sam))
plotTry("mf_partial_F",                 fpart.plot(NSH3f.sam))
plotTry("mf_cor_params",                cor.plot(NSH3f.sam),width=2000,height=2000)
plotTry("mf_observation_var_by_source", obsvar.plot(NSH3f.sam))
plotTry("mf_retrospective_stock",       plot(NSH3f.retro))
plotTry("mf_retrospective_SAM_parameters", retroParams(NSH3f.retro),width=2000,height=1600)

#- Mohn's rho as calculated by this version
mohn <- lapply(list(sf=NSH.retro,mf=NSH3f.retro),function(r)
          lapply(c(ssb="ssb",fbar="fbar",rec="rec"),function(tp)
            mohns.rho(r,ref.year=max(as.numeric(names(r))),span=nRetro,type=tp)))

#===============================================================================
# Numbers for the comparison
#===============================================================================
fq  <- function(x) as.data.frame(x)[,c("age","year","area","data")]
ext <- function(s) list(nlogl=s@nlogl,nopar=s@nopar,params=s@params,stock.n=fq(s@stock.n),
                        harvest=fq(s@harvest),residuals=s@residuals,ssb=ssb(s),fbar=fbar(s),
                        rec=rec(s),catchabilities=catchabilities(s),obs.var=obs.var(s))
res <- list(version=as.character(packageVersion("FLSAM",lib.loc=lib)),timing=unlist(timing),
            sf=ext(NSH.sam),mf=ext(NSH3f.sam),
            sf.retro=lapply(NSH.retro,function(x) list(nlogl=x@nlogl,ssb=ssb(x),fbar=fbar(x),rec=rec(x))),
            mf.retro=lapply(NSH3f.retro,function(x) list(nlogl=x@nlogl,ssb=ssb(x),fbar=fbar(x),rec=rec(x))),
            mohn=mohn)
if(reuse) res$timing <- readRDS(file.path(outDir,"results.rds"))$timing
saveRDS(res,file.path(outDir,"results.rds"))
if(!reuse) save(NSH,NSH.tun,NSH.ctrl,NSH.sam,NSH.retro,NSH3f,NSH3.ctrl,NSH3f.sam,NSH3f.retro,file=file.path(outDir,"assessment.RData"))
cat("done",label,"\n"); print(unlist(timing))
