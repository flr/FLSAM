#-------------------------------------------------------------------------------
# Test of tag-recapture data in FLSAM with the WGWIDE 2026 NEA mackerel assessment
# (C:/Repositories/2026_mac.27.nea_assessment, model 0.01).
# The ICES input files are converted to FLR objects, the configuration of the
# stockassessment fit is translated with conf2ctrl(), and the FLSAM fit is
# compared with the original sam.fit.
#-------------------------------------------------------------------------------
lib.new  <- "C:/Users/NielsHintzen/Repositories/FLSAM_libs/new"
assDir   <- "C:/Repositories/2026_mac.27.nea_assessment"
dataDir  <- file.path(assDir,"boot/data/input_data_2026WGWIDE")
outDir   <- "C:/Users/NielsHintzen/Repositories/FLSAM/dev/mac27/output"
dir.create(outDir,showWarnings=FALSE,recursive=TRUE)

suppressMessages({library(FLSAM,lib.loc=lib.new); library(ggplot2)})

#- Original stockassessment fit of model 0.01 ------------------------------------
load(file.path(assDir,"model","fit_model_0.01.RData"))      # fit
fit.orig <- fit

#- ICES input files (as data0.01.r) --------------------------------------------------
rd <- function(f) read.ices(file.path(dataDir,f))
cn <- rd("cn.dat"); cw <- rd("cw.dat"); dw <- rd("dw.dat"); lf <- rd("lf.dat"); lw <- rd("lw.dat")
mo <- rd("mo.dat"); nm <- rd("nm.dat"); pf <- rd("pf.dat"); pm <- rd("pm.dat"); sw <- rd("sw.dat")
surveys <- rd("survey.dat")
recap <- read.table(file.path(dataDir,"tag_steel.dat"),header=TRUE)
recap <- recap[recap$Type==1 & recap$RecaptureY<=2006,]
recap <- recap[recap[,1]>=min(as.numeric(rownames(sw))),]
age   <- recap$ReleaseY - recap$Yearclass
recap <- recap[age>1 & age<12,]
recap <- recap[recap$ReleaseY!=recap$RecaptureY,]
recap <- recap[recap$RecaptureY!=2025,]
recap$R <- as.numeric(recap$R)

#- FLStock (years of the biology, catches NA where not available) ------------------------
yrs  <- rownames(sw); ages <- colnames(sw)
fq   <- function(m){
  a <- array(NA,dim=c(length(ages),length(yrs)),dimnames=list(age=ages,year=yrs))
  mm <- t(m); a[rownames(mm),colnames(mm)] <- mm
  FLQuant(a,dimnames=list(age=ages,year=yrs))
}
MAC <- FLStock(catch.n=fq(cn))
MAC@catch.wt <- fq(cw); MAC@landings.wt <- fq(lw); MAC@discards.wt <- fq(dw)
MAC@landings.n <- MAC@catch.n * fq(lf); MAC@discards.n <- MAC@catch.n * (1-fq(lf))
MAC@stock.wt <- fq(sw); MAC@mat <- fq(mo); MAC@m <- fq(nm)
MAC@harvest.spwn <- fq(pf); MAC@m.spwn <- fq(pm)
MAC@range[c("plusgroup","minfbar","maxfbar")] <- c(max(as.numeric(ages)),4,8)
name(MAC) <- "NEA mackerel"
#- catch weights and landing fractions of the last (survey-only) year as in the year before, as FLSAM does when extending a stock
for(sl in c("catch.wt","landings.wt","discards.wt")) slot(MAC,sl)[,ac(2026)] <- slot(MAC,sl)[,ac(2025)]

#- Catch weights (fixed variances before 2000, as data0.01.r) -------------------------------
W <- MAC@catch.n; W[] <- NA; W[,ac(1998:1999)] <- 10
catch.vars <- FLQuants(residual=W)

#- FLIndices --------------------------------------------------------------------------------
MAC.tun <- FLIndices(lapply(names(surveys),function(nm){
  m   <- surveys[[nm]]
  bio <- all(colnames(m)=="-1")
  idx <- FLQuant(t(m),dimnames=list(age=if(bio) "all" else colnames(m),year=rownames(m)))
  x   <- FLIndex(index=idx,name=nm)
  x@type <- ifelse(bio,"biomass","number")
  x@range[c("startf","endf")] <- attr(m,"time")
  x
}))
names(MAC.tun) <- names(surveys)

#- SAM data from the FLR objects, compared with the original data --------------------------
dataS <- FLSAM2SAM(FLStocks(residual=MAC),MAC.tun,NULL,catch.vars,recap)
cat("\n== Data comparison FLSAM vs original sam.fit\n")
d0 <- fit.orig$data
o0 <- order(d0$aux[,"year"],d0$aux[,"fleet"],d0$aux[,"age"])
key  <- function(d) paste(d$aux[,"year"],d$aux[,"fleet"],d$aux[,"age"],if("RecaptureY" %in% colnames(d$aux)) d$aux[,"RecaptureY"] else "",sep="_")
k0 <- key(d0); k1 <- key(dataS)
cat("observations original:",length(k0)," FLSAM:",length(k1)," (FLSAM includes",sum(is.na(dataS$logobs)),"missing)\n")
m  <- match(k0,k1)
cat("all original observations found:",!any(is.na(m)),"\n")
cat("logobs identical:",isTRUE(all.equal(d0$logobs,dataS$logobs[m])),"\n")
cat("weights identical:",isTRUE(all.equal(d0$weight,dataS$weight[m])),"\n")
for(nn in c("stockMeanWeight","catchMeanWeight","natMor","propMat","propF","propM","landFrac")){
  a0 <- matrix(d0[[nn]],nrow=dim(d0[[nn]])[1],dimnames=list(dimnames(d0[[nn]])[[1]],NULL))
  a1 <- matrix(dataS[[nn]],nrow=dim(dataS[[nn]])[1],dimnames=list(dimnames(dataS[[nn]])[[1]],NULL))
  yy <- intersect(rownames(a0),rownames(a1))
  cat(sprintf("%-16s identical on common years (%s-%s): %s\n",nn,min(yy),max(yy),isTRUE(all.equal(unname(a0[yy,]),unname(a1[yy,])))))
}
for(nn in c("sampleTimes","fleetTypes"))
  cat(sprintf("%-16s identical: %s\n",nn,isTRUE(all.equal(unname(c(d0[[nn]])),unname(c(dataS[[nn]]))))))

#- Configuration from the original fit, and the FLSAM fit -------------------------------------
ctrl <- conf2ctrl(fit.orig$conf,dataS)
ctrl@name <- "NEA mackerel"
cat("\nTag records in the control:",nrow(ctrl@recapture),"\n")
ctrl@residuals <- FALSE
t1 <- system.time(sam <- FLSAM(MAC,MAC.tun,ctrl,catch.vars=catch.vars,silent=TRUE))["elapsed"]
cat("\n== Fit comparison (",round(t1),"s)\n")
cat("nlogl original:",fit.orig$opt$objective,"  FLSAM:",nlogl(sam),"\n")
ssb0 <- ssbtable(fit.orig)[,1]; ssb1 <- ssb(sam)$value
cat("max rel. difference SSB:",max(abs(ssb1/ssb0-1)),"\n")
fb0  <- fbartable(fit.orig)[,1]; fb1 <- fbar(sam)$value
cat("max rel. difference Fbar:",max(abs(fb1/fb0-1)),"\n")
print(tag.pars(sam))

#- Tag plot and residuals (discrete OSA residuals for the recaptures) --------------------------
png(file.path(outDir,"mac_tags.png"),2400,1600,res=150); print(tags.plot(sam)); dev.off()
png(file.path(outDir,"mac_stock.png"),1600,1200,res=150); print(plot(sam)); dev.off()
ctrl@residuals <- TRUE
t2 <- system.time(sam.res <- FLSAM(MAC,MAC.tun,ctrl,catch.vars=catch.vars,starting.values=sam,silent=TRUE))["elapsed"]
cat("\nfit with residuals:",round(t2),"s; residuals by fleet:\n")
print(aggregate(std.res~fleet,data=residuals(sam.res),FUN=function(x) c(n=length(x),mean=mean(x),sd=sd(x))))

#- Retrospective with tag data (2 peels) -------------------------------------------------------
ctrl@residuals <- FALSE
r <- retro(MAC,MAC.tun,ctrl,retro=2,catch.vars=catch.vars)
cat("\nretro peels:",names(r),"\n")
#- the same peels with stockassessment::retro. Its year matrix gives per fleet the FIRST
#  year that is removed (catch, egg, IESSNS, tagAsIndex, rec-IESSNS, recaptures by release year)
Ryear  <- rbind(c(2026,2026,2026,2025,2026,2025),
                c(2025,2025,2025,2025,2025,2024))
ret.sa <- stockassessment::retro(fit.orig,year=Ryear,ncores=1)
print(data.frame(peel=names(r),
                 nlogl.FLSAM=sapply(r,nlogl),
                 nlogl.stockassessment=c(fit.orig$opt$objective,sapply(ret.sa,function(x) x$opt$objective)),
                 tags=sapply(r,function(x) nrow(x@control@recapture))))
png(file.path(outDir,"mac_retro.png"),1600,1200,res=150); print(plot(r,ci=TRUE)); dev.off()
save(MAC,MAC.tun,ctrl,catch.vars,sam,sam.res,r,file=file.path(outDir,"mac_flsam.RData"))
cat("\nTEST COMPLETED\n")
