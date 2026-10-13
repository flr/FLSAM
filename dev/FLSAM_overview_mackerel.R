#===============================================================================
# FLSAM 3.1 - overview of the functions, using the built-in mackerel objects
#
# Data (WGWIDE 2026 NEA mackerel assessment, model 0.01; 1998 onwards):
#   MAC       FLStock, catches 1998-2025, biology to 2026
#   MAC.tun   FLIndices: Egg-survey (SSB index), IESSNS, tagAsIndex, rec-IESSNS
#   MAC.ctrl  FLSAM.control incl. the steel-tag recaptures (MAC.ctrl@recapture)
#   MAC.cvars FLQuants: fixed catch variances (weights) for 1998-1999
#   MAC.sam   FLSAM, saved without the stockassessment fit
#
# This example shows what the herring example (FLSAM_overview.R) cannot:
# tag-recapture data (fleet type 5), fixed catch variances (catch.vars),
# AR(1) correlated observations, a biomass (SSB) index, AR(1) correlated F,
# and a retrospective with tag data and catch variances.
# One fit takes 1-2 minutes; running everything about 20-30 minutes.
# NEW = new in FLSAM 3.1, CHANGED = changed behaviour in 3.1.
#===============================================================================

library(FLSAM)
library(ggplot2)

data(MAC)             # MAC, MAC.tun, MAC.ctrl, MAC.cvars
data(MAC.sam)         # MAC.sam

#- The stored assessment can be used directly for the FLR-based output
MAC.sam
plot(MAC.sam)
tag.pars(MAC.sam)                          # NEW: release survival and recapture overdispersion

#===============================================================================
# 1. Data and configuration
#===============================================================================
summary(MAC)
names(MAC.tun)
sapply(MAC.tun,function(x) x@type)         # "biomass": Egg-survey is an SSB index (age "all")
MAC.cvars$residual[,ac(1998:2001)]         # weights used as fixed variances (fixVarToWeight)

#- NEW: tag-recapture data (SAM fleet type 5) is part of the control object
head(MAC.ctrl@recapture)                   # ReleaseY, Yearclass, Nscan, R, Type, RecaptureY, ...
nrow(MAC.ctrl@recapture)
MAC.ctrl@fleets                            # Recaptures = fleet type 5, always the last fleet

#- Building the control from scratch: the recaptures and catch variances are arguments
ctrl <- FLSAM.control(MAC,MAC.tun,catch.vars=MAC.cvars,recapture=MAC.ctrl@recapture)
ctrl@fleets

#- The configuration of model 0.01 (see data-raw/MAC.R)
MAC.ctrl@states                            # F states, plus group bound to age 11
MAC.ctrl@cor.F                             # 2 = AR(1) correlation of the F increments
MAC.ctrl@obs.vars
MAC.ctrl@cor.obs.Flag                      # AR(1) correlated ages for IESSNS and tagAsIndex
MAC.ctrl@cor.obs
MAC.ctrl@fixVarToWeight                    # NEW: per fleet
MAC.ctrl@biomassTreat                      # Egg-survey: 0 = SSB index

#- Translation between FLR and SAM formats (the recaptures are added as fleet type 5)
dat  <- FLSAM2SAM(FLStocks(residual=MAC),MAC.tun,MAC.ctrl@sumFleets,MAC.cvars,MAC.ctrl@recapture)
table(dat$fleetTypes[dat$aux[,"fleet"]])   # observations per fleet type
conf <- ctrl2conf(MAC.ctrl,dat)
conf$keyCorObs
ctrl.back <- conf2ctrl(conf,dat)           # CHANGED: also rebuilds the recapture data
nrow(ctrl.back@recapture)

#- Dropping a fleet: the recapture data is dropped with the Recaptures fleet
ctrl.notag <- drop.from.control(MAC.ctrl,fleets="Recaptures")
ctrl.notag@fleets

#===============================================================================
# 2. Running the assessment
#===============================================================================
#- NEW: starting values from the stored assessment -> fast convergence
#  (catch.vars are needed for every fit with fixed catch variances)
sam <- FLSAM(MAC,MAC.tun,MAC.ctrl,catch.vars=MAC.cvars,starting.values=MAC.sam)
c(stored=nlogl(MAC.sam),refit=nlogl(sam))  # 585.6876 as in the WGWIDE assessment
fit <- samfit(sam)                         # NEW: the stockassessment fit is kept
stockassessment::modeltable(fit)

#- Without tag data, for comparison (different data: no likelihood-ratio test)
sam.notag <- FLSAM(MAC,MAC.tun,ctrl.notag,catch.vars=MAC.cvars)
plot(FLSAMs(with.tags=sam,without.tags=sam.notag),ci=TRUE)

#- Fix parameters, e.g. no correlation between the F increments of the ages
#  (cor.F = SAM itrans_rho; 0 on the transformed scale is a correlation of 0)
setParsTable()
subset(parameter.table(sam),name=="itrans_rho")
sp <- data.frame(par="cor.F",number=0,value=0)
sam.rho0 <- FLSAM(MAC,MAC.tun,MAC.ctrl,catch.vars=MAC.cvars,set.pars=sp,starting.values=sam)
lr.test(sam.rho0,sam)

#===============================================================================
# 3. Estimates
#===============================================================================
head(ssb(sam)); head(fbar(sam)); head(rec(sam))
catchabilities(sam)                        # Egg-survey: one catchability for the SSB index
obs.var(sam)
cor.obs(sam)                               # AR(1) correlation parameters
obscov(sam)[["IESSNS"]]                    # resulting observation covariance between ages
cor.F(sam)
tag.pars(sam)                              # NEW
parameter.table(sam)                       # NEW: all parameters with readable names
MAC.stk <- MAC + sam                       # estimates into the FLStock

#===============================================================================
# 4. Plots and diagnostics
#===============================================================================
plot(sam)
plot(MAC,sam)
tags.plot(sam)                             # NEW: observed and predicted recaptures
params.plot(sam)                           # NEW
cor.plot(sam)                              # parameter correlations (readable names)
ssbfbar.ellipse(sam,years=2021:2025)       # NEW: analytical SSB-Fbar ellipses
otolith(sam,year=2025,n=1000)
obsvar.plot(sam)
obscv.plot(sam)
obscor.plot(sam)                           # correlated observations (IESSNS, tagAsIndex)
obscor.empirical.plot(sam)
obssd.plot(sam)
data.plot(sam)
fselectivity.plot(sam,period=5)
procerr.plot(MAC+sam,weight="stock.wt",type="n",rel=TRUE)

#- One-step-ahead residuals; the recaptures get discrete OSA residuals (NEW)
res <- residuals(sam)
aggregate(std.res~fleet,data=res,FUN=function(x) round(c(n=length(x),mean=mean(x),sd=sd(x)),2))
head(subset(res,fleet=="Recaptures"))      # RecaptureY identifies the recapture year
rt <- residual.tests(sam)
subset(rt,p<0.05)
residual.diagnostics(sam)
model.description(sam)

jt <- jitter.test(sam,nojit=3,seed=1)
attr(jt,"table")

#===============================================================================
# 5. Retrospective with tag data and fixed catch variances
#===============================================================================
#  NEW: catch.vars are windowed per peel; recaptures are trimmed by recapture year
#  and release year; each peel starts from the previous one (use.start=TRUE)
MAC.ctrl.r <- MAC.ctrl; MAC.ctrl.r@residuals <- FALSE    # faster
MAC.retro <- retro(MAC,MAC.tun,MAC.ctrl.r,retro=3,catch.vars=MAC.cvars)
sapply(MAC.retro,function(x) nrow(x@control@recapture))  # tag records per peel
plot(MAC.retro)
plot(MAC.retro,ci=TRUE)
mohns.rho(MAC.retro,span=3,type="ssb")
attr(mohns.rho(MAC.retro,span=3,type="fbar"),"mohn")
retroParams(MAC.retro)
params.plot(MAC.retro)
retroSelectivity(MAC.retro,2019:2025)
ssbfbar.ellipse(MAC.retro,years=2022)

#- Leave one survey out (CHANGED: catch.vars are passed on)
MAC.loo <- loo(MAC,MAC.tun,MAC.ctrl.r,catch.vars=MAC.cvars)
plot(MAC.loo)

#===============================================================================
# 6. Forecasts and reference points
#===============================================================================
fc <- FLSAM.forecast(sam,fscale=c(1,1,1),ave.years=2023:2025,rec.years=2010:2025)
fc
FLSAM.refpts(sam,c("0.1dYPR","0.35SPR"))   # Fmax is not defined for mackerel
yp <- FLSAM.ypr(sam)
attr(yp,"refpts")
FLSAM.catchbyfleet(sam,obs.show=TRUE)[1:5,]

#===============================================================================
# 7. Output table
#===============================================================================
out <- FLSAM.out(MAC,MAC.tun,sam,format="TABLE %i.")
head(out,40)
