#===============================================================================
# FLSAM 3.1 - overview of all functions, using the built-in herring objects
#
# Data: NSH (FLStock), NSH.tun (FLIndices), NSH.ctrl (FLSAM.control),
#       NSH.sam (FLSAM, saved without the stockassessment fit) and
#       HERAS.sams (FLSAMs, 4 configurations of the HERAS catchabilities)
#
# Run section by section. Running everything takes about 10-15 minutes
# (retrospective, leave-one-out, jitter test and Monte Carlo runs).
# NEW = new in FLSAM 3.1, CHANGED = changed behaviour in 3.1.
#===============================================================================

library(FLSAM)        # loads FLCore and stockassessment
library(ggplot2)      # all ggplot-based FLSAM plots return ggplot objects

data(NSH)             # NSH, NSH.tun, NSH.ctrl
data(NSH.sam)         # NSH.sam
data(HERAS.sams)      # HERAS.sams

#- NEW: colour palette used by all FLSAM plots (PFA colours first, then 7 more)
FLSAM.colours(4)
FLSAM.colours(11)
# options(FLSAM.palette=c("black","red","blue","darkgreen"))   # use your own palette

#===============================================================================
# 1. Configuration: FLSAM.control
#===============================================================================

#- Default configuration from the data (first guess, based on stockassessment::defcon)
ctrl <- FLSAM.control(NSH,NSH.tun)
ctrl                                       # printing an FLSAM.control
slotNames(ctrl)

#- The configuration used for the North Sea herring example
NSH.ctrl@states                            # coupling of the F states
NSH.ctrl@catchabilities                    # coupling of the survey catchabilities
NSH.ctrl@obs.vars                          # coupling of the observation variances
NSH.ctrl@logN.vars                         # coupling of the log(N) process variances
NSH.ctrl@f.vars                            # coupling of the log(F) process variances
NSH.ctrl@cor.F                             # 0 independent, 1 compound symmetry, 2 AR(1)
NSH.ctrl@cor.obs.Flag                      # observation correlation: ID, AR or US
NSH.ctrl@srr                               # stock-recruitment model (0 = random walk)
NSH.ctrl@biomassTreat                      # treatment of biomass indices (SCAI: 0 = SSB index)

#- Bind parameters and renumber the couplings with update()
ctrl <- NSH.ctrl
ctrl@catchabilities["HERAS",ac(5:8)] <- 101
ctrl <- update(ctrl)                       # CHANGED: now also renumbers weight/maturity/M couplings
ctrl@catchabilities

#- Drop fleets or ages
drop.from.control(NSH.ctrl,fleets="IBTS0")
drop.from.control(NSH.ctrl,ages=8)

#- Validity check                          # CHANGED: crashed before
validObject(NSH.ctrl)

#- NEW: fixVarToWeight per fleet (index.var / catch.vars used as fixed variances)
NSH.ctrl@fixVarToWeight

#- NEW: hockeyStickCurve fixes the smoothness of the bent hyperbola (srr 63)
NSH.ctrl@hockeyStickCurve

#- NEW: objects saved with older FLSAM versions are upgraded automatically
#  when used (FLSAM, retro, printing, ...); or explicitly:
NSH.ctrl <- upgradeFLSAM(NSH.ctrl)
NSH.sam  <- upgradeFLSAM(NSH.sam)
HERAS.sams <- upgradeFLSAM(HERAS.sams)

#- Translation between FLR and SAM formats
dat  <- FLSAM2SAM(FLStocks(residual=NSH),NSH.tun,NSH.ctrl@sumFleets)   # SAM data list
conf <- ctrl2conf(NSH.ctrl,dat)                                       # SAM configuration list
str(conf,max.level=1)
ctrl.back <- conf2ctrl(conf,dat)                                      # and back

#===============================================================================
# 2. Running the assessment
#===============================================================================

#- Fit; the stockassessment fit is stored in the object (NEW: keep.fit=TRUE)
sam <- FLSAM(NSH,NSH.tun,NSH.ctrl)
sam                                        # printing an FLSAM
name(sam) <- "North Sea herring"

#- NEW: access the stockassessment fit and use any stockassessment function
fit <- samfit(sam)
class(fit)
stockassessment::modeltable(fit)
stockassessment::fbarplot(fit)

#- Return only the raw sam fit, and convert it afterwards
fit2 <- FLSAM(NSH,NSH.tun,NSH.ctrl,return.fit=TRUE)
sam2 <- SAM2FLR(fit2,NSH.ctrl)

#- Without storing the fit (smaller objects)
sam.small <- FLSAM(NSH,NSH.tun,NSH.ctrl,keep.fit=FALSE)

#- Arguments of stockassessment::sam.fit can be passed on
sam.quick <- FLSAM(NSH,NSH.tun,NSH.ctrl,newtonsteps=0,silent=TRUE)

#- NEW: fix parameters with set.pars (all fixed-effect parameters)
setParsTable()                             # FLSAM names that can be used; SAM names work as well
sp <- data.frame(par=c("obs.vars","cor.F"),number=c(1,0),value=c(log(0.3),2))
sam.fixed <- FLSAM(NSH,NSH.tun,NSH.ctrl,set.pars=sp)
c(free=nopar(sam),fixed=nopar(sam.fixed))

#- Fast fit for simulations, results in the FLStock   # CHANGED: works when surveys extend beyond the catches
NSH.fast <- FLSAM.fast(NSH,NSH.tun,NSH.ctrl)
NSH.fast[["residual"]]@stock.n[,ac(2009:2011)]

#===============================================================================
# 3. Slots of the FLSAM object
#===============================================================================
control(sam)
nopar(sam)
n.states(sam)
states(sam)
nlogl(sam)
vcov(sam)[1:5,1:5]
rescov(sam)[1:5,1:5]
obscov(sam)                                # observation covariances (only for AR/US structures)
head(params(sam))
head(residuals(sam))
info(sam)
stock.n(sam)[,ac(2009:2012)]
harvest(sam)[,ac(2009:2012)]
components(sam)                            # only for partial surveys (fleet type 6)

#- Replacement methods
sam.copy <- sam
nopar(sam.copy) <- 19
nlogl(sam.copy) <- 167.4
info(sam.copy)  <- info(sam)               # CHANGED: no longer corrupts the values

#===============================================================================
# 4. Estimates with confidence intervals
#===============================================================================
head(ssb(sam))
head(tsb(sam))
head(fbar(sam))
head(catch(sam))
head(rec(sam))
head(n(sam))
head(f(sam))
head(catch.n(sam))
catchabilities(sam)
obs.var(sam)
n.var(sam)
f.var(sam)
cor.F(sam)
try(power.law.exps(sam))                   # error: no power-law catchabilities in this model
try(cor.obs(sam))                          # error: no correlated observations in this model
try(catch.scale(sam))                      # error: no catch multipliers in this model

#- Stock-recruitment parameters             # CHANGED: no global 'ctrl', no duplicated rows
ctrl.bh <- NSH.ctrl; ctrl.bh@srr <- 2L
sam.bh  <- FLSAM(NSH,NSH.tun,ctrl.bh)
rec.pars(sam.bh)

#- Parameters
coef(sam)
coefficients(sam)
head(params(sam))

#- Model comparison
AIC(sam)
AIC(HERAS.sams)
lr.test(HERAS.sams)
lr.test(sam.fixed,sam)

#- FLSAMs: a list of FLSAM objects
sams <- FLSAMs(free=sam,fixed=sam.fixed)
ssb(sams)[1:5,]

#- Put the estimates into the FLStock
NSH.sam.stk <- NSH + sam                   # same as merge(NSH, sam)
NSH.sam.stk <- merge(NSH,sam)
NSH.stks    <- NSH + HERAS.sams            # FLStocks, one per assessment

#===============================================================================
# 5. Plots (ggplot2; CHANGED: were lattice)
#===============================================================================
plot(sam)                                  # SSB, F and recruitment with confidence band
plot(sam,sam.fixed)                        # two assessments
plot(HERAS.sams)                           # FLSAMs
plot(HERAS.sams,NSH)                       # with the FLStock
plot(NSH,sam)
cor.plot(sam)                              # parameter correlations
fpart.plot(sam)                            # partial Fbar per catch fleet
procerr.plot(NSH+sam,weight="stock.wt",type="n",rel=TRUE)   # process error in N
procerr.plot(NSH+sam,weight="stock.wt",type="mort",rel=TRUE)
procerr.plot(NSH+sam,weight="stock.wt",type="tsb")

#- Base-graphics plots
residual.diagnostics(sam)                  # 6 panels per fleet and age (CHANGED: loess smoother drawn)
obsvar.plot(sam)                           # observation variances by data source
obscv.plot(sam)                            # observation variance vs CV
try(obscor.plot(sam))                      # error: no correlated observations in this model
ctrl.ar <- NSH.ctrl; ctrl.ar@cor.obs.Flag[names(ctrl.ar@fleets)=="HERAS"] <- "AR"   # factor without names: use an index
ctrl.ar@cor.obs["HERAS",] <- c(-1,rep(0,7)); ctrl.ar <- update(ctrl.ar)
sam.ar  <- FLSAM(NSH,NSH.tun,ctrl.ar)      # AR(1) correlation between HERAS ages
obscor.plot(sam.ar)
cor.obs(sam.ar)
otolith(sam,year=2011,n=1000)              # joint uncertainty of SSB and Fbar
try(comp.plot(sam))                        # error: no partial surveys in this model

#- Plots are ggplot objects and can be modified
plot(sam) + theme_minimal() + ggtitle("North Sea herring")

#===============================================================================
# 6. NEW diagnostics (as in stockassessment)
#===============================================================================
#- Tests on the one-step-ahead residuals: bias, variance, autocorrelation in
#  time and age, mean-variance relationship and normality
rt <- residual.tests(sam)
subset(rt,p<0.05)

obscor.empirical.plot(sam)                 # empirical residual correlation between ages
obssd.plot(sam)                            # observation sd by fleet and age
obssd.plot(sam,show.rel.w=TRUE)            # relative weights
data.plot(sam)                             # data overview
fselectivity.plot(sam)                     # F / max(F) by age, coloured by period
fselectivity.plot(sam,period=5)
model.description(sam)                     # text description of the model

#- Jitter test: refit from perturbed starting values
jt <- jitter.test(sam,nojit=5,seed=1)
attr(jt,"table")
plot(jt)

#- Prediction-variance link: estimate a power relation between prediction and
#  observation variance for HERAS, and plot it
ctrl.pv <- NSH.ctrl
ctrl.pv@predVarObsLink["HERAS",ac(1:8)] <- 0
sam.pv  <- FLSAM(NSH,NSH.tun,ctrl.pv)
predsd.plot(sam.pv,fleet="HERAS")
lr.test(sam,sam.pv)

#===============================================================================
# 7. Retrospective and leave-one-out analyses
#===============================================================================
NSH.retro <- retro(NSH,NSH.tun,NSH.ctrl,retro=3)
plot(NSH.retro)
mohns.rho(NSH.retro,type="ssb")            # CHANGED: defaults work, sign as Mohn (1999)
attr(mohns.rho(NSH.retro,type="fbar"),"mohn")   # NEW: mean Mohn's rho
retroResiduals(NSH.retro,"HERAS",2005:2011)
retroResiduals(NSH.retro,"HERAS",2005:2011,cohort=TRUE)
retroSelectivity(NSH.retro,2006:2011)
retroParams(NSH.retro)

NSH.loo <- loo(NSH,NSH.tun,NSH.ctrl)       # leave one fleet out
plot(NSH.loo)
NSH.loi <- loi(NSH,NSH.tun,NSH.ctrl)       # keep one fleet in
ssb(NSH.loo)[1:5,]                         # (lr.test is not valid here: different data)
# NSH.looi <- looi(NSH,NSH.tun,NSH.ctrl)   # all combinations (CHANGED: works without 'type'; slow)

#===============================================================================
# 8. NEW forecasts, reference points and per-recruit analyses
#===============================================================================
fc <- FLSAM.forecast(sam,fscale=c(1,1,1),ave.years=2009:2011,rec.years=2000:2011)
fc
plot(attr(fc,"samforecast"))               # the stockassessment forecast plot
ggplot(fc,aes(year,value)) +
  geom_ribbon(aes(ymin=lbnd,ymax=ubnd),fill=FLSAM.colours(1),alpha=0.15) +
  geom_line(colour=FLSAM.colours(1)) +
  facet_wrap(~quantity,scales="free_y") + theme_bw()

mf <- FLSAM.modelforecast(sam,fscale=c(1,1,1))
head(mf)

FLSAM.refpts(sam,c("Max","0.1dYPR","0.35SPR"))   # per-recruit reference points (MSY needs srr != 0)
yp <- FLSAM.ypr(sam)
attr(yp,"refpts")
FLSAM.catchbyfleet(sam,obs.show=TRUE)[1:5,]

#===============================================================================
# 9. Simulation and MSE
#===============================================================================
#- Draw assessments from the estimated uncertainty
sims <- simulate(NSH,NSH.tun,NSH.ctrl,n=10,set.seed=1)   # CHANGED: set.pars and set.seed are optional now
#  the simulated objects hold N and F only (no derived SSB/Fbar): compute these via the FLStock
sims.stk <- NSH + sims                     # FLStocks with simulated stock.n and harvest
sims.ssb <- do.call(rbind,lapply(seq_along(sims.stk),function(i) data.frame(as.data.frame(ssb(sims.stk[[i]])),sim=i)))
ggplot(sims.ssb,aes(year,data,group=sim)) + geom_line(colour=FLSAM.colours(1),alpha=0.4) +
  labs(y="SSB",title="Simulated SSB") + theme_bw()

#- Monte Carlo: simulate new data, refit, return an FLStock with iterations
NSH.mc <- monteCarloStock(NSH,NSH.tun,sam,realisations=3,ncores=2)
NSH.mc@stock.n[,"2011"]
# monteCarloStockMscan(NSH,NSH.tun,sam,realisations=3,mscan=seq(-0.1,0.1,0.05),ncores=2)  # profile on M

#- One assessment per iteration (MSE)       # CHANGED: works again, optional parallel
NSH.iter <- propagate(NSH,3)
mse      <- FLSAM.MSE(NSH.iter,NSH.tun,NSH.ctrl,ncores=1)
mse[["residual"]]@stock.n[,"2011"]

#- Kobe plot data (plotkobe requires the 'kobe' package)
kb <- kobeFLSAM(sam,what=c("trks","pts"))
head(kb$trks)
# plotkobe(sam)

#===============================================================================
# 10. Output tables
#===============================================================================
out <- FLSAM.out(NSH,NSH.tun,sam,format="TABLE %i.")
head(out,40)
# write(out,file="NSH_sam.out")
out2 <- sam.out(NSH,NSH.tun,sam)            # alias
formatQuant(NSH@catch.n[,ac(2009:2011)])[1:5]
formatMatrix(NSH.ctrl@states)

#===============================================================================
# 11. Internal / legacy helpers (exported, rarely needed)
#===============================================================================
af(c("a","b"))                             # as.factor
anf(factor(c("1","2")))                    # factor -> numeric
diags(sam)[1:3,]                           # residual data for custom plots (needs plyr)
# SIM2FLR(), sam.fitfast(), clean.void.catches(), getLowerBounds(), getUpperBounds():
#   internal helpers; FLSAM2par() and FLSAM(..., starting.values=): legacy pin
#   functionality, not maintained
