#-------------------------------------------------------------------------------
# MAC.R - North East Atlantic mackerel example data
#
# Built from the input data of the WGWIDE 2026 NEA mackerel assessment
# (ICES, mac.27.nea, model 0.01; repository 2026_mac.27.nea_assessment):
# catch-at-age and biology 1998-2025 (biology to 2026), egg survey (SSB index),
# IESSNS, tag-based index, IESSNS recruitment index, steel-tag recaptures, and
# fixed catch variances before 2000. The configuration is that of model 0.01.
#
# Creates data/MAC.RData (MAC, MAC.tun, MAC.ctrl, MAC.cvars) and
# data/MAC.sam.RData (MAC.sam).
#-------------------------------------------------------------------------------
library(FLSAM)

assDir  <- "C:/Repositories/2026_mac.27.nea_assessment"
dataDir <- file.path(assDir,"boot/data/input_data_2026WGWIDE")

#- ICES input files (as data0.01.r of the assessment) -----------------------------
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
rownames(recap) <- NULL

#- FLStock: biology to 2026 (survey year), catches to 2025 --------------------------
yrs  <- rownames(sw); ages <- colnames(sw)
fq   <- function(m){
  a  <- array(NA,dim=c(length(ages),length(yrs)),dimnames=list(age=ages,year=yrs))
  mm <- t(m); a[rownames(mm),colnames(mm)] <- mm
  FLQuant(a,dimnames=list(age=ages,year=yrs))
}
MAC <- FLStock(catch.n=fq(cn))
MAC@catch.wt <- fq(cw); MAC@landings.wt <- fq(lw); MAC@discards.wt <- fq(dw)
MAC@landings.n <- MAC@catch.n * fq(lf); MAC@discards.n <- MAC@catch.n * (1-fq(lf))
MAC@stock.wt <- fq(sw); MAC@mat <- fq(mo); MAC@m <- fq(nm)
MAC@harvest.spwn <- fq(pf); MAC@m.spwn <- fq(pm)
for(sl in c("catch.wt","landings.wt","discards.wt")) slot(MAC,sl)[,"2026"] <- slot(MAC,sl)[,"2025"]
MAC@catch     <- computeCatch(MAC)
MAC@landings  <- computeLandings(MAC)
MAC@discards  <- computeDiscards(MAC)
MAC@range[c("plusgroup","minfbar","maxfbar")] <- c(12,4,8)
name(MAC) <- "NEA mackerel"
desc(MAC) <- "WGWIDE 2026 NEA mackerel assessment input data (model 0.01)"

#- Catch weights: fixed variance (weight 10) before 2000 --------------------------------
W <- MAC@catch.n; W[] <- NA; W[,ac(1998:1999)] <- 10
MAC.cvars <- FLQuants(residual=W)

#- FLIndices -----------------------------------------------------------------------------
MAC.tun <- FLIndices(lapply(names(surveys),function(nm){
  m   <- surveys[[nm]]
  bio <- all(colnames(m)=="-1")
  idx <- FLQuant(t(m),dimnames=list(age=if(bio) "all" else colnames(m),year=rownames(m)))
  x   <- FLIndex(index=idx,name=nm)
  x@type <- ifelse(bio,"biomass","number")
  x@range[c("startf","endf")] <- attr(m,"time")
  x@range["plusgroup"] <- NA                                # no plus groups in the surveys
  x
}))
names(MAC.tun) <- names(surveys)

#- Configuration of model 0.01 ---------------------------------------------------------------
MAC.ctrl <- FLSAM.control(MAC,MAC.tun,catch.vars=MAC.cvars,recapture=recap)
MAC.ctrl@name <- "NEA mackerel"
MAC.ctrl@fixVarToWeight[] <- TRUE                          # catch weights are fixed variances
MAC.ctrl@states["catch unique",]   <- c(0:9,9)             # F states, plus group bound to age 11
MAC.ctrl@f.vars["catch unique",]   <- c(0,0,1,1,1,1,1,1,1,2,2)
MAC.ctrl@cor.F                     <- 2                    # AR(1) correlation of F increments
MAC.ctrl@obs.vars["catch unique",] <- c(0,0,1,1,1,1,1,1,1,2,2)
MAC.ctrl@obs.vars["Egg-survey",1]  <- 3
MAC.ctrl@obs.vars["IESSNS",]       <- c(-1,4,5,5,5,5,5,5,6,6,-1)
MAC.ctrl@obs.vars["tagAsIndex",]   <- c(7,8,9,9,9,9,9,9,9,9,-1)
MAC.ctrl@obs.vars["rec-IESSNS",1]  <- 10
MAC.ctrl@cor.obs.Flag[names(MAC.ctrl@fleets) %in% c("IESSNS","tagAsIndex")] <- "AR"
MAC.ctrl@cor.obs["IESSNS",]        <- c(-1,0,0,0,0,0,0,0,0,-1)
MAC.ctrl@cor.obs["tagAsIndex",]    <- c(1,1,1,1,1,1,1,1,1,-1)
MAC.ctrl@residuals <- TRUE
MAC.ctrl <- update(MAC.ctrl)

save(MAC,MAC.tun,MAC.ctrl,MAC.cvars,file="../data/MAC.RData",compress="xz")

#- Assessment (stored without the stockassessment fit to keep the package small) --------------
MAC.sam <- FLSAM(MAC,MAC.tun,MAC.ctrl,catch.vars=MAC.cvars,keep.fit=FALSE)
save(MAC.sam,file="../data/MAC.sam.RData",compress="xz")
