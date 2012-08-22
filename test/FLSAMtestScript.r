#-------------------------------------------------------------------------------
# In this document we will show the functionality and reporting possibilities of
# the FLR rapper to SAM: FLSAM. Please note that the FLSAM package is no more than just
# a simple wrapper to the ADMB stock assessment model SAM, developed by Anders
# Nielsen (DTU-Aqua) and allows you to use the benefits of FLR and R in a convenient way.
#-------------------------------------------------------------------------------

### ============================================================================
### Setup assessment
### ============================================================================


# First, load the FLSAM package and datasets. The datasets contain an already complete
# FLR stock object and an FLR tuning indices object.

library(FLSAM)
data(NSH)
data(NSH.sam)

# To be able to run the assessment, a control file, that specifies the assessment settings
# needs to be created. However, while loading data in the first step, a defaul control
# object is created already. To show you how it's done though, an example is given below,
# using the FLSAM.control function.

if(!exists("NSH.ctrl")) NSH.ctrl <- FLSAM.control(NSH,NSH.tun)

### ============================================================================
### Run the assessment
### ============================================================================

# Finally, all datasets and control objects are ready to perform the assessment, hence
# call the FLSAM function.

NSH.sam <- FLSAM(NSH,NSH.tun,NSH.ctrl)

# Give the object a name and update the stock object
name(NSH.sam) <- "North Sea Herring"
NSH     <- NSH + NSH.sam

### ============================================================================
### Investigate the results
### ============================================================================

pdf(file.path(tempdir(),paste(name(NSH.sam),".pdf",sep="")))

# Standard outputs
ssb(NSH.sam)
fbar(NSH.sam)
tsb(NSH.sam)
rec(NSH.sam)
f(NSH.sam)
n(NSH.sam)

# The residual diagnostic plots
residual.diagnostics(NSH.sam)

#- The sam object
print(plot(NSH.sam))

#- Plot the catchabilities of the surveys
catch <- catchabilities(NSH.sam)
print(xyplot(value+ubnd+lbnd ~ age | fleet,catch,
          scale=list(alternating=FALSE,y=list(relation="free")),as.table=TRUE,
          type="l",lwd=c(2,1,1),col=c("black","grey","grey"),
          subset=fleet %in% c("HERAS"),
          main="Survey catchability parameters",ylab="Catchability",xlab="Age"))

#- Plot the observation variances
obsvar.plot(NSH.sam)

#- Plot the observation CVs
obscv.plot(NSH.sam)

#- Plot the model uncertainty (with very small number of re-draws to make it a bit quicker)
otolith(NSH.sam,year=2011,plot=T,n=100)

#- Plot the correlation of the parameters
cor.plot(NSH.sam)

dev.off()

### ============================================================================
### Perform retro analyses and drop some surveys from the model to test performance
### ============================================================================

retro.NSH <- retro(NSH,NSH.tun,NSH.ctrl,5)

looi.NSH  <- looi(NSH,NSH.tun,NSH.ctrl,type="LOO")

# Perform AIC to see which one performs best
AIC(looi.NSH)

### ============================================================================
### Sample from the model uncertainty and construct new stock realisations
### ============================================================================

mc.NSH    <- monteCarloStock(NSH,NSH.sam,100)