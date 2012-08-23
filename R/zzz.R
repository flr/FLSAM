.onLoad <- function(lib,pkg) {
packageStartupMessage(sprintf("\nFLSAM %s (%s)",packageDescription("FLSAM")$Version,
           packageDescription("FLSAM")$Date))
dis.msg <- '--------------------------------------------------------------------------------
* This package is an FLR wrapper to SAM, the State-space Assessment Model
    developed by Anders Nielsen, DTU-Aqua. For details regarding SAM, 
    please consult <http://www.stockassessment.org>.

* SAM and FLSAM come with ABSOLUTELY NO WARRANTY.

* Please report all bugs, questions and feature requestions to the FLR
    mailing list, flr-list@flr-project.org. We´re here to help! 

* For more information on learning to use FLSAM, run vignette("FLSAM"). For 
    the full manual as a PDF file, run vignette("FLSAM-manual"). Full
    documentation is also availabe through the standard R help system - use
    help(package="FLSAM") to get started there.
--------------------------------------------------------------------------------
'
packageStartupMessage(dis.msg)
}


