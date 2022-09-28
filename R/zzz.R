.onAttach <- function(lib,pkg) {

packageStartupMessage(sprintf("\nFLSAM %s (%s)",packageDescription("FLSAM")$Version,
           packageDescription("FLSAM")$Date), "\n",
"------------------------------------------------------------------------------

* NOTE: This version of FLSAM requires a particular version of the
    'stockassessment' package available as package 'stockassessmentComp'.
    Both 'stockassessment' and 'stockassessmentComp' can be installed in
    the same machine, but should not be loaded in the same R session.

* This package is an FLR wrapper to SAM, the State-space Assessment Model
    developed by Anders Nielsen, DTU-Aqua. For details regarding SAM, 
    please consult <http://www.stockassessment.org>.

* SAM and FLSAM come with ABSOLUTELY NO WARRANTY.

* Please report all bugs, questions and feature requestions to the FLR
    mailing list, flr-list@flr-project.org. We're here to help! 

* For more information on learning to use FLSAM, run vignette('FLSAM'). For 
    the full manual as a PDF file, run vignette('FLSAM-manual'). Full
    documentation is also availabe through the standard R help system - use
    help(package='FLSAM') to get started there.
------------------------------------------------------------------------------")


if("stockassessment" %in% (.packages())){
  warning("'stockassessment' already attached, functions will be masked.")
}

if(packageVersion("stockassessmentComp") != "0.5.4")
  stop("Version of stockassessment is not compatible with FLSAM, see above.")
}
