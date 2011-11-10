.onLoad <- function(lib,pkg) {
   require(methods)
   cat(sprintf("FLSAM %s development version (%s)\n",
          packageDescription("FLSAM")$Version,
          packageDescription("FLSAM")$Date))
   cat("
FLSAM v 0.40 (2011-11-09)
--------------------------------------------------------------------------------
* This package is an FLR wrapper to SAM, the State-space Assessment Model
    developed by Anders Nielsen <an@aqua.dtu.dk>. It is intended to be used
    only in close collaboration with experts in stock assessment modelling.
    Correct model specification, thorough validation of model results, and model
    diagnostics are required to obtain valid output.\n

* SAM and FLSAM come with ABSOLUTELY NO WARRANTY.
    For details regarding SAM, please consult <http://www.stockassessment.org>.
    For details and support regarding the FLSAM package, please consult the HAWG
    repository webpage, <http://hawg.googlecode.com>, where FLSAM is hosted, or
    the FLSAM users mailing list <http://groups.google.com/group/FLSAM>.\n

* If a 'bug' occurs, that is suspected to be caused by the assessment model
    itself (and not the FLR-interface), then please verify that the current
    official version of SAM available at <http://www.stockassessment.org>
    contains the bug before contacting the author with a bug report.\n

* The State-space Assessment Model (SAM) is constantly developing, so not
    all features that are available in SAM are available through FLSAM.
    Please consult http://www.stockassessment.org or the authors if your
    fish stock require features not covered.
--------------------------------------------------------------------------------")
#  cat("---------------------------------------------------------\n")
#  cat("* Note that FLR packages are under constant development,\n") 
#  cat("    please report bugs to flr-team@flr-project.org.\n")
#  cat("* New documentation can be found in vignetes,\n") 
#  cat("    run vignette(package=\"FLCore\") to get a list.\n")
#  cat("* For more information go to http:\\\\flr-project.org or\n")
#  cat("    subscribe the mailing list flr-list@flr-project.org.\n")
#  cat("---------------------------------------------------------\n")
}


