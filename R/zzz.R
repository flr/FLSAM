.onLoad <- function(lib,pkg) {
packageStartupMessage(sprintf("\nFLSAM %s (%s)",packageDescription("FLSAM")$Version,
           packageDescription("FLSAM")$Date))
dis.msg <- '--------------------------------------------------------------------------------
* This package is an FLR wrapper to SAM, the State-space Assessment Model
    developed by Anders Nielsen <an@aqua.dtu.dk>. It is intended to be used
    only in close collaboration with experts in stock assessment modelling.
    Correct model specification, thorough validation of model results, and model
    diagnostics are required to obtain valid output.

* SAM and FLSAM come with ABSOLUTELY NO WARRANTY.

* For details regarding SAM, please consult <http://www.stockassessment.org>.
    For details and support regarding the FLSAM package, please consult the FLSAM
    wiki, <http://code.google.com/p/hawg/wiki/FLSAM>, where FLSAM is hosted, or
    the associated FLSAM users mailing list <http://groups.google.com/group/FLSAM>.

* If a "bug" occurs, that is suspected to be caused by the assessment model
    itself (and not the FLR-interface), then please verify that the current
    official version of SAM available at <http://www.stockassessment.org>
    contains the bug before contacting the author with a bug report.

* The State-space Assessment Model (SAM) is constantly developing, so not
    all features that are available in SAM are available through FLSAM.
    Please consult <http://www.stockassessment.org> or the authors if your
    fish stock requires features not covered.

* For help on FLSAM, try help(package="FLSAM"). For the full manual as a PDF file, 
  use vignette("FLSAM")
--------------------------------------------------------------------------------
'
packageStartupMessage(dis.msg)
}


