.onAttach <- function(lib,pkg) {

packageStartupMessage(sprintf("\nFLSAM %s (%s)",packageDescription("FLSAM")$Version,
           packageDescription("FLSAM")$Date), "\n",
"------------------------------------------------------------------------------

* This package is an FLR wrapper to SAM, the State-space Assessment Model
    developed by Anders Nielsen, DTU-Aqua, and uses the 'stockassessment'
    package (>= 0.12.0). Install or update it using
      remotes::install_github('fishfollower/SAM/stockassessment')
    For details regarding SAM, please consult <https://www.stockassessment.org>.

* SAM and FLSAM come with ABSOLUTELY NO WARRANTY.

* Please report bugs, questions and feature requests at
    <https://github.com/flr/FLSAM/issues>.

* Documentation is available through the standard R help system - use
    help(package='FLSAM') to get started.
------------------------------------------------------------------------------")

if(packageVersion("stockassessment") < "0.12.0")
  stop("Version of stockassessment is not compatible with FLSAM (>= 0.12.0 required), see above.")
}
