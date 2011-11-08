.onLoad <- function(lib,pkg) {
   require(methods)
   cat(sprintf("FLSAM %s development version (%s)\n",
          packageDescription("FLSAM")$Version,
          packageDescription("FLSAM")$Date))
   cat("<Insert Health Warning Here>\n")
#  cat("---------------------------------------------------------\n")
#  cat("* Note that FLR packages are under constant development,\n") 
#  cat("    please report bugs to flr-team@flr-project.org.\n")
#  cat("* New documentation can be found in vignetes,\n") 
#  cat("    run vignette(package=\"FLCore\") to get a list.\n")
#  cat("* For more information go to http:\\\\flr-project.org or\n")
#  cat("    subscribe the mailing list flr-list@flr-project.org.\n")
#  cat("---------------------------------------------------------\n")
}


