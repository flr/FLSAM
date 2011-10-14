FLSAM <-function(stck,tun,ctrl,temp.dir=tempdir(),batch.mode=FALSE) {
  #---------------------------------------------------
  # Output FLR objects into a format for SAM to read
  #---------------------------------------------------
  #General Setup
  admb.stem <- "ssass" 
  
  #Write output files
  FLR2SAM(stck,tun,ctrl,temp.dir)

  #---------------------------------------------------
  # We're ready! Run the executable
  #---------------------------------------------------
  admb.args <-  "-nr 2 -noinit -iprint 1"
  #Platform specific issues
  if (.Platform$OS.type=="unix") {
    admb.exec <- file.path(system.file("bin", "linux", package="FLSAM",
                   mustWork=TRUE), admb.stem)
    file.copy(admb.exec, temp.dir)
  } else if (.Platform$OS.type == "windows") {
    admb.exec <- file.path(system.file("bin", "windows", package="FLSAM", mustWork=TRUE),
      sprintf("%s.exe",admb.stem))
    file.copy(admb.exec, temp.dir)
  } else {
    stop(sprintf("Platform type, %s, is not currently supported.",R.version$os))
  }

  #Run!
  cmd <- sprintf("./%s -nr 2 -noinit -iprint 1" ,basename(admb.exec))
  olddir <- setwd(temp.dir)
  rtn <- system(cmd)
  setwd(olddir)
  if(rtn!=0) {
    if(batch.mode) {
      return(NULL)
    } else {
      stop(sprintf("An error occurred while running ADMB. Return code %s.",rtn))
    }
  }

  #---------------------------------------------------
  # Now read the results from the assessment
  #---------------------------------------------------
  res <- SAM2FLR(ctrl,temp.dir,admb.stem)

  return(res)
}

