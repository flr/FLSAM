monteCarloStock <- function(stck,sam,realisations,run.dir=tempdir()){

  require(MASS)
  ctrl              <- sam@control
  
  #-Create new stock object with nr of realisations in iter slot
  mcstck            <- propagate(stck,iter=realisations)
  mcstck            <- window(mcstck,start=range(sam)["minyear"],end=range(sam)["maxyear"])
  mcstck@stock.n[]  <- NA
  mcstck@harvest[]  <- NA

  #-Generate new parameter realisations from vcov
  random.param      <- mvrnorm(realisations,sam@params$value,sam@vcov)
  save(random.param,file=file.path(run.dir,"random.param.RData"))
  
  #-Control settings
  n.states          <- length(unique(ctrl@states[names(which(ctrl@fleets==0)),]))
  yrs               <- dims(sam)$minyear:dims(sam)$maxyear
  ages              <- dims(sam)$age
  
  #-Extract the state variables
  u                 <- random.param[,which(colnames(random.param)=="U")]
  ca                <- random.param[,which(colnames(random.param)=="logCatch")]
  idxNs             <- c(mapply(seq,from=seq(1,ncol(u),ages+n.states),
                                    to  =seq(1,ncol(u),ages+n.states)+ages-1,
                                    by  =1))
  idxFs             <- c(mapply(seq,from=seq(1,ncol(u),ages+n.states)+ages,
                                    to  =seq(1,ncol(u),ages+n.states)+n.states+ages-1,
                                    by  =1))
  mcstck@stock.n[]  <- exp(aperm(array(u[,idxNs],dim=c(realisations,ages,    length(yrs))),perm=c(2,3,1)))
  mcstck@harvest[]  <- exp(aperm(array(u[,idxFs],dim=c(realisations,n.states,length(yrs))),perm=c(2,3,1)))[ctrl@states[names(which(NSH.ctrl@fleets==0)),],,]
  mcstck@catch[]    <- exp(aperm(array(ca       ,dim=c(realisations,         length(yrs))),perm=c(2,1)))

  return(mcstck)}

