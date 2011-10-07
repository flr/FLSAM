setClass("FLSAM",
  	representation(
      "FLComp",
      nopar    = "integer",
      nlogl    = "numeric",
      maxgrad  = "numeric",
      logDetHess = "numeric",
      params   = "data.frame",
      stock.n  = "FLQuant",
      harvest  ="FLQuant",
      fit      = "data.frame"),
  	prototype=prototype(),
  	validity=function(object){
                	# Everything is fine
                	return(TRUE)}
)

