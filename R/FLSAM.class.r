setClass("FLSAM",
  	representation(
      "FLComp",
      control  = "FLSAM.control",
      nopar    = "integer",
      nlogl    = "numeric",
      params   = "data.frame",
      stock.n  = "FLQuant",
      harvest  = "FLQuant",
      residuals = "data.frame"),
  	prototype=prototype(),
  	validity=function(object){
                	# Everything is fine
                	return(TRUE)}
)

