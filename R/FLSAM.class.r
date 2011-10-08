setClass("FLSAM",
  	representation(
      "FLComp",
      control  = "FLSAM.control",
      nopar    = "integer",
      nlogl    = "numeric",
      params   = "data.frame",
      stock.n  = "FLQuantPoint",
      harvest  ="FLQuantPoint",
      residuals = "data.frame"),
  	prototype=prototype(),
  	validity=function(object){
                	# Everything is fine
                	return(TRUE)}
)

