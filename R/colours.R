#' Colour palette of the FLSAM plots
#'
#' Returns \code{n} distinct (qualitative) colours for the ggplot2 plots of
#' FLSAM. The palette starts with dark blue (\code{"#004868"}), orange
#' (\code{"#FF800A"}), light blue (\code{"#00BADC"}) and green
#' (\code{"#31C381"}), followed by seven contrasting colours. When more colours
#' are needed than the palette holds, the colours are recycled and the plots
#' distinguish the series by line type. Set your own palette with
#' \code{options(FLSAM.palette=c(...))}.
#'
#' @param n Number of colours.
#' @return A character vector of colours.
#' @examples
#' FLSAM.colours(4)
#' FLSAM.colours(11)
#' \dontrun{options(FLSAM.palette=c("black","red","blue","darkgreen"))}
#' @export
FLSAM.colours <- function(n=4){
  pal <- getOption("FLSAM.palette",c(darkblue="#004868",orange="#FF800A",lightblue="#00BADC",green="#31C381",
                                     pink="#CC79A7",purple="#7F3C8D",red="#D55E00",yellow="#E6C200",
                                     black="#000000",grey="#8C8C8C",brown="#8C510A"))
  pal <- unname(pal)
  return(rep(pal,length.out=n))
}

#- Colour of a single role in the palette
.col <- function(role=c("main","accent","light","extra")){
  FLSAM.colours(4)[match(match.arg(role),c("main","accent","light","extra"))]
}

#- Matching colour and linetype scales for a discrete variable: colours are
#  recycled when there are more levels than colours, then the linetype changes
.discreteScales <- function(n,pal=NULL,name=waiver()){
  if(is.null(pal)) pal <- unname(getOption("FLSAM.palette",FLSAM.colours(11)))
  nc  <- length(pal)
  list(scale_colour_manual(values=rep(pal,length.out=n),name=name),
       scale_linetype_manual(values=rep(c("solid","dashed","dotted","dotdash"),each=nc,length.out=n),name=name))
}

#- Colours of the retrospective selectivity plot (as in FLSAM <= 3.0, FLCore lattice theme)
.retroColours <- c("#0072B2","#E69F00","#009E73","#D55E00","#56B4E9","#F0E442","#CC79A7")
