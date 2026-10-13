#' Colour palette of the FLSAM plots
#'
#' Returns \code{n} colours for the ggplot2 plots of FLSAM. The default palette
#' is dark blue (\code{"#004868"}), orange (\code{"#FF800A"}), light blue
#' (\code{"#00BADC"}) and green (\code{"#31C381"}); more colours are
#' interpolated between these. Set your own palette with
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
  pal <- getOption("FLSAM.palette",c(darkblue="#004868",orange="#FF800A",lightblue="#00BADC",green="#31C381"))
  pal <- unname(pal)
  if(n<=length(pal)) return(pal[seq_len(n)])
  return(grDevices::colorRampPalette(pal)(n))
}

#- Colour of a single role in the palette
.col <- function(role=c("main","accent","light","extra")){
  FLSAM.colours(4)[match(match.arg(role),c("main","accent","light","extra"))]
}
