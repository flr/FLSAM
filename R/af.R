#' Conversion helpers
#'
#' \code{af} converts to a factor, \code{anf} converts a factor to numeric.
#'
#' @param x A vector.
#' @return A factor or numeric vector.
#' @keywords internal
#' @rdname af
#' @export
`af` <-
function(x){return(as.factor(x))}

#hello world
