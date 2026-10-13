#' Example datasets
#'
#' FLR objects to run example assessments.
#' \itemize{
#'   \item \code{data(NSH)}: \code{NSH} (\code{FLStock}), \code{NSH.tun}
#'     (\code{FLIndices}) and \code{NSH.ctrl} (\code{\link{FLSAM.control-class}})
#'     for North Sea autumn spawning herring, truncated to 1992 onwards so
#'     that the model converges quickly. The indices are the larval index (SCAI),
#'     the acoustic survey (HERAS), IBTS-Q1 and IBTS0.
#'   \item \code{data(NSH.sam)}: \code{NSH.sam}, the \code{\link{FLSAM-class}}
#'     result of the assessment (without a stored fit, see \code{\link{samfit}}).
#'   \item \code{data(HERAS.sams)}: an \code{\link{FLSAMs}} object with four
#'     assessments that differ in the binding of the HERAS catchabilities
#'     ("all.free", "four.lvls", "three.lvls", "two.lvls"), e.g. for
#'     \code{\link{lr.test}}.
#'   \item \code{data(SOL)}, \code{data(SOL.sam)}: \code{SOL}, \code{SOL.tun},
#'     \code{SOL.ctrl} and \code{SOL.sam}, a sole example with the
#'     corresponding assessment.
#' }
#' Objects created with FLSAM versions before 3.1 can be updated with
#' \code{\link{upgradeFLSAM}}.
#'
#' @source ICES. 2012. Report of the Herring Assessment Working Group for the
#'   Area South of 62 N (HAWG). ICES CM 2012/ACOM:06.
#' @examples
#' data(NSH)
#' data(NSH); data(NSH.sam)
#' plot(NSH.sam)
#' NSH <- NSH + NSH.sam
#' data(HERAS.sams)
#' lr.test(HERAS.sams)
#' @docType data
#' @keywords datasets
#' @name FLSAM-data
#' @aliases NSH NSH.tun NSH.ctrl NSH.sam HERAS.sams SOL SOL.tun SOL.ctrl SOL.sam
NULL
