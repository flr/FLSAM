#' FLSAM: an FLR wrapper to the State-space Assessment Model (SAM)
#'
#' FLSAM runs the state-space assessment model SAM (package
#' \pkg{stockassessment}, \url{https://github.com/fishfollower/SAM}) from FLR
#' objects. The main steps are:
#' \itemize{
#'   \item set up a configuration with \code{\link{FLSAM.control}};
#'   \item fit the model with \code{\link{FLSAM}}, which returns an
#'     \code{\link{FLSAM-class}} object that also stores the underlying
#'     \code{sam} fit (\code{\link{samfit}});
#'   \item extract results with the \code{\link{accessors}} (\code{ssb},
#'     \code{fbar}, \code{rec}, \code{catchabilities}, \code{obs.var}, ...) or merge them into an \code{FLStock}
#'     (\code{stock + sam});
#'   \item diagnose the fit: \code{\link{residual.diagnostics}},
#'     \code{\link{residual.tests}}, \code{\link{procerr.plot}},
#'     \code{\link{cor.plot}}, \code{\link{obscor.empirical.plot}},
#'     \code{\link{obssd.plot}}, \code{\link{jitter.test}},
#'     \code{\link{retro}} and \code{\link{mohns.rho}}, \code{\link{looi}},
#'     \code{\link{lr.test}}, \code{\link{AIC}};
#'   \item forecasts and reference points: \code{\link{FLSAM.forecast}},
#'     \code{\link{FLSAM.modelforecast}}, \code{\link{FLSAM.refpts}},
#'     \code{\link{FLSAM.ypr}}.
#' }
#' All plots are made with \pkg{ggplot2} (except a few base-graphics plots:
#' \code{residual.diagnostics}, \code{otolith}, \code{comp.plot},
#' \code{obsvar.plot}, \code{obscv.plot}, \code{obscor.plot}).
#'
#' @name FLSAM-package
#' @aliases FLSAM-package
#' @rawNamespace exportPattern("^[^\\.]")
#' @import FLCore
#' @import stockassessment
#' @importFrom utils packageDescription packageVersion capture.output relist write.table
#' @importFrom parallel detectCores makeCluster clusterEvalQ stopCluster
#' @importFrom TMB MakeADFun oneStepPredict sdreport
#' @importFrom ellipse plotcorr ellipse
#' @importFrom MASS kde2d
#' @importFrom graphics abline axis barplot contour grid layout lcm legend lines par points polygon text title
#' @importFrom grDevices colorRampPalette
#' @importFrom stats acf aggregate approxfun as.ts setNames cov2cor density loess.smooth na.omit na.pass nlminb pchisq qqline qqnorm plogis xtabs cor ave
#' @importFrom methods as is new slot slot<- slotNames validObject .hasSlot setOldClass
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon geom_col geom_tile geom_text geom_path geom_hline facet_wrap labs theme_bw theme element_text expand_limits scale_fill_gradient2 scale_fill_gradient scale_fill_brewer scale_colour_manual scale_colour_viridis_c scale_size_area scale_fill_manual scale_colour_gradient scale_linetype_manual waiver
"_PACKAGE"

utils::globalVariables(c("year","value","lbnd","ubnd","name","age","fleet","std.res","retro","cohort",
                         "data","bdiff","x","Var1","Var2","cor","age1","age2","y","sel","p","obs",
                         "quantity","harvest","stock","nameOrig","param","sign","idx","period","%dopar%"))
