#-------------------------------------------------------------------------------
# Documentation of S4 generics and methods that are defined across files
#-------------------------------------------------------------------------------

#' Estimates with confidence intervals from FLSAM objects
#'
#' Extract model estimates in an easy to use format, corrected for the
#' binding of parameters and back-transformed from log space, together with
#' their 95\% confidence intervals. All functions work on \code{\link{FLSAM-class}}
#' and \code{\link{FLSAMs}} objects (the latter add a column \code{name}).
#' \itemize{
#'   \item \code{ssb}, \code{tsb}, \code{fbar}, \code{catch}: spawning stock
#'     biomass, total stock biomass, mean fishing mortality and total catch (in
#'     weight) by year;
#'   \item \code{n}, \code{rec}: numbers at age and recruitment;
#'   \item \code{f}: fishing mortality at age (by catch fleet);
#'   \item \code{catch.n}: observed and modelled catch numbers (from the residuals);
#'   \item \code{components}: proportions of partial-survey components;
#'   \item \code{catchabilities}, \code{power.law.exps}, \code{obs.var}:
#'     catchability, density-dependent catchability power and observation
#'     standard deviation by fleet and age;
#'   \item \code{n.var}, \code{f.var}: process standard deviations of log(N)
#'     and log(F);
#'   \item \code{cor.F}, \code{cor.obs}: correlation parameters of the F
#'     process and of the observations;
#'   \item \code{catch.scale}: catch multipliers;
#'   \item \code{rec.pars}: parameters of the stock-recruitment model
#'     (\code{srr a}, \code{srr b}, ...).
#' }
#'
#' @param object An \code{FLSAM} or \code{FLSAMs} object.
#' @param ... Not used.
#' @return A \code{data.frame} with meta-data columns (\code{name},
#'   \code{fleet}, \code{year}, \code{age}, as applicable) and \code{value},
#'   \code{CV} (standard deviation on log scale), \code{lbnd} and \code{ubnd}
#'   (95\% confidence limits).
#' @examples
#' data(NSH); data(NSH.sam)
#' catchabilities(NSH.sam)
#' obs.var(NSH.sam)
#' head(ssb(NSH.sam))
#' head(fbar(NSH.sam))
#' head(rec(NSH.sam))
#' data(HERAS.sams)
#' head(ssb(HERAS.sams))
#' @name accessors
NULL

#' Parameter estimates of FLSAM objects
#'
#' \code{coef} and \code{coefficients} return the estimated (fixed-effect)
#' parameters, \code{params} all parameters, random effects and derived
#' quantities, with their standard deviations. Parameter names are those of
#' SAM, e.g. \code{logFpar} (catchabilities), \code{logSdLogFsta},
#' \code{logSdLogN}, \code{logSdLogObs}, \code{transfIRARdist},
#' \code{itrans_rho}, \code{rec_pars}, \code{logN}, \code{logF},
#' \code{logssb}, \code{logfbar}.
#'
#' @param object An \code{\link{FLSAM-class}} or \code{\link{FLSAMs}} object.
#' @param ... Not used.
#' @return A \code{data.frame} with columns \code{name}, \code{value} and
#'   \code{std.dev} (and \code{object.name} for \code{FLSAMs}).
#' @seealso \code{\link{accessors}}, \code{\link{FLSAM-slots}}
#' @examples
#' data(NSH); data(NSH.sam)
#' coef(NSH.sam)
#' @name parameter.accessors
NULL

#' Merge FLSAM results into FLStock objects
#'
#' \code{merge(stock, sam)} and \code{stock + sam} copy the estimated numbers
#' at age (\code{stock.n}) and fishing mortality (\code{harvest}) of an
#' \code{\link{FLSAM-class}} object into an \code{FLStock}, trimming to the
#' common year range. With an \code{\link{FLSAMs}} object an \code{FLStocks}
#' is returned; with an \code{FLStocks} (multi-fleet) the results are copied
#' into each stock.
#'
#' @param x,e1,e2 An \code{FLStock}, \code{FLStocks}, \code{FLSAM} or \code{FLSAMs} object.
#' @param y An \code{FLSAM} object.
#' @param ... Not used.
#' @return An \code{FLStock} or \code{FLStocks} object.
#' @examples
#' data(NSH); data(NSH.sam)
#' NSH <- NSH + NSH.sam
#' @name FLSAM-merge
NULL

#' Akaike information criterion of FLSAM objects
#'
#' \code{AIC = 2 * nopar + 2 * nlogl}.
#'
#' @param object An \code{\link{FLSAM-class}} or \code{\link{FLSAMs}} object.
#' @param ... Not used.
#' @param k Not used.
#' @return A numeric value (vector for \code{FLSAMs}).
#' @examples
#' data(HERAS.sams)
#' AIC(HERAS.sams)
#' @name AIC-FLSAM
NULL

#' Simulate assessments from the estimated uncertainty
#'
#' Fits the assessment and draws \code{n} parameter vectors (fixed and random
#' effects) from a multivariate normal distribution with the joint precision
#' matrix of the fit; each draw is converted into an \code{FLSAM} object.
#'
#' @param x An \code{FLStock}.
#' @param y An \code{FLIndices}.
#' @param z An \code{\link{FLSAM.control-class}}.
#' @param n Number of draws.
#' @param set.pars Optional \code{set.pars} (see \code{\link{FLSAM}}).
#' @param set.seed Optional random seed.
#' @return An \code{\link{FLSAMs}} object with \code{n} elements. These hold the
#'   simulated parameters, \code{stock.n} and \code{harvest}, but no derived
#'   quantities: calculate SSB, Fbar etc. via the \code{FLStock}
#'   (\code{stock + sims}).
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' sims <- simulate(NSH, NSH.tun, NSH.ctrl, n=10, set.seed=1)
#' }
#' @name simulate-FLSAM
NULL

#' Likelihood ratio test between FLSAM assessments
#'
#' Compares nested models with a likelihood ratio test: the difference in
#' negative log-likelihood, degrees of freedom and the p-value from a
#' chi-squared distribution. Models must be ordered by increasing or
#' decreasing number of parameters.
#'
#' @param object An \code{\link{FLSAMs}} object, or an \code{FLSAM} object
#'   followed by other \code{FLSAM} objects in \code{...}.
#' @param ... Further \code{FLSAM} objects.
#' @param type \code{"sequential"} compares model 1 vs 2, 2 vs 3, ...;
#'   \code{"first"} compares every model against model 1.
#' @return A table with the test statistics.
#' @section Warning: Only valid for nested models fitted to the same data.
#' @examples
#' data(HERAS.sams)
#' lr.test(HERAS.sams)
#' lr.test(HERAS.sams[[1]], HERAS.sams[[2]], HERAS.sams[[3]])
#' @name lr.test
NULL

#' Leave-one-in and leave-one-out analyses
#'
#' Runs a series of assessments where the survey fleets are left out of or
#' kept in the assessment. \code{loo} leaves one fleet out at a time,
#' \code{loi} keeps only one fleet in, and \code{looi} with
#' \code{type="full"} (the default) runs all combinations. Partial-survey
#' fleets are treated as one fleet.
#'
#' @param stck An \code{FLStock} (or \code{FLStocks}).
#' @param tun An \code{FLIndices} object.
#' @param ctrl An \code{\link{FLSAM.control-class}} object.
#' @param type \code{"full"}, \code{"LOO"} or \code{"LOI"}.
#' @param base.run Optional \code{FLSAM} with all fleets, to avoid refitting it.
#' @param set.pars Optional \code{set.pars} (see \code{\link{FLSAM}}).
#' @param ... Not used.
#' @return An \code{\link{FLSAMs}} object named by the fleet combination.
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' one.out <- loo(NSH, NSH.tun, NSH.ctrl)
#' plot(one.out)
#' }
#' @name looi
NULL

#' Retrospective analysis with FLSAM
#'
#' Repeats the assessment removing one year of data at a time. Survey years
#' beyond the last catch year are kept consistently in each peel.
#'
#' @param stock An \code{FLStock} or \code{FLStocks} object.
#' @param indices An \code{FLIndices} object.
#' @param control An \code{\link{FLSAM.control-class}} object.
#' @param retro Number of retrospective years (peels).
#' @param year.range Alternatively, the assessment years to run.
#' @param set.pars Optional \code{set.pars} (see \code{\link{FLSAM}}).
#' @param catch.vars Optional catch weights (see \code{\link{FLSAM}}), cut to
#'   the years of every peel.
#' @param use.start Use the estimates of the previous (one year longer) peel
#'   as starting values for the next peel (default \code{TRUE}, see
#'   \code{\link{startingValues}}). Set to \code{FALSE} to start every peel
#'   from the default starting values.
#' @details Tag-recapture data in the control object (see
#'   \code{\link{FLSAM.control}}) are cut per peel such that the lag between the
#'   last recapture year and the last year of the assessment (at most one year)
#'   is kept.
#' @return An \code{\link{FLSAMs}} object named by assessment year; failed
#'   peels are dropped with a warning.
#' @seealso \code{\link{mohns.rho}}, \code{\link{retroResiduals}}
#' @examples
#' \donttest{
#' data(NSH); data(NSH.sam)
#' res <- retro(NSH, NSH.tun, NSH.ctrl, retro=3)
#' plot(res)
#' mohns.rho(res, type="ssb")
#' }
#' @name retro
NULL

#' Residual data for ggplot-style diagnostics
#'
#' Returns the residuals of an \code{\link{FLSAM-class}} object with lagged
#' residuals, normal quantiles and the estimated numbers at age added, for
#' use in custom diagnostic plots (requires \pkg{plyr}).
#'
#' @param object An \code{FLSAM} object.
#' @param ... Not used.
#' @return A \code{data.frame}.
#' @name diags
NULL

#' FLSAMs class
#'
#' A list of \code{\link{FLSAM-class}} objects, e.g. the result of
#' \code{\link{retro}} or \code{\link{looi}}. Create with \code{FLSAMs(...)}
#' from \code{FLSAM} objects or a list of them.
#'
#' @param object An \code{FLSAM} object, a list of \code{FLSAM} objects, or missing.
#' @param ... Further \code{FLSAM} objects.
#' @return An \code{FLSAMs} object.
#' @examples
#' data(HERAS.sams)
#' names(HERAS.sams)
#' plot(HERAS.sams)
#' @name FLSAMs
#' @aliases FLSAMs-class
NULL
