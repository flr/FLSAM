# accessors.R - DESC
# /home/mosquia/Projects/FLR/pkgs/stable/FLSAM/R/accessors.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


#' Slot accessors for FLSAM objects
#'
#' Get and set the slots of an \code{\link{FLSAM-class}} object, e.g.
#' \code{control(sam)}, \code{nopar(sam)}, \code{stock.n(sam)},
#' \code{harvest(sam)}, \code{params(sam)}, \code{residuals(sam)},
#' \code{vcov(sam)}, \code{rescov(sam)}, \code{obscov(sam)} and
#' \code{info(sam)}. For estimates with confidence intervals, see
#' \code{\link{accessors}}; for the stored stockassessment fit see
#' \code{\link{samfit}}.
#'
#' @param object An \code{FLSAM} object.
#' @param value Replacement value.
#' @param ... Not used.
#' @return The slot content (or the modified object for replacement methods).
#' @examples
#' data(NSH); data(NSH.sam)
#' nopar(NSH.sam)
#' head(residuals(NSH.sam))
#' stock.n(NSH.sam)[,"2011"]
#' @name FLSAM-slots
NULL

# control {{{

#' @rdname FLSAM-slots
setMethod("control", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "control"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("control", signature(object="FLSAM", value="FLSAM.control"),
  function(object, value) {
    slot(object, "control") <- value
    return(object)
  }
)
# }}}

# nopar {{{
#' @rdname FLSAM-slots
setGeneric("nopar", function(object, ...) standardGeneric("nopar"))

#' @rdname FLSAM-slots
setMethod("nopar", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "nopar"))
  }
)

#' @rdname FLSAM-slots
setGeneric("nopar<-", function(object, ..., value) standardGeneric("nopar<-"))

#' @rdname FLSAM-slots
setReplaceMethod("nopar", signature(object="FLSAM", value="numeric"),
  function(object, value) {
    slot(object, "nopar") <- as.integer(value)
    return(object)
  }
)
# }}}

# n.states {{{
#' @rdname FLSAM-slots
setGeneric("n.states", function(object, ...) standardGeneric("n.states"))

#' @rdname FLSAM-slots
setMethod("n.states", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "n.states"))
  }
)

#' @rdname FLSAM-slots
setGeneric("n.states<-", function(object, ..., value) standardGeneric("n.states<-"))

#' @rdname FLSAM-slots
setReplaceMethod("n.states", signature(object="FLSAM", value="numeric"),
  function(object, value) {
    slot(object, "n.states") <- as.integer(value)
    return(object)
  }
)
# }}}

# states {{{
#' @rdname FLSAM-slots
setGeneric("states", function(object, ...) standardGeneric("states"))

#' @rdname FLSAM-slots
setMethod("states", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "states"))
  }
)

#' @rdname FLSAM-slots
setGeneric("states<-", function(object, ..., value) standardGeneric("states<-"))

#' @rdname FLSAM-slots
setReplaceMethod("states", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "states") <- as.integer(value)
    return(object)
  }
)
# }}}

# components {{{
#' @noRd
setGeneric("components", function(object, ...) standardGeneric("components"))

#' @noRd
setMethod("components", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "components"))
  }
)

#' @rdname FLSAM-slots
setGeneric("components<-", function(object, ..., value) standardGeneric("components<-"))

#' @rdname FLSAM-slots
setReplaceMethod("components", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "components") <- value
    return(object)
  }
)
# }}}

# stock.wt {{{

#' @rdname FLSAM-slots
setMethod("stock.wt", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "stock.wt"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("stock.wt", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "stock.wt") <- value
    return(object)
  }
)
# }}}

# catch.wt {{{

#' @rdname FLSAM-slots
setMethod("catch.wt", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "catch.wt"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("catch.wt", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "catch.wt") <- value
    return(object)
  }
)
# }}}

# mat {{{

#' @rdname FLSAM-slots
setMethod("mat", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "mat"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("mat", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "mat") <- value
    return(object)
  }
)
# }}}

# m {{{

#' @rdname FLSAM-slots
setMethod("m", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "m"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("m", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "m") <- value
    return(object)
  }
)
# }}}

# stock.n {{{

#' @rdname FLSAM-slots
setMethod("stock.n", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "stock.n"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("stock.n", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "stock.n") <- value
    return(object)
  }
)
# }}}

# harvest {{{

#' @rdname FLSAM-slots
setMethod("harvest", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "harvest"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("harvest", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "harvest") <- value
    return(object)
  }
)
# }}}

# nlogl {{{
#' @rdname FLSAM-slots
setGeneric("nlogl", function(object, ...) standardGeneric("nlogl"))

#' @rdname FLSAM-slots
setMethod("nlogl", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "nlogl"))
  }
)

#' @rdname FLSAM-slots
setGeneric("nlogl<-", function(object, ..., value) standardGeneric("nlogl<-"))

#' @rdname FLSAM-slots
setReplaceMethod("nlogl", signature(object="FLSAM", value="numeric"),
  function(object, value) {
    slot(object, "nlogl") <- value
    return(object)
  }
)
# }}}

# vcov {{{

#' @rdname FLSAM-slots
setMethod("vcov", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "vcov"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("vcov", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "vcov") <- value
    return(object)
  }
)
# }}}

# rescov {{{
#' @rdname FLSAM-slots
setGeneric("rescov", function(object, ...) standardGeneric("rescov"))

#' @rdname FLSAM-slots
setMethod("rescov", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "rescov"))
  }
)

#' @rdname FLSAM-slots
setGeneric("rescov<-", function(object, ..., value) standardGeneric("rescov<-"))

#' @rdname FLSAM-slots
setReplaceMethod("rescov", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "rescov") <- value
    return(object)
  }
)
# }}}

# obscov {{{
#' @rdname FLSAM-slots
setGeneric("obscov", function(object, ...) standardGeneric("obscov"))

#' @rdname FLSAM-slots
setMethod("obscov", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "obscov"))
  }
)

#' @rdname FLSAM-slots
setGeneric("obscov<-", function(object, ..., value) standardGeneric("obscov<-"))

#' @rdname FLSAM-slots
setReplaceMethod("obscov", signature(object="FLSAM", value="list"),
  function(object, value) {
    slot(object, "obscov") <- value
    return(object)
  }
)
# }}}

# params {{{

#' @rdname FLSAM-slots
setMethod("params", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "params"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("params", signature(object="FLSAM", value="data.frame"),
  function(object, value) {
    slot(object, "params") <- value
    return(object)
  }
)
# }}}

# residuals {{{

#' @rdname FLSAM-slots
setMethod("residuals", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "residuals"))
  }
)

#' @rdname FLSAM-slots
setReplaceMethod("residuals", signature(object="FLSAM", value="data.frame"),
  function(object, value) {
    slot(object, "residuals") <- value
    return(object)
  }
)
# }}}

# info {{{
#' @rdname FLSAM-slots
setGeneric("info", function(object, ...) standardGeneric("info"))

#' @rdname FLSAM-slots
setMethod("info", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "info"))
  }
)

#' @rdname FLSAM-slots
setGeneric("info<-", function(object, ..., value) standardGeneric("info<-"))

#' @rdname FLSAM-slots
setReplaceMethod("info", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "info") <- value
    return(object)
  }
)
# }}}
