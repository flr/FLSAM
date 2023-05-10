# accessors.R - DESC
# /home/mosquia/Projects/FLR/pkgs/stable/FLSAM/R/accessors.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# control {{{

setMethod("control", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "control"))
  }
)

setReplaceMethod("control", signature(object="FLSAM", value="FLSAM.control"),
  function(object, value) {
    slot(object, "control") <- value
    return(object)
  }
)
# }}}

# nopar {{{
setGeneric("nopar", function(object, ...) standardGeneric("nopar"))

setMethod("nopar", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "nopar"))
  }
)

setGeneric("nopar<-", function(object, ..., value) standardGeneric("nopar<-"))

setReplaceMethod("nopar", signature(object="FLSAM", value="numeric"),
  function(object, value) {
    slot(object, "nopar") <- as.integer(value)
    return(object)
  }
)
# }}}

# n.states {{{
setGeneric("n.states", function(object, ...) standardGeneric("n.states"))

setMethod("n.states", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "n.states"))
  }
)

setGeneric("n.states<-", function(object, ..., value) standardGeneric("n.states<-"))

setReplaceMethod("n.states", signature(object="FLSAM", value="numeric"),
  function(object, value) {
    slot(object, "n.states") <- as.integer(value)
    return(object)
  }
)
# }}}

# states {{{
setGeneric("states", function(object, ...) standardGeneric("states"))

setMethod("states", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "states"))
  }
)

setGeneric("states<-", function(object, ..., value) standardGeneric("states<-"))

setReplaceMethod("states", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "states") <- as.integer(value)
    return(object)
  }
)
# }}}

# components {{{
setGeneric("components", function(object, ...) standardGeneric("components"))

setMethod("components", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "components"))
  }
)

setGeneric("components<-", function(object, ..., value) standardGeneric("components<-"))

setReplaceMethod("components", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "components") <- as.integer(value)
    return(object)
  }
)
# }}}

# stock.wt {{{

setMethod("stock.wt", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "stock.wt"))
  }
)

setReplaceMethod("stock.wt", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "stock.wt") <- value
    return(object)
  }
)
# }}}

# catch.wt {{{

setMethod("catch.wt", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "catch.wt"))
  }
)

setReplaceMethod("catch.wt", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "catch.wt") <- value
    return(object)
  }
)
# }}}

# mat {{{

setMethod("mat", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "mat"))
  }
)

setReplaceMethod("mat", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "mat") <- value
    return(object)
  }
)
# }}}

# m {{{

setMethod("m", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "m"))
  }
)

setReplaceMethod("m", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "m") <- value
    return(object)
  }
)
# }}}

# stock.n {{{

setMethod("stock.n", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "stock.n"))
  }
)

setReplaceMethod("stock.n", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "stock.n") <- value
    return(object)
  }
)
# }}}

# harvest {{{

setMethod("harvest", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "harvest"))
  }
)

setReplaceMethod("harvest", signature(object="FLSAM", value="FLQuant"),
  function(object, value) {
    slot(object, "harvest") <- value
    return(object)
  }
)
# }}}

# nlogl {{{
setGeneric("nlogl", function(object, ...) standardGeneric("nlogl"))

setMethod("nlogl", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "nlogl"))
  }
)

setGeneric("nlogl<-", function(object, ..., value) standardGeneric("nlogl<-"))

setReplaceMethod("nlogl", signature(object="FLSAM", value="numeric"),
  function(object, value) {
    slot(object, "nlogl") <- value
    return(object)
  }
)
# }}}

# vcov {{{

setMethod("vcov", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "vcov"))
  }
)

setReplaceMethod("vcov", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "vcov") <- value
    return(object)
  }
)
# }}}

# rescov {{{
setGeneric("rescov", function(object, ...) standardGeneric("rescov"))

setMethod("rescov", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "rescov"))
  }
)

setGeneric("rescov<-", function(object, ..., value) standardGeneric("rescov<-"))

setReplaceMethod("rescov", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "rescov") <- value
    return(object)
  }
)
# }}}

# obscov {{{
setGeneric("obscov", function(object, ...) standardGeneric("obscov"))

setMethod("obscov", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "obscov"))
  }
)

setGeneric("obscov<-", function(object, ..., value) standardGeneric("obscov<-"))

setReplaceMethod("obscov", signature(object="FLSAM", value="list"),
  function(object, value) {
    slot(object, "obscov") <- value
    return(object)
  }
)
# }}}

# params {{{

setMethod("params", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "params"))
  }
)

setReplaceMethod("params", signature(object="FLSAM", value="data.frame"),
  function(object, value) {
    slot(object, "params") <- value
    return(object)
  }
)
# }}}

# residuals {{{

setMethod("residuals", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "residuals"))
  }
)

setReplaceMethod("residuals", signature(object="FLSAM", value="data.frame"),
  function(object, value) {
    slot(object, "residuals") <- value
    return(object)
  }
)
# }}}

# info {{{
setGeneric("info", function(object, ...) standardGeneric("info"))

setMethod("info", signature(object="FLSAM"),
  function(object) {
    return(slot(object, "info"))
  }
)

setGeneric("info<-", function(object, ..., value) standardGeneric("info<-"))

setReplaceMethod("info", signature(object="FLSAM", value="matrix"),
  function(object, value) {
    slot(object, "info") <- as.integer(value)
    return(object)
  }
)
# }}}
