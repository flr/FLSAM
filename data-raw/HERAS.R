# HERAS.R - DESC
# /HERAS.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

library(FLSAM)

data(NSH)

ctrl.list <- list(all.free=3:10,
                  four.lvls=c(3,3,4,4,5,5,6,6),
                  three.lvls=c(3,3,4,4,5,5,5,5),
                  two.lvls=c(3,3,3,4,4,4,4,4))

ctrls <- lapply(ctrl.list,
  function(x) {
    tmp.ctrl <- NSH.ctrl
    tmp.ctrl@catchabilities["HERAS",ac(1:8)] <- x
    return(update(tmp.ctrl))
})

sam.list <- lapply(ctrls,FLSAM,stck=NSH,tun=NSH.tun)                  
HERAS.sams <- FLSAMs(sam.list)

save(HERAS.sams, file="../data/HERAS.sams.RData", compress="xz")

