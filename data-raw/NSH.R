# NSH.R - DESC
# /NSH.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


library(FLSAM)

load("NSH.RData")

NSH.ctrl <- FLSAM.control(NSH, NSH.tun)

save(NSH, NSH.tun, NSH.ctrl, file="../data/NSH.RData", compress="xz")

NSH.sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)

save(NSH.sam, file="../data/NSH.sam.RData", compress="xz")
