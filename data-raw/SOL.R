# SOL.R - DESC
# /SOL.R

# Copyright European Union, 2018
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.


library(FLSAM)

load("SOL.RData")

SOL.ctrl <- FLSAM.control(SOL, SOL.tun)

save(SOL, SOL.tun, SOL.ctrl, file="../data/SOL.RData", compress="xz")
