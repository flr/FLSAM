# test-MSE.R - DESC
# /test-MSE.R

# Copyright Iago MOSQUEIRA (WMR), 2021
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

library(FLSAM)
data(NSH)
data(NSH.sam)

NSH.iter <- expand(NSH,iter=1:5)
# res <- FLSAM.MSE(NSH.iter,NSH.tun,NSH.ctrl,starting.sam=NSH.sam)
