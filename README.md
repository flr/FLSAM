# FLSAM

- Version: 3.1.0
- Date: 2026-10-13
- Author and maintainer: N.T. Hintzen <nhintzen AT pelagicfish.eu>
- Original co-author: M.R. Payne
- Repository: <https://github.com/flr/FLSAM/>
- Bug reports: <https://github.com/flr/FLSAM/issues>
- Changes: see [NEWS.md](NEWS.md)

## Overview

FLSAM runs SAM, the State-space Assessment Model developed by Anders Nielsen
(DTU Aqua), from [FLR](https://flr-project.org) objects. The model itself is
fitted by the [stockassessment](https://github.com/fishfollower/SAM) package;
FLSAM translates `FLStock`, `FLIndices` and an `FLSAM.control` configuration
into SAM input, and the results back into an `FLSAM` object.

Main features:

- configuration of all SAM options of stockassessment 0.12 through `FLSAM.control`,
  including fleet-specific weighting, AR/US observation correlation,
  prediction-variance links, biological process models and stock-recruitment models;
- the `stockassessment` fit is stored in the `FLSAM` object (`samfit()`), so all
  stockassessment functions can be used as well;
- fixing any parameter (`set.pars`) and starting values from an earlier
  assessment (`starting.values`), also with one year more or less data;
- diagnostics: one-step-ahead residuals and residual tests, process error,
  observation variances and correlations, parameter correlations, jitter test,
  retrospective analysis with Mohn's rho, leave-one-in/leave-one-out analyses,
  likelihood-ratio tests;
- forecasts, reference points and per-recruit analyses through stockassessment;
- simulation, Monte Carlo and MSE helpers;
- ggplot2 plots with a configurable colour palette (`FLSAM.colours()`).

## Installation

FLSAM needs the `stockassessment` package (>= 0.12.0):

    install.packages("stockassessment", repos = "https://fishfollower.r-universe.dev")

or

    remotes::install_github("fishfollower/SAM/stockassessment")

Then install FLSAM from the FLR repository:

    install.packages("FLSAM", repos = "https://flr.r-universe.dev")

or the development version from GitHub:

    remotes::install_github("flr/FLSAM", ref = "dev-sam-align")

## Getting started

    library(FLSAM)
    data(NSH)                              # North Sea herring: NSH, NSH.tun, NSH.ctrl
    sam <- FLSAM(NSH, NSH.tun, NSH.ctrl)   # run the assessment
    plot(sam)                              # SSB, F and recruitment
    residual.diagnostics(sam)              # residual diagnostics by fleet and age
    residual.tests(sam)                    # statistical tests on the residuals
    NSH.retro <- retro(NSH, NSH.tun, NSH.ctrl, retro = 5)
    mohns.rho(NSH.retro, type = "ssb")

A complete overview of all functions, using the built-in herring data, is
given in [dev/FLSAM_overview.R](dev/FLSAM_overview.R).

## Documentation

- Help pages: `help(package = "FLSAM")`, start with `?FLSAM`
- FLR: <https://flr-project.org>
- SAM: <https://www.stockassessment.org>

## License

Copyright (c) 2004-2026 The FLR Team. Released under the
[GPL](http://www.gnu.org/licenses/gpl-2.0.html).

## Contact

You are welcome to submit suggestions and bug reports at
<https://github.com/flr/FLSAM/issues>, to send a pull request on
<https://github.com/flr/FLSAM/>, or to e-mail the maintainer at
nhintzen AT pelagicfish.eu.
