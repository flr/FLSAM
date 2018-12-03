
# FLSAM
- Version: 2.1.0
- Date: 2018-11-27
- Author: M.R. Payne <mpa ΑΤ aqua.dtu.dk>, N.T. Hintzen <niels.hintzen ΑΤ wur.nl>
- Maintainer: M.R. Payne <mpa ΑΤ aqua.dtu.dk>, N.T. Hintzen <niels.hintzen ΑΤ wur.nl>
- Repository: <https://github.com/flr/FLCore/>
- Bug reports: <https://github.com/flr/FLCore/issues>

## Overview
The FLSAM package provides an FLR version of SAM, the State-space Assessment Model developed by Anders Nielsen, DTU-Aqua. For details regarding SAM, please consult <http://www.stockassessment.org>. 

To install this package, start R and enter:

	install.packages("FLSAM", repos="http://flr-project.org/R")

or download from the [FLSAM releases page](https://github.com/flr/FLSAM/releases/latest)

**PLEASE NOTE** that this version of FLSAM requires a particular version of the 'stockassessment' package. Please install using

  devtools::install_github('fishwollower/SAM/stockassessment', ref='component')

## Documentation
- [Help pages](http://flr-project.org/FLSAM)

## Build Status
[![Travis Build Status](https://travis-ci.org/flr/FLSAM.svg?branch=master)](https://travis-ci.org/flr/FLSAM)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/flr/FLSAM?branch=master&svg=true)](https://ci.appveyor.com/project/flr/FLSAM)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/FLSAM)](https://cran.r-project.org/package=FLSAM)

## Releases
- [All release](https://github.com/flr/FLSAM/releases/)

## License
Copyright (c) 2004-2015 The FLR Team. Released under the [GPL](http://www.gnu.org/licenses/gpl-2.0.html).

## Contact
You are welcome to:

- Submit suggestions and bug-reports at: <https://github.com/flr/FLSAM/issues>
- Send a pull request on: <https://github.com/flr/FLSAM/>
- Compose a friendly e-mail to: <flrteam AT flr-project.org>
