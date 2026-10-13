# FLSAM 3.1.0

Alignment with stockassessment (SAM) 0.12.0 (master of fishfollower/SAM).

## New features
- `FLSAM` objects store the stockassessment fit in the new slot `fit`
  (`FLSAM(..., keep.fit=TRUE)`); access it with `samfit()` to use any
  stockassessment function.
- `set.pars` accepts all fixed-effect SAM parameters (FLSAM or SAM names, see
  `setParsTable()`).
- `FLSAM.control`: `fixVarToWeight` per fleet, new slot `hockeyStickCurve`
  (srr 63), all `biomassTreat` codes of SAM documented and validated.
- Wrappers: `FLSAM.forecast()`, `FLSAM.modelforecast()`, `FLSAM.refpts()`,
  `FLSAM.ypr()`, `FLSAM.catchbyfleet()`, `model.description()`.
- Diagnostics (ggplot2): `residual.tests()`, `jitter.test()`,
  `obscor.empirical.plot()`, `obssd.plot()`, `predsd.plot()`, `data.plot()`,
  `fselectivity.plot()`.
- Objects saved with older versions (missing the new slots) keep working:
  `FLSAM()`, `ctrl2conf()`, `update()`, `drop.from.control()`, `SAM2FLR()`,
  `FLSAM.out()` and printing upgrade them automatically; `upgradeFLSAM()`
  upgrades `FLSAM.control`, `FLSAM` and `FLSAMs` objects explicitly.

## Changes
- All lattice plots are now ggplot2 plots and return `ggplot` objects: the
  `plot()` methods, `fpart.plot()`, `cor.plot()`, `procerr.plot()`,
  `retroResiduals()`, `retroSelectivity()`, `retroParams()`.
- `mohns.rho()` follows Mohn (1999) and `stockassessment::mohn`:
  `100*(peel-ref)/ref` (sign reversed compared to earlier versions); the mean
  is returned in attribute `"mohn"`.
- Documentation generated with roxygen2; all help pages rewritten.
- `.onAttach` requires stockassessment >= 0.12.0 (was exactly 0.12.0).

## Bug fixes
- `ctrl2conf()` ignored `fixVarToWeight`.
- `FLSAM2SAM()` used `catch.wt` for landing and discard weights.
- `update()` on `FLSAM.control` did not renumber stock/catch weight, maturity
  and natural mortality couplings.
- `conf2ctrl()` returned `logN.vars` as a matrix.
- `validObject()` on `FLSAM.control` failed.
- `mohns.rho()` failed without `ref.year`.
- `looi()` without `type` failed; `retro()` on `FLStocks` ignored `set.pars`
  in the first peel.
- `rec.pars()` used a global `ctrl` and returned duplicated rows.
- `FLSAM.MSE()` did not run; `FLSAM.fast()` failed when surveys extend beyond
  the catch data.
- `monteCarloStock()` / `monteCarloStockMscan()` failed with `set.pars` and
  when realisations did not converge; `ncores` default needed `parallel`.
- `residual.diagnostics()` never drew the smoother; `components<-` and
  `info<-` corrupted the values; `plotkobe()` used a global `NSH.sam`;
  `SIM2FLR()` removed the wrong parameters when several were fixed.
