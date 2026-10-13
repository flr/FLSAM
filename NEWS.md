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
- Plot colours from `FLSAM.colours()`: a qualitative palette starting with
  dark blue #004868, orange #FF800A, light blue #00BADC and green #31C381;
  series beyond the palette get another line type. Change with
  `options(FLSAM.palette=...)`. `retroSelectivity()` keeps the colours of
  earlier versions; `retroParams()` restarts the colours in every facet;
  `fselectivity.plot()` colours years by period (`period=10`).
- `FLSAM(..., starting.values=)` works again (`startingValues()`): an earlier
  assessment (`FLSAM` with or without stored fit, or a `sam` fit), also with
  one year less or more data; states are matched by year. `retro()` uses the
  previous peel as starting values for the next (`use.start=TRUE`, default).
- `pin.R` (`FLSAM2par()`, legacy pin functionality) removed.
- `parameter.table()`: all estimated parameters with readable names (e.g.
  "catchability HERAS 5-8"), type, back-transformed value and 95% CI;
  `params.plot()` compares them between assessments or retrospective peels;
  `cor.plot()` and `retroParams()` use the readable names.
- `ssbfbar.ellipse()`: analytical confidence ellipse of SSB and Fbar for one or
  more years and assessments.
- Tag-recapture data (SAM fleet type 5): `FLSAM.control(..., recapture=)`,
  discrete one-step-ahead residuals for the recaptures, `tag.pars()`,
  `tags.plot()`; `retro()` cuts the tag data per peel. Tested against the
  WGWIDE 2026 NEA mackerel assessment (identical fit and retrospective peels).
- `retro()` gains `catch.vars` (catch weights were not passed to the peels).
- `plot(FLSAMs, ci=TRUE)`: confidence bands for retrospective peels or several
  assessments.
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
- `simulate()` could not be called without `set.pars`, and failed for more
  than one simulation with `set.pars`.
- `sam.out()` called itself (infinite recursion) instead of `FLSAM.out()`.
- `obscor.plot()` gave an unclear error without correlated observations.
- `conf2ctrl()` did not set `logP.vars` and set a wrong plus group, so
  `SAM2FLR()` failed for fits made outside FLSAM; `SAM2FLR()` identifies the
  catch fleets by fleet type instead of by the word "catch" in their name.
