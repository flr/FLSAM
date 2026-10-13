# FLSAM 3.0.0.9005 → 3.1.0: overzicht van de wijzigingen

Branch: `dev-sam-align` (vanaf tag `v3.0.0.9005` = oude master `89bb528`).
Doelversie van SAM: `stockassessment` 0.12.0 (master van fishfollower/SAM), **niet** `F_TAC`.

Testen: `dev/compare_versions.R` (deel A: oud vs nieuw, deel B: nieuwe functionaliteit).
Resultaat van deel A: 136 onderdelen bit-identiek, de rest alleen afrondingsverschillen of
verklaarde verschillen (zie onderaan). `R CMD check`: 0 errors, 0 warnings, 2 notes.

## Per bestand

Legenda: **B** = bugfix, **U** = uitbreiding (afstemming op SAM), **G** = lattice → ggplot2,
**D** = alleen documentatie (roxygen), **N** = nieuw.

| Bestand | Functie / onderdeel | Wat | Soort |
|---|---|---|---|
| `R/FLSAM.control.R` | klasse `FLSAM.control` | nieuw slot `hockeyStickCurve`; alle slots gedocumenteerd | U, D |
| | `FLSAM.control()` | `fixVarToWeight` per vloot (met vlootnamen), `hockeyStickCurve` | U |
| | `drop.from.control()` | `fixVarToWeight` per vloot meenemen; oude objecten upgraden | U |
| | `update()` | slotnamen met hoofdletter → koppelingen van gewichten/maturity/M werden niet hernummerd | B |
| | validity | herschreven: crashte (`plus.group` is vector, vloottypes 6/7 afgekeurd, foute dimensiecheck) | B |
| `R/FLR2SAM.R` | `ctrl2conf()` | `fixVarToWeight` werd altijd 0; `hockeyStickCurve` doorgeven | B, U |
| | `conf2ctrl()` | `logN.vars` werd matrix; `plus.group`, `hockeyStickCurve` gezet | B, U |
| | `FLSAM2SAM()` | landings-/discardgewicht uit `landings.wt`/`discards.wt` (val terug op `catch.wt` bij NA); argument `fleets=` | B |
| `R/FLSAM.R` | klasse `FLSAM` | nieuw slot `fit` (klasse `sam` of `NULL`) | U |
| | `FLSAM()` | `keep.fit=TRUE`; `set.pars` voor alle vaste parameters; `hockeyStickCurve` | U |
| | `setParsTable()`, `.setPars()`, `.setParsSAM()` | vertaling FLSAM- ↔ SAM-namen; deterministische TMB-map | N |
| | `FLSAMs`, `sam.fitfast()` e.d. | alleen documentatie | D |
| `R/SAM2FLR.R` | `SAM2FLR()` | `keep.fit`; `plogis` i.p.v. `boot::inv.logit`; upgrade oude control | U |
| | `SIM2FLR()` | haalde verkeerde parameters weg bij meerdere vaste waarden in één parameter | B |
| `R/FLSAM.R` | `startingValues()`, `FLSAM(..., starting.values=)` | startwaarden uit een eerder assessment (ook met een jaar minder/meer data; logN/logF op jaar uitgelijnd; ook uit objecten zonder fit); vervangt de pin-functionaliteit | N |
| `R/pin.R` | `FLSAM2par()` | code ongewijzigd (pin-functionaliteit wordt niet meer gebruikt); alleen een documentatiestub (intern) | D |
| `R/accessors.R` | `components<-`, `info<-` | `as.integer()` vernielde de waarden | B |
| | overige accessors | documentatie | D |
| `R/methods.R` | `rec.pars()` (via `.extract.fleet.parameters`) | gebruikte globale `ctrl`, dubbele rijen, max. 2 parameters | B |
| | `mohns.rho()` | defaults werkten niet, teken omgedraaid t.o.v. Mohn (1999), gemiddelde als attribuut | B |
| | overige methods | documentatie | D |
| `R/methods_plots.R` | `plot()` (7 methodes), `fpart.plot()`, `cor.plot()`, `procerr.plot()`, `retroResiduals()`, `retroSelectivity()`, `retroParams()` | lattice → ggplot2 (geven `ggplot` terug) | G |
| | `retroParams()` | parameternamen bijgewerkt (oude ADMB-namen) | B |
| | `plot(FLSAMs, FLStock)` | `futureYrs` had geen default | B |
| | `plotkobe()` | gebruikte globale `NSH.sam` | B |
| | `obscor.plot()` | overbodige `require(ellipse)` weg | – |
| | `comp.plot()`, `otolith()`, `obsvar.plot()`, `obscv.plot()` | ongewijzigd (base graphics), alleen documentatie | D |
| `R/looi_retro.R` | `looi()` zonder `type` | faalde | B |
| | `retro()` voor `FLStocks` | eerste peel negeerde `set.pars` | B |
| `R/FLSAM.MSE.R` | `FLSAM.MSE()` | herschreven: draaide niet (`parallel`, `it`, `%dofuture%` ongedefinieerd); nu per iteratie, optioneel parallel (`ncores`) | B |
| | `FLSAM.fast()` | faalde als surveys verder lopen dan vangsten | B |
| | `.fillStock()` | hulpfunctie | N |
| `R/monteCarloStock.R` | `monteCarloStock()`, `monteCarloStockMscan()` | `set.pars`-tak (ongedefinieerde `stcks`/`catch.vars`), niet-geconvergeerde realisaties, typfout `lengt(`, `saveNLLs=NULL` | B |
| `R/diagnostics.R` | `residual.diagnostics()` | loess-smoother werd nooit getekend | B |
| `R/kobe-FLSAM.R` | `kobeFLSamFn()` | overbodige `require(ellipse)` weg; documentatie | D |
| `R/zzz.R` | `.onAttach()` | eiste exact 0.12.0, verwees naar oude SAM-branch | B |
| `R/FLSAM.out.R` | `FLSAM.out()` | upgrade oud FLSAM-object vóór de `validObject`-check | U |
| `R/af.R`, `R/anf.r` | – | alleen documentatie | D |
| `R/upgrade.R` | `upgradeFLSAM()`, `samfit()`, `.upgradeControl()`, `show`/`print`-methodes | oude objecten (zonder nieuwe slots) upgraden, ook automatisch bij printen; toegang tot de SAM-fit | N |
| `R/sam_tools.R` | `FLSAM.forecast()`, `FLSAM.modelforecast()`, `FLSAM.refpts()`, `FLSAM.ypr()`, `FLSAM.catchbyfleet()`, `model.description()` | wrappers rond stockassessment, geven data.frames terug | N |
| `R/diagnostics_sam.R` | `residual.tests()`, `jitter.test()`, `obscor.empirical.plot()`, `obssd.plot()`, `predsd.plot()`, `data.plot()`, `fselectivity.plot()` | diagnostics zoals in SAM, in ggplot2 | N |
| `R/colours.R` | `FLSAM.colours()`, `.col()` | PFA-kleuren als standaardpalet van alle ggplot2-plots (`options(FLSAM.palette=)`) | N |
| `R/methods_plots.R`, `R/diagnostics_sam.R` | ggplot2-plots | kleuren uit `FLSAM.colours()` | U |
| `dev/her27/` | `run_version.R`, `run_extras.R`, `compare_her27.R` | vergelijking op de HAWG 2026 haring-assessment (sf + mf, retro 10): 116/116 identiek | N |
| `R/FLSAM-package.R`, `R/docs.R`, `R/data.R` | – | package-, groeps- en datadocumentatie; imports | N |
| `DESCRIPTION` | – | versie 3.1.0, Imports (+ggplot2, stats, graphics, grDevices, parallel), Suggests, Collate | U |
| `NAMESPACE` | – | nu gegenereerd door roxygen2 (`exportPattern` behouden) | D |
| `man/` | – | 19 oude handgeschreven Rd's (deels uit het ADMB-tijdperk) vervangen door 55 gegenereerde | D |
| `data/*.RData` | `NSH.ctrl`, `NSH.sam`, `SOL.ctrl`, `SOL.sam`, `HERAS.sams` | geüpgraded naar de nieuwe klassedefinities (inhoud ongewijzigd) | U |
| `NEWS.md`, `README.md` | – | changelog; versie en installatie | D |
| `dev/compare_versions.R`, `dev/overzicht_wijzigingen.md` | – | testscript en dit overzicht (uitgesloten van de build) | N |

## Bugs (stap 6)

Alle onderstaande bugs zijn gefixt, conform je opdracht.

Al gemeld in de inventarisatie: `fixVarToWeight` genegeerd; `catch.wt` voor landings/discards;
`update()`-slotnamen; `conf2ctrl` `logN.vars`; validity; `mohns.rho` (default + teken); `looi()`
zonder `type`; `FLSAM.MSE`; `FLSAM.fast`; `rec.pars`; `monteCarloStock`; `residual.diagnostics`.

Tijdens het werk extra gevonden:
- `components<-` en `info<-` vernielden de waarden (`as.integer()` op een matrix).
- `plotkobe()` gebruikte een globale `NSH.sam`.
- `retro()` voor `FLStocks` negeerde `set.pars` in de eerste peel.
- `SIM2FLR()` verwijderde verkeerde parameters als er meerdere in één parametervector vastgezet waren.
- `monteCarloStockMscan()`: typfout `lengt(`; de default `ncores=detectCores()-1` faalde als `parallel` niet geladen was.
- `plot(FLSAMs, FLStock)` had geen default voor `futureYrs`.
- `simulate()` kon niet zonder `set.pars` worden aangeroepen; met `set.pars` faalde het bij meer dan één simulatie (`SIM2FLR` haalde per simulatie opnieuw parameters weg).
- `sam.out()` riep zichzelf aan (oneindige recursie) in plaats van `FLSAM.out()`.
- `obscor.plot()` gaf een onduidelijke foutmelding zonder gecorreleerde observaties.
- `.onAttach` stopte bij elke stockassessment-versie ≠ 0.12.0.
- In de oude `set.pars` werd de volgorde van de parameters in de TMB-map bepaald door `jitter()`, dus willekeurig.

Rectificatie: in een tussenbericht noemde ik een ontbrekende import van `mvrnorm` in `simulate()`.
Dat klopt niet: die loopt via FLCore en werkte gewoon.

## Oude objecten (zonder de nieuwe slots `hockeyStickCurve` en `fit`)
- Werken zonder aanpassing in: `FLSAM()`, `retro()`, `looi()`/`loo()`/`loi()`, `drop.from.control()`,
  `update()`, `ctrl2conf()`, `SAM2FLR()`, `FLSAM.out()`, alle accessors en plots, en bij printen.
  Dit is getest met de oorspronkelijke `NSH`, `NSH.sam` en `HERAS.sams` uit versie 3.0.0.9005
  (deel B8 van het testscript).
- `upgradeFLSAM(object)` upgradet expliciet (`FLSAM.control`, `FLSAM`, `FLSAMs`); daarna opnieuw `save()`n.
- Alleen een directe `validObject(oud_object)` faalt nog (interne S4-check); na `upgradeFLSAM()` werkt die wel.
- De data-objecten in het package zelf zijn al geüpgraded.

## Bevindingen in SAM zelf (niet in FLSAM opgelost)
- `stockassessment::hcr()` faalt in 0.12.0 met "Forecast type not implemented".
- `stochasticReferencepoints()` wordt in 0.12.0 niet geëxporteerd en `referencepoints()` is defunct.
  Daarom zijn er geen wrappers voor `hcr` en stochastische referentiepunten; die werken pas in `F_TAC`.
- `modelDescription()` herhaalt de tekst per vloot als `maxAgePlusGroup` vlootspecifiek is.
  De wrapper `model.description()` vangt dat af.

## Verschillen in de vergelijking (deel A), alle verklaard
- `disMeanWeight` (NSH): nu `discards.wt` i.p.v. `catch.wt` (bugfix). Dit heeft geen effect op de fit.
- `mohns.rho`: zelfde waarden, omgekeerd teken (bugfix).
- `set.pars`-run:
  - schattingen verschillen ≤ 1e-10 relatief;
  - `vcov` bevat dezelfde waarden in een andere volgorde (de willekeurige volgorde in de oude versie);
  - de one-step-ahead residuen van het eerste datajaar (1992) verschillen. Die zijn numeriek instabiel en hangen af van de parametervolgorde. Alle andere jaren zijn gelijk.

## Niet gedaan / open punten
- `pin.R` en `FLSAM(..., starting.values=)` zijn op jouw verzoek ongewijzigd gelaten (worden niet meer
  gebruikt). `starting.values` werkt daardoor nog steeds niet (`updateStart()` bestaat niet); verwijderen kan desgewenst.
- Base-graphics plots (`residual.diagnostics`, `otolith`, `comp.plot`, `obsvar.plot`, `obscv.plot`,
  `obscor.plot`) zijn niet omgezet: de opdracht was lattice → ggplot2.
- `R CMD check`-notes: `require()` van Suggests-pakketten (`plyr`, `kobe`, `doParallel`, `RColorBrewer`,
  `reshape2`) in bestaande functies; ongewijzigd gelaten.
- `sam.fitfast()` wordt nergens gebruikt en is niet getest tegen stockassessment 0.12.
- De pkgdown-site in `docs/` is niet opnieuw gegenereerd.
