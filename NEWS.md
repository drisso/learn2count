# Version 0.2.0 (2022-03-08)

- Major restructuring of the package
- New function `simdata` now substitutes the three original simulation functions
  (`nbinom.Simdata`, `pois.simdata`, `zinb.simdata`)
- Function `simdata` does not set seed internally unlike its predecessor
