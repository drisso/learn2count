# Version 0.3.0 (2023-03-23)

- Added methods for guided structure learning of DAGs

# Version 0.2.1 (2023-03-10)

- Use `rowPair` slot to store graph rather than metadata

# Version 0.2.0 (2022-03-14)

- Major restructuring of the package
- New function `simdata` now substitutes the three original simulation functions
  (`nbinom.Simdata`, `pois.simdata`, `zinb.simdata`)
- Function `simdata` does not set seed internally unlike its predecessor
- The older simulation functions have been retained (for now) to test the new function
- Added unit tests
- Improved documentation
- Addded S4 methods for a wrapper function, `PCzinb`, which should be the main function used by users
- Added SummarizedExperiment support
- Added vignette
- Added new function, `QPtransform` to pre-process real data.
- Added new function, `prediction_scores` to evaluate the estimated graph.
