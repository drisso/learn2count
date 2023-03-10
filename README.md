# The `learn2count` package

This package implements algorithms for structure learning of graphical models for count data.

The function `PCzinb` implements three algorithms to estimate the structure of a graph from the input data.

The function `simdata` can be used to simulate data.

## Installation

The preferred way to install the package is
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("drisso/learn2count")
```

## Usage

Please, see the [vignette](vignettes/intro.Rmd) for detailed examples of the package usage.

## Versions of this package

The analyses and figures of the Nguyen et al. (2023) paper were done with package version `0.1.3`, which can be found [here](https://github.com/drisso/learn2count/releases/tag/v0.1.3). Please use this version to reproduce the results of the paper.

For virtually all other uses, we recommend using the latest stable version of the package (corresponding to the `master` branch).

## References

Nguyen, Van den Berge, Chiogna, Risso (2023). Structure learning for zero- inflated counts, with an application to single-cell RNA sequencing data. _Annals of Applied Statistics_. [In Print.](https://www.e-publications.org/ims/submission/AOAS/user/submissionFile/55336?confirm=5ee2c2d1)
