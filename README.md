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

The analyses and figures of the Nguyen et al. (2022) paper were done with package version `0.3.0`, which can be found [here](https://github.com/drisso/learn2count/releases/tag/v0.3.0). Please use this version to reproduce the results of the paper.

For virtually all other uses, we recommend using the latest stable version of the package (corresponding to the `master` branch).

## References

Nguyen, Van den Berge, Chiogna, Risso (2023). [Structure learning for zero- inflated counts, with an application to single-cell RNA sequencing data](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-17/issue-3/Structure-learning-for-zero-inflated-counts-with-an-application-to/10.1214/23-AOAS1732.full). _Annals of Applied Statistics_.

Nguyen, Chiogna, Risso, Banzato (2024). Guided structure learning of DAGs for count data. _Statistical Modelling. In print_. [Preprint](https://doi.org/10.48550/arXiv.2206.09754).
