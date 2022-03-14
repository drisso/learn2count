This package implements algorithms for structure learning of graphical models for count data.

The function `PCzinb` implements three algorithms to estimate the structure of a graph from the input data.

The function `simdata` can be used to simulate data.

Please, see the [vignette](vignettes/intro.Rmd) for detailed examples of the package usage.

The preferred way to install the package is
```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("learn2count")
```
