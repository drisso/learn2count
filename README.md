This package implements algorithms for structure learning of graphical models for count data.

For estimating a structure from given data, use function:
- `pois.wald` for Poisson node conditional models; 
- `nbscale.noT` (or `nb.wald` using glm approach) for Negative Binomial node conditional models;
-  `zinb1.noT` and `zinb0.noT` for zero-inflated negative binomial node conditional models.

For simulating a data set, see functions `pois.simdata`, `nbinom.Simdata`, and `zinb.simdata` for Poisson node conditional models, Negative Binomial node conditional models, and zero-inflated negative binomial node conditional models, respectively.

To install the package, you can run the following in R:

```{r}
library(remotes)
install_github("drisso/learn2count")
```
