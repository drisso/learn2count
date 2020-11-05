This package implements algorithms for structure learning of graphical models for count data. For estimating a structure from give data, use function PCPoisson.R for Poisson node conditional models; function PCnb-optim.R (or PCnbinom.R using glm approach) for Negative Binomial node conditional models; and function PCzinb1noT.R and PCnb0noT.R for zero-inflated negative binomial node conditional models.

For simulating a data set, see functions datasimPois.R, datasinNB.R, datasimzinb.R for Poisson node conditional models, Negative Binomial node conditional models, and zero-inflated negative binomial node conditional models respectively.

optimnb.R and optimzinb.R include all neccessary functions to estimate parameters in each model.
