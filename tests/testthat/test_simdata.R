context("Test simulation function.")

test_that("simulation function works as expected", {

    set.seed(13124)
    p <- 4
    n <- 10
    B <- matrix(c(0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0), nrow=p)

    set.seed(123)
    expect_silent(zinb <- simdata(n, p, B, mu = 5, mu_noise=1, family="ZINB", theta=1, pi=0.2))
    set.seed(123)
    expect_silent(nb <- simdata(n, p, B, mu = 5, mu_noise=1, family="NB", theta=1))
    set.seed(123)
    expect_silent(poi <- simdata(n, p, B, mu = 5, mu_noise=1))

    set.seed(123)
    tmp_nb <- nbinom.Simdata(n, p, B, mu = 5, mu.nois=1,theta=1)
    set.seed(123)
    tmp_zinb <- zinb.simdata(n, p, B, mu = 5, mu.nois=1,theta=1,pi=0.2)
    set.seed(123)
    tmp_poi <- pois.simdata (n, p, B, lambda = 5, lambda.c=1)

    expect_equal(nb, tmp_nb)
    expect_equal(zinb, tmp_zinb)
    expect_equal(poi, tmp_poi)

})