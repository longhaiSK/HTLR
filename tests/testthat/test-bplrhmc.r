load("bplrhmc_expect.rda")
set.seed(SEED)

data("colon")
dat <- split_data(colon$X, colon$y, p.train = 0.9)

expect_equal_htlr <- function(actual, expected)
{
  expect_equal(actual$mcdeltas, expected$mcdeltas)
  expect_equal(as.vector(actual$mclogw), expected$mclogw)
  expect_equal(as.vector(actual$mcloglike), expected$mcloglike)
  expect_equal(as.vector(actual$mcuvar), expected$mcuvar)
  expect_equal(as.vector(actual$mchmcrej), expected$mchmcrej)
  expect_equal(actual$mcsigmasbt, expected$mcsigmasbt)
  expect_equal(actual$mcvardeltas, expected$mcvardeltas)
  expect_equal(unname(actual$probs_pred), unname(expected$probs_pred))
}

test_that("t(1, -10), eta = 10, lasso, bcbc, bi", {
  set.seed(SEED)
  fit <- htlr_fit(X_tr = dat$x.tr, y_tr = dat$y.tr, X_ts = dat$x.te,
                  fsel = 1L:500, ptype = "t", alpha = 1, s = -10, eta = 10,
                  initial_state = "bcbcsfrda", iters_h = 10, iters_rmc = 10, thin = 1)
  expect_equal_htlr(fit, colon500_t1_s10_eta10_bcbc_bi)
})



