skip_on_travis()
skip_on_cran()

load("bplrhmc_expect.rda")

set.seed(SEED)
data("colon")
dat <- split_data(colon$X, colon$y, p.train = 0.9)

set.seed(SEED)
mlr <- gendata_MLR(n = 50, p = 10)
dat3 <- split_data(mlr$X, mlr$y, p.train = 0.8)

expect_equal_predict <- function(actual, expected) {
  expect_equal(unname(actual$probs_pred), unname(expected$probs_pred))
}

expect_equal_htlr <- function(actual, expected)
{
  expect_equal(actual$mcdeltas, expected$mcdeltas)
  expect_equal(as.vector(actual$mclogw), expected$mclogw)
  expect_equal(as.vector(actual$mcloglike), expected$mcloglike)
  expect_equal(as.vector(actual$mcuvar), expected$mcuvar)
  expect_equal(as.vector(actual$mchmcrej), expected$mchmcrej)
  expect_equal(actual$mcsigmasbt, expected$mcsigmasbt)
  expect_equal(actual$mcvardeltas, expected$mcvardeltas)
  expect_equal_predict(actual, expected)
}

test_that("colon500_t1_s10_eta10_bcbc_bi", {
  set.seed(SEED)
  fit <- htlr_fit(X_tr = dat$x.tr, y_tr = dat$y.tr, X_ts = dat$x.te,
                  fsel = 1L:500, ptype = "t", alpha = 1, s = -10, eta = 10,
                  initial_state = "bcbcsfrda", iters_h = 5, iters_rmc = 5, thin = 1)
  expect_equal_htlr(fit, colon500_t1_s10_eta10_bcbc_bi)
})

test_that("sim_t1_s10_eta10_bcbc_mul", {
  set.seed(SEED)
  fit <- htlr_fit(X_tr = dat3$x.tr, y_tr = dat3$y.tr, X_ts = dat3$x.te,
                  ptype = "t", alpha = 1, s = -10, eta = 10,
                  initial_state = "bcbcsfrda", iters_h = 5, iters_rmc = 5, thin = 1)
  expect_equal_htlr(fit, sim_t1_s10_eta10_bcbc_mul)
})

test_that("colon500_neg1_s10_eta0_bcbc_bi", {
  set.seed(SEED)
  fit <- htlr_fit(X_tr = dat$x.tr, y_tr = dat$y.tr, X_ts = dat$x.te,
                  fsel = 1L:500, ptype = "neg", alpha = 1, s = -10,
                  initial_state = "bcbcsfrda", iters_h = 5, iters_rmc = 5, thin = 1)
  expect_equal_htlr(fit, colon500_neg1_s10_eta0_bcbc_bi)
})

test_that("colon500_ghs1_s10_eta0_bcbc_bi", {
  set.seed(SEED)
  fit <- htlr_fit(X_tr = dat$x.tr, y_tr = dat$y.tr, X_ts = dat$x.te,
                  fsel = 1L:500, ptype = "ghs", alpha = 1, s = -10,
                  initial_state = "bcbcsfrda", iters_h = 5, iters_rmc = 5, thin = 1)
  expect_equal_htlr(fit, colon500_ghs1_s10_eta0_bcbc_bi)
})

# test_that("colon500_t1_s10_eta10_lasso_bi", {
#   set.seed(SEED)
#   fit <- htlr_fit(X_tr = dat$x.tr, y_tr = dat$y.tr, X_ts = dat$x.te,
#                   fsel = 1L:500, ptype = "t", alpha = 1, s = -10, eta = 10,
#                   initial_state = "lasso", iters_h = 0, iters_rmc = 1, thin = 1)
#   expect_equal_htlr(fit, colon500_t1_s10_eta10_lasso_bi)
# })

# HoSC <- list("x.train" = read.csv("../../data-raw/HoCS_train_data.csv", header = F),
#              "x.test" = read.csv("../../data-raw/HoCS_test_data.csv", header = F),
#              "y.train" = rep(1:3, each = 10),
#              "y.test" = c(rep(1, 50), rep(2, 27), rep(3, 52)))
# 
# test_that("data_with_singular_col", {
#   fit.1 <- htlr(X = HoSC$x.train, y = HoSC$y.train, init = "bcbc", iter = 10)
#   fit.2 <- htlr_fit(X_tr = HoSC$x.train, y_tr = HoSC$y.train, ptype = "t",
#                     initial_state = "lasso", iters_h = 5, iters_rmc = 5, thin = 1)
#   
#   expect_equal(length(fit.1$featur$fsel), 86)
#   expect_equal(length(fit.2$featur$fsel), 86)
#   
#   pred <- predict(fit.1, HoSC$x.test)
# })
