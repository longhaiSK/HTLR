SEED <- 1234

set.seed(SEED)
dat <- gendata_MLR(n = 100, p = 10)
dat <- split_data(dat$X, dat$y, p.train = 0.8)

set.seed(SEED)
fit.with.warmup <- htlr(X = dat$x.tr, y = dat$y.tr, iter = 10, keep.warmup.hist = T)

set.seed(SEED)
fit.wout.warmup <- htlr(X = dat$x.tr, y = dat$y.tr, iter = 10, keep.warmup.hist = F) 

test_that("predict() can handle fit with warmup records", {
  expect_equal(predict(fit.wout.warmup, dat$x.te), predict(fit.with.warmup, dat$x.te))
})

test_that("as.matrix() can handle fit with warmup records", {
  # include.warmup wouldn't have effect if warmup record is not available
  expect_equal(as.matrix(fit.with.warmup, include.warmup = F, k = 1), as.matrix(fit.wout.warmup, include.warmup = T, k = 1))
  expect_equal(as.matrix(fit.with.warmup, include.warmup = F, k = 1), as.matrix(fit.wout.warmup, include.warmup = F, k = 1))
})
