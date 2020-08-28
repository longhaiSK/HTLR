skip_on_cran()

load("gendata_expect.rda")

test_that("sim MLR works", {
  set.seed(SEED)
  actual <- gendata_MLR(n = 50, p = 10)
  expect_equal(unname(actual$X), dat.mlr$X)
  expect_equal(actual$y, dat.mlr$y)
  expect_equal(actual$deltas, dat.mlr$deltas)
})

test_that("sim FAM works", {
  n <- 100
  p <- 10
  
  means <- rbind(
    c(0, 1, 0),
    c(0, 0, 0),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1)
  ) * 2
  means <- rbind(means, matrix(0, p - 5, 3))
  
  A <- diag(1, p)
  A[1:5, 1:3] <- rbind(
    c(1, 0, 0),
    c(2, 1, 0),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1)
  )
  set.seed(SEED)
  actual <- gendata_FAM(n, means, A, sd_g = 0.5, stdx = TRUE)
  
  expect_equal(unname(actual$X), dat.fam$X)
  expect_equal(actual$y, dat.fam$y)
  expect_equal(actual$muj, dat.fam$muj)
  expect_equal(actual$SGM, dat.fam$SGM)
})
