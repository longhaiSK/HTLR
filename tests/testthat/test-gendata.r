skip_on_travis()
skip_on_cran()

SEED <- 1001

test_that("sim MLR works", {
  set.seed(SEED)
  expect <- HTLR.old::htlr_gendata(n = 50, p = 10)
  
  set.seed(SEED)
  actual <- gendata_MLR(n = 50, p = 10)
  
  expect_equal(unname(actual$X), expect$X)
  expect_equal(actual$y, expect$y)
  expect_equal(actual$deltas, expect$deltas)
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
  expect <- HTLR.old::gendata_fam(n, means, A, sd_g = 0.5, stdx = TRUE)
  
  set.seed(SEED)
  actual <- gendata_FAM(n, means, A, sd_g = 0.5, stdx = TRUE)
  
  expect_equal(unname(actual$X), expect$X)
  expect_equal(actual$y, expect$y)
  expect_equal(actual$muj, expect$muj)
  expect_equal(actual$SGM, expect$SGM)
})
