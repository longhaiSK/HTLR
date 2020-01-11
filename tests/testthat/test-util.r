set.seed(1234)

mat_80x90 <- matrix(rnorm(n = 80 * 90, mean = 100, sd = 50), 80, 90)
mat_800x900 <- matrix(rnorm(n = 800 * 900, mean = 100, sd = 50), 800, 900)

test_that("std works", {
  mat_80x90_nuj <- apply(mat_80x90, 2, median)
  mat_80x90_sdj <- apply(mat_80x90, 2, sd)
  mat_80x90_std_exp <- sweep(mat_80x90, 2, mat_80x90_nuj, "-")
  mat_80x90_std_exp <- sweep(mat_80x90_std_exp, 2, mat_80x90_sdj, "/")  
  
  mat_80x90_std_act <- std(mat_80x90)
  expect_equal(c(mat_80x90_std_act), c(mat_80x90_std_exp))
  expect_equal(attr(mat_80x90_std_act, "center"), mat_80x90_nuj)
  expect_equal(attr(mat_80x90_std_act, "scale"), mat_80x90_sdj)
  expect_equal(attr(mat_80x90_std_act, "nonsingular"), 1L:ncol(mat_80x90))
})

## compute V (delta)
comp_vardeltas_r <- function(deltas)
{
  K <- ncol(deltas)
  SUMdeltas <- rowSums(deltas)
  SUMsqdeltas <- rowSums(deltas^2)
  SUMsqdeltas - SUMdeltas^2 / (K + 1)
}

test_that("comp_vardeltas works", {
  expect_equal(c(comp_vardeltas(mat_80x90)), comp_vardeltas_r(mat_80x90))
})

comp_lsl_r <- function(lv)
{
  log_sum_exp(cbind(0, lv))
}

test_that("comp_lsl works", {
  expect_equal(comp_lsl(mat_80x90), comp_lsl_r(mat_80x90))
})

log_normcons_r <- function(lv)
{
  sum(comp_lsl_r(lv))
}

test_that("log_normcons works", {
  expect_equal(log_normcons(mat_80x90), log_normcons_r(mat_80x90))
})
