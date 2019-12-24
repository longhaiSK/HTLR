set.seed(1234)

mat_80x90 <- matrix(rnorm(n = 80 * 90, mean = 100, sd = 50), 80, 90)

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
