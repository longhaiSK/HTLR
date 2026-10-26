test_that("binary vectors and one-column matrices produce identical fits", {
  X <- cbind(a = seq_len(12), b = rep(c(-1, 1), 6))
  y <- rep(1:2, each = 6)
  deltas <- c(0.1, -0.2, 0.3)
  fit <- function(init) htlr(X, y, init = init, iter = 2, warmup = 1,
                             leap = 1, leap.warm = 1)
  set.seed(42)
  expect_message(vector_fit <- fit(deltas), "skip standardizing")
  set.seed(42)
  expect_message(matrix_fit <- fit(matrix(deltas, ncol = 1)), "skip standardizing")
  expect_identical(vector_fit$mcdeltas, matrix_fit$mcdeltas)
  expect_identical(vector_fit$feature, matrix_fit$feature)
  expect_identical(vector_fit$mc.param$init, "customized")
})

test_that("invalid numeric vectors fail with actionable errors", {
  X <- matrix(seq_len(24), 12, 2)
  binary <- rep(1:2, each = 6)
  expect_error(htlr_fit(X, binary, initial_state = c(0, 0)),
               "Expected: nrow=3, ncol=1; Actual: nrow=2, ncol=1", fixed = TRUE)
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(htlr_fit(X, binary, initial_state = c(0, bad, 0)),
                 "finite coefficients")
  }
  expect_error(htlr_fit(X, rep(1:3, each = 4), initial_state = rep(0, 6)),
               "only for binary responses")
})

test_that("custom matrix validation preserves coefficients and metadata", {
  for (K in 1:2) {
    deltas <- matrix(seq_len(3 * K) / 10, 3, K,
                     dimnames = list(c("Intercept", "b", "a"),
                                     paste0("class", seq_len(K))))
    expect_identical(validate_initial_matrix(deltas, 2, K), deltas)
  }
  expect_identical(validate_initial_matrix(matrix(0L, 3, 1), 2, 1),
                   matrix(0L, 3, 1))
})

test_that("custom matrices require numeric finite coefficients and exact shape", {
  for (bad in list(matrix("0", 3, 1), matrix(TRUE, 3, 1),
                   matrix(1i, 3, 1))) {
    expect_error(validate_initial_matrix(bad, 2, 1), "numeric matrix")
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    deltas <- matrix(0, 3, 1)
    deltas[2, 1] <- bad
    expect_error(validate_initial_matrix(deltas, 2, 1), "finite coefficients")
  }
  expect_error(validate_initial_matrix(matrix(0, 2, 1), 2, 1),
               "Expected: nrow=3, ncol=1; Actual: nrow=2, ncol=1", fixed = TRUE)
  expect_error(validate_initial_matrix(matrix(0, 3, 1), 2, 2),
               "Expected: nrow=3, ncol=2; Actual: nrow=3, ncol=1", fixed = TRUE)
})

test_that("htlr_fit validates custom matrices before standardization", {
  X <- matrix(seq_len(24), 8, 3)
  y <- rep(1:2, each = 4)
  # A character matrix previously entered the standardization path.
  expect_error(htlr_fit(X, y, initial_state = matrix("0", 4, 1)),
               "numeric matrix")
  expect_error(htlr_fit(X, y, initial_state = matrix(Inf, 4, 1)),
               "finite coefficients")
  expect_error(htlr_fit(X, y, fsel = c(3, 1),
                        initial_state = matrix(0, 4, 1)),
               "Expected: nrow=3, ncol=1", fixed = TRUE)
})

test_that("binary and multiclass matrices retain positional fitting behavior", {
  X <- cbind(a = seq_len(12), b = rep(c(-1, 1), 6),
             c = rep(1:3, 4))
  for (K in 1:2) {
    y <- rep(seq_len(K + 1), length.out = nrow(X))
    deltas <- matrix(seq_len(3 * K) / 100, 3, K)
    fit <- function(stdzx) htlr_fit(
      X, y, fsel = c(3, 1), initial_state = deltas, stdzx = stdzx, ptype = "t",
      iters_h = 0, iters_rmc = 1, leap_L = 1, leap_L_h = 1)
    set.seed(4)
    expect_warning(expect_message(standardized <- fit(TRUE), "skip standardizing"),
                   "positional compatibility")
    set.seed(4)
    expect_warning(unstandardized <- fit(FALSE), "positional compatibility")
    expect_equal(standardized$mcdeltas, unstandardized$mcdeltas)
    expect_identical(unname(standardized$feature$fsel), c(3, 1))
    expect_equal(standardized$feature$X, cbind(Intercept = 1, X[, c(3, 1)]))
    expect_equal(standardized$feature$nuj, c(0, 0))
    expect_equal(standardized$feature$sdj, c(1, 1))
    expect_identical(standardized$mc.param$init, "customized")
    expect_equal(dim(standardized$mcdeltas)[1:2], c(3, K))
  }
})
