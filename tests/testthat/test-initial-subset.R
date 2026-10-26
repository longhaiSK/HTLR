test_that("subset initializers zero-fill in final feature and class order", {
  X <- matrix(seq_len(36), 12, 3, dimnames = list(NULL, c("a", "b", "c")))
  for (K in 1:2) {
    classes <- as.character(seq_len(K) + 1)
    d <- matrix(seq_len(2 * K) / 10, 2, K)
    expected <- matrix(0, 4, K)
    expected[1, ] <- d[1, rev(seq_len(K))]
    expected[3, ] <- d[2, rev(seq_len(K))]
    for (features in list("c", 3L)) {
      init <- htlr_init(d, features, rev(classes))
      actual <- align_initial_state(init, X, c(2, 3, 1), classes)
      expect_identical(actual, expected)
      expect_identical(init$deltas, d)
      expect_identical(init$features, features)
      expect_identical(init$classes, rev(classes))
    }
  }
})

test_that("named subset rows and class metadata survive expansion", {
  X <- matrix(1:24, 8, 3, dimnames = list(NULL, c("a", "b", "c")))
  d <- matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2,
              dimnames = list(c("(Intercept)", "b"), c("3", "2")))
  actual <- align_initial_state(htlr_init(d, "b"), X, 3:1, c("2", "3"))
  expected <- matrix(c(0.3, 0, 0.4, 0, 0.1, 0, 0.2, 0), 4, 2,
                     dimnames = list(c("(Intercept)", "c", "b", "a"), c("2", "3")))
  expect_identical(actual, expected)
})

test_that("subset metadata rejects features outside the final model", {
  X <- matrix(1:24, 8, 3, dimnames = list(NULL, c("a", "b", "c")))
  for (features in list("c", 3L, "missing", 4L)) {
    init <- htlr_init(matrix(c(0.1, 0.2), 2, 1), features)
    expect_error(align_initial_state(init, X, 1:2, "2"),
                 "Initializer features outside final 'fsel'")
  }
  expect_error(align_initial_state(htlr_init(matrix(0, 2, 2), "a"), X, 1:3, "2"),
               "Expected: nrow=2, ncol=1", fixed = TRUE)
  expect_error(align_initial_state(htlr_init(matrix(0, 2, 1), "a", "wrong"), X, 1:3, "2"),
               "Initializer classes mismatch")
  # Direct matrices cannot silently acquire subset semantics.
  expect_error(align_initial_state(matrix(0, 2, 1), X, 1:3, "2"),
               "Expected: nrow=4, ncol=1", fixed = TRUE)
})

test_that("intercept-only metadata initializes all slopes to zero", {
  X <- matrix(1:24, 8, 3, dimnames = list(NULL, c("a", "b", "c")))
  for (features in list(character(), integer())) {
    init <- htlr_init(matrix(c(0.1, 0.2), 1, 2), features)
    expected <- matrix(c(0.1, 0, 0, 0, 0.2, 0, 0, 0), 4, 2)
    expect_identical(align_initial_state(init, X, 1:3, c("2", "3")), expected)
  }
})

test_that("subset and explicitly zero-filled matrices produce identical draws", {
  X <- cbind(a = seq_len(12), b = rep(c(-1, 1), 6), c = rep(1:3, 4))
  for (K in 1:2) {
    y <- rep(seq_len(K + 1), length.out = 12)
    d <- matrix(seq_len(2 * K) / 100, 2, K)
    expanded <- matrix(0, 3, K)
    expanded[c(1, 3), ] <- d
    for (features in list("b", 2L)) {
      for (stdx in c(FALSE, TRUE)) {
        fit <- function(init) suppressMessages(htlr(
          X, y, fsel = c(3, 2), init = init, stdx = stdx,
          iter = 2, warmup = 1, leap = 1, leap.warm = 1))
        set.seed(42)
        actual <- fit(htlr_init(d, features))
        set.seed(42)
        expect_warning(expected <- fit(expanded), "positional compatibility")
        expect_identical(actual$mcdeltas, expected$mcdeltas)
        expect_identical(actual$feature, expected$feature)
        expect_identical(unname(actual$feature$fsel), c(3, 2))
        expect_identical(actual$feature$X, cbind(Intercept = 1, X[, c(3, 2)]))
        expect_identical(actual$mc.param$init, "customized")
      }
    }
  }
})
