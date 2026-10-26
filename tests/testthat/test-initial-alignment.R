test_that("named coefficients align both axes including the intercept", {
  X <- matrix(1:24, 8, 3, dimnames = list(NULL, c("a", "b", "c")))
  deltas <- matrix(1:8 / 10, 4, 2,
                    dimnames = list(c("b", "(Intercept)", "c", "a"), c("3", "2")))
  aligned <- align_initial_state(deltas, X, 1:3, c("2", "3"))
  expect_identical(aligned, deltas[c(2, 4, 1, 3), c(2, 1)])
  expect_identical(align_initial_state(deltas, X[, c(3, 1, 2)], 1:3, c("2", "3")),
                   deltas[c(2, 3, 4, 1), c(2, 1)])
})

test_that("explicit names and indices align coefficients without inferring identity", {
  X <- matrix(1:24, 8, 3, dimnames = list(NULL, c("a", "b", "c")))
  deltas <- matrix(1:8 / 10, 4, 2)
  init <- htlr_init(deltas, c("c", "a", "b"), c("3", "2"))
  expect_s3_class(init, "htlr.init")
  expect_identical(init$features, c("c", "a", "b"))
  expect_identical(init$classes, c("3", "2"))
  expect_identical(init$intercept, 1L)
  expect_identical(align_initial_state(init, X, 1:3, c("2", "3")),
                   deltas[c(1, 3, 4, 2), c(2, 1)])
  indices <- htlr_init(deltas, c(3, 1, 2), c("3", "2"))
  expect_identical(align_initial_state(indices, unname(X), c(2, 3, 1), c("2", "3")),
                   deltas[c(1, 4, 2, 3), c(2, 1)])
  init$intercept <- 2L
  expect_error(align_initial_state(init, X, 1:3, c("2", "3")), "first row")
})

test_that("ambiguous feature, intercept and class metadata fail", {
  X <- matrix(1:16, 8, 2, dimnames = list(NULL, c("a", "b")))
  d <- matrix(0, 3, 1, dimnames = list(c("Intercept", "a", "b"), "2"))
  for (rn in list(c("Intercept", "a", "a"), c("Intercept", "a", ""),
                  c("Intercept", "a", NA_character_))) {
    rownames(d) <- rn
    expect_error(align_initial_state(d, X, 1:2, "2"), "unique, nonmissing, nonempty")
  }
  rownames(d) <- c("offset", "a", "b")
  expect_error(align_initial_state(d, X, 1:2, "2"), "exactly one")
  rownames(d) <- c("Intercept", "(Intercept)", "b")
  expect_error(align_initial_state(d, X, 1:2, "2"), "exactly one")
  rownames(d) <- c("Intercept", "a", "wrong")
  expect_error(align_initial_state(d, X, 1:2, "2"),
               "Expected: a, b; Actual: a, wrong", fixed = TRUE)
  rownames(d) <- c("Intercept", "a", "b")
  colnames(d) <- "wrong"
  expect_error(align_initial_state(d, X, 1:2, "2"), "Initializer classes mismatch")
  colnames(d) <- "2"
  expect_error(align_initial_state(d, unname(X), 1:2, "2"), "X column names")
  colnames(X) <- c("a", "a")
  expect_error(align_initial_state(d, X, 1:2, "2"), "X column names")
  for (features in list(c("a", "a"), c("a", NA), c(1, 1), c(0, 2), c(1.5, 2)))
    expect_error(htlr_init(unname(d), features), "unique nonempty names")
  expect_error(htlr_init(unname(d), c("a", "b"), c("2", "3")), "one name per")
  expect_error(htlr_init(d, c("b", "a")), "row names must")
  expect_error(htlr_init(d, c("a", "b"), "wrong"), "column names must")
  expect_error(htlr_init(matrix(Inf, 3, 1), 1:2), "finite coefficients")
  expect_error(htlr_init(matrix(0, 2, 1), 1:2), "mismatch data")
  bad <- matrix(0, 3, 2, dimnames = list(NULL, c("2", "2")))
  expect_error(align_initial_state(bad, X, 1:2, c("2", "3")), "unique, nonmissing")
})

test_that("unnamed rows retain compatibility with an explicit warning for fsel", {
  X <- matrix(1:24, 8, 3)
  d <- matrix(0, 3, 1)
  expect_warning(actual <- align_initial_state(d, X, c(3, 1), "2"),
                 "positional compatibility")
  expect_identical(actual, d)
  expect_silent(align_initial_state(matrix(0, 4, 1), X, 1:3, "2"))
})

test_that("named and explicit initializers reach the sampler in aligned order", {
  X <- cbind(a = seq_len(12), b = rep(c(-1, 1), 6))
  for (K in 1:2) {
    y <- factor(rep(seq_len(K + 1), length.out = 12),
                labels = letters[seq_len(K + 1)])
    d <- matrix(seq_len(3 * K) / 100, 3, K,
                dimnames = list(c("Intercept", "a", "b"), letters[2:(K + 1)]))
    shuffled <- d[c(3, 1, 2), rev(seq_len(K)), drop = FALSE]
    explicit <- htlr_init(unname(d[c(1, 3, 2), rev(seq_len(K)), drop = FALSE]),
                          c("b", "a"), rev(letters[2:(K + 1)]))
    fit <- function(init) htlr(X, y, init = init, stdx = FALSE,
                               iter = 2, warmup = 1, leap = 1, leap.warm = 1)
    set.seed(23)
    reference <- fit(unname(d))
    for (init in list(shuffled, explicit)) {
      set.seed(23)
      actual <- fit(init)
      expect_identical(actual$mcdeltas, reference$mcdeltas)
      expect_identical(actual$feature, reference$feature)
      expect_identical(actual$mc.param$init, "customized")
    }
  }
})
