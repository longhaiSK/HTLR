test_that("rankings are complete permutations and leave screening to fsel", {
  X <- cbind(noise = rep(c(-1, 1), 20), signal = rep(c(0, 3), each = 20) +
               rep(seq_len(20) / 100, 2), weak = seq_len(40))
  y <- rep(1:2, each = 20)
  expect_identical(order_plain(X, y), 1:3)
  expect_identical(order_plain(X[, FALSE, drop = FALSE], y), integer())
  for (rank_fn in list(order_ftest, order_kruskal)) {
    ranking <- rank_fn(X, y)
    expect_identical(sort(ranking), 1:3)
    expect_identical(ranking[1], 2L)
    set.seed(1)
    fit <- htlr(X, y, fsel = ranking[1:2], init = "random",
                iter = 2, warmup = 1, leap = 1, leap.warm = 1)
    expect_identical(unname(fit$feature$fsel), ranking[1:2])
  }
})

test_that("lasso initialization keeps every input feature and uses real metadata", {
  set.seed(2)
  X <- matrix(rnorm(80 * 5), 80, 5, dimnames = list(NULL, letters[1:5]))
  for (C in 2:3) {
    y <- rep(seq_len(C), length.out = 80)
    set.seed(3)
    d <- lasso_deltas(X, y, lambda = c(0.1, 0.05))
    expect_identical(rownames(d), c("(Intercept)", colnames(X)))
    expect_identical(colnames(d), as.character(seq_len(C))[-1])
    expect_equal(dim(d), c(6, C - 1))
    expect_true(all(is.finite(d)))
    expect_identical(align_initial_state(d, X, 1:5, colnames(d)), d)
    set.seed(3)
    unnamed <- lasso_deltas(unname(X), y, lambda = c(0.1, 0.05))
    expect_null(rownames(unnamed))
    expect_equal(unname(unnamed), unname(d))
    expect_silent(align_initial_state(unnamed, unname(X), 1:5, colnames(d)))
  }
})

test_that("legacy lasso screening is explicit and zero-fills omitted initial slopes", {
  set.seed(4)
  X <- matrix(rnorm(80 * 5), 80, 5, dimnames = list(NULL, letters[1:5]))
  y <- rep(1:2, each = 40)
  selected <- c(4L, 2L)
  ranking <- function(X, y) c(selected, 5L, 3L, 1L)
  set.seed(5)
  expect_warning(legacy <- lasso_deltas(X, y, lambda = c(0.1, 0.05),
                                        rank_fn = ranking, k = 2), "deprecated")
  set.seed(5)
  subset <- lasso_deltas(X[, selected], y, lambda = c(0.1, 0.05))
  expect_equal(legacy[c(1, selected + 1), , drop = FALSE], subset)
  expect_equal(unname(legacy[c(2, 4, 6), , drop = FALSE]), matrix(0, 3, 1))
  expect_identical(rownames(legacy), c("(Intercept)", colnames(X)))
  init <- htlr_init(subset, colnames(X)[selected])
  expect_equal(align_initial_state(init, X, 1:5, "2"), legacy)
  fit <- function(init) htlr(X, y, init = init, stdx = FALSE,
                            iter = 2, warmup = 1, leap = 1, leap.warm = 1)
  set.seed(6)
  full_fit <- fit(legacy)
  set.seed(6)
  subset_fit <- fit(init)
  expect_identical(subset_fit$mcdeltas, full_fit$mcdeltas)
  expect_identical(unname(subset_fit$feature$fsel), 1:5)
  for (k in list(0, -1, 6, 1.5, NA_real_))
    expect_error(suppressWarnings(lasso_deltas(X, y, k = k)), "'k' must be")
  expect_error(suppressWarnings(lasso_deltas(X, y, k = 2,
                                             rank_fn = function(X, y) c(1, 1))),
               "'fsel' must select")
  expect_error(suppressWarnings(lasso_deltas(X, y, k = 2,
                                             rank_fn = function(X, y) 1:2)),
               "rank every input feature")
})

test_that("BCBCSF internal feature selection does not select final coefficient rows", {
  X <- matrix(seq_len(24) / 100, 8, 3, dimnames = list(NULL, c("a", "b", "c")))
  y <- rep(1:2, each = 4)
  testthat::local_mocked_bindings(
    bcbcsf_fitpred = function(...) list(fit_bcbcsf = list(fsel = c(3L, 1L))),
    bcbcsf_sumfit = function(...) list(muj = matrix(c(0.1, 0.2, 0.3, 0.4), 2, 2)))
  d <- bcbcsf_deltas(X, y)
  expect_identical(rownames(d), c("(Intercept)", colnames(X)))
  expect_identical(colnames(d), "2")
  expect_equal(dim(d), c(4, 1))
  expect_equal(unname(d[-1, 1]), c(0.2, 0, 0.2))
  expect_identical(align_initial_state(d, X, 1:3, "2"), d)
})
