test_that("fsel resolves source columns without changing supplied order", {
  X <- matrix(1:24, 8, 3, dimnames = list(NULL, c("a", "b", "c")))
  expect_identical(resolve_fsel(NULL, X), 1:3)
  expect_identical(resolve_fsel(c("c", "a"), X), c(3L, 1L))
  expect_identical(resolve_fsel(c(3, 1), X), c(3, 1))
  expect_identical(resolve_fsel(c(first = 3L, second = 1L), X), c(3L, 1L))
  expect_identical(resolve_fsel(2L, unname(X)), 2L)
  for (bad in list(integer(), c(1, 1), 0, -1, 4, 1.5, NA_real_, Inf,
                   TRUE, 1i, matrix(1, 1, 1), list(1)))
    expect_error(resolve_fsel(bad, X), "'fsel' must select")
  for (bad in list(c("a", "a"), c("a", NA), ""))
    expect_error(resolve_fsel(bad, X), "unique, nonmissing, nonempty")
  expect_error(resolve_fsel("missing", X), "Unknown 'fsel' feature names: missing", fixed = TRUE)
  expect_error(resolve_fsel("a", unname(X)), "X column names")
  colnames(X) <- c("a", "a", "c")
  expect_error(resolve_fsel("a", X), "X column names")
})

test_that("fit validates fsel before preprocessing and initialization", {
  X <- matrix(1:24, 8, 3)
  y <- rep(1:2, each = 4)
  expect_error(htlr(X, y, fsel = c(2, 2), init = "random"), "'fsel' must select")
  expect_error(htlr_fit(X, y, fsel = 1.5, initial_state = matrix(0, 2, 1)),
               "'fsel' must select")
  expect_error(htlr(X, y, stdx = NA), "'stdzx' must be TRUE or FALSE", fixed = TRUE)
  expect_error(htlr(X, y, fsel = "missing"), "X column names")
})

test_that("name selection and numeric selection give identical fits and predictions", {
  X <- cbind(a = seq_len(12), b = rep(c(-1, 1), 6), c = rep(1:3, 4))
  y <- rep(1:2, each = 6)
  fit <- function(fsel, init) htlr(X, y, fsel = fsel, init = init,
                                  iter = 4, warmup = 2, leap = 1, leap.warm = 1)
  for (init in list("random", htlr_init(matrix(c(0.1, 0.2), 2, 1), "a"))) {
    set.seed(21)
    by_name <- suppressMessages(fit(c("c", "a"), init))
    set.seed(21)
    by_index <- suppressMessages(fit(c(3L, 1L), init))
    expect_identical(by_name$mcdeltas, by_index$mcdeltas)
    expect_identical(by_name$feature, by_index$feature)
    expect_identical(unname(by_name$feature$fsel), c(3L, 1L))
    expect_equal(predict(by_name, X), predict(by_index, X))
  }
  set.seed(21)
  all_default <- htlr(X, y, init = "random", iter = 4, warmup = 2, leap = 1, leap.warm = 1)
  set.seed(21)
  all_null <- fit(NULL, "random")
  expect_identical(all_default$mcdeltas, all_null$mcdeltas)
  expect_identical(all_default$feature, all_null$feature)
})

test_that("standardization preserves selected order and removes singular columns", {
  X <- cbind(a = c(1:11, 100), constant = rep(4, 12), b = rep(c(-3, 7), 6))
  y <- rep(1:2, each = 6)
  set.seed(21)
  expect_message(fit <- htlr(X, y, fsel = c("b", "constant", "a"), init = "random",
                             iter = 2, warmup = 1, leap = 1, leap.warm = 1),
                 "dropped 1 singular")
  selected <- X[, c(3, 1)]
  centers <- apply(selected, 2, median)
  scales <- apply(selected, 2, sd)
  standardized <- sweep(sweep(selected, 2, centers, "-"), 2, scales, "/")
  expect_identical(unname(fit$feature$fsel), c(3L, 1L))
  expect_equal(fit$feature$nuj, unname(centers))
  expect_equal(fit$feature$sdj, unname(scales))
  expect_equal(fit$feature$X, cbind(Intercept = 1, standardized))
  expect_true(fit$feature$stdx)
  manual_X <- X
  manual_X[, c(3, 1)] <- standardized
  manual_fit <- fit
  manual_fit$feature$nuj <- c(0, 0)
  manual_fit$feature$sdj <- c(1, 1)
  expect_equal(predict(fit, X), predict(manual_fit, manual_X))
  expect_error(suppressMessages(htlr(X, y, fsel = 2, init = "random")),
               "No nonsingular features remain")
})

test_that("custom coefficients retain raw scale and accurate standardization metadata", {
  X <- cbind(a = c(1:11, 100), constant = rep(4, 12), b = rep(c(-3, 7), 6))
  y <- rep(1:2, each = 6)
  d <- matrix(c(0.01, 0.02, 0.03, 0.04), 4, 1,
              dimnames = list(c("Intercept", colnames(X)), NULL))
  for (init in list(d, htlr_init(matrix(c(0.01, 0.03), 2, 1), "constant"))) {
    fit <- function(stdx, X) htlr(X, y, init = init, stdx = stdx,
                                  iter = 4, warmup = 2, leap = 1, leap.warm = 1)
    set.seed(21)
    expect_message(actual <- fit(TRUE, as.data.frame(X)), "skip standardizing")
    set.seed(21)
    expected <- fit(FALSE, X)
    expect_identical(actual$mcdeltas, expected$mcdeltas)
    expect_identical(actual$feature, expected$feature)
    expect_false(actual$feature$stdx)
    expect_equal(actual$feature$nuj, rep(0, 3))
    expect_equal(actual$feature$sdj, rep(1, 3))
    expect_identical(actual$feature$X, cbind(Intercept = 1, X))
    expect_equal(predict(actual, X), predict(expected, X))
  }
})
