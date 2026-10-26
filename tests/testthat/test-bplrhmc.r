legacy_fixture <- readRDS(test_path("fixtures", "legacy-3.1-1.rds"))

for (case_name in names(legacy_fixture$cases)) {
  test_that(paste("sampler matches legacy reference for", case_name), {
    fixture <- legacy_fixture
    case <- fixture$cases[[case_name]]
    dat <- fixture$inputs[[case$input]]
    old_kind <- RNGkind()
    on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
    do.call(RNGkind, as.list(fixture$metadata$rng_kind))
    set.seed(fixture$metadata$seed)
    fit_case <- function() do.call(htlr_fit,
      c(list(X_tr = dat$x.tr, y_tr = dat$y.tr, X_ts = dat$x.te), case$args))
    if (case$args$ptype == "t") {
      fit <- fit_case()
    } else {
      expect_warning(fit <- fit_case(), "random logw currently only supports t prior",
                     fixed = TRUE)
    }
    for (field in names(case$expected)) {
      actual <- fit[[field]]
      expected <- case$expected[[field]]
      if (field %in% c("mclogw", "mcloglike", "mcuvar", "mchmcrej")) {
        actual <- as.vector(actual)
        expected <- as.vector(expected)
      }
      if (field == "probs_pred") {
        actual <- unname(actual)
        expected <- unname(expected)
      }
      # Allow small differences in floating-point results between platforms.
      expect_equal(actual, expected, tolerance = 1e-7, info = field)
    }
  })
}

for (ptype in c("t", "neg", "ghs")) {
  for (classes in c(2L, 3L)) {
    test_that(paste(ptype, "sampler has valid traces for", classes, "classes"), {
      X <- cbind(a = seq(-1, 1, length.out = 18), b = rep(c(-1, 0, 1), 6))
      y <- rep(seq_len(classes), length.out = nrow(X))
      set.seed(1001)
      fit_case <- function() htlr_fit(X, y, ptype = ptype, eta = 0,
        initial_state = "random", iters_h = 2, iters_rmc = 4,
        leap_L = 2, leap_L_h = 2, rep.legacy = FALSE)
      if (ptype == "t") {
        fit <- fit_case()
      } else {
        expect_warning(fit <- fit_case(), "random logw currently only supports t prior",
                       fixed = TRUE)
      }
      expect_equal(dim(fit$mcdeltas), c(3L, classes - 1L, 5L))
      for (field in c("mcdeltas", "mclogw", "mcloglike", "mcsigmasbt",
                      "mcvardeltas", "mcuvar", "mchmcrej"))
        expect_true(all(is.finite(fit[[field]])), info = field)
      expect_true(all(fit$mcsigmasbt > 0))
      expect_true(all(fit$mcvardeltas >= 0))
      expect_true(all(fit$mchmcrej >= 0 & fit$mchmcrej <= 1))
      probs <- predict(fit, X)
      expect_equal(dim(probs), c(nrow(X), classes))
      expect_true(all(is.finite(probs) & probs >= 0 & probs <= 1))
      expect_equal(rowSums(probs), rep(1, nrow(X)))

      # Check the saved log likelihood against the multinomial model.
      for (i in seq_len(dim(fit$mcdeltas)[3])) {
        deltas <- matrix(fit$mcdeltas[, , i], nrow = ncol(X) + 1)
        logits <- cbind(0, fit$feature$X %*% deltas)
        weights <- exp(logits - apply(logits, 1, max))
        probabilities <- weights / rowSums(weights)
        expected <- sum(log(probabilities[cbind(seq_along(y), y)]))
        expect_equal(as.vector(fit$mcloglike)[i], expected, tolerance = 1e-10)
      }
    })
  }
}
