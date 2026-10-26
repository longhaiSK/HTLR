library(HTLR.old)

# Run each callback with a fixed seed. The library must use a normal stack.
sample_sigma <- function(symbol, K, vardeltas) {
  .C(symbol, p = as.integer(length(vardeltas)),
     sigmasbt = rep(1, length(vardeltas)), vardeltas = vardeltas,
     K = as.integer(K), alpha = 1, log_aw = -3,
     PACKAGE = "HTLR.old")$sigmasbt
}

for (symbol in c("spl_sgm_neg", "spl_sgm_ghs")) {
  for (K in c(1L, 2L)) {
    set.seed(1001)
    actual <- sample_sigma(symbol, K, c(0.2, 1, 4))
    set.seed(1001)
    expected <- sample_sigma(symbol, K, c(0.2, 1, 4))
    stopifnot(identical(actual, expected),
              all(is.finite(actual)), all(actual > 0))
  }
}

sample_scale <- function(eta) {
  .C("samplew", logw = -3, K = 2L, p = 3L,
     vardeltas = c(0.2, 1, 4), nu = 1, s = -3, eta = eta,
     PACKAGE = "HTLR.old")$logw
}

set.seed(1001)
actual <- sample_scale(1)
set.seed(1001)
expected <- sample_scale(1)
stopifnot(identical(actual, expected), is.finite(actual), actual != -3,
          identical(sample_scale(0), -3))
