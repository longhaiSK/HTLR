# Run from the repository root with the patched HTLR.old package installed.
# See tests/testthat/fixtures/README.md before you replace the fixtures.
devtools::load_all()
stopifnot(packageVersion("HTLR.old") == package_version("3.1-1"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
seed <- 1001L
fields <- c("mcdeltas", "mclogw", "mcloglike", "mcuvar", "mchmcrej",
            "mcsigmasbt", "mcvardeltas", "probs_pred")

set.seed(seed)
data("colon", package = "HTLR")
binary <- split_data(colon$X, colon$y, p.train = 0.9)
set.seed(seed)
mlr <- gendata_MLR(n = 50, p = 10)
multiclass <- split_data(mlr$X, mlr$y, p.train = 0.8)
inputs <- list(binary = binary, multiclass = multiclass)
cases <- list(
  colon_t = list(input = "binary", prior = list(ptype = "t", eta = 10)),
  sim_t = list(input = "multiclass", prior = list(ptype = "t", eta = 10)),
  colon_neg = list(input = "binary", prior = list(ptype = "neg", eta = 0)),
  colon_ghs = list(input = "binary", prior = list(ptype = "ghs", eta = 0))
)
for (name in names(cases)) {
  case <- cases[[name]]
  dat <- inputs[[case$input]]
  args <- c(list(alpha = 1, s = -10, initial_state = "bcbcsfrda",
                 iters_h = 5, iters_rmc = 5, thin = 1), case$prior)
  if (case$input == "binary") args$fsel <- 1L:500
  set.seed(seed)
  fit <- do.call(HTLR.old::htlr_fit,
                 c(list(X_tr = dat$x.tr, y_tr = dat$y.tr, X_ts = dat$x.te), args))
  cases[[name]] <- list(input = case$input, args = args, expected = fit[fields])
}

set.seed(seed)
mlr_expected <- HTLR.old::htlr_gendata(n = 50, p = 10)
means <- rbind(c(0, 1, 0), c(0, 0, 0), c(0, 0, 1),
               c(0, 0, 1), c(0, 0, 1), matrix(0, 5, 3)) * 2
A <- diag(1, 10)
A[1:5, 1:3] <- rbind(c(1, 0, 0), c(2, 1, 0), c(0, 0, 1),
                     c(0, 0, 1), c(0, 0, 1))
fam_args <- list(n = 100, muj = means, A = A, sd_g = 0.5, stdx = TRUE)
set.seed(seed)
fam_expected <- do.call(HTLR.old::gendata_fam, fam_args)
fixture <- list(
  metadata = list(
    package = "HTLR.old", version = as.character(packageVersion("HTLR.old")),
    source_commit = "ac8e7b34bf3a039ae2042cedcd67693ec7379c62",
    seed = seed, rng_kind = RNGkind(), R = R.version.string,
    platform = R.version$platform,
    dependencies = vapply(c("HTLR", "BCBCSF", "glmnet", "Rcpp", "RcppArmadillo"),
                          function(pkg) as.character(packageVersion(pkg)), ""),
    session = sessionInfo()
  ),
  inputs = inputs, cases = cases,
  generators = list(
    mlr = list(args = list(n = 50, p = 10), expected = mlr_expected[c("X", "y", "deltas")]),
    fam = list(args = fam_args, expected = fam_expected[c("X", "y", "muj", "SGM")])
  )
)
saveRDS(fixture, "tests/testthat/fixtures/legacy-3.1-1.rds", version = 2, compress = "xz")
