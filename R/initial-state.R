as.htlr.init <- function(x) UseMethod("as.htlr.init")

#' @exportS3Method
as.htlr.init.cv.glmnet <- function(x)
{
  coefs <- coef(x, s = "lambda.min") |> Reduce(f = cbind) |> as.matrix()
  deltas <- coefs[, -1, drop = FALSE] - coefs[, 1]
  return (deltas)
}

resolve_fsel <- function(fsel, X)
{
  if (is.null(fsel)) fsel <- seq_len(ncol(X))
  if (is.character(fsel) && is.null(dim(fsel))) {
    validate_init_names(colnames(X), "X column names")
    validate_init_names(fsel, "'fsel'")
    missing <- fsel[!fsel %in% colnames(X)]
    if (length(missing))
      stop(sprintf("Unknown 'fsel' feature names: %s.", paste(missing, collapse = ", ")),
           call. = FALSE)
    fsel <- match(fsel, colnames(X))
  }
  if (!is.numeric(fsel) || is.complex(fsel) || !is.null(dim(fsel)) ||
      !length(fsel) || any(!is.finite(fsel)) ||
      any(fsel != floor(fsel)) || any(fsel < 1 | fsel > ncol(X)) ||
      anyDuplicated(fsel))
    stop("'fsel' must select at least one feature using unique positive integer indices within ncol(X), or valid column names.",
         call. = FALSE)
  unname(fsel)
}

#' Describe External Initial Coefficients
#'
#' @param deltas A finite numeric matrix. The first row is the intercept and
#' remaining rows correspond to \code{features}, in order.
#' @param features Unique feature names or positive original column indices.
#' @param classes Optional names of the non-reference response classes, in
#' coefficient-column order. The reference class is the first level of
#' \code{factor(y)} and has no coefficient column.
#' @return An \code{htlr.init} object with coefficients, feature and class
#' metadata, and \code{intercept = 1L} identifying the intercept row.
#' @details Features may be a subset of the final fitted features specified by
#' \code{fsel}. Omitted feature coefficients start at zero in every class;
#' the supplied intercept is retained. Initializer features outside \code{fsel}
#' are rejected. Direct coefficient matrices still require all fitted features.
#' Coefficients use the scale of the supplied predictors and skip feature standardization.
#' Metadata is authoritative; matrix dimnames, when supplied, must agree with it.
#' @export
#' @examples
#' init <- htlr_init(matrix(0, 3, 1), features = c("gene_b", "gene_a"),
#'                   classes = "case")
#' # Initialize only gene_b while fitting both features.
#' X <- cbind(gene_a = seq_len(8), gene_b = rep(c(-1, 1), 4))
#' y <- factor(rep(c("control", "case"), each = 4),
#'             levels = c("control", "case"))
#' init <- htlr_init(matrix(c(0.1, 0.2), 2, 1), "gene_b", "case")
#' set.seed(1)
#' fit <- htlr(X, y, init = init, stdx = FALSE,
#'             iter = 2, warmup = 1, leap = 1, leap.warm = 1)
htlr_init <- function(deltas, features, classes = NULL)
{
  valid_features <- (is.character(features) &&
    !anyNA(features) && all(nzchar(features))) ||
    (is.numeric(features) && !is.complex(features) &&
     all(is.finite(features)) && all(features > 0 & features == floor(features)))
  if (!valid_features || anyDuplicated(features))
    stop("'features' must be unique nonempty names or positive integer column indices.",
         call. = FALSE)
  if (!is.matrix(deltas))
    stop("'deltas' must be a numeric matrix.", call. = FALSE)
  deltas <- validate_initial_matrix(deltas, length(features), ncol(deltas))
  if (ncol(deltas) < 1)
    stop("'deltas' must have at least one class column.", call. = FALSE)
  if (!is.null(colnames(deltas)))
    validate_init_names(colnames(deltas), "Initializer column names")
  if (!is.null(classes)) {
    validate_init_names(classes, "Initializer classes")
    if (length(classes) != ncol(deltas))
      stop("'classes' must have one name per non-reference coefficient column.", call. = FALSE)
    if (!is.null(colnames(deltas)) && !identical(colnames(deltas), classes))
      stop("Coefficient column names must agree with 'classes' in order.", call. = FALSE)
  }
  if (!is.null(rownames(deltas))) {
    rn <- rownames(deltas)
    validate_init_names(rn, "Initializer row names")
    if (!rn[1] %in% c("Intercept", "(Intercept)") ||
        !identical(rn[-1], as.character(features)))
      stop("Coefficient row names must identify the intercept first, followed by 'features' in order.",
           call. = FALSE)
  }
  structure(list(deltas = deltas, features = features, classes = classes,
                 intercept = 1L), class = "htlr.init")
}

validate_init_names <- function(x, label)
{
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
    stop(paste0(label, " must be unique, nonmissing, nonempty names."), call. = FALSE)
}

match_init_names <- function(actual, expected, label)
{
  validate_init_names(actual, paste(label, "names"))
  validate_init_names(expected, paste("Expected", label, "names"))
  if (length(actual) != length(expected) || !setequal(actual, expected))
    stop(sprintf("%s mismatch. Expected: %s; Actual: %s.", label,
                 paste(expected, collapse = ", "), paste(actual, collapse = ", ")),
         call. = FALSE)
  match(expected, actual)
}

align_initial_state <- function(initial_state, X, fsel, classes)
{
  explicit <- inherits(initial_state, "htlr.init")
  if (explicit) {
    if (!identical(initial_state$intercept, 1L))
      stop("Initializer 'intercept' must identify the first row (1L).", call. = FALSE)
    # Revalidate objects in case their metadata has been edited after construction.
    initial_state <- htlr_init(initial_state$deltas, initial_state$features,
                               initial_state$classes)
    deltas <- initial_state$deltas
    features <- initial_state$features
    class_names <- initial_state$classes
    if (is.null(class_names)) class_names <- colnames(deltas)
  } else {
    deltas <- initial_state
    features <- NULL
    class_names <- colnames(deltas)
  }
  deltas <- validate_initial_matrix(deltas,
                                    if (explicit) length(features) else length(fsel),
                                    length(classes))
  if (explicit && is.numeric(features)) {
    selected <- if (is.character(fsel)) match(fsel, colnames(X)) else unname(fsel)
    deltas <- expand_initial_features(deltas, features, selected)
  } else if (explicit || !is.null(rownames(deltas))) {
    validate_init_names(colnames(X), "X column names")
    expected <- colnames(X)[if (is.character(fsel)) match(fsel, colnames(X)) else fsel]
    if (!explicit) {
      rn <- rownames(deltas)
      validate_init_names(rn, "Initializer row names")
      intercept <- which(rn %in% c("Intercept", "(Intercept)"))
      if (length(intercept) != 1L)
        stop("Named coefficients must have exactly one 'Intercept' or '(Intercept)' row.", call. = FALSE)
      features <- rn[-intercept]
      deltas <- deltas[c(intercept, setdiff(seq_len(nrow(deltas)), intercept)), , drop = FALSE]
    }
    if (explicit) {
      validate_init_names(expected, "Selected feature names")
      deltas <- expand_initial_features(deltas, features, expected)
    } else {
      order <- match_init_names(features, expected, "Initializer features")
      deltas <- deltas[c(1L, order + 1L), , drop = FALSE]
    }
  } else if (!identical(unname(fsel), seq_len(ncol(X))) &&
             !identical(unname(fsel), as.numeric(seq_len(ncol(X))))) {
    warning("Unnamed initial coefficients use positional compatibility with non-default 'fsel'; use named coefficients or htlr_init() to identify features.",
            call. = FALSE)
  }
  if (!is.null(class_names)) {
    order <- match_init_names(class_names, classes, "Initializer classes")
    deltas <- deltas[, order, drop = FALSE]
  }
  deltas
}

expand_initial_features <- function(deltas, features, selected)
{
  if (anyNA(selected) || anyDuplicated(selected))
    stop("Final selected features must be nonmissing and unique.", call. = FALSE)
  outside <- features[!features %in% selected]
  if (length(outside))
    stop(sprintf("Initializer features outside final 'fsel': %s. Selected features: %s.",
                 paste(outside, collapse = ", "), paste(selected, collapse = ", ")),
         call. = FALSE)
  rn <- if (!is.null(rownames(deltas))) c(rownames(deltas)[1], as.character(selected))
  expanded <- matrix(0, length(selected) + 1L, ncol(deltas),
                     dimnames = if (!is.null(dimnames(deltas))) list(rn, colnames(deltas)))
  expanded[1, ] <- deltas[1, ]
  expanded[match(features, selected) + 1L, ] <- deltas[-1, , drop = FALSE]
  expanded
}
