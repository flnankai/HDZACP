.hdzacp_check_matrix <- function(x, min_n = 4L) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix with observations in rows.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < min_n || ncol(x) < 1L) {
    stop("`x` has too few rows or no columns.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("`x` must contain only finite values.", call. = FALSE)
  }
  x
}

.hdzacp_mn <- function(n, eta) {
  if (length(eta) != 1L || !is.finite(eta) || eta <= 0 || eta >= 0.5) {
    stop("`eta` must be a finite number strictly between 0 and 0.5.", call. = FALSE)
  }
  mn <- as.integer(ceiling(n * eta))
  if (mn < 1L || 2L * mn >= n) {
    stop("`eta` leaves no admissible split point.", call. = FALSE)
  }
  mn
}

.hdzacp_seed <- function(seed, i = 0L, j = 0L, offset = 0L) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max - 1L, 1L)
  modulus <- 2147483646
  as.integer((as.double(seed) + 104729 * as.double(i) +
    1000003 * as.double(j) + 1009 * as.double(offset)) %% modulus + 1)
}

.hdzacp_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.hdzacp_rank_matrix <- function(x, ties = c("error", "random"), seed = NULL) {
  ties <- match.arg(ties)
  duplicated_columns <- vapply(seq_len(ncol(x)), function(j) anyDuplicated(x[, j]) > 0L,
    logical(1L))
  if (any(duplicated_columns) && ties == "error") {
    stop(
      "Ties were found in ", sum(duplicated_columns), " coordinate(s). ",
      "The exact rank calibration assumes continuous margins. Use ",
      "`ties = \"random\"` only with permutation calibration, or preprocess the data.",
      call. = FALSE
    )
  }
  make <- function() {
    ans <- vapply(seq_len(ncol(x)), function(j) {
      as.integer(rank(x[, j], ties.method = if (ties == "random") "random" else "first"))
    }, integer(nrow(x)))
    if (!is.matrix(ans)) ans <- matrix(ans, ncol = ncol(x))
    ans
  }
  if (ties == "random") {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max - 1L, 1L)
    ranks <- .hdzacp_with_seed(seed, make())
  } else {
    ranks <- make()
  }
  list(ranks = ranks, ties = any(duplicated_columns), tie_seed = seed)
}

.hdzacp_clip_p <- function(p) pmin(1 - 1e-15, pmax(1e-15, p))

.hdzacp_cauchy_p <- function(p_sum, p_max, weight = 0.5) {
  if (length(weight) != 1L || !is.finite(weight) || weight <= 0 || weight >= 1) {
    stop("`cauchy_weight` must lie strictly between zero and one.", call. = FALSE)
  }
  score <- weight / tan(pi * .hdzacp_clip_p(p_sum)) +
    (1 - weight) / tan(pi * .hdzacp_clip_p(p_max))
  stats::pcauchy(score, lower.tail = FALSE)
}

.hdzacp_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.hdzacp_parallel_lapply <- function(X, FUN, cores = 1L, ...) {
  cores <- max(1L, as.integer(cores))
  if (cores == 1L || length(X) <= 1L) return(lapply(X, FUN, ...))
  cores <- min(cores, length(X))
  cluster <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, suppressPackageStartupMessages(library(HDZACP)))
  parallel::parLapply(cluster, X, FUN, ...)
}

.hdzacp_require <- function(package, source = NULL) {
  if (requireNamespace(package, quietly = TRUE)) return(invisible(TRUE))
  extra <- if (is.null(source)) "" else paste0(" Install it from ", source, ".")
  stop("Package `", package, "` is required for this method.", extra, call. = FALSE)
}

.hdzacp_method_order <- c("ZAC", "ZAM", "ZAS", "KDist", "HDD", "E-Divisive", "gSeg")

.hdzacp_result_row <- function(method, statistic, p.value, alpha, estimate,
                               calibration, selected_component = NA_character_) {
  data.frame(
    method = method,
    statistic = as.double(statistic),
    p.value = as.double(p.value),
    reject = is.finite(p.value) && p.value <= alpha,
    estimate = as.integer(estimate),
    estimate_fraction = as.double(estimate),
    calibration = calibration,
    selected_component = selected_component,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.hdzacp_segments <- function(n, changepoints) {
  cp <- sort(unique(as.integer(changepoints)))
  cp <- cp[is.finite(cp) & cp >= 1L & cp < n]
  ends <- c(cp, n)
  starts <- c(1L, cp + 1L)
  data.frame(start = starts, end = ends, length = ends - starts + 1L)
}

#' Official empirical-data links
#'
#' The package contains no empirical observations. This function returns the
#' landing pages and direct download links used by the article reproduction.
#'
#' @return A data frame of data-set names, landing pages, and download URLs.
#' @export
hdzacp_data_urls <- function() {
  data.frame(
    data = c("BrainCloud (GSE30272)", "Gas Sensor Array Drift"),
    landing_page = c(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30272",
      "https://archive.ics.uci.edu/dataset/270/gas+sensor+array+drift+dataset+at+different+concentrations"
    ),
    download_url = c(
      paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30272/",
        "matrix/GSE30272_series_matrix.txt.gz"),
      paste0("https://archive.ics.uci.edu/static/public/270/",
        "gas+sensor+array+drift+dataset+at+different+concentrations.zip")
    ),
    stringsAsFactors = FALSE
  )
}


