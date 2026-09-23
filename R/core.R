.hdzacp_cache <- new.env(parent = emptyenv())

#' Exact null moments of the coordinatewise rank statistic
#'
#' @param n Number of observations.
#' @param eta Boundary fraction; splits are scanned from `ceiling(n * eta)` to
#'   `n - ceiling(n * eta)`.
#' @return A data frame with split index `k`, exact null mean `mu`, and exact
#'   null standard deviation `sd`. Full vectors used internally are attributes.
#' @export
hdzacp_null_moments <- function(n, eta = 0.1) {
  n <- as.integer(n)
  mn <- .hdzacp_mn(n, eta)
  raw <- hdzacp_null_moments_cpp(n, mn)
  answer <- data.frame(k = raw$k, mu = raw$mu, sd = raw$sd)
  attr(answer, "mu_full") <- raw$mu_full
  attr(answer, "sd_full") <- raw$sd_full
  attr(answer, "table_bytes") <- raw$table_bytes
  answer
}

.hdzacp_moment_vectors <- function(n, eta) {
  key <- paste(n, format(eta, digits = 16), sep = ":")
  if (!exists(key, envir = .hdzacp_cache, inherits = FALSE)) {
    tab <- hdzacp_null_moments(n, eta)
    assign(key, list(
      table = tab,
      mu_full = attr(tab, "mu_full"),
      sd_full = attr(tab, "sd_full")
    ), envir = .hdzacp_cache)
  }
  get(key, envir = .hdzacp_cache, inherits = FALSE)
}

#' Compute the ZAS and ZAM scans without calibration
#'
#' @param x Numeric data matrix, with observations in rows and coordinates in
#'   columns.
#' @param eta Boundary fraction.
#' @param ties How to handle within-coordinate ties. The default errors because
#'   the exact moment theory assumes continuous margins. Random tie breaking is
#'   intended for permutation analyses only.
#' @param tie_seed Seed used when `ties = "random"`.
#' @return An object of class `hdzacp_scan` containing ZAS/ZAM statistics and
#'   split-after estimates.
#' @export
hdzacp_scan <- function(x, eta = 0.1, ties = c("error", "random"), tie_seed = NULL) {
  x <- .hdzacp_check_matrix(x)
  ties <- match.arg(ties)
  n <- nrow(x)
  mn <- .hdzacp_mn(n, eta)
  ranked <- .hdzacp_rank_matrix(x, ties, tie_seed)
  moments <- .hdzacp_moment_vectors(n, eta)
  raw <- hdzacp_scan_rank_cpp(
    ranked$ranks, seq_len(n), moments$mu_full, moments$sd_full, mn
  )
  answer <- list(
    statistics = c(ZAS = unname(raw[["S_sum"]]), ZAM = unname(raw[["S_max"]])),
    estimates = c(ZAS = as.integer(raw[["cp_sum"]]), ZAM = as.integer(raw[["cp_max"]])),
    estimate_fraction = c(ZAS = raw[["cp_sum"]] / n, ZAM = raw[["cp_max"]] / n),
    diagnostics = raw[c("dense_mid", "coord_mid_deficit", "coord_mid_raw")],
    n = n,
    p = ncol(x),
    eta = eta,
    m = mn,
    ties = ranked$ties,
    tie_seed = ranked$tie_seed,
    call = match.call()
  )
  class(answer) <- "hdzacp_scan"
  answer
}

#' Control numerical asymptotic calibration
#'
#' @param ou_paths Number of independently simulated OU paths.
#' @param ou_grid Positive even number of OU grid subintervals.
#' @param threads OpenMP threads used only for the OU simulation.
#' @param ou_seed OU-reference seed.
#' @param max_control Named arguments passed to [hdzacp_max_constants()].
#' @return A control list.
#' @export
hdzacp_asymptotic_control <- function(ou_paths = 20000L, ou_grid = 8192L,
                                       threads = 1L, ou_seed = 20260909L,
                                       max_control = list()) {
  list(
    ou_paths = as.integer(ou_paths), ou_grid = as.integer(ou_grid),
    threads = as.integer(threads), ou_seed = as.integer(ou_seed),
    max_control = max_control
  )
}

#' Simulate the Ornstein--Uhlenbeck reference for ZAS
#'
#' The reference distribution is the supremum of a stationary Gaussian OU
#' process on the interval determined by `eta`. Coarse and fine grids share the
#' same paths, which permits a direct discretization diagnostic.
#'
#' @param eta Boundary fraction.
#' @param paths Number of paths.
#' @param grid Positive even number of fine-grid subintervals.
#' @param threads OpenMP threads.
#' @param seed Integer seed.
#' @return An object of class `hdzacp_ou_reference` with coarse/fine maxima.
#' @export
hdzacp_ou_reference <- function(eta = 0.1, paths = 200000L, grid = 32768L,
                                 threads = 1L, seed = 20260909L) {
  .hdzacp_mn(100L, eta)
  paths <- as.integer(paths)
  grid <- as.integer(grid)
  if (paths < 1L || grid < 2L || grid %% 2L != 0L) {
    stop("`paths` must be positive and `grid` must be a positive even integer.", call. = FALSE)
  }
  maxima <- hdzacp_ou_max_cpp(
    paths, grid, log((1 - eta) / eta), as.integer(threads), as.integer(seed)
  )
  answer <- list(
    coarse = maxima[, "coarse"], fine = maxima[, "fine"], eta = eta,
    paths = paths, grid = grid, threads = as.integer(threads), seed = as.integer(seed)
  )
  class(answer) <- "hdzacp_ou_reference"
  answer
}

.hdzacp_ou_tail <- function(statistic, reference) {
  if (inherits(reference, "hdzacp_ou_reference")) reference <- reference$fine
  reference <- sort(as.double(reference))
  if (!length(reference) || any(!is.finite(reference))) {
    stop("The OU reference must contain finite simulated maxima.", call. = FALSE)
  }
  (1 + length(reference) - findInterval(statistic, reference)) /
    (length(reference) + 1)
}

#' Test for a high-dimensional distributional change point
#'
#' `ZAS` is the dense sum aggregation, `ZAM` is the sparse maximum aggregation,
#' and `ZAC` adaptively combines their evidence. The default uses a shared
#' whole-vector permutation orbit. Consequently each permutation keeps all
#' coordinates of an observation together and preserves cross-sectional
#' dependence. ZAC is ranked symmetrically within that same orbit; it is not
#' obtained by simply applying a Cauchy CDF to two observed permutation p-values.
#'
#' @param x Numeric data matrix with observations in rows.
#' @param methods Any subset of `"ZAC"`, `"ZAM"`, and `"ZAS"`.
#' @param calibration `"permutation"` (recommended), `"asymptotic"`, or its
#'   alias `"analytic"`.
#' @param eta Boundary fraction.
#' @param alpha Test level.
#' @param permutations Number of random permutations; the orbit has one
#'   additional observed ordering.
#' @param seed Integer seed passed explicitly to the C++ permutation generator.
#' @param cauchy_weight ZAS weight for asymptotic ZAC. Permutation ZAC follows
#'   the article's equal-weight symmetric orbit and therefore requires `0.5`.
#' @param ties Tie handling. Random tie breaking is permitted only with
#'   permutation calibration.
#' @param tie_seed Seed for random tie breaking.
#' @param ou_reference Optional object from [hdzacp_ou_reference()] or a numeric
#'   vector of fine-grid maxima.
#' @param asymptotic_control Control from [hdzacp_asymptotic_control()].
#' @param keep_permutations Retain component and Cauchy orbit statistics.
#' @param ... Additional arguments passed by the ZAC, ZAM, and ZAS convenience
#'   wrappers to `hdzacp_test()`.
#' @return An object of class `hdzacp_test`; `$results` is a tidy method table.
#' @export
hdzacp_test <- function(x, methods = c("ZAC", "ZAM", "ZAS"),
                         calibration = c("permutation", "asymptotic", "analytic"),
                         eta = 0.1, alpha = 0.05, permutations = 999L,
                         seed = 20260821L, cauchy_weight = 0.5,
                         ties = c("error", "random"), tie_seed = seed,
                         ou_reference = NULL,
                         asymptotic_control = hdzacp_asymptotic_control(),
                         keep_permutations = FALSE) {
  x <- .hdzacp_check_matrix(x)
  methods <- unique(toupper(methods))
  if (!length(methods) || any(!methods %in% c("ZAC", "ZAM", "ZAS"))) {
    stop("`methods` must be a nonempty subset of ZAC, ZAM, and ZAS.", call. = FALSE)
  }
  calibration <- match.arg(calibration)
  if (calibration == "analytic") calibration <- "asymptotic"
  ties <- match.arg(ties)
  if (ties == "random" && calibration != "permutation") {
    stop("Random tie breaking is supported only with permutation calibration.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must lie strictly between zero and one.", call. = FALSE)
  }
  n <- nrow(x)
  mn <- .hdzacp_mn(n, eta)
  ranked <- .hdzacp_rank_matrix(x, ties, tie_seed)
  moments <- .hdzacp_moment_vectors(n, eta)

  if (calibration == "permutation") {
    permutations <- as.integer(permutations)
    if (permutations < 1L) stop("`permutations` must be positive.", call. = FALSE)
    if (!isTRUE(all.equal(cauchy_weight, 0.5))) {
      stop("Permutation ZAC uses the article's equal weight, `cauchy_weight = 0.5`.", call. = FALSE)
    }
    raw <- hdzacp_test_rank_cpp(
      ranked$ranks, mn, permutations, as.integer(seed), moments$mu_full,
      moments$sd_full, isTRUE(keep_permutations)
    )
    p_sum <- raw$p_sum
    p_max <- raw$p_max
    p_cauchy <- raw$p_cauchy
    statistic_cauchy <- raw$T_cauchy
    calibration_label <- "whole-vector symmetric-orbit permutation"
    reference <- NULL
    constants <- NULL
  } else {
    if (ranked$ties) stop("Asymptotic calibration requires tie-free coordinates.", call. = FALSE)
    raw <- hdzacp_scan_rank_cpp(
      ranked$ranks, seq_len(n), moments$mu_full, moments$sd_full, mn
    )
    control <- asymptotic_control
    if (is.null(ou_reference)) {
      ou_reference <- do.call(hdzacp_ou_reference, c(list(eta = eta), list(
        paths = control$ou_paths, grid = control$ou_grid,
        threads = control$threads, seed = control$ou_seed
      )))
    }
    p_sum <- .hdzacp_ou_tail(raw[["S_sum"]], ou_reference)
    constants <- do.call(hdzacp_max_constants, c(
      list(n = n, p = ncol(x), eta = eta), control$max_control
    ))
    normalized_max <- (raw[["S_max"]] - constants$b) / constants$a
    p_max <- -expm1(-exp(-normalized_max))
    p_cauchy <- .hdzacp_cauchy_p(p_sum, p_max, cauchy_weight)
    statistic_cauchy <- cauchy_weight / tan(pi * .hdzacp_clip_p(p_sum)) +
      (1 - cauchy_weight) / tan(pi * .hdzacp_clip_p(p_max))
    calibration_label <- "asymptotic OU/Gumbel/Cauchy"
    reference <- ou_reference
  }

  cp_sum <- as.integer(raw[["cp_sum"]])
  cp_max <- as.integer(raw[["cp_max"]])
  selected <- if (p_sum <= p_max) "ZAS" else "ZAM"
  cp_cauchy <- if (selected == "ZAS") cp_sum else cp_max
  all_rows <- rbind(
    data.frame(method = "ZAC", statistic = statistic_cauchy, p.value = p_cauchy,
      reject = p_cauchy <= alpha, estimate = cp_cauchy,
      estimate_fraction = cp_cauchy / n, calibration = calibration_label,
      selected_component = selected),
    data.frame(method = "ZAM", statistic = raw[["S_max"]], p.value = p_max,
      reject = p_max <= alpha, estimate = cp_max, estimate_fraction = cp_max / n,
      calibration = calibration_label, selected_component = NA_character_),
    data.frame(method = "ZAS", statistic = raw[["S_sum"]], p.value = p_sum,
      reject = p_sum <= alpha, estimate = cp_sum, estimate_fraction = cp_sum / n,
      calibration = calibration_label, selected_component = NA_character_)
  )
  results <- all_rows[match(methods, all_rows$method), , drop = FALSE]
  rownames(results) <- NULL
  answer <- list(
    results = results,
    components = list(
      ZAS = c(statistic = raw[["S_sum"]], p.value = p_sum, estimate = cp_sum),
      ZAM = c(statistic = raw[["S_max"]], p.value = p_max, estimate = cp_max)
    ),
    settings = list(
      n = n, p = ncol(x), eta = eta, m = mn, alpha = alpha,
      calibration = calibration, permutations = if (calibration == "permutation") permutations else NA_integer_,
      seed = as.integer(seed), cauchy_weight = cauchy_weight,
      ties = ties, ties_found = ranked$ties, tie_seed = ranked$tie_seed
    ),
    asymptotic = if (calibration == "asymptotic") list(
      ou_reference = reference, max_constants = constants
    ) else NULL,
    permutation_statistics = if (calibration == "permutation" && keep_permutations) raw$orbit else NULL,
    call = match.call()
  )
  class(answer) <- "hdzacp_test"
  answer
}

#' @rdname hdzacp_test
#' @export
zac_test <- function(x, ...) hdzacp_test(x, methods = "ZAC", ...)

#' @rdname hdzacp_test
#' @export
zam_test <- function(x, ...) hdzacp_test(x, methods = "ZAM", ...)

#' @rdname hdzacp_test
#' @export
zas_test <- function(x, ...) hdzacp_test(x, methods = "ZAS", ...)

#' @export
print.hdzacp_scan <- function(x, ...) {
  cat("HDZACP uncalibrated rank scan\n")
  cat("n =", x$n, ", p =", x$p, ", eta =", format(x$eta), "\n")
  print(data.frame(method = names(x$statistics), statistic = unname(x$statistics),
    estimate = unname(x$estimates), estimate_fraction = unname(x$estimate_fraction)),
    row.names = FALSE)
  invisible(x)
}

#' @export
print.hdzacp_test <- function(x, ...) {
  cat("HDZACP high-dimensional distributional change-point test\n")
  cat("n =", x$settings$n, ", p =", x$settings$p,
      ", calibration =", x$settings$calibration, "\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
summary.hdzacp_test <- function(object, ...) object$results

#' @export
plot.hdzacp_test <- function(x, ...) {
  values <- -log10(pmax(x$results$p.value, .Machine$double.xmin))
  graphics::barplot(values, names.arg = x$results$method,
    ylab = expression(-log[10](italic(p))), col = "#2B6F9C", ...)
  graphics::abline(h = -log10(x$settings$alpha), lty = 2, col = "#B53A3A")
  invisible(x)
}



