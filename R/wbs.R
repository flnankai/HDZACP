.hdzacp_wbs_intervals <- function(n, n_intervals, min_interval, include_full, seed) {
  n_intervals <- as.integer(n_intervals)
  min_interval <- as.integer(min_interval)
  if (n_intervals < 0L || min_interval < 4L || min_interval > n) {
    stop("Invalid WBS interval controls.", call. = FALSE)
  }
  max_length <- if (include_full) n - 1L else n
  if (n_intervals > 0L && min_interval > max_length) {
    stop("No random interval satisfies `min_interval` after reserving the full interval.",
      call. = FALSE)
  }
  lengths <- if (n_intervals > 0L) seq.int(min_interval, max_length) else integer()
  counts <- n - lengths + 1L
  total <- sum(as.double(counts))
  if (n_intervals > total) {
    stop("More random intervals were requested than are available.", call. = FALSE)
  }
  random <- if (n_intervals) .hdzacp_with_seed(seed, {
    ids <- sort(sample.int(total, n_intervals, replace = FALSE))
    cumulative <- cumsum(counts)
    length_index <- findInterval(ids - 1, cumulative) + 1L
    previous <- c(0, cumulative)[length_index]
    start <- as.integer(ids - previous - 1)
    len <- lengths[length_index]
    cbind(start = start, end = start + len)
  }) else matrix(integer(), 0L, 2L, dimnames = list(NULL, c("start", "end")))
  if (include_full) random <- rbind(random, c(start = 0L, end = n))
  random <- unique(random)
  random[order(random[, "end"] - random[, "start"], random[, "start"]), , drop = FALSE]
}

.hdzacp_validate_intervals <- function(intervals, n, min_interval) {
  intervals <- as.matrix(intervals)
  if (ncol(intervals) != 2L || !is.numeric(intervals) || any(!is.finite(intervals))) {
    stop("`intervals` must be a finite two-column numeric matrix.", call. = FALSE)
  }
  storage.mode(intervals) <- "integer"
  colnames(intervals) <- c("start", "end")
  if (any(intervals[, 1L] < 0L) || any(intervals[, 2L] > n) ||
      any(intervals[, 2L] - intervals[, 1L] < min_interval)) {
    stop("Every interval must use split-after boundaries in [0,n] and meet `min_interval`.",
      call. = FALSE)
  }
  unique(intervals)
}

.hdzacp_wbs_select <- function(evidence, method, n, threshold) {
  p_name <- switch(method, ZAS = "p_sum", ZAM = "p_max", ZAC = "p_cauchy")
  cp_name <- switch(method, ZAS = "cp_sum", ZAM = "cp_max", ZAC = "cp_cauchy")
  selected <- integer()
  recurse <- function(left, right) {
    eligible <- which(
      evidence$start >= left & evidence$end <= right &
      evidence[[cp_name]] > left & evidence[[cp_name]] < right
    )
    if (!length(eligible)) return(invisible(NULL))
    ordered <- eligible[order(
      evidence[[p_name]][eligible],
      evidence$end[eligible] - evidence$start[eligible],
      evidence$start[eligible]
    )]
    best <- ordered[1L]
    if (!is.finite(evidence[[p_name]][best]) || evidence[[p_name]][best] > threshold) {
      return(invisible(NULL))
    }
    cp <- as.integer(evidence[[cp_name]][best])
    if (cp %in% selected) return(invisible(NULL))
    selected <<- c(selected, cp)
    recurse(left, cp)
    recurse(cp, right)
    invisible(NULL)
  }
  recurse(0L, n)
  sort(unique(selected))
}

#' Wild binary segmentation with ZAC, ZAM, and ZAS
#'
#' Each random interval is sampled uniformly from all eligible intervals, and
#' the full interval can be added separately. A change point is the index of the
#' final observation before a distributional change. `stopping = "unadjusted"`
#' reproduces the article's simulation and empirical WBS convention. The safer
#' general default, `"bonferroni"`, adjusts by the number of tested intervals.
#' `"paper"` uses the article's analytic interval adjustment and truncated
#' component p-values for asymptotic ZAC.
#'
#' @param x Numeric observation-by-coordinate matrix.
#' @param methods Any subset of ZAC, ZAM, and ZAS.
#' @param calibration Calibration passed to [hdzacp_test()].
#' @param eta,alpha Boundary fraction and level.
#' @param n_intervals Number of random intervals, excluding the optional full
#'   interval.
#' @param min_interval Minimum interval length.
#' @param include_full Include `[0,n]`.
#' @param intervals Optional user-supplied split-after interval matrix. When
#'   supplied, random interval controls are ignored.
#' @param permutations Local permutation count.
#' @param stopping One of `"bonferroni"`, `"paper"`, or `"unadjusted"`.
#' @param seed Seed for interval generation and local tests.
#' @param ties,tie_seed Tie controls passed to [hdzacp_test()].
#' @param ou_reference,asymptotic_control Asymptotic controls.
#' @param keep_evidence Retain the interval-level table.
#' @param ... Additional arguments passed by the ZAC, ZAM, and ZAS convenience
#'   wrappers to `hdzacp_wbs()`.
#' @return An object of class `hdzacp_wbs`.
#' @export
hdzacp_wbs <- function(x, methods = c("ZAC", "ZAM", "ZAS"),
                        calibration = c("permutation", "asymptotic", "analytic"),
                        eta = 0.1, alpha = 0.05, n_intervals = 50L,
                        min_interval = 40L, include_full = TRUE,
                        intervals = NULL, permutations = 99L,
                        stopping = c("bonferroni", "paper", "unadjusted"),
                        seed = 20260824L, ties = c("error", "random"),
                        tie_seed = seed, ou_reference = NULL,
                        asymptotic_control = hdzacp_asymptotic_control(),
                        keep_evidence = TRUE) {
  x <- .hdzacp_check_matrix(x, min_n = 5L)
  methods <- unique(toupper(methods))
  if (!length(methods) || any(!methods %in% c("ZAC", "ZAM", "ZAS"))) {
    stop("`methods` must be a nonempty subset of ZAC, ZAM, and ZAS.", call. = FALSE)
  }
  calibration <- match.arg(calibration)
  if (calibration == "analytic") calibration <- "asymptotic"
  stopping <- match.arg(stopping)
  ties <- match.arg(ties)
  n <- nrow(x)
  if (is.null(intervals)) {
    intervals <- .hdzacp_wbs_intervals(
      n, n_intervals, min_interval, include_full, as.integer(seed)
    )
  } else {
    intervals <- .hdzacp_validate_intervals(intervals, n, min_interval)
  }
  if (!nrow(intervals)) stop("No intervals are available.", call. = FALSE)
  if (calibration == "asymptotic" && is.null(ou_reference)) {
    control <- asymptotic_control
    ou_reference <- hdzacp_ou_reference(
      eta = eta, paths = control$ou_paths, grid = control$ou_grid,
      threads = control$threads, seed = control$ou_seed
    )
  }
  rows <- vector("list", nrow(intervals))
  for (i in seq_len(nrow(intervals))) {
    left <- intervals[i, "start"]
    right <- intervals[i, "end"]
    local_x <- x[seq.int(left + 1L, right), , drop = FALSE]
    local_seed <- .hdzacp_seed(seed, i, right - left, 1L)
    fit <- hdzacp_test(
      local_x, methods = c("ZAC", "ZAM", "ZAS"), calibration = calibration,
      eta = eta, alpha = alpha, permutations = permutations, seed = local_seed,
      ties = ties, tie_seed = .hdzacp_seed(tie_seed, i, right - left, 2L),
      ou_reference = ou_reference, asymptotic_control = asymptotic_control
    )
    by_method <- fit$results[match(c("ZAS", "ZAM", "ZAC"), fit$results$method), ]
    rows[[i]] <- data.frame(
      start = left, end = right, length = right - left,
      p_sum = by_method$p.value[1L], p_max = by_method$p.value[2L],
      p_cauchy = by_method$p.value[3L],
      cp_sum = left + by_method$estimate[1L],
      cp_max = left + by_method$estimate[2L],
      cp_cauchy = left + by_method$estimate[3L],
      selected_component = by_method$selected_component[3L],
      stringsAsFactors = FALSE
    )
  }
  evidence <- do.call(rbind, rows)
  m_eff <- nrow(evidence)
  threshold <- if (stopping == "unadjusted") alpha else alpha / m_eff
  if (calibration == "permutation" && stopping != "unadjusted" &&
      threshold < 1 / (as.integer(permutations) + 1)) {
    warning(
      "The adjusted threshold is below the smallest attainable permutation p-value; ",
      "increase `permutations` or use a different stopping rule.", call. = FALSE
    )
  }
  if (stopping == "paper" && calibration == "asymptotic") {
    q_c <- alpha / (4 * m_eff)
    evidence$p_cauchy <- .hdzacp_cauchy_p(
      pmin(evidence$p_sum, 1 - q_c), pmin(evidence$p_max, 1 - q_c), 0.5
    )
  }
  changepoints <- setNames(lapply(methods, function(method) {
    .hdzacp_wbs_select(evidence, method, n, threshold)
  }), methods)
  segments <- lapply(changepoints, function(cp) .hdzacp_segments(n, cp))
  answer <- list(
    changepoints = changepoints, segments = segments,
    evidence = if (keep_evidence) evidence else NULL,
    intervals = intervals,
    thresholds = setNames(rep(threshold, length(methods)), methods),
    settings = list(
      n = n, p = ncol(x), eta = eta, alpha = alpha,
      calibration = calibration, permutations = as.integer(permutations),
      stopping = stopping, n_intervals = nrow(intervals), seed = as.integer(seed)
    ),
    call = match.call()
  )
  class(answer) <- "hdzacp_wbs"
  answer
}

#' @rdname hdzacp_wbs
#' @export
zac_wbs <- function(x, ...) hdzacp_wbs(x, methods = "ZAC", ...)

#' @rdname hdzacp_wbs
#' @export
zam_wbs <- function(x, ...) hdzacp_wbs(x, methods = "ZAM", ...)

#' @rdname hdzacp_wbs
#' @export
zas_wbs <- function(x, ...) hdzacp_wbs(x, methods = "ZAS", ...)

#' @export
print.hdzacp_wbs <- function(x, ...) {
  cat("HDZACP wild binary segmentation\n")
  cat("n =", x$settings$n, ", p =", x$settings$p,
      ", stopping =", x$settings$stopping, "\n")
  for (method in names(x$changepoints)) {
    locations <- x$changepoints[[method]]
    cat(method, ": ", if (length(locations)) paste(locations, collapse = ", ") else "none", "\n", sep = "")
  }
  invisible(x)
}

#' @export
plot.hdzacp_wbs <- function(x, ...) {
  methods <- names(x$changepoints)
  graphics::plot(c(1, x$settings$n), c(0.5, length(methods) + 0.5),
    type = "n", yaxt = "n", xlab = "Observation index (split after k)",
    ylab = "Method", ...)
  graphics::axis(2, at = seq_along(methods), labels = methods, las = 1)
  for (i in seq_along(methods)) {
    graphics::segments(1, i, x$settings$n, i, col = "grey80")
    if (length(x$changepoints[[i]])) {
      graphics::points(x$changepoints[[i]], rep(i, length(x$changepoints[[i]])),
        pch = 19, col = "#B53A3A")
    }
  }
  invisible(x)
}


