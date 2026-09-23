#' KDist single-change comparison
#'
#' Calls `KDist::kcpd_single(type = "e-dist", method = "permutation")` from the
#' KDist package (Zhang et al.; GitHub repository
#' `zhangxiany-tamu/KDist`).
#'
#' @inheritParams hdzacp_test
#' @return A one-row standardized result data frame.
#' @export
compare_kdist <- function(x, alpha = 0.05, permutations = 999L,
                           seed = 20260821L) {
  .hdzacp_require("KDist", "GitHub repository `zhangxiany-tamu/KDist`")
  x <- .hdzacp_check_matrix(x)
  fit <- KDist::kcpd_single(
    data = x, type = "e-dist", B = as.integer(permutations), alpha = alpha,
    seeds = as.integer(seed), method = "permutation", num_cores = 1L
  )
  estimate <- as.integer(fit$locations[1L])
  data.frame(
    method = "KDist", statistic = as.double(fit$statistic[1L]),
    p.value = as.double(fit$pvalue[1L]), p_value = as.double(fit$pvalue[1L]),
    reject = is.finite(fit$pvalue[1L]) && fit$pvalue[1L] <= alpha,
    estimate = estimate, estimate_fraction = estimate / nrow(x),
    calibration = "KDist permutation", package = "KDist",
    stringsAsFactors = FALSE
  )
}

.hdzacp_hdd_pvalue <- function(distance_matrix, observed, permutations, seed) {
  n <- nrow(distance_matrix)
  cp <- as.integer(observed[["changepoint"]])
  fixed_index <- cp - 1L
  permutation_statistics <- numeric(permutations)
  .hdzacp_with_seed(seed, {
    if (fixed_index >= 1L && fixed_index <= n) {
      remaining <- setdiff(seq_len(n), fixed_index)
      for (b in seq_len(permutations)) {
        order <- append(sample(remaining), fixed_index, after = fixed_index - 1L)
        permutation_statistics[b] <- hdzacp_hdd_stat_cpp(distance_matrix, order)[["statistic"]]
      }
    } else {
      for (b in seq_len(permutations)) {
        permutation_statistics[b] <- hdzacp_hdd_stat_cpp(
          distance_matrix, sample.int(n)
        )[["statistic"]]
      }
    }
  })
  (1 + sum(permutation_statistics >= observed[["statistic"]])) /
    (permutations + 1)
}

#' HDD single-change comparison
#'
#' Reproduces the `dist1`/`test_statistic` method of the HDDchangepoint package
#' (GitHub repository `rezadrikvandi/HDDchangepoint`) with an algebraically
#' equivalent compiled distance scan. The article's conditional permutation
#' rule and a plus-one Monte Carlo p-value are retained. This calibration is
#' therefore labelled `article-fast` rather than represented as an unmodified
#' call to the package's default testing wrapper.
#'
#' @inheritParams compare_kdist
#' @return A one-row standardized result data frame.
#' @export
compare_hdd <- function(x, alpha = 0.05, permutations = 999L,
                         seed = 20260821L) {
  .hdzacp_require("HDDchangepoint", "GitHub repository `rezadrikvandi/HDDchangepoint`")
  x <- .hdzacp_check_matrix(x)
  distance <- hdzacp_hdd_distance_cpp(x)
  observed <- hdzacp_hdd_stat_cpp(distance, seq_len(nrow(x)))
  p_value <- .hdzacp_hdd_pvalue(distance, observed, as.integer(permutations), seed)
  estimate <- as.integer(observed[["changepoint"]])
  data.frame(
    method = "HDD", statistic = as.double(observed[["statistic"]]),
    p.value = p_value, p_value = p_value, reject = p_value <= alpha,
    estimate = estimate, estimate_fraction = estimate / nrow(x),
    calibration = "article-fast conditional permutation plus-one",
    package = "HDDchangepoint", stringsAsFactors = FALSE
  )
}

#' E-Divisive single-change comparison
#'
#' Calls `ecp::e.divisive()` with the article's energy exponent and sequential
#' permutation settings.
#'
#' @inheritParams compare_kdist
#' @param eta Minimum segment fraction.
#' @return A one-row standardized result data frame.
#' @export
compare_edivisive <- function(x, alpha = 0.05, permutations = 999L,
                               seed = 20260821L, eta = 0.1) {
  .hdzacp_require("ecp")
  x <- .hdzacp_check_matrix(x)
  fit <- .hdzacp_with_seed(seed, ecp::e.divisive(
    X = x, sig.lvl = alpha, R = as.integer(permutations), k = NULL,
    min.size = .hdzacp_mn(nrow(x), eta), alpha = 1
  ))
  interior <- setdiff(fit$estimates, c(1L, nrow(x) + 1L))
  estimate <- if (length(interior)) as.integer(interior[1L]) else NA_integer_
  p_value <- if (length(fit$p.values)) as.double(fit$p.values[1L]) else NA_real_
  data.frame(
    method = "E-Divisive", statistic = NA_real_, p.value = p_value,
    p_value = p_value, reject = isTRUE(fit$k.hat > 1L), estimate = estimate,
    estimate_fraction = estimate / nrow(x), calibration = "sequential permutation",
    package = "ecp", stringsAsFactors = FALSE
  )
}

#' Graph-based gSeg single-change comparison
#'
#' Constructs an MST with `ade4::mstree(stats::dist(x))` and calls
#' `gSeg::gseg1(statistics = "m")` with permutation calibration.
#'
#' @inheritParams compare_edivisive
#' @return A one-row standardized result data frame.
#' @export
compare_gseg <- function(x, alpha = 0.05, permutations = 999L,
                          seed = 20260821L, eta = 0.1) {
  .hdzacp_require("gSeg")
  .hdzacp_require("ade4")
  x <- .hdzacp_check_matrix(x)
  mn <- .hdzacp_mn(nrow(x), eta)
  graph <- ade4::mstree(stats::dist(x))
  fit <- .hdzacp_with_seed(seed, {
    invisible(utils::capture.output(result <- gSeg::gseg1(
      n = nrow(x), E = graph, statistics = "m", n0 = mn,
      n1 = nrow(x) - mn, pval.appr = FALSE, skew.corr = FALSE,
      pval.perm = TRUE, B = as.integer(permutations)
    )))
    result
  })
  observed <- fit$scanZ$max.type$Zmax
  permuted <- fit$pval.perm$max.type$maxZs
  p_value <- (1 + sum(permuted >= observed)) / (length(permuted) + 1)
  estimate <- as.integer(fit$scanZ$max.type$tauhat)
  data.frame(
    method = "gSeg", statistic = as.double(observed), p.value = p_value,
    p_value = p_value, reject = p_value <= alpha, estimate = estimate,
    estimate_fraction = estimate / nrow(x), calibration = "permutation plus-one",
    package = "gSeg + ade4", stringsAsFactors = FALSE
  )
}

#' Compare all seven single-change methods used in the article
#'
#' @param x Numeric observation-by-coordinate matrix.
#' @param methods Any subset of ZAC, ZAM, ZAS, KDist, HDD, E-Divisive, and gSeg.
#' @inheritParams hdzacp_test
#' @return A standardized result data frame. The `package` column identifies
#'   every external implementation.
#' @export
hdzacp_compare <- function(x, methods = .hdzacp_method_order, eta = 0.1,
                            alpha = 0.05, permutations = 999L,
                            seed = 20260821L, ties = c("error", "random")) {
  x <- .hdzacp_check_matrix(x)
  ties <- match.arg(ties)
  methods <- unique(methods)
  aliases <- c("E-Div." = "E-Divisive", EDIVISIVE = "E-Divisive")
  methods[methods %in% names(aliases)] <- aliases[methods[methods %in% names(aliases)]]
  if (any(!methods %in% .hdzacp_method_order)) {
    stop("Unknown comparison method.", call. = FALSE)
  }
  rows <- list()
  z_methods <- intersect(c("ZAC", "ZAM", "ZAS"), methods)
  if (length(z_methods)) {
    zfit <- hdzacp_test(
      x, methods = z_methods, calibration = "permutation", eta = eta,
      alpha = alpha, permutations = permutations, seed = .hdzacp_seed(seed, 1L),
      ties = ties, tie_seed = .hdzacp_seed(seed, 1L, 1L)
    )$results
    zfit$p_value <- zfit$p.value
    zfit$package <- "HDZACP"
    rows[["Z"]] <- zfit[, c("method", "statistic", "p.value", "p_value", "reject",
      "estimate", "estimate_fraction", "calibration", "package")]
  }
  if ("KDist" %in% methods) rows[["KDist"]] <- compare_kdist(
    x, alpha, permutations, .hdzacp_seed(seed, 2L)
  )
  if ("HDD" %in% methods) rows[["HDD"]] <- compare_hdd(
    x, alpha, permutations, .hdzacp_seed(seed, 3L)
  )
  if ("E-Divisive" %in% methods) rows[["E-Divisive"]] <- compare_edivisive(
    x, alpha, permutations, .hdzacp_seed(seed, 4L), eta
  )
  if ("gSeg" %in% methods) rows[["gSeg"]] <- compare_gseg(
    x, alpha, permutations, .hdzacp_seed(seed, 5L), eta
  )
  answer <- do.call(rbind, rows)
  answer <- answer[match(methods, answer$method), , drop = FALSE]
  rownames(answer) <- NULL
  answer
}

.hdzacp_hdd_fixed_p <- function(distance, cp_estimate, permutations, seed) {
  n <- nrow(distance)
  split_after <- as.integer(cp_estimate) - 1L
  fixed_index <- split_after
  if (split_after < 2L || split_after >= n) return(1)
  distance_squared <- distance^2
  total_sum <- rowSums(distance)
  total_squared_sum <- rowSums(distance_squared)
  fixed_statistic <- function(order) {
    left <- order[seq_len(split_after)]
    left_sum <- rowSums(distance[, left, drop = FALSE])
    left_squared_sum <- rowSums(distance_squared[, left, drop = FALSE])
    right_size <- n - split_after
    mean(left_squared_sum / split_after +
      (total_squared_sum - left_squared_sum) / right_size -
      2 * (left_sum / split_after) * ((total_sum - left_sum) / right_size))
  }
  observed <- fixed_statistic(seq_len(n))
  permuted <- .hdzacp_with_seed(seed, {
    remaining <- setdiff(seq_len(n), fixed_index)
    vapply(seq_len(permutations), function(b) {
      order <- append(sample(remaining), fixed_index, after = fixed_index - 1L)
      fixed_statistic(order)
    }, numeric(1L))
  })
  (1 + sum(permuted >= observed)) / (permutations + 1)
}

.hdzacp_hdd_multiple <- function(x, alpha, permutations, n_intervals,
                                  min_segment, seed) {
  selected <- integer()
  recurse <- function(local_x, offset, recursion_index) {
    n <- nrow(local_x)
    if (n < 2L * min_segment + 1L) return(invisible(NULL))
    intervals <- .hdzacp_wbs_intervals(
      n, n_intervals, 2L * min_segment, TRUE,
      .hdzacp_seed(seed, recursion_index, n, 1L)
    )
    candidates <- lapply(seq_len(nrow(intervals)), function(i) {
      left <- intervals[i, "start"]
      right <- intervals[i, "end"]
      subset <- local_x[seq.int(left + 1L, right), , drop = FALSE]
      distance <- hdzacp_hdd_distance_cpp(subset)
      stat <- hdzacp_hdd_stat_cpp(distance, seq_len(nrow(subset)))
      c(cp = left + as.integer(stat[["changepoint"]]),
        statistic = stat[["statistic"]])
    })
    candidates <- as.data.frame(do.call(rbind, candidates))
    candidates <- candidates[candidates$cp > min_segment &
      candidates$cp <= n - min_segment, , drop = FALSE]
    if (!nrow(candidates)) return(invisible(NULL))
    cp <- as.integer(candidates$cp[which.max(candidates$statistic)[1L]])
    distance <- hdzacp_hdd_distance_cpp(local_x)
    p_value <- .hdzacp_hdd_fixed_p(
      distance, cp, permutations, .hdzacp_seed(seed, recursion_index, n, 2L)
    )
    if (!is.finite(p_value) || p_value > alpha) return(invisible(NULL))
    split_after <- cp - 1L
    if (split_after < min_segment || n - split_after < min_segment) return(invisible(NULL))
    global_cp <- offset + split_after
    selected <<- c(selected, global_cp)
    recurse(local_x[seq_len(split_after), , drop = FALSE], offset, 2L * recursion_index)
    recurse(local_x[seq.int(split_after + 1L, n), , drop = FALSE],
      global_cp, 2L * recursion_index + 1L)
    invisible(NULL)
  }
  recurse(x, 0L, 1L)
  sort(unique(selected))
}

.hdzacp_gseg_multiple <- function(x, alpha, permutations, min_segment, seed) {
  selected <- integer()
  recurse <- function(local_x, offset, node) {
    n <- nrow(local_x)
    if (n < 2L * min_segment + 1L) return(invisible(NULL))
    graph <- ade4::mstree(stats::dist(local_x))
    fit <- .hdzacp_with_seed(.hdzacp_seed(seed, node, n), {
      invisible(utils::capture.output(result <- gSeg::gseg1(
        n = n, E = graph, statistics = "m", n0 = min_segment,
        n1 = n - min_segment, pval.appr = FALSE, skew.corr = FALSE,
        pval.perm = TRUE, B = permutations
      )))
      result
    })
    observed <- fit$scanZ$max.type$Zmax
    permuted <- fit$pval.perm$max.type$maxZs
    p_value <- (1 + sum(permuted >= observed)) / (length(permuted) + 1)
    if (!is.finite(p_value) || p_value > alpha) return(invisible(NULL))
    cp <- as.integer(fit$scanZ$max.type$tauhat)
    if (cp < min_segment || n - cp < min_segment) return(invisible(NULL))
    global <- offset + cp
    selected <<- c(selected, global)
    recurse(local_x[seq_len(cp), , drop = FALSE], offset, 2L * node)
    recurse(local_x[seq.int(cp + 1L, n), , drop = FALSE], global, 2L * node + 1L)
    invisible(NULL)
  }
  recurse(x, 0L, 1L)
  sort(unique(selected))
}

#' Compare all seven multiple-change estimators used in the article
#'
#' @param x Numeric observation-by-coordinate matrix.
#' @param methods Any subset of ZAC, ZAM, ZAS, KDist, HDD, E-Divisive, and gSeg.
#' @param eta,alpha Boundary fraction and level.
#' @param permutations Local permutation count.
#' @param n_intervals Random WBS intervals (the full interval is additional).
#' @param min_interval,min_segment Minimum interval and segment lengths.
#' @param seed Integer seed.
#' @return A named list of estimated split-after locations, with package and
#'   calibration metadata in attributes.
#' @export
hdzacp_compare_multiple <- function(x, methods = .hdzacp_method_order,
                                     eta = 0.1, alpha = 0.05,
                                     permutations = 99L, n_intervals = 50L,
                                     min_interval = 40L, min_segment = 20L,
                                     seed = 20260824L) {
  x <- .hdzacp_check_matrix(x)
  methods <- unique(methods)
  aliases <- c("E-Div." = "E-Divisive", EDIVISIVE = "E-Divisive")
  methods[methods %in% names(aliases)] <- aliases[methods[methods %in% names(aliases)]]
  if (!length(methods) || any(!methods %in% .hdzacp_method_order)) {
    stop("Unknown comparison method.", call. = FALSE)
  }
  z <- hdzacp_wbs(
    x, eta = eta, alpha = alpha, n_intervals = n_intervals,
    min_interval = min_interval, permutations = permutations,
    stopping = "unadjusted", seed = .hdzacp_seed(seed, 1L)
  )$changepoints
  if ("KDist" %in% methods) .hdzacp_require("KDist", "GitHub repository `zhangxiany-tamu/KDist`")
  if ("HDD" %in% methods) .hdzacp_require("HDDchangepoint", "GitHub repository `rezadrikvandi/HDDchangepoint`")
  if ("E-Divisive" %in% methods) .hdzacp_require("ecp")
  if ("gSeg" %in% methods) .hdzacp_require("gSeg")
  if ("gSeg" %in% methods) .hdzacp_require("ade4")
  kdist <- if ("KDist" %in% methods) {
    .hdzacp_with_seed(.hdzacp_seed(seed, 2L), KDist::kcpd_wbs(
      data = x, type = "e-dist", M = as.integer(n_intervals),
      B = as.integer(permutations), alpha = alpha,
      seeds = .hdzacp_seed(seed, 2L), num_cores = 1L
    ))$locations
  } else integer()
  hdd <- if ("HDD" %in% methods) {
    .hdzacp_hdd_multiple(
      x, alpha, as.integer(permutations), as.integer(n_intervals),
      as.integer(min_segment), .hdzacp_seed(seed, 3L)
    )
  } else integer()
  ediv_cp <- integer()
  if ("E-Divisive" %in% methods) {
    ediv <- .hdzacp_with_seed(.hdzacp_seed(seed, 4L), ecp::e.divisive(
      X = x, sig.lvl = alpha, R = as.integer(permutations), k = NULL,
      min.size = as.integer(min_segment), alpha = 1
    ))
    ediv_cp <- setdiff(ediv$estimates, c(1L, nrow(x) + 1L)) - 1L
  }
  gseg <- if ("gSeg" %in% methods) {
    .hdzacp_gseg_multiple(
      x, alpha, as.integer(permutations), as.integer(min_segment),
      .hdzacp_seed(seed, 5L)
    )
  } else integer()
  answer <- list(
    ZAC = sort(unique(as.integer(z$ZAC))),
    ZAM = sort(unique(as.integer(z$ZAM))),
    ZAS = sort(unique(as.integer(z$ZAS))),
    KDist = sort(unique(as.integer(kdist))),
    HDD = sort(unique(as.integer(hdd))),
    `E-Divisive` = sort(unique(as.integer(ediv_cp))),
    gSeg = sort(unique(as.integer(gseg)))
  )
  answer <- answer[methods]
  attr(answer, "packages") <- c(
    ZAC = "HDZACP", ZAM = "HDZACP", ZAS = "HDZACP", KDist = "KDist",
    HDD = "HDDchangepoint (article-fast)", `E-Divisive` = "ecp",
    gSeg = "gSeg + ade4"
  )[methods]
  answer
}

#' Multiple-change evaluation metrics
#'
#' @param truth,estimate True and estimated split-after locations.
#' @param n Series length.
#' @return Number recovery, scaled Hausdorff distance, directional distances,
#'   and adjusted Rand index.
#' @export
hdzacp_segmentation_metrics <- function(truth, estimate, n) {
  sanitize <- function(cp) {
    cp <- sort(unique(as.integer(round(cp))))
    cp[is.finite(cp) & cp >= 1L & cp < n]
  }
  truth <- sanitize(truth)
  estimate <- sanitize(estimate)
  largest_true_segment <- max(diff(c(0L, truth, n)))
  if (!length(estimate)) {
    true_to_estimated <- n / largest_true_segment
    estimated_to_true <- 0
  } else if (!length(truth)) {
    true_to_estimated <- 0
    estimated_to_true <- n / largest_true_segment
  } else {
    distances <- abs(outer(truth, estimate, "-"))
    true_to_estimated <- max(apply(distances, 1L, min)) / largest_true_segment
    estimated_to_true <- max(apply(distances, 2L, min)) / largest_true_segment
  }
  labels <- function(cp) findInterval(seq_len(n), cp, left.open = TRUE)
  contingency <- table(labels(truth), labels(estimate))
  choose2 <- function(x) x * (x - 1) / 2
  total_pairs <- choose2(sum(contingency))
  index <- sum(choose2(contingency))
  row_pairs <- sum(choose2(rowSums(contingency)))
  col_pairs <- sum(choose2(colSums(contingency)))
  expected <- if (total_pairs) row_pairs * col_pairs / total_pairs else 0
  maximum <- 0.5 * (row_pairs + col_pairs)
  ari <- if (abs(maximum - expected) < .Machine$double.eps) {
    as.double(identical(labels(truth), labels(estimate)))
  } else (index - expected) / (maximum - expected)
  c(
    n_true = length(truth), n_estimated = length(estimate),
    count_error = length(estimate) - length(truth),
    exact_count = as.double(length(estimate) == length(truth)),
    scaled_hausdorff = max(true_to_estimated, estimated_to_true),
    true_to_estimated = true_to_estimated,
    estimated_to_true = estimated_to_true,
    ARI = ari
  )
}

