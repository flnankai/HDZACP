#' Simulate the covariance-preserving alternatives used by HDZACP
#'
#' Generate an independent-component sequence and mix its coordinates to give
#' independent, AR(1), or compound-symmetry dependence.  Changed innovations
#' are standardized to have mean zero and variance one, so the alternatives
#' change marginal shape while preserving the population mean and covariance.
#'
#' @param n Number of observations.
#' @param p Number of coordinates.
#' @param change_points Integer split indices.  `pattern = "single"` accepts
#'   zero or one split; `pattern = "epidemic"` requires two splits.
#' @param s Number of leading independent components whose marginal law changes.
#' @param dependence Cross-coordinate dependence structure.
#' @param rho Correlation parameter.  AR(1) requires `abs(rho) < 1`; compound
#'   symmetry requires `-1 / (p - 1) < rho < 1` when `p > 1`.
#' @param alternative Changed innovation law.  `"normal"` gives the null.
#' @param df Degrees of freedom for the standardized Student alternative; must
#'   exceed two.
#' @param gamma_shape Shape parameter for the standardized Gamma alternative.
#' @param pattern A permanent single change or an epidemic change that returns
#'   to the standard-normal law after the second split.
#' @param seed Optional random-number seed.
#'
#' @return An `n` by `p` numeric matrix.  Design information is stored in the
#'   `"hdzacp_design"` attribute.
#' @export
simulate_hdzacp <- function(
    n,
    p,
    change_points = numeric(),
    s = 0,
    dependence = c("AR1", "CS", "independent"),
    rho = 0.5,
    alternative = c("student", "exponential", "gamma", "normal"),
    df = 3,
    gamma_shape = 4,
    pattern = c("single", "epidemic"),
    seed = NULL) {
  dependence <- match.arg(dependence)
  alternative <- match.arg(alternative)
  pattern <- match.arg(pattern)

  n <- .hdzacp_scalar_integer(n, "n", lower = 2L)
  p <- .hdzacp_scalar_integer(p, "p", lower = 1L)
  s <- .hdzacp_scalar_integer(s, "s", lower = 0L, upper = p)
  if (length(change_points)) {
    if (any(!is.finite(change_points)) ||
        any(abs(change_points - round(change_points)) > 0)) {
      stop("`change_points` must contain finite integer split indices.",
           call. = FALSE)
    }
    change_points <- sort(unique(as.integer(round(change_points))))
    if (any(change_points < 1L | change_points >= n)) {
      stop("Every change point must lie in 1, ..., n - 1.", call. = FALSE)
    }
  } else {
    change_points <- integer()
  }
  if (pattern == "single" && length(change_points) > 1L) {
    stop("`pattern = \"single\"` accepts at most one change point.",
         call. = FALSE)
  }
  if (pattern == "epidemic" && length(change_points) != 2L) {
    stop("`pattern = \"epidemic\"` requires exactly two change points.",
         call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho)) {
    stop("`rho` must be one finite number.", call. = FALSE)
  }
  if (dependence == "AR1" && abs(rho) >= 1) {
    stop("AR(1) dependence requires `abs(rho) < 1`.", call. = FALSE)
  }
  if (dependence == "CS" && p > 1L &&
      (rho <= -1 / (p - 1) || rho >= 1)) {
    stop("Compound symmetry requires -1/(p-1) < rho < 1.",
         call. = FALSE)
  }
  if (alternative == "student" &&
      (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 2)) {
    stop("The standardized Student alternative requires `df > 2`.",
         call. = FALSE)
  }
  if (alternative == "gamma" &&
      (!is.numeric(gamma_shape) || length(gamma_shape) != 1L ||
       !is.finite(gamma_shape) || gamma_shape <= 0)) {
    stop("`gamma_shape` must be positive.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .hdzacp_scalar_integer(seed, "seed", lower = 0L)
  }

  changed_rows <- integer()
  if (length(change_points) && s > 0L && alternative != "normal") {
    changed_rows <- if (pattern == "single") {
      seq.int(change_points[[1L]] + 1L, n)
    } else {
      seq.int(change_points[[1L]] + 1L, change_points[[2L]])
    }
  }

  generate <- function() {
    innovations <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    if (length(changed_rows)) {
      count <- length(changed_rows) * s
      changed <- switch(
        alternative,
        student = sqrt((df - 2) / df) * stats::rt(count, df = df),
        exponential = stats::rexp(count) - 1,
        gamma = (stats::rgamma(count, shape = gamma_shape) - gamma_shape) /
          sqrt(gamma_shape),
        normal = stats::rnorm(count)
      )
      innovations[changed_rows, seq_len(s)] <- matrix(
        changed, nrow = length(changed_rows), ncol = s
      )
    }

    if (dependence == "independent" || p == 1L) {
      x <- innovations
    } else if (dependence == "AR1") {
      x <- innovations
      innovation_scale <- sqrt(1 - rho^2)
      for (j in 2:p) {
        x[, j] <- rho * x[, j - 1L] + innovation_scale * innovations[, j]
      }
    } else {
      sigma <- matrix(rho, nrow = p, ncol = p)
      diag(sigma) <- 1
      lower <- t(chol(sigma))
      x <- innovations %*% t(lower)
    }
    storage.mode(x) <- "double"
    x
  }

  x <- if (is.null(seed)) generate() else .hdzacp_with_seed(seed, generate())
  attr(x, "hdzacp_design") <- list(
    n = n,
    p = p,
    change_points = change_points,
    s = s,
    dependence = dependence,
    rho = rho,
    alternative = alternative,
    df = if (alternative == "student") df else NA_real_,
    gamma_shape = if (alternative == "gamma") gamma_shape else NA_real_,
    pattern = pattern,
    seed = seed
  )
  x
}

#' Degrees-of-freedom schedule used in the article simulations
#'
#' The AR(1) schedule weakens a Student-tail change as its sparsity index grows.
#' The compound-symmetry schedule agrees with the AR(1) rule through
#' `transition_start` and blends into the dense-end inverse-log schedule by
#' `transition_end`.
#'
#' @param s Number(s) of changed independent components.
#' @param p Dimension of the observation vector.
#' @param dependence Dependence design determining the schedule.
#' @param df_at_s1 Degrees of freedom at the sparse endpoint.
#' @param df_slope Sparse-side slope multiplying `log2(s)`.
#' @param df_dense_slope Dense-side slope multiplying `log2(p / s)`.
#' @param transition_start,transition_end Compound-symmetry transition limits.
#' @param df_cap Upper cap for the degrees of freedom.
#' @param ... Reserved for forward compatibility; currently must be empty.
#'
#' @return A numeric vector with the same length as `s`.
#' @export
hdzacp_df <- function(
    s,
    p = 200,
    dependence = c("AR1", "CS"),
    df_at_s1 = 2.0001,
    df_slope = 0.10,
    df_dense_slope = 0.05,
    transition_start = 20,
    transition_end = 50,
    df_cap = 30,
    ...) {
  dependence <- match.arg(dependence)
  dots <- list(...)
  if (length(dots)) {
    stop("Unused argument(s): ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) ||
      p < 1 || abs(p - round(p)) > 0) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }
  p <- as.integer(round(p))
  if (!is.numeric(s) || !length(s) || any(!is.finite(s)) ||
      any(s < 1 | s > p) || any(abs(s - round(s)) > 0)) {
    stop("`s` must contain integers between 1 and `p`.", call. = FALSE)
  }
  s <- as.double(s)
  base_inputs <- c(df_at_s1 = df_at_s1, df_slope = df_slope, df_cap = df_cap)
  if (any(!is.finite(base_inputs)) || df_at_s1 <= 2 || df_slope < 0 ||
      df_cap < df_at_s1) {
    stop("Invalid degrees-of-freedom schedule parameters.", call. = FALSE)
  }
  sparse_df <- df_at_s1 + df_slope * log2(s)
  if (dependence == "AR1") return(pmin(df_cap, sparse_df))

  cs_inputs <- c(df_dense_slope = df_dense_slope,
    transition_start = transition_start, transition_end = transition_end)
  if (any(!is.finite(cs_inputs)) || df_dense_slope < 0 || transition_start < 1 ||
      transition_start >= transition_end) {
    stop("Invalid compound-symmetry schedule parameters.", call. = FALSE)
  }

  dense_df <- df_at_s1 + df_dense_slope * log2(pmax(1, p / s))
  weight <- ifelse(
    s <= transition_start,
    0,
    ifelse(
      s >= transition_end,
      1,
      log2(s / transition_start) / log2(transition_end / transition_start)
    )
  )
  pmin(df_cap, (1 - weight) * sparse_df + weight * dense_df)
}

.hdzacp_scalar_integer <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      abs(x - round(x)) > 0 || x < lower || x > upper) {
    stop(sprintf("`%s` must be an integer in [%s, %s].", name, lower, upper),
         call. = FALSE)
  }
  as.integer(round(x))
}



