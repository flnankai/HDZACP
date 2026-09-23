.hdzacp_reproduction_methods <- c(
  "ZAC", "ZAM", "ZAS", "KDist", "HDD", "E-Divisive", "gSeg"
)

.hdzacp_reproduction_colours <- c(
  ZAC = "#D55E00", ZAM = "#CC79A7", ZAS = "#0072B2",
  KDist = "#009E73", HDD = "#E69F00", `E-Divisive` = "#4D4D4D",
  gSeg = "#56B4E9"
)

.hdzacp_normalize_method <- function(x) {
  x <- as.character(x)
  x[x %in% c("KDist_e_dist", "KDist-e-dist")] <- "KDist"
  x[x %in% c("HDDchangepoint_dist1", "HDD-dist1")] <- "HDD"
  x[x %in% c("ecp_E_divisive", "E-Div.", "E-Divisive")] <- "E-Divisive"
  x[x %in% c("gSeg_MST_max", "gSeg-MST")] <- "gSeg"
  x
}

.hdzacp_reproduction_controls <- function(repetitions, permutations, cores,
                                           quick, analytic = FALSE) {
  repetitions <- .hdzacp_scalar_integer(repetitions, "repetitions", 1L)
  cores <- .hdzacp_scalar_integer(cores, "cores", 1L)
  if (!analytic) {
    permutations <- .hdzacp_scalar_integer(permutations, "permutations", 1L)
  }
  if (isTRUE(quick)) {
    repetitions <- min(repetitions, if (analytic) 20L else 10L)
    if (!analytic) permutations <- min(permutations, 19L)
    cores <- min(cores, 2L)
  }
  list(repetitions = repetitions, permutations = permutations, cores = cores)
}

.hdzacp_invoke <- function(fun, args) {
  if (is.character(fun)) {
    if (!exists(fun, mode = "function", inherits = TRUE)) {
      stop("Required HDZACP function `", fun, "` is not available.", call. = FALSE)
    }
    fun <- get(fun, mode = "function", inherits = TRUE)
  }
  fml <- names(formals(fun))
  if (!"..." %in% fml) args <- args[names(args) %in% fml]
  do.call(fun, args)
}

.hdzacp_result_frame <- function(result, methods = NULL) {
  if (is.list(result) && !is.data.frame(result) &&
      is.data.frame(result$results)) result <- result$results
  if (!is.data.frame(result)) {
    result <- tryCatch(as.data.frame(result), error = function(e) NULL)
  }
  if (is.null(result) || !nrow(result)) {
    stop("The method call returned no tabular results.", call. = FALSE)
  }
  names(result) <- sub("^p_value$", "p.value", names(result))
  names(result) <- sub("^pvalue$", "p.value", names(result))
  if (!"method" %in% names(result)) stop("Results have no `method` column.", call. = FALSE)
  result$method <- .hdzacp_normalize_method(result$method)
  if (!"p.value" %in% names(result)) result$p.value <- NA_real_
  if (!"statistic" %in% names(result)) result$statistic <- NA_real_
  if (!"estimate" %in% names(result)) {
    candidate <- intersect(
      c("split_after_index", "estimated_change", "raw_cp_index"), names(result)
    )
    result$estimate <- if (length(candidate)) result[[candidate[[1L]]]] else NA_real_
  }
  if (!"reject" %in% names(result)) result$reject <- NA
  result$status <- "ok"
  result$error <- ""
  keep <- c("method", "statistic", "p.value", "reject", "estimate", "status", "error")
  result <- result[, keep, drop = FALSE]
  if (!is.null(methods)) {
    missing <- setdiff(methods, result$method)
    if (length(missing)) {
      result <- rbind(result, data.frame(
        method = missing, statistic = NA_real_, p.value = NA_real_,
        reject = NA, estimate = NA_real_, status = "error",
        error = "Method did not return a result.", stringsAsFactors = FALSE
      ))
    }
    result <- result[match(methods, result$method), , drop = FALSE]
  }
  rownames(result) <- NULL
  result
}

.hdzacp_error_frame <- function(methods, error) {
  data.frame(
    method = methods, statistic = NA_real_, p.value = NA_real_, reject = NA,
    estimate = NA_real_, status = "error", error = conditionMessage(error),
    stringsAsFactors = FALSE
  )
}

.hdzacp_compare_once <- function(x, methods, eta, alpha, permutations, seed,
                                  force_location = FALSE) {
  tryCatch({
    result <- .hdzacp_invoke("hdzacp_compare", list(
      x = x, methods = methods, eta = eta, alpha = alpha,
      permutations = permutations, seed = seed, cores = 1L,
      force_location = force_location, return_location = force_location
    ))
    .hdzacp_result_frame(result, methods)
  }, error = function(e) .hdzacp_error_frame(methods, e))
}

.hdzacp_wilson <- function(successes, total, level = 0.95) {
  if (!total) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  phat <- successes / total
  denominator <- 1 + z^2 / total
  center <- (phat + z^2 / (2 * total)) / denominator
  half <- z / denominator * sqrt(
    phat * (1 - phat) / total + z^2 / (4 * total^2)
  )
  c(lower = max(0, center - half), upper = min(1, center + half))
}

.hdzacp_binary_summary <- function(raw, keys, outcome, value_name) {
  groups <- split(raw, interaction(raw[keys], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(group) {
    valid <- group$status == "ok" & !is.na(group[[outcome]])
    total <- sum(valid)
    successes <- if (total) sum(as.logical(group[[outcome]][valid])) else 0L
    ci <- .hdzacp_wilson(successes, total)
    answer <- group[1L, keys, drop = FALSE]
    answer$replications <- nrow(group)
    answer$valid_replications <- total
    answer$failures <- nrow(group) - total
    answer$successes <- successes
    answer[[value_name]] <- if (total) successes / total else NA_real_
    answer$mc_se <- if (total) sqrt(answer[[value_name]] *
      (1 - answer[[value_name]]) / total) else NA_real_
    answer$ci95_lower <- ci[["lower"]]
    answer$ci95_upper <- ci[["upper"]]
    answer
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer
}

.hdzacp_write_config <- function(config, path) {
  values <- vapply(config, function(x) paste(x, collapse = ","), character(1L))
  utils::write.csv(
    data.frame(setting = names(values), value = unname(values),
               stringsAsFactors = FALSE),
    path, row.names = FALSE
  )
}

.hdzacp_plot_method_curves <- function(summary, value, output, ylab,
                                        ylim = c(0, 1), nominal = NULL) {
  methods <- intersect(.hdzacp_reproduction_methods, unique(summary$method))
  s_values <- sort(unique(summary$s))
  grDevices::pdf(output, width = 8.2, height = 5.8, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5.1, 5.2, 1.2, 1.0), las = 1)
  graphics::plot(
    range(log2(s_values)), ylim, type = "n", xaxt = "n",
    xlab = "Number of changed independent components, s", ylab = ylab,
    bty = "l"
  )
  graphics::axis(1, at = log2(s_values), labels = s_values)
  graphics::abline(h = pretty(ylim), col = "#ECECEC", lwd = 0.7)
  if (!is.null(nominal)) graphics::abline(h = nominal, lty = 3, col = "#666666")
  for (i in seq_along(methods)) {
    d <- summary[summary$method == methods[[i]], , drop = FALSE]
    d <- d[match(s_values, d$s), , drop = FALSE]
    graphics::lines(
      log2(s_values), d[[value]], type = "o", pch = 14L + i,
      lwd = 1.6, cex = 0.85, col = .hdzacp_reproduction_colours[[methods[[i]]]]
    )
  }
  graphics::legend(
    "bottom", inset = c(0, -0.28), xpd = NA, horiz = TRUE, bty = "n",
    legend = methods, col = .hdzacp_reproduction_colours[methods],
    pch = 15L + seq_along(methods), lty = 1, lwd = 1.6, cex = 0.82
  )
  invisible(output)
}

.hdzacp_size_task <- function(index, tasks, settings) {
  task <- tasks[index, , drop = FALSE]
  x <- simulate_hdzacp(
    settings$n_values[[task$n_index]], settings$p_values[[task$p_index]],
    dependence = task$dependence, rho = settings$rho, alternative = "normal",
    seed = .hdzacp_seed(settings$seed, task$replication,
                        task$scenario_index, 0L)
  )
  out <- .hdzacp_compare_once(
    x, settings$methods, settings$eta, settings$alpha,
    settings$permutations,
    .hdzacp_seed(settings$seed, task$replication, task$scenario_index, 1L)
  )
  out$n <- nrow(x)
  out$p <- ncol(x)
  out$dependence <- task$dependence
  out$rho <- settings$rho
  out$replication <- task$replication
  out$reject <- is.finite(out$p.value) & out$p.value <= settings$alpha
  out
}

#' Reproduce the seven-method Gaussian-null size study
#'
#' This is the paper preset: AR(1) and compound-symmetry correlation `0.5`,
#' trimming `0.1`, nominal level `0.05`, 1,000 replications, and 199
#' whole-vector permutations in each of the twelve `(n,p,dependence)` cells.
#'
#' @param output_dir Directory in which CSV files are written.
#' @param n_values,p_values Sample-size and dimension grids.
#' @param dependence Dependence structures to evaluate.
#' @param rho Cross-coordinate correlation.
#' @param repetitions Monte Carlo replications per configuration.
#' @param permutations Whole-vector permutations per test.
#' @param cores Number of outer parallel workers.
#' @param eta Trimming fraction.
#' @param alpha Nominal level.
#' @param methods Methods passed to [hdzacp_compare()].
#' @param seed Reproducibility seed.
#' @param quick If `TRUE`, retain every design cell but reduce repetitions,
#'   permutations, and workers for a smoke test.
#' @param ... Reserved for forward compatibility; currently must be empty.
#'
#' @return A list containing `raw`, `summary`, `config`, and output paths.
#' @export
reproduce_hdzacp_size <- function(
    output_dir,
    ...,
    n_values = c(100L, 200L),
    p_values = c(100L, 200L, 400L),
    dependence = c("AR1", "CS"),
    rho = 0.5,
    repetitions = 1000L,
    permutations = 199L,
    cores = 25L,
    eta = 0.10,
    alpha = 0.05,
    methods = .hdzacp_reproduction_methods,
    seed = 20260822L,
    quick = FALSE) {
  if (length(list(...))) stop("Unused arguments supplied through `...`.", call. = FALSE)
  output_dir <- .hdzacp_dir(output_dir)
  n_values <- sort(unique(as.integer(n_values)))
  p_values <- sort(unique(as.integer(p_values)))
  dependence <- unique(match.arg(dependence, c("AR1", "CS", "independent"),
                                 several.ok = TRUE))
  controls <- .hdzacp_reproduction_controls(
    repetitions, permutations, cores, quick
  )
  settings <- list(
    n_values = n_values, p_values = p_values, dependence = dependence,
    rho = rho, repetitions = controls$repetitions,
    permutations = controls$permutations, cores = controls$cores,
    eta = eta, alpha = alpha, methods = methods, seed = seed, quick = quick
  )
  scenarios <- expand.grid(
    n_index = seq_along(n_values), p_index = seq_along(p_values),
    dependence = dependence, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scenarios$scenario_index <- seq_len(nrow(scenarios))
  tasks <- scenarios[rep(seq_len(nrow(scenarios)), each = controls$repetitions), ]
  tasks$replication <- rep(seq_len(controls$repetitions), nrow(scenarios))
  pieces <- .hdzacp_parallel_lapply(
    seq_len(nrow(tasks)), .hdzacp_size_task, tasks = tasks,
    settings = settings, cores = controls$cores
  )
  raw <- do.call(rbind, pieces)
  raw <- raw[, c("replication", "n", "p", "dependence", "rho", "method",
                 "statistic", "p.value", "reject", "estimate", "status", "error")]
  summary <- .hdzacp_binary_summary(
    raw, c("n", "p", "dependence", "rho", "method"), "reject",
    "empirical_size"
  )
  summary$alpha <- alpha
  summary$permutations <- controls$permutations
  summary <- summary[order(summary$dependence, summary$n, summary$p,
                           match(summary$method, methods)), , drop = FALSE]
  paths <- list(
    raw = file.path(output_dir, "hdzacp_size_replications.csv"),
    summary = file.path(output_dir, "hdzacp_size_summary.csv"),
    config = file.path(output_dir, "hdzacp_size_config.csv")
  )
  utils::write.csv(raw, paths$raw, row.names = FALSE)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  .hdzacp_write_config(settings, paths$config)
  invisible(list(raw = raw, summary = summary, config = settings, paths = paths))
}

.hdzacp_power_scenarios <- function() {
  data.frame(
    scenario = c("student_ar1", "student_cs", "exponential_ar1", "gamma_ar1"),
    dependence = c("AR1", "CS", "AR1", "AR1"),
    alternative = c("student", "student", "exponential", "gamma"),
    figure = c(
      "fig_power_ar_student_t.pdf", "fig_power_cs_student_t.pdf",
      "fig_power_ar_exp_skew.pdf", "fig_power_ar_gamma_skew.pdf"
    ),
    stringsAsFactors = FALSE
  )
}

.hdzacp_power_task <- function(index, tasks, settings) {
  task <- tasks[index, , drop = FALSE]
  scenario <- settings$scenarios[task$scenario_index, , drop = FALSE]
  s <- settings$s_values[[task$s_index]]
  df <- if (scenario$alternative == "student") {
    hdzacp_df(s, settings$p, scenario$dependence)
  } else NA_real_
  x <- simulate_hdzacp(
    settings$n, settings$p, change_points = floor(settings$n / 2), s = s,
    dependence = scenario$dependence, rho = settings$rho,
    alternative = scenario$alternative,
    df = if (is.finite(df)) df else 3,
    gamma_shape = settings$gamma_shape, pattern = "single",
    seed = .hdzacp_seed(settings$seed, task$replication,
                        100L * task$scenario_index + task$s_index, 0L)
  )
  out <- .hdzacp_compare_once(
    x, settings$methods, settings$eta, settings$alpha,
    settings$permutations,
    .hdzacp_seed(settings$seed, task$replication,
                 100L * task$scenario_index + task$s_index, 1L)
  )
  out$scenario <- scenario$scenario
  out$dependence <- scenario$dependence
  out$alternative <- scenario$alternative
  out$s <- s
  out$df <- df
  out$replication <- task$replication
  out$reject <- is.finite(out$p.value) & out$p.value <= settings$alpha
  out
}

#' Reproduce all four paper power experiments
#'
#' Runs the AR(1) and compound-symmetry Student-tail experiments and the AR(1)
#' centered-exponential and standardized-Gamma experiments.  The four returned
#' PDF files have the names used by the article source.
#'
#' @param output_dir Output directory.
#' @param scenarios Any subset of `"student_ar1"`, `"student_cs"`,
#'   `"exponential_ar1"`, and `"gamma_ar1"`.
#' @param n,p Sample size and dimension.
#' @param s_values Sparsity grid.
#' @param rho Correlation parameter.
#' @param gamma_shape Shape for the Gamma alternative.
#' @param repetitions,permutations,cores Monte Carlo controls.
#' @param eta,alpha Scan trimming and nominal level.
#' @param methods Methods passed to [hdzacp_compare()].
#' @param seed Reproducibility seed.
#' @param quick Run a downscaled smoke test while preserving every scenario.
#' @param ... Reserved for forward compatibility; currently must be empty.
#'
#' @return A list containing raw and summarized power and generated paths.
#' @export
reproduce_hdzacp_power <- function(
    output_dir,
    ...,
    scenarios = c("student_ar1", "student_cs", "exponential_ar1", "gamma_ar1"),
    n = 200L,
    p = 200L,
    s_values = c(1L, 2L, 5L, 10L, 20L, 30L, 40L, 50L, 60L, 80L, 120L, 200L),
    rho = 0.5,
    gamma_shape = 4,
    repetitions = 1000L,
    permutations = 199L,
    cores = 25L,
    eta = 0.10,
    alpha = 0.05,
    methods = .hdzacp_reproduction_methods,
    seed = 20260821L,
    quick = FALSE) {
  if (length(list(...))) stop("Unused arguments supplied through `...`.", call. = FALSE)
  output_dir <- .hdzacp_dir(output_dir)
  available <- .hdzacp_power_scenarios()
  if (any(!scenarios %in% available$scenario)) stop("Unknown power scenario.", call. = FALSE)
  scenarios <- available[match(unique(scenarios), available$scenario), , drop = FALSE]
  n <- .hdzacp_scalar_integer(n, "n", 8L)
  p <- .hdzacp_scalar_integer(p, "p", 2L)
  s_values <- sort(unique(as.integer(s_values)))
  if (any(s_values < 1L | s_values > p)) stop("`s_values` must lie in 1, ..., p.", call. = FALSE)
  controls <- .hdzacp_reproduction_controls(repetitions, permutations, cores, quick)
  settings <- list(
    scenarios = scenarios, n = n, p = p, s_values = s_values, rho = rho,
    gamma_shape = gamma_shape, repetitions = controls$repetitions,
    permutations = controls$permutations, cores = controls$cores,
    eta = eta, alpha = alpha, methods = methods, seed = seed, quick = quick
  )
  cells <- expand.grid(
    scenario_index = seq_len(nrow(scenarios)), s_index = seq_along(s_values),
    KEEP.OUT.ATTRS = FALSE
  )
  tasks <- cells[rep(seq_len(nrow(cells)), each = controls$repetitions), ]
  tasks$replication <- rep(seq_len(controls$repetitions), nrow(cells))
  pieces <- .hdzacp_parallel_lapply(
    seq_len(nrow(tasks)), .hdzacp_power_task, tasks = tasks,
    settings = settings, cores = controls$cores
  )
  raw <- do.call(rbind, pieces)
  raw <- raw[, c("replication", "scenario", "dependence", "alternative", "s", "df",
                 "method", "statistic", "p.value", "reject", "estimate",
                 "status", "error")]
  summary <- .hdzacp_binary_summary(
    raw, c("scenario", "dependence", "alternative", "s", "method"),
    "reject", "empirical_power"
  )
  df_map <- unique(raw[, c("scenario", "s", "df"), drop = FALSE])
  summary <- merge(summary, df_map, by = c("scenario", "s"), all.x = TRUE, sort = FALSE)
  summary$alpha <- alpha
  summary$permutations <- controls$permutations
  summary <- summary[order(match(summary$scenario, scenarios$scenario), summary$s,
                           match(summary$method, methods)), , drop = FALSE]
  paths <- list(
    raw = file.path(output_dir, "hdzacp_power_replications.csv"),
    summary = file.path(output_dir, "hdzacp_power_summary.csv"),
    config = file.path(output_dir, "hdzacp_power_config.csv")
  )
  utils::write.csv(raw, paths$raw, row.names = FALSE)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  .hdzacp_write_config(settings[names(settings) != "scenarios"], paths$config)
  figure_paths <- setNames(character(nrow(scenarios)), scenarios$scenario)
  for (i in seq_len(nrow(scenarios))) {
    figure_paths[[i]] <- file.path(output_dir, scenarios$figure[[i]])
    .hdzacp_plot_method_curves(
      summary[summary$scenario == scenarios$scenario[[i]], , drop = FALSE],
      "empirical_power", figure_paths[[i]], "Empirical power", c(0, 1), alpha
    )
  }
  paths$figures <- figure_paths
  invisible(list(raw = raw, summary = summary, config = settings, paths = paths))
}

.hdzacp_location_task <- function(index, tasks, settings) {
  task <- tasks[index, , drop = FALSE]
  dependence <- settings$dependence[[task$dependence_index]]
  s <- settings$s_values[[task$s_index]]
  df <- hdzacp_df(s, settings$p, dependence)
  tau <- floor(settings$n / 2)
  x <- simulate_hdzacp(
    settings$n, settings$p, tau, s, dependence, settings$rho,
    "student", df = df, pattern = "single",
    seed = .hdzacp_seed(settings$seed, task$replication,
                        100L * task$dependence_index + task$s_index, 0L)
  )
  out <- .hdzacp_compare_once(
    x, settings$methods, settings$eta, settings$alpha,
    settings$permutations,
    .hdzacp_seed(settings$seed, task$replication,
                 100L * task$dependence_index + task$s_index, 1L),
    force_location = TRUE
  )
  e_row <- match("E-Divisive", out$method)
  if (!is.na(e_row) && !is.finite(out$estimate[[e_row]]) &&
      requireNamespace("ecp", quietly = TRUE)) {
    forced <- tryCatch(.hdzacp_with_seed(
      .hdzacp_seed(settings$seed, task$replication,
        100L * task$dependence_index + task$s_index, 2L),
      ecp::e.divisive(
        X = x, sig.lvl = settings$alpha, R = settings$permutations, k = 1L,
        min.size = ceiling(settings$n * settings$eta), alpha = 1
      )
    ), error = function(e) NULL)
    if (!is.null(forced)) {
      interior <- setdiff(forced$estimates, c(1L, settings$n + 1L))
      if (length(interior)) out$estimate[[e_row]] <- as.integer(interior[[1L]])
    }
  }
  out$dependence <- dependence
  out$s <- s
  out$df <- df
  out$replication <- task$replication
  out$tau <- tau
  out$absolute_error <- abs(as.double(out$estimate) - tau)
  out$scaled_absolute_error <- out$absolute_error / settings$n
  out
}

.hdzacp_summarize_location <- function(raw) {
  keys <- c("dependence", "s", "df", "method")
  groups <- split(raw, interaction(raw[keys], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(group) {
    valid <- group$status == "ok" & is.finite(group$scaled_absolute_error)
    answer <- group[1L, keys, drop = FALSE]
    answer$replications <- nrow(group)
    answer$valid_replications <- sum(valid)
    answer$failures <- sum(!valid)
    answer$mean_abs_error <- if (any(valid)) mean(group$absolute_error[valid]) else NA_real_
    answer$median_abs_error <- if (any(valid)) stats::median(group$absolute_error[valid]) else NA_real_
    answer$mean_scaled_abs_error <- if (any(valid)) mean(group$scaled_absolute_error[valid]) else NA_real_
    answer$median_scaled_abs_error <- if (any(valid)) stats::median(group$scaled_absolute_error[valid]) else NA_real_
    answer
  })
  do.call(rbind, rows)
}

#' Reproduce the single-change localization experiment
#'
#' @inheritParams reproduce_hdzacp_power
#' @param dependence Dependence structures; the paper uses both AR(1) and CS.
#'
#' @return A list containing raw and summarized errors and the PDF path.
#' @export
reproduce_hdzacp_location <- function(
    output_dir,
    ...,
    n = 200L,
    p = 200L,
    s_values = c(1L, 2L, 5L, 10L, 20L, 30L, 40L, 50L, 60L, 80L, 120L, 200L),
    dependence = c("AR1", "CS"),
    rho = 0.5,
    repetitions = 1000L,
    permutations = 199L,
    cores = 25L,
    eta = 0.10,
    alpha = 0.05,
    methods = .hdzacp_reproduction_methods,
    seed = 20260821L,
    quick = FALSE) {
  if (length(list(...))) stop("Unused arguments supplied through `...`.", call. = FALSE)
  output_dir <- .hdzacp_dir(output_dir)
  n <- .hdzacp_scalar_integer(n, "n", 8L)
  p <- .hdzacp_scalar_integer(p, "p", 2L)
  s_values <- sort(unique(as.integer(s_values)))
  dependence <- unique(match.arg(dependence, c("AR1", "CS"), several.ok = TRUE))
  controls <- .hdzacp_reproduction_controls(repetitions, permutations, cores, quick)
  settings <- list(
    n = n, p = p, s_values = s_values, dependence = dependence, rho = rho,
    repetitions = controls$repetitions, permutations = controls$permutations,
    cores = controls$cores, eta = eta, alpha = alpha, methods = methods,
    seed = seed, quick = quick
  )
  cells <- expand.grid(
    dependence_index = seq_along(dependence), s_index = seq_along(s_values),
    KEEP.OUT.ATTRS = FALSE
  )
  tasks <- cells[rep(seq_len(nrow(cells)), each = controls$repetitions), ]
  tasks$replication <- rep(seq_len(controls$repetitions), nrow(cells))
  pieces <- .hdzacp_parallel_lapply(
    seq_len(nrow(tasks)), .hdzacp_location_task, tasks = tasks,
    settings = settings, cores = controls$cores
  )
  raw <- do.call(rbind, pieces)
  summary <- .hdzacp_summarize_location(raw)
  summary <- summary[order(match(summary$dependence, dependence), summary$s,
                           match(summary$method, methods)), , drop = FALSE]
  paths <- list(
    raw = file.path(output_dir, "hdzacp_location_replications.csv"),
    summary = file.path(output_dir, "hdzacp_location_summary.csv"),
    config = file.path(output_dir, "hdzacp_location_config.csv"),
    figure = file.path(output_dir, "fig_location_t_scaled_error.pdf")
  )
  utils::write.csv(raw, paths$raw, row.names = FALSE)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  .hdzacp_write_config(settings, paths$config)
  grDevices::pdf(paths$figure, width = 8.2, height = 8.6, useDingbats = FALSE)
  old <- graphics::par(mfrow = c(length(dependence), 1L), mar = c(4.4, 5.2, 2.0, 1.0))
  on.exit({graphics::par(old); grDevices::dev.off()}, add = TRUE)
  for (dep in dependence) {
    d <- summary[summary$dependence == dep, , drop = FALSE]
    ymax <- max(d$mean_scaled_abs_error, na.rm = TRUE)
    .hdzacp_plot_curves_current_device(
      d, "mean_scaled_abs_error", "Mean scaled absolute error",
      c(0, max(0.05, 1.05 * ymax)), paste0(dep, " dependence")
    )
  }
  invisible(list(raw = raw, summary = summary, config = settings, paths = paths))
}

.hdzacp_plot_curves_current_device <- function(summary, value, ylab, ylim, title = "") {
  methods <- intersect(.hdzacp_reproduction_methods, unique(summary$method))
  s_values <- sort(unique(summary$s))
  graphics::plot(range(log2(s_values)), ylim, type = "n", xaxt = "n",
                 xlab = "Number of changed independent components, s",
                 ylab = ylab, main = title, bty = "l")
  graphics::axis(1, at = log2(s_values), labels = s_values)
  graphics::abline(h = pretty(ylim), col = "#ECECEC", lwd = 0.7)
  for (i in seq_along(methods)) {
    d <- summary[summary$method == methods[[i]], , drop = FALSE]
    d <- d[match(s_values, d$s), , drop = FALSE]
    graphics::lines(log2(s_values), d[[value]], type = "o", pch = 14L + i,
                    lwd = 1.4, cex = 0.72,
                    col = .hdzacp_reproduction_colours[[methods[[i]]]])
  }
  graphics::legend("topright", methods,
                   col = .hdzacp_reproduction_colours[methods], lty = 1,
                   pch = 15L + seq_along(methods), bty = "n", cex = 0.72,
                   ncol = 2L)
}

.hdzacp_extract_changepoints <- function(result, methods) {
  changepoints <- if (is.list(result) && !is.null(result$changepoints)) {
    result$changepoints
  } else result
  if (!is.list(changepoints)) {
    if (length(methods) != 1L) stop("WBS result has no named change-point list.")
    changepoints <- setNames(list(changepoints), methods)
  }
  if (is.null(names(changepoints))) names(changepoints) <- methods[seq_along(changepoints)]
  names(changepoints) <- .hdzacp_normalize_method(names(changepoints))
  out <- setNames(vector("list", length(methods)), methods)
  for (method in methods) {
    value <- changepoints[[method]]
    if (is.null(value)) value <- integer()
    out[[method]] <- sort(unique(as.integer(value[is.finite(value)])))
  }
  out
}

.hdzacp_multiple_task <- function(index, tasks, settings) {
  task <- tasks[index, , drop = FALSE]
  dependence <- settings$dependence[[task$dependence_index]]
  s <- settings$s_values[[task$s_index]]
  df <- hdzacp_df(s, settings$p, dependence)
  truth <- floor(settings$n * c(0.3, 0.7))
  x <- simulate_hdzacp(
    settings$n, settings$p, truth, s, dependence, settings$rho,
    "student", df = df, pattern = "epidemic",
    seed = .hdzacp_seed(settings$seed, task$replication,
                        100L * task$dependence_index + task$s_index, 0L)
  )
  result <- tryCatch(
    .hdzacp_invoke("hdzacp_compare_multiple", list(
      x = x, eta = settings$eta, alpha = settings$alpha,
      methods = settings$methods,
      permutations = settings$permutations,
      n_intervals = settings$n_intervals,
      min_interval = settings$min_interval,
      min_segment = settings$min_segment,
      seed = .hdzacp_seed(settings$seed, task$replication,
        100L * task$dependence_index + task$s_index, 1L)
    )),
    error = function(e) e
  )
  if (inherits(result, "condition")) {
    rows <- lapply(settings$methods, function(method) data.frame(
      method = method, estimated_locations = "", estimated_count = NA_integer_,
      count_error = NA_integer_, exact_count = NA,
      scaled_hausdorff = NA_real_, directed_true_to_estimated = NA_real_,
      directed_estimated_to_true = NA_real_, adjusted_rand = NA_real_,
      status = "error", error = conditionMessage(result), stringsAsFactors = FALSE
    ))
  } else {
    estimates <- .hdzacp_extract_changepoints(result, settings$methods)
    rows <- lapply(settings$methods, function(method) {
      estimate <- estimates[[method]]
      metric <- .hdzacp_invoke("hdzacp_segmentation_metrics", list(
        truth = truth, estimate = estimate, n = settings$n
      ))
      metric <- unlist(metric)
      data.frame(
        method = method,
        estimated_locations = paste(estimate, collapse = ";"),
        estimated_count = length(estimate),
        count_error = length(estimate) - length(truth),
        exact_count = length(estimate) == length(truth),
        scaled_hausdorff = as.double(metric[["scaled_hausdorff"]]),
        directed_true_to_estimated = as.double(metric[["true_to_estimated"]]),
        directed_estimated_to_true = as.double(metric[["estimated_to_true"]]),
        adjusted_rand = as.double(metric[["ARI"]]),
        status = "ok", error = "", stringsAsFactors = FALSE
      )
    })
  }
  out <- do.call(rbind, rows)
  out$replication <- task$replication
  out$dependence <- dependence
  out$s <- s
  out$df <- df
  out$truth <- paste(truth, collapse = ";")
  out
}

.hdzacp_summarize_multiple <- function(raw) {
  keys <- c("dependence", "s", "df", "method")
  groups <- split(raw, interaction(raw[keys], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(group) {
    valid <- group$status == "ok" & is.finite(group$estimated_count)
    answer <- group[1L, keys, drop = FALSE]
    answer$replications <- nrow(group)
    answer$valid_replications <- sum(valid)
    answer$failures <- sum(!valid)
    for (name in c("exact_count", "count_error", "scaled_hausdorff", "adjusted_rand")) {
      answer[[paste0("mean_", name)]] <- if (any(valid)) {
        mean(as.double(group[[name]][valid]))
      } else NA_real_
    }
    answer$mean_abs_count_error <- if (any(valid)) {
      mean(abs(group$count_error[valid]))
    } else NA_real_
    answer
  })
  do.call(rbind, rows)
}

#' Reproduce the Gaussian--Student--Gaussian multiple-change experiment
#'
#' The paper preset uses AR(1), changes at `0.3 n` and `0.7 n`, 50 WBS
#' intervals plus the full interval, minimum interval length 40, minimum
#' segment length 20, 100 replications, and 99 local permutations.  CS may be
#' added through `dependence = c("AR1", "CS")`.
#'
#' @inheritParams reproduce_hdzacp_location
#' @param n_intervals Number of random WBS intervals.
#' @param min_interval Minimum WBS interval length.
#' @param min_segment Minimum terminal segment length.
#'
#' @return A list containing raw results, summaries, count frequencies, and
#'   one PDF for each requested dependence structure.
#' @export
reproduce_hdzacp_multiple <- function(
    output_dir,
    ...,
    n = 200L,
    p = 200L,
    s_values = c(1L, 2L, 5L, 10L, 20L, 30L, 40L, 50L, 60L, 80L, 120L, 200L),
    dependence = "AR1",
    rho = 0.5,
    repetitions = 100L,
    permutations = 99L,
    cores = 25L,
    eta = 0.10,
    alpha = 0.05,
    methods = .hdzacp_reproduction_methods,
    n_intervals = 50L,
    min_interval = 40L,
    min_segment = 20L,
    seed = 20260824L,
    quick = FALSE) {
  if (length(list(...))) stop("Unused arguments supplied through `...`.", call. = FALSE)
  output_dir <- .hdzacp_dir(output_dir)
  n <- .hdzacp_scalar_integer(n, "n", 8L)
  p <- .hdzacp_scalar_integer(p, "p", 2L)
  s_values <- sort(unique(as.integer(s_values)))
  dependence <- unique(match.arg(dependence, c("AR1", "CS"), several.ok = TRUE))
  controls <- .hdzacp_reproduction_controls(repetitions, permutations, cores, quick)
  if (isTRUE(quick)) n_intervals <- min(as.integer(n_intervals), 10L)
  settings <- list(
    n = n, p = p, s_values = s_values, dependence = dependence, rho = rho,
    repetitions = controls$repetitions, permutations = controls$permutations,
    cores = controls$cores, eta = eta, alpha = alpha, methods = methods,
    n_intervals = as.integer(n_intervals), min_interval = as.integer(min_interval),
    min_segment = as.integer(min_segment), seed = seed, quick = quick
  )
  cells <- expand.grid(
    dependence_index = seq_along(dependence), s_index = seq_along(s_values),
    KEEP.OUT.ATTRS = FALSE
  )
  tasks <- cells[rep(seq_len(nrow(cells)), each = controls$repetitions), ]
  tasks$replication <- rep(seq_len(controls$repetitions), nrow(cells))
  pieces <- .hdzacp_parallel_lapply(
    seq_len(nrow(tasks)), .hdzacp_multiple_task, tasks = tasks,
    settings = settings, cores = controls$cores
  )
  raw <- do.call(rbind, pieces)
  raw <- raw[, c("replication", "dependence", "s", "df", "truth", "method",
                 "estimated_locations", "estimated_count", "count_error",
                 "exact_count", "scaled_hausdorff", "adjusted_rand",
                 "directed_true_to_estimated", "directed_estimated_to_true",
                 "status", "error")]
  summary <- .hdzacp_summarize_multiple(raw)
  summary <- summary[order(match(summary$dependence, dependence), summary$s,
                           match(summary$method, methods)), , drop = FALSE]
  valid <- raw[raw$status == "ok" & is.finite(raw$estimated_count), , drop = FALSE]
  count_frequency <- as.data.frame(table(
    dependence = valid$dependence, s = valid$s, method = valid$method,
    estimated_count = valid$estimated_count
  ), stringsAsFactors = FALSE)
  count_frequency <- count_frequency[count_frequency$Freq > 0, , drop = FALSE]
  names(count_frequency)[names(count_frequency) == "Freq"] <- "frequency"
  paths <- list(
    raw = file.path(output_dir, "multiple_cp_raw.csv"),
    summary = file.path(output_dir, "multiple_cp_summary.csv"),
    count_frequency = file.path(output_dir, "multiple_cp_count_frequency.csv"),
    config = file.path(output_dir, "multiple_cp_config.csv")
  )
  utils::write.csv(raw, paths$raw, row.names = FALSE)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  utils::write.csv(count_frequency, paths$count_frequency, row.names = FALSE)
  .hdzacp_write_config(settings, paths$config)
  figure_paths <- setNames(character(length(dependence)), dependence)
  for (dep in dependence) {
    filename <- if (dep == "AR1") "fig_multiple_cp_normal_t_normal_ar.pdf" else
      "fig_multiple_cp_normal_t_normal_cs.pdf"
    figure_paths[[dep]] <- file.path(output_dir, filename)
    d <- summary[summary$dependence == dep, , drop = FALSE]
    grDevices::pdf(figure_paths[[dep]], width = 8.4, height = 10.2,
                   useDingbats = FALSE)
    old <- graphics::par(mfrow = c(3, 1), mar = c(4.2, 5.0, 1.2, 1.0))
    .hdzacp_plot_curves_current_device(d, "mean_exact_count",
      "Pr(estimated number = 2)", c(0, 1), "")
    ymax <- max(d$mean_scaled_hausdorff, na.rm = TRUE)
    .hdzacp_plot_curves_current_device(d, "mean_scaled_hausdorff",
      "Mean scaled Hausdorff loss", c(0, max(0.05, 1.05 * ymax)), "")
    .hdzacp_plot_curves_current_device(d, "mean_adjusted_rand",
      "Mean adjusted Rand index", c(0, 1), "")
    graphics::par(old)
    grDevices::dev.off()
  }
  paths$figures <- figure_paths
  invisible(list(raw = raw, summary = summary, count_frequency = count_frequency,
                 config = settings, paths = paths))
}

.hdzacp_simulate_null_design <- function(n, p, design, rho, block_size, seed) {
  if (design != "block") {
    return(simulate_hdzacp(
      n, p, dependence = if (design == "ar1") "AR1" else "independent",
      rho = rho, alternative = "normal", seed = seed
    ))
  }
  z <- .hdzacp_with_seed(seed, matrix(stats::rnorm(n * p), n, p))
  x <- z
  for (start in seq.int(1L, p, by = block_size)) {
    index <- start:min(p, start + block_size - 1L)
    q <- length(index)
    if (q > 1L) {
      sigma <- matrix(rho, q, q); diag(sigma) <- 1
      x[, index] <- z[, index, drop = FALSE] %*% chol(sigma)
    }
  }
  x
}

.hdzacp_analytic_task <- function(index, tasks, settings, ou_reference) {
  task <- tasks[index, , drop = FALSE]
  n <- settings$pairs$n[[task$pair_index]]
  p <- settings$pairs$p[[task$pair_index]]
  design <- settings$designs[[task$design_index]]
  x <- .hdzacp_simulate_null_design(
    n, p, design, settings$rho, settings$block_size,
    .hdzacp_seed(settings$seed, task$replication,
                 100L * task$design_index + task$pair_index, 0L)
  )
  result <- tryCatch(
    .hdzacp_invoke("hdzacp_test", list(
      x = x, methods = c("ZAC", "ZAM", "ZAS"),
      calibration = "asymptotic", eta = settings$eta, alpha = 0.05,
      seed = .hdzacp_seed(settings$seed, task$replication,
                          100L * task$design_index + task$pair_index, 1L),
      ou_reference = ou_reference,
      asymptotic_control = hdzacp_asymptotic_control(
        ou_paths = settings$ou_paths, ou_grid = settings$ou_grid,
        threads = 1L, ou_seed = settings$seed,
        max_control = list(cache_dir = settings$cache_dir, verbose = FALSE)
      )
    )),
    error = function(e) e
  )
  out <- if (inherits(result, "condition")) {
    .hdzacp_error_frame(c("ZAC", "ZAM", "ZAS"), result)
  } else .hdzacp_result_frame(result, c("ZAC", "ZAM", "ZAS"))
  out$n <- n; out$p <- p; out$design <- design
  out$replication <- task$replication
  out
}

.hdzacp_ou_draws <- function(reference) {
  candidates <- list(reference)
  if (is.list(reference)) candidates <- c(
    list(reference$fine, reference$draws, reference$reference, reference$maxima), candidates
  )
  for (value in candidates) {
    if (is.matrix(value)) {
      column <- if ("fine" %in% colnames(value)) "fine" else ncol(value)
      return(as.double(value[, column]))
    }
    if (is.numeric(value) && length(value) > 10L) return(as.double(value))
  }
  stop("Could not extract OU reference draws.", call. = FALSE)
}

.hdzacp_draw_analytic_size <- function(summary, output, pairs, designs) {
  methods <- c("ZAS", "ZAM", "ZAC")
  labels <- c(ZAS = "SUM", ZAM = "MAX", ZAC = "Cauchy")
  design_labels <- c(independent = "Independent", block = "Blocks", ar1 = "AR(1)")
  grDevices::pdf(output, width = 6.5, height = 6.4, useDingbats = FALSE)
  old <- graphics::par(mfrow = c(3, 3), mar = c(2.4, 3.2, 2.1, 0.7),
                       oma = c(1.8, 0.5, 0.2, 0.2), las = 1)
  on.exit({graphics::par(old); grDevices::dev.off()}, add = TRUE)
  panel <- 0L
  for (method in methods) for (design in designs) {
    panel <- panel + 1L
    d <- summary[summary$method == method & summary$design == design &
                   abs(summary$alpha - 0.05) < 1e-12, , drop = FALSE]
    d <- d[match(paste(pairs$n, pairs$p), paste(d$n, d$p)), , drop = FALSE]
    graphics::plot(NA, NA, xlim = c(0.55, nrow(pairs) + 0.45), ylim = c(0, 0.10),
                   axes = FALSE, xlab = "", ylab = "", bty = "n")
    graphics::rect(0.55, 0.03, nrow(pairs) + 0.45, 0.07,
                   col = "#EFEFEF", border = NA)
    graphics::abline(h = 0.05, lty = 2, col = "#777777")
    graphics::segments(seq_len(nrow(pairs)), d$ci95_lower,
                       seq_len(nrow(pairs)), d$ci95_upper)
    graphics::points(seq_len(nrow(pairs)), d$empirical_size, pch = 19, cex = 0.7)
    graphics::axis(1, at = seq_len(nrow(pairs)), labels = seq_len(nrow(pairs)))
    graphics::axis(2, at = c(0, 0.05, 0.10),
                   labels = if (design == designs[[1L]]) c("0", ".05", ".10") else FALSE)
    graphics::box()
    graphics::title(main = sprintf("(%s) %s | %s", letters[[panel]],
      labels[[method]], design_labels[[design]]), cex.main = 0.82)
  }
  graphics::mtext("Setting", side = 1, outer = TRUE, line = 0.4)
}

.hdzacp_draw_analytic_cdf <- function(raw, constants, output, pairs) {
  wanted <- unique(rbind(pairs[1L, , drop = FALSE],
                        pairs[min(2L, nrow(pairs)), , drop = FALSE],
                        pairs[max(1L, nrow(pairs) - 1L), , drop = FALSE],
                        pairs[nrow(pairs), , drop = FALSE]))
  grDevices::pdf(output, width = 6.5, height = 5.5, useDingbats = FALSE)
  old <- graphics::par(mfrow = c(2, 2), mar = c(3.4, 3.8, 2.1, 0.8),
                       oma = c(1.0, 0.3, 0.2, 0.2), las = 1)
  on.exit({graphics::par(old); grDevices::dev.off()}, add = TRUE)
  xx <- seq(-3, 7, length.out = 2001L)
  for (i in seq_len(4L)) {
    pair <- wanted[min(i, nrow(wanted)), , drop = FALSE]
    key <- paste0("n", pair$n, "_p", pair$p)
    cst <- constants[[key]]
    d <- raw[raw$design == "independent" & raw$method == "ZAM" &
               raw$n == pair$n & raw$p == pair$p & raw$status == "ok", ]
    z <- (d$statistic - cst$b) / cst$a
    graphics::plot(stats::ecdf(z), xlim = c(-3, 7), ylim = c(0, 1),
                   verticals = TRUE, do.points = FALSE, lwd = 1.3,
                   xlab = "", ylab = "Cumulative probability", bty = "l",
                   main = sprintf("(%s) n = %d, p = %d", letters[[i]], pair$n, pair$p))
    graphics::lines(xx, exp(-exp(-xx)), lty = 2, col = "#A34B13", lwd = 1.4)
    graphics::abline(h = 0.95, lty = 3, col = "#999999")
    graphics::legend("bottomright", c("Empirical", "Standard Gumbel"),
                     lty = c(1, 2), col = c("black", "#A34B13"), bty = "n", cex = 0.75)
  }
  graphics::mtext("Standardized MAX statistic", side = 1, outer = TRUE, line = 0.1)
}

#' Reproduce the direct asymptotic-calibration appendix
#'
#' The full preset evaluates six `(n,p)` pairs under independent coordinates,
#' independent blocks of size three with correlation `0.3`, and AR(1)
#' correlation `0.3`.  It uses 2,000 null samples per configuration and an
#' independent 200,000-path, 32,768-step OU reference distribution.  No
#' permutation or empirical critical-value fitting is used.
#'
#' @param output_dir Output directory.
#' @param pairs Two-column data frame or matrix of `(n,p)` pairs.
#' @param designs Any subset of `"independent"`, `"block"`, and `"ar1"`.
#' @param rho Correlation for block and AR(1) designs.
#' @param block_size Block size.
#' @param repetitions Monte Carlo replications per configuration.
#' @param alphas Nominal levels summarized from common p-values.
#' @param eta Trimming fraction.
#' @param ou_paths,ou_grid OU reference resolution.
#' @param cores Number of outer workers.
#' @param seed Reproducibility seed.
#' @param quick Downscale replications and OU resolution without dropping cells.
#' @param ... Reserved for forward compatibility; currently must be empty.
#'
#' @return A list containing raw p-values, size summaries, MAX quantiles,
#'   calibration objects, and output paths.
#' @export
reproduce_hdzacp_analytic_size <- function(
    output_dir,
    ...,
    pairs = data.frame(
      n = c(100L, 200L, 200L, 500L, 500L, 1000L),
      p = c(100L, 200L, 500L, 200L, 500L, 500L)
    ),
    designs = c("independent", "block", "ar1"),
    rho = 0.30,
    block_size = 3L,
    repetitions = 2000L,
    alphas = c(0.01, 0.05, 0.10),
    eta = 0.10,
    ou_paths = 200000L,
    ou_grid = 32768L,
    cores = 16L,
    seed = 20260909L,
    quick = FALSE) {
  if (length(list(...))) stop("Unused arguments supplied through `...`.", call. = FALSE)
  output_dir <- .hdzacp_dir(output_dir)
  pairs <- as.data.frame(pairs)
  if (!all(c("n", "p") %in% names(pairs))) stop("`pairs` needs columns `n` and `p`.")
  pairs$n <- as.integer(pairs$n); pairs$p <- as.integer(pairs$p)
  designs <- unique(match.arg(designs, c("independent", "block", "ar1"),
                              several.ok = TRUE))
  controls <- .hdzacp_reproduction_controls(repetitions, NA_integer_, cores,
                                             quick, analytic = TRUE)
  if (isTRUE(quick)) {
    ou_paths <- min(as.integer(ou_paths), 2000L)
    ou_grid <- min(as.integer(ou_grid), 1024L)
  }
  cache_dir <- .hdzacp_dir(file.path(output_dir, "calibration_cache"))
  settings <- list(
    pairs = pairs, designs = designs, rho = rho, block_size = block_size,
    repetitions = controls$repetitions, alphas = alphas, eta = eta,
    ou_paths = ou_paths, ou_grid = ou_grid, cores = controls$cores,
    seed = seed, quick = quick, cache_dir = cache_dir
  )
  ou_reference <- .hdzacp_invoke("hdzacp_ou_reference", list(
    eta = eta, paths = ou_paths, grid = ou_grid,
    seed = .hdzacp_seed(seed, 0L, 0L, 101L), threads = controls$cores,
    cache = TRUE, cache_dir = cache_dir
  ))
  saveRDS(ou_reference, file.path(output_dir, "ou_reference.rds"))
  constants <- setNames(vector("list", nrow(pairs)),
                        paste0("n", pairs$n, "_p", pairs$p))
  for (i in seq_len(nrow(pairs))) {
    constants[[i]] <- .hdzacp_invoke("hdzacp_max_constants", list(
      n = pairs$n[[i]], p = pairs$p[[i]], eta = eta,
      cache_dir = cache_dir, t_nodes = 17L, endpoint_M = 1024L,
      moment_M = 4096L, ode_tol = 1e-10, verbose = FALSE
    ))
    saveRDS(constants[[i]], file.path(
      output_dir, sprintf("max_constants_n%d_p%d.rds", pairs$n[[i]], pairs$p[[i]])
    ))
  }
  cells <- expand.grid(
    pair_index = seq_len(nrow(pairs)), design_index = seq_along(designs),
    KEEP.OUT.ATTRS = FALSE
  )
  tasks <- cells[rep(seq_len(nrow(cells)), each = controls$repetitions), ]
  tasks$replication <- rep(seq_len(controls$repetitions), nrow(cells))
  pieces <- .hdzacp_parallel_lapply(
    seq_len(nrow(tasks)), .hdzacp_analytic_task, tasks = tasks,
    settings = settings, ou_reference = ou_reference, cores = controls$cores
  )
  raw <- do.call(rbind, pieces)
  raw <- raw[, c("replication", "design", "n", "p", "method", "statistic",
                 "p.value", "status", "error")]
  expanded <- do.call(rbind, lapply(alphas, function(alpha) {
    d <- raw
    d$alpha <- alpha
    d$reject <- is.finite(d$p.value) & d$p.value <= alpha
    d
  }))
  summary <- .hdzacp_binary_summary(
    expanded, c("design", "n", "p", "method", "alpha"),
    "reject", "empirical_size"
  )
  summary <- summary[order(match(summary$design, designs), summary$n, summary$p,
                           match(summary$method, c("ZAS", "ZAM", "ZAC")),
                           summary$alpha), , drop = FALSE]
  quantiles <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    n <- pairs$n[[i]]; p <- pairs$p[[i]]
    d <- raw[raw$design == "independent" & raw$n == n & raw$p == p &
               raw$method == "ZAM" & raw$status == "ok", ]
    cst <- constants[[paste0("n", n, "_p", p)]]
    critical <- cst$b + cst$a * (-log(-log(0.95)))
    empirical <- stats::quantile(d$statistic, 0.95, names = FALSE, na.rm = TRUE)
    data.frame(n = n, p = p, a = cst$a, b = cst$b,
               analytic_critical_095 = critical,
               empirical_quantile_095 = empirical,
               difference = critical - empirical)
  }))
  paths <- list(
    raw = file.path(output_dir, "analytic_size_replications.csv"),
    summary = file.path(output_dir, "analytic_size_summary.csv"),
    quantiles = file.path(output_dir, "analytic_max_quantiles.csv"),
    config = file.path(output_dir, "analytic_size_config.rds"),
    size_figure = file.path(output_dir, "fig_analytic_size.pdf"),
    cdf_figure = file.path(output_dir, "fig_analytic_max_cdf.pdf")
  )
  utils::write.csv(raw, paths$raw, row.names = FALSE)
  utils::write.csv(summary, paths$summary, row.names = FALSE)
  utils::write.csv(quantiles, paths$quantiles, row.names = FALSE)
  saveRDS(settings, paths$config)
  .hdzacp_draw_analytic_size(summary, paths$size_figure, pairs, designs)
  .hdzacp_draw_analytic_cdf(raw, constants, paths$cdf_figure, pairs)
  invisible(list(raw = raw, summary = summary, max_quantiles = quantiles,
                 ou_reference = ou_reference, max_constants = constants,
                 config = settings, paths = paths))
}

#' Reproduce the HDZACP article
#'
#' Runs all simulation sections and, when requested, both empirical analyses.
#' Empirical observations are downloaded into `data_dir`; they are never stored
#' in the installed package or source repository.
#'
#' @param output_dir Parent output directory.
#' @param sections Any subset of `"size"`, `"power"`, `"location"`,
#'   `"multiple"`, `"analytic-size"`, `"braincloud"`, and `"gas-sensor"`.
#' @param quick Use the smoke-test preset in every selected section.
#' @param cores Optional common worker count. If `NULL`, each simulation uses
#'   its paper default (25 for the main simulations and 16 for analytic size).
#' @param data_dir Directory for externally downloaded empirical files.
#' @param download Download missing empirical data from the official source.
#' @param seed Seed passed to both empirical analyses.
#'
#' @return A named list of reproduction results.
#' @export
reproduce_hdzacp_article <- function(
    output_dir,
    sections = c("size", "power", "location", "multiple", "analytic-size",
                 "braincloud", "gas-sensor"),
    quick = FALSE,
    cores = NULL,
    data_dir = NULL,
    download = TRUE,
    seed = 20260821L) {
  output_dir <- .hdzacp_dir(output_dir)
  if (is.null(data_dir)) data_dir <- file.path(output_dir, "external-data")
  allowed <- c("size", "power", "location", "multiple", "analytic-size",
               "braincloud", "gas-sensor")
  sections <- match.arg(sections, allowed, several.ok = TRUE)
  answer <- list()
  common <- function(default) if (is.null(cores)) default else cores
  if ("size" %in% sections) answer$size <- reproduce_hdzacp_size(
    file.path(output_dir, "size"), cores = common(25L), quick = quick
  )
  if ("power" %in% sections) answer$power <- reproduce_hdzacp_power(
    file.path(output_dir, "power"), cores = common(25L), quick = quick
  )
  if ("location" %in% sections) answer$location <- reproduce_hdzacp_location(
    file.path(output_dir, "location"), cores = common(25L), quick = quick
  )
  if ("multiple" %in% sections) answer$multiple <- reproduce_hdzacp_multiple(
    file.path(output_dir, "multiple"), cores = common(25L), quick = quick
  )
  if ("analytic-size" %in% sections) answer$analytic_size <-
    reproduce_hdzacp_analytic_size(
      file.path(output_dir, "analytic-size"), cores = common(16L), quick = quick
    )
  if ("braincloud" %in% sections) answer$braincloud <-
    reproduce_hdzacp_braincloud(
      data_dir = file.path(data_dir, "braincloud"),
      output_dir = file.path(output_dir, "braincloud"), download = download,
      seed = seed, quick = quick
    )
  if ("gas-sensor" %in% sections) answer$gas_sensor <-
    reproduce_hdzacp_gas_sensor(
      zip_file = file.path(data_dir, "gas-sensor", "gas_sensor_array_drift.zip"),
      raw_dir = file.path(data_dir, "gas-sensor", "raw"),
      output_dir = file.path(output_dir, "gas-sensor"), download = download,
      seed = seed, quick = quick
    )
  invisible(answer)
}

