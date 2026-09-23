#' Download the BrainCloud GEO Series Matrix
#'
#' @param path Destination `.txt.gz` path outside the installed package.
#' @param overwrite Replace an existing file.
#' @return The normalized destination path, invisibly.
#' @export
download_braincloud <- function(path, overwrite = FALSE) {
  url <- hdzacp_data_urls()$download_url[1L]
  if (file.exists(path) && file.info(path)$size > 0 && !overwrite) {
    return(invisible(normalizePath(path, winslash = "/")))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::download.file(url, path, mode = "wb", quiet = FALSE)
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("BrainCloud download did not create a valid file.", call. = FALSE)
  }
  invisible(normalizePath(path, winslash = "/"))
}

.hdzacp_strip_geo_quotes <- function(x) gsub('^"|"$', "", x)

.hdzacp_geo_header <- function(x) {
  .hdzacp_strip_geo_quotes(strsplit(x, "\t", fixed = TRUE)[[1L]])
}

.hdzacp_characteristic <- function(characteristics, key) {
  prefix <- paste0(tolower(key), ":")
  hit <- characteristics[startsWith(tolower(characteristics), prefix)]
  if (!length(hit)) return(NA_character_)
  trimws(sub("^[^:]+:", "", hit[[1L]]))
}

.hdzacp_read_braincloud <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  header <- character()
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) stop("GEO table-begin marker was not found.", call. = FALSE)
    if (identical(line, "!series_matrix_table_begin")) break
    header <- c(header, line)
  }
  sample_rows <- lapply(header[grepl("^!Sample_", header)], .hdzacp_geo_header)
  accession_index <- which(vapply(sample_rows, function(z) {
    identical(z[[1L]], "!Sample_geo_accession")
  }, logical(1L)))
  if (!length(accession_index)) stop("Sample accessions were not found.", call. = FALSE)
  accessions <- sample_rows[[accession_index[[1L]]]][-1L]
  characteristic_rows <- sample_rows[vapply(sample_rows, function(z) {
    identical(z[[1L]], "!Sample_characteristics_ch1")
  }, logical(1L))]
  characteristics <- lapply(seq_along(accessions), function(i) {
    values <- vapply(characteristic_rows, function(z) z[i + 1L], character(1L))
    values[nzchar(values)]
  })
  numeric_value <- function(key) suppressWarnings(as.numeric(vapply(
    characteristics, .hdzacp_characteristic, character(1L), key = key
  )))
  character_value <- function(key) vapply(
    characteristics, .hdzacp_characteristic, character(1L), key = key
  )
  metadata <- data.frame(
    accession = accessions, age = numeric_value("age"),
    sv1 = numeric_value("surrogate variable 1 (sv1)"),
    sv2 = numeric_value("surrogate variable 2 (sv2)"),
    sex = character_value("sex"), race = character_value("race"),
    stringsAsFactors = FALSE
  )
  expression_table <- utils::read.delim(
    connection, header = TRUE, sep = "\t", quote = "\"", comment.char = "",
    check.names = FALSE, fill = TRUE, na.strings = c("NA", "NULL", "")
  )
  expression_table <- expression_table[
    !is.na(expression_table[[1L]]) &
      expression_table[[1L]] != "!series_matrix_table_end", , drop = FALSE
  ]
  probe_id <- as.character(expression_table[[1L]])
  expression <- as.matrix(expression_table[-1L])
  storage.mode(expression) <- "double"
  if (!all(accessions %in% colnames(expression))) {
    stop("The GEO expression table is missing sample columns.", call. = FALSE)
  }
  expression <- expression[, accessions, drop = FALSE]
  rownames(expression) <- make.unique(probe_id)
  list(metadata = metadata, expression = expression)
}

.hdzacp_robust_scale <- function(x) {
  centers <- apply(x, 2L, stats::median, na.rm = TRUE)
  scales <- apply(x, 2L, stats::mad, constant = 1.4826, na.rm = TRUE)
  fallback <- !is.finite(scales) | scales <= sqrt(.Machine$double.eps)
  if (any(fallback)) {
    scales[fallback] <- apply(x[, fallback, drop = FALSE], 2L, stats::sd, na.rm = TRUE)
  }
  scales[!is.finite(scales) | scales <= sqrt(.Machine$double.eps)] <- 1
  sweep(sweep(x, 2L, centers, "-"), 2L, scales, "/")
}

#' Prepare the BrainCloud data used in the article
#'
#' Adults aged at least 20 are ordered by age and accession. For each probe,
#' the intercept, SV1, and SV2 are removed (age is not regressed out); finite
#' nonconstant probes are ranked by residual variance, the top features are
#' retained, and columns are median/MAD standardized.
#'
#' @param path Downloaded GEO Series Matrix path.
#' @param age_min Minimum age.
#' @param max_features Maximum number of residual-variance-ranked probes.
#' @param adjust_sv Regress out SV1 and SV2.
#' @param robust_standardize Apply median/MAD standardization.
#' @return A list with `data`, ordered `metadata`, and `selected_features`.
#' @export
prepare_braincloud <- function(path, age_min = 20, max_features = 2000L,
                                adjust_sv = TRUE, robust_standardize = TRUE) {
  raw <- .hdzacp_read_braincloud(path)
  metadata <- raw$metadata
  eligible <- is.finite(metadata$age) & metadata$age >= age_min
  if (adjust_sv) eligible <- eligible & is.finite(metadata$sv1) & is.finite(metadata$sv2)
  metadata <- metadata[eligible, , drop = FALSE]
  metadata <- metadata[order(metadata$age, metadata$accession), , drop = FALSE]
  expression <- t(raw$expression[, metadata$accession, drop = FALSE])
  expression <- expression[, colSums(!is.finite(expression)) == 0L, drop = FALSE]
  if (adjust_sv) {
    design <- cbind(intercept = 1, sv1 = metadata$sv1, sv2 = metadata$sv2)
    expression <- qr.resid(qr(design), expression)
  }
  feature_variance <- apply(expression, 2L, stats::var)
  usable <- is.finite(feature_variance) & feature_variance > 0
  expression <- expression[, usable, drop = FALSE]
  feature_variance <- feature_variance[usable]
  selected <- order(feature_variance, decreasing = TRUE)[seq_len(min(
    as.integer(max_features), ncol(expression)
  ))]
  expression <- expression[, selected, drop = FALSE]
  selected_variance <- feature_variance[selected]
  if (robust_standardize) expression <- .hdzacp_robust_scale(expression)
  if (any(!is.finite(expression))) stop("Non-finite values remain after preprocessing.", call. = FALSE)
  rownames(expression) <- metadata$accession
  list(
    data = expression, metadata = metadata,
    selected_features = data.frame(
      rank = seq_along(selected), probe_id = colnames(expression),
      residual_variance_before_scaling = as.double(selected_variance),
      stringsAsFactors = FALSE
    )
  )
}

#' Reproduce the BrainCloud seven-method analysis
#'
#' @param data_file Existing Series Matrix path. If `NULL`, a path below
#'   `data_dir` is used.
#' @param data_dir Download directory; no data are written into the package.
#' @param output_dir Result directory.
#' @param download Download a missing file from NCBI GEO.
#' @param max_features,eta,alpha,permutations,seed Article controls.
#' @param quick Use a small feature/permutation smoke configuration.
#' @return A list containing the seven-method table and prepared-data audit.
#' @export
reproduce_hdzacp_braincloud <- function(
    data_file = NULL, data_dir = file.path("data", "braincloud"),
    output_dir = file.path("output", "braincloud"), download = TRUE,
    max_features = 2000L, eta = 0.1, alpha = 0.05,
    permutations = 999L, seed = 20260821L, quick = FALSE) {
  output_dir <- .hdzacp_dir(output_dir)
  if (is.null(data_file)) data_file <- file.path(data_dir, "GSE30272_series_matrix.txt.gz")
  if (!file.exists(data_file)) {
    if (!download) stop("BrainCloud file is missing and `download = FALSE`.", call. = FALSE)
    download_braincloud(data_file)
  }
  if (quick) {
    max_features <- min(as.integer(max_features), 100L)
    permutations <- min(as.integer(permutations), 19L)
  }
  prepared <- prepare_braincloud(data_file, max_features = max_features)
  results <- hdzacp_compare(
    prepared$data, eta = eta, alpha = alpha, permutations = permutations, seed = seed
  )
  utils::write.csv(results, file.path(output_dir, "braincloud_seven_methods.csv"), row.names = FALSE)
  utils::write.csv(prepared$selected_features,
    file.path(output_dir, "braincloud_selected_features.csv"), row.names = FALSE)
  utils::write.csv(prepared$metadata,
    file.path(output_dir, "braincloud_ordered_samples.csv"), row.names = FALSE)
  audit <- data.frame(
    source_url = hdzacp_data_urls()$download_url[1L], n = nrow(prepared$data),
    p = ncol(prepared$data), age_min = 20, max_features = max_features,
    eta = eta, permutations = permutations, seed = seed
  )
  utils::write.csv(audit, file.path(output_dir, "braincloud_audit.csv"), row.names = FALSE)
  saveRDS(list(results = results, audit = audit),
    file.path(output_dir, "braincloud_results.rds"))
  list(results = results, audit = audit, metadata = prepared$metadata,
    selected_features = prepared$selected_features, data_file = normalizePath(data_file, winslash = "/"))
}

#' Download and extract the UCI Gas Sensor Array Drift data
#'
#' @param zip_file Destination archive path.
#' @param raw_dir Extraction directory.
#' @param overwrite Redownload the archive.
#' @return Paths to `batch1.dat`, ..., `batch10.dat`.
#' @export
download_gas_sensor <- function(zip_file, raw_dir, overwrite = FALSE) {
  expected <- file.path(raw_dir, paste0("batch", 1:10, ".dat"))
  if (all(file.exists(expected)) && !overwrite) return(invisible(expected))
  dir.create(dirname(zip_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(zip_file) || file.info(zip_file)$size <= 0 || overwrite) {
    utils::download.file(hdzacp_data_urls()$download_url[2L], zip_file,
      mode = "wb", quiet = FALSE)
  }
  utils::unzip(zip_file, exdir = raw_dir)
  nested <- list.files(raw_dir, pattern = "^batch[0-9]+\\.dat$",
    recursive = TRUE, full.names = TRUE)
  for (path in nested) {
    target <- file.path(raw_dir, basename(path))
    if (!identical(normalizePath(path, winslash = "/", mustWork = FALSE),
                   normalizePath(target, winslash = "/", mustWork = FALSE))) {
      file.copy(path, target, overwrite = TRUE)
    }
  }
  if (!all(file.exists(expected))) stop("Could not locate batch1.dat through batch10.dat.", call. = FALSE)
  invisible(expected)
}

.hdzacp_read_gas_batch <- function(path, batch) {
  tokens <- scan(path, what = character(), quiet = TRUE)
  if (length(tokens) %% 129L != 0L) stop("Unexpected row width in ", path, call. = FALSE)
  raw <- matrix(tokens, ncol = 129L, byrow = TRUE)
  header <- raw[, 1L]
  gas <- suppressWarnings(as.integer(sub(";.*$", "", header)))
  concentration <- suppressWarnings(as.double(sub("^[^;]*;", "", header)))
  feature_tokens <- raw[, -1L, drop = FALSE]
  feature_indices <- as.integer(sub(":.*$", "", feature_tokens[1L, ]))
  if (!identical(feature_indices, seq_len(128L))) {
    stop("Features are not the dense sequence 1,...,128.", call. = FALSE)
  }
  x <- matrix(as.double(sub("^[0-9]+:", "", feature_tokens)),
    nrow = nrow(feature_tokens), ncol = ncol(feature_tokens))
  colnames(x) <- sprintf("sensor_feature_%03d", seq_len(ncol(x)))
  list(
    metadata = data.frame(batch = as.integer(batch), source_row = seq_len(nrow(x)),
      gas_class = gas, concentration = concentration),
    x = x
  )
}

.hdzacp_adjust_concentration <- function(x, concentration, df, adjust_scale) {
  if (any(!is.finite(x)) || any(!is.finite(concentration)) || any(concentration <= 0)) {
    stop("Concentration adjustment requires finite data and positive concentrations.", call. = FALSE)
  }
  z <- .hdzacp_robust_scale(x)
  design <- cbind(intercept = 1, splines::ns(log(concentration), df = df))
  design_qr <- qr(design)
  residual <- qr.resid(design_qr, z)
  if (adjust_scale) {
    floor_value <- pmax(apply(residual^2, 2L, stats::median) * 0.01, 1e-8)
    log_squared <- log(sweep(residual^2, 2L, floor_value, "+"))
    fitted_log_variance <- design %*% qr.coef(design_qr, log_squared)
    fitted_log_variance <- pmax(pmin(fitted_log_variance, 12), -12)
    residual <- residual / sqrt(exp(fitted_log_variance))
  }
  .hdzacp_robust_scale(residual)
}

#' Prepare the Gas Sensor Array Drift analysis matrix
#'
#' @param raw_dir Directory containing the ten UCI batch files.
#' @param gas_class Target class; the article uses class 1 (ethanol).
#' @param balance_per_batch Maximum observations retained from each batch.
#' @param concentration_df Natural-spline degrees of freedom.
#' @param adjust_conditional_scale Adjust both conditional location and scale.
#' @param seed Sampling and within-batch ordering seed.
#' @return A list with adjusted data, ordered metadata, and a batch audit.
#' @export
prepare_gas_sensor <- function(raw_dir, gas_class = 1L, balance_per_batch = 50L,
                               concentration_df = 4L,
                               adjust_conditional_scale = TRUE,
                               seed = 20260821L) {
  paths <- file.path(raw_dir, paste0("batch", 1:10, ".dat"))
  if (!all(file.exists(paths))) stop("The ten UCI batch files are required.", call. = FALSE)
  pieces <- lapply(seq_along(paths), function(i) .hdzacp_read_gas_batch(paths[[i]], i))
  metadata_all <- do.call(rbind, lapply(pieces, `[[`, "metadata"))
  x_all <- do.call(rbind, lapply(pieces, `[[`, "x"))
  keep <- metadata_all$gas_class == gas_class
  metadata <- metadata_all[keep, , drop = FALSE]
  x <- x_all[keep, , drop = FALSE]
  finite <- colSums(!is.finite(x)) == 0L
  variances <- apply(x[, finite, drop = FALSE], 2L, stats::var)
  x <- x[, which(finite)[is.finite(variances) & variances > 0], drop = FALSE]
  adjusted <- .hdzacp_adjust_concentration(
    x, metadata$concentration, as.integer(concentration_df), adjust_conditional_scale
  )
  batches <- sort(unique(metadata$batch))
  audit <- do.call(rbind, lapply(batches, function(batch) {
    ii <- which(metadata$batch == batch)
    data.frame(
      batch = batch, n_all_gases = sum(metadata_all$batch == batch),
      n_target_gas = length(ii), n_unique_concentrations = length(unique(metadata$concentration[ii])),
      concentration_min = min(metadata$concentration[ii]),
      concentration_max = max(metadata$concentration[ii])
    )
  }))
  audit$n_selected <- pmin(audit$n_target_gas, as.integer(balance_per_batch))
  selected <- integer()
  for (batch in batches) {
    candidates <- which(metadata$batch == batch)
    selected <- c(selected, .hdzacp_with_seed(seed + 1000L + batch,
      sample(candidates, min(balance_per_batch, length(candidates)), replace = FALSE)))
  }
  metadata <- metadata[selected, , drop = FALSE]
  adjusted <- adjusted[selected, , drop = FALSE]
  ordering <- .hdzacp_with_seed(seed + 2000L, unlist(
    lapply(split(seq_len(nrow(metadata)), metadata$batch), sample), use.names = FALSE
  ))
  metadata <- metadata[ordering, , drop = FALSE]
  adjusted <- adjusted[ordering, , drop = FALSE]
  metadata$analysis_index <- seq_len(nrow(metadata))
  metadata$sample_id <- sprintf("batch%02d_row%04d", metadata$batch, metadata$source_row)
  rownames(adjusted) <- metadata$sample_id
  if (any(!is.finite(adjusted))) stop("Non-finite values remain after preprocessing.", call. = FALSE)
  list(data = adjusted, metadata = metadata, audit = audit)
}

.hdzacp_gas_figure <- function(changepoints, metadata, path) {
  methods <- names(changepoints)
  boundaries <- cumsum(as.integer(table(factor(metadata$batch, levels = sort(unique(metadata$batch))))))
  boundaries <- head(boundaries, -1L)
  grDevices::pdf(path, width = 11, height = 7, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(6.5, 7, 1, 1), xpd = NA, cex.axis = 1.15, cex.lab = 1.25)
  graphics::plot(c(1, nrow(metadata)), c(0.5, length(methods) + 0.5), type = "n",
    yaxt = "n", xlab = "Analysis index", ylab = "", bty = "l")
  graphics::axis(2, at = seq_along(methods), labels = methods, las = 1, cex.axis = 1.1)
  graphics::mtext("Method", side = 2, line = 5, cex = 1.25)
  for (b in boundaries) graphics::abline(v = b, col = "grey70", lty = 3)
  colors <- grDevices::hcl.colors(length(methods), "Dark 3")
  for (i in seq_along(methods)) {
    graphics::segments(1, i, nrow(metadata), i, col = "grey85")
    if (length(changepoints[[i]])) graphics::points(
      changepoints[[i]], rep(i, length(changepoints[[i]])), pch = 19,
      col = colors[i], cex = 1.25
    )
  }
  graphics::legend("bottom", inset = c(0, -0.23), horiz = TRUE, bty = "n",
    legend = c("estimated change point", "batch boundary"), pch = c(19, NA),
    lty = c(NA, 3), col = c("#3B6FB6", "grey60"), cex = 1.05)
}

#' Reproduce the Gas Sensor global and multiple-change analyses
#'
#' The analysis is descriptive rather than causal: concentration support and
#' batch are intertwined in the source data even after the documented robust
#' conditional location/scale adjustment.
#'
#' @param zip_file,raw_dir Archive and extraction paths outside the package.
#' @param output_dir Result directory.
#' @param download Download missing data from UCI.
#' @param balance_per_batch Maximum target-gas observations per batch.
#' @param eta,alpha,permutations,permutations_multiple,n_intervals,min_interval,min_segment,seed
#'   Article controls.
#' @param quick Reduce permutation/interval counts for a smoke run.
#' @return Global results, multiple-change estimates, metadata, and audit.
#' @export
reproduce_hdzacp_gas_sensor <- function(
    zip_file = file.path("data", "gas_sensor_drift", "gas_sensor_array_drift.zip"),
    raw_dir = file.path("data", "gas_sensor_drift", "raw"),
    output_dir = file.path("output", "gas_sensor_drift"), download = TRUE,
    balance_per_batch = 50L, eta = 0.1, alpha = 0.05,
    permutations = 999L, permutations_multiple = 199L,
    n_intervals = 100L, min_interval = 28L, min_segment = 14L,
    seed = 20260821L, quick = FALSE) {
  output_dir <- .hdzacp_dir(output_dir)
  expected <- file.path(raw_dir, paste0("batch", 1:10, ".dat"))
  if (!all(file.exists(expected))) {
    if (!download) stop("Gas Sensor batch files are missing and `download = FALSE`.", call. = FALSE)
    download_gas_sensor(zip_file, raw_dir)
  }
  if (quick) {
    permutations <- min(as.integer(permutations), 19L)
    permutations_multiple <- min(as.integer(permutations_multiple), 19L)
    n_intervals <- min(as.integer(n_intervals), 10L)
  }
  prepared <- prepare_gas_sensor(
    raw_dir, gas_class = 1L, balance_per_batch = balance_per_batch,
    concentration_df = 4L, adjust_conditional_scale = TRUE, seed = seed
  )
  global <- hdzacp_compare(
    prepared$data, eta = eta, alpha = alpha, permutations = permutations, seed = seed
  )
  multiple <- hdzacp_compare_multiple(
    prepared$data, eta = eta, alpha = alpha, permutations = permutations_multiple,
    n_intervals = n_intervals, min_interval = min_interval,
    min_segment = min_segment, seed = seed
  )
  multi_table <- do.call(rbind, lapply(names(multiple), function(method) {
    cp <- multiple[[method]]
    if (!length(cp)) return(data.frame(method = method, change_point = NA_integer_))
    data.frame(method = method, change_point = cp)
  }))
  utils::write.csv(global, file.path(output_dir, "gas_sensor_global_results.csv"), row.names = FALSE)
  utils::write.csv(multi_table, file.path(output_dir, "gas_sensor_multiple_results.csv"), row.names = FALSE)
  utils::write.csv(prepared$audit, file.path(output_dir, "gas_sensor_batch_audit.csv"), row.names = FALSE)
  utils::write.csv(prepared$metadata, file.path(output_dir, "gas_sensor_analysis_samples.csv"), row.names = FALSE)
  figure <- file.path(output_dir, "fig_gas_sensor_multiple_cp.pdf")
  .hdzacp_gas_figure(multiple, prepared$metadata, figure)
  audit <- list(
    source_url = hdzacp_data_urls()$download_url[2L], n = nrow(prepared$data),
    p = ncol(prepared$data), batch_sizes = table(prepared$metadata$batch),
    settings = list(eta = eta, alpha = alpha, permutations = permutations,
      permutations_multiple = permutations_multiple, n_intervals = n_intervals,
      min_interval = min_interval, min_segment = min_segment, seed = seed)
  )
  saveRDS(list(global = global, multiple = multiple, audit = audit),
    file.path(output_dir, "gas_sensor_results.rds"))
  list(global = global, multiple = multiple, metadata = prepared$metadata,
    batch_audit = prepared$audit, audit = audit, figure = figure)
}
