#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
option_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

if (!requireNamespace("HDZACP", quietly = TRUE)) {
  stop("Install HDZACP before running this script.", call. = FALSE)
}

output_dir <- option_value("output-dir", file.path("results", "hdzacp_article"))
cores <- option_value("cores")
seed <- option_value("seed")

call_args <- list(quick = quick, output_dir = output_dir)
if (!is.null(cores)) call_args$cores <- as.integer(cores)
if (!is.null(seed)) call_args$seed <- as.integer(seed)

runner <- getExportedValue("HDZACP", "reproduce_hdzacp_article")
formal_names <- names(formals(runner))
if (!"..." %in% formal_names) {
  call_args <- call_args[names(call_args) %in% formal_names]
}

message(if (quick) "Running HDZACP smoke reproduction." else
  "Running the full HDZACP article reproduction.")
do.call(runner, call_args)
