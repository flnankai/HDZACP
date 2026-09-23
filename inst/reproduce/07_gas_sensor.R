#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
value <- function(name, default = NULL) {
  z <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(z)) sub(paste0("^--", name, "="), "", z[[length(z)]]) else default
}
if (!requireNamespace("HDZACP", quietly = TRUE)) stop("Install HDZACP first.")
fun <- getExportedValue("HDZACP", "reproduce_hdzacp_gas_sensor")
a <- list(quick = quick, output_dir = value("output-dir", "results/gas_sensor"))
if (!is.null(value("data-zip"))) a$data_zip <- value("data-zip")
if (!is.null(value("cores"))) a$cores <- as.integer(value("cores"))
if (!is.null(value("seed"))) a$seed <- as.integer(value("seed"))
if (!"..." %in% names(formals(fun))) a <- a[names(a) %in% names(formals(fun))]
do.call(fun, a)
