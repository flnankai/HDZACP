#' HDZACP: high-dimensional distributional change points
#'
#' HDZACP implements the ZAS (sum), ZAM (maximum), and ZAC (adaptive Cauchy)
#' rank-scan procedures developed in *Nonparametric Change-Point Detection and
#' Inference for High-Dimensional Distributions*. The package provides
#' permutation-calibrated single- and multiple-change procedures, optional
#' asymptotic calibration, comparison-method wrappers, and complete
#' reproduction entry points.
#'
#' @keywords internal
#' @useDynLib HDZACP, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

utils::globalVariables(c(
  "method", "p.value", "reject", "estimate", "statistic", "s",
  "power", "dependence", "metric", "value"
))
