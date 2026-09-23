#' Deterministic asymptotic constants for ZAM
#'
#' Computes the current manuscript's deterministic subcritical normalization
#' for the ZAM maximum scan. The calculation uses the Hardy ODE, Bellman
#' terminal weights, the exact Bernoulli endpoint recursion implemented in C++,
#' and deterministic quadrature. It never fits a Gumbel law to simulated null
#' statistics.
#'
#' The present numerical implementation deliberately stops for critical or
#' boundary regimes instead of substituting a surrogate formula. Permutation
#' calibration is consequently the package default and is also preferable when
#' the asymptotic weak-dependence conditions are doubtful.
#'
#' @param n,p Sample size and dimension.
#' @param eta Boundary fraction.
#' @param cache_dir Directory for deterministic intermediate calculations.
#' @param t_nodes Number of integration nodes on the half interval.
#' @param endpoint_M Endpoint-recursion truncation.
#' @param moment_M Moment-series truncation.
#' @param ode_tol ODE tolerance.
#' @param verbose Emit numerical progress messages.
#' @return A list containing scale `a`, center `b`, branch information, numerical
#'   diagnostics, and the constants grid.
#' @export
hdzacp_max_constants <- function(n, p, eta = 0.1,
                                  cache_dir = file.path(
                                    tools::R_user_dir("HDZACP", "cache"),
                                    "max-calibration"
                                  ),
                                  t_nodes = 17L, endpoint_M = 1024L,
                                  moment_M = 4096L, ode_tol = 1e-10,
                                  verbose = interactive()) {
  .hdzacp_require("deSolve")
  latest_max_constants(
    n = as.integer(n), p = as.integer(p), eta = eta,
    cache_dir = cache_dir, t_nodes = as.integer(t_nodes),
    endpoint_M = as.integer(endpoint_M), moment_M = as.integer(moment_M),
    ode_tol = ode_tol, verbose = verbose
  )
}

#' Inspect the asymptotic ZAM branch
#'
#' @inheritParams hdzacp_max_constants
#' @return A one-row data frame identifying the subcritical, critical, or
#'   boundary regime.
#' @export
hdzacp_max_branch <- function(n, p, eta = 0.1) {
  .hdzacp_require("deSolve")
  latest_max_branch(as.integer(n), as.integer(p), eta)
}
