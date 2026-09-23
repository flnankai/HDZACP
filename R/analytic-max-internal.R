# Numerical calibration utilities for the corrected 2026-09-09 manuscript.
# Deterministic calculations only; no fitted simulation quantiles.
# Implements the subcritical normalization of body.tex, with the exact
# Bernoulli endpoint recursion and Bellman terminal weights. Critical and
# boundary inputs fail explicitly rather than substituting another law.
# Requires Rcpp and deSolve. Example:
# source("hdcp_latest_max_calibration.R"); latest_max_constants(200, 200)
# The cache stores deterministic constants only. The reported diagnostics
# compare truncation and integration resolutions; they are empirical
# numerical sensitivity checks, not rigorous uniform error bounds.

latest_max_hardy_b <- function(t = 0.1, q_max = 24, tol = 1e-11) {
  if (!requireNamespace("deSolve", quietly = TRUE)) stop("deSolve is required")
  stopifnot(length(t) == 1L, t > 0, t < 1)
  shoot <- function(theta, tt) {
    aa <- (1 - sqrt(1 - 4 * theta)) / 2
    q0 <- 1e-5
    x0 <- plogis(qlogis(tt) + q0)
    vv <- tt * (1 - tt)
    b2 <- aa * (1 - aa) * (1 - 2 * tt) * (aa - .5) /
      (vv * (2 - 3 * aa))
    u0 <- tt + aa * (x0 - tt) + b2 * (x0 - tt)^2
    w0 <- qlogis(u0) - qlogis(tt)
    okay <- TRUE
    rhs <- function(q, w, par) {
      x <- plogis(qlogis(tt) + q)
      # Stable evaluation of x-u using survival probabilities near x=1.
      gap <- plogis(-qlogis(tt) - w[1]) - plogis(-qlogis(tt) - q)
      if (!is.finite(gap) || gap <= 0) stop("Hardy branch crossed u=x")
      list(theta * q * x * (1 - x) / gap)
    }
    ans <- tryCatch(suppressWarnings(deSolve::ode(w0, c(q0, q_max), rhs,
      parms = NULL, method = "ode45", rtol = tol, atol = tol, maxsteps = 200000)),
      error = function(e) NULL)
    !is.null(ans) && nrow(ans) == 2L && ans[2, 2] < q_max
  }
  one <- function(tt) {
    if (shoot(.25 - 1e-9, tt)) return(.25)
    lo <- tt * (1 - tt) / 2
    hi <- .25 - 1e-9
    for (ii in seq_len(38L)) {
      mid <- (lo + hi) / 2
      if (shoot(mid, tt)) lo <- mid else hi <- mid
    }
    (lo + hi) / 2
  }
  bp <- one(min(t, 1 - t))
  list(b = bp, q_max = q_max, ode_tolerance = tol,
       method = "central-branch deterministic Hardy ODE shooting")
}

latest_max_branch <- function(n, p, eta = .1, hardy = latest_max_hardy_b(eta)) {
  b <- hardy$b
  N <- log(n)
  gamma <- log(p) / N
  qc <- if (b == .25) Inf else 1 / sqrt(1 - 4 * b)
  gammac <- if (is.infinite(qc)) Inf else (qc - 1)^2 / (4 * qc)
  epsN <- N^(-.25)
  branch <- ifelse(gamma < gammac - epsN, "subcritical",
    ifelse(gamma > gammac + epsN, "boundary", "critical"))
  ac <- (1 - sqrt(1 - 4 * b)) / 2
  individual <- .5 - ac / 2 - b * c(-.5, .5) * log(1 / eta)
  data.frame(n = n, p = p, eta = eta, b_eta = b, q_c = qc,
    gamma_c = gammac, gamma_n = gamma, epsilon_N = epsN,
    branch = branch, alpha_eta = sum(pmax(individual, 0)),
    beta_eta = sum(.5 * (individual > 0) + 1.5 * (individual == 0)))
}


latest_max_kl <- function(x, t) {
  z <- numeric(length(x))
  nz <- x > 0; no <- x < 1
  z[nz] <- z[nz] + x[nz] * log(x[nz] / t)
  z[no] <- z[no] + (1 - x[no]) * log((1 - x[no]) / (1 - t))
  z
}

latest_max_bellman <- function(t, theta, ode_tol = 1e-10, q_max = 28) {
  stopifnot(theta > 0, theta < .25, t > 0, t < 1)
  aa <- (1 - sqrt(1 - 4 * theta)) / 2
  kap <- aa / 2
  vv <- t * (1 - t)
  b2 <- aa * (1 - aa) * (1 - 2 * t) * (aa - .5) /
    (vv * (2 - 3 * aa))
  dlogA0 <- b2 / (1 - aa)
  rhs <- function(q, state, parms) {
    x <- plogis(qlogis(t) + q)
    u <- plogis(qlogis(t) + state[1])
    gap <- if (q > 0) plogis(-qlogis(t) - state[1]) -
      plogis(-qlogis(t) - q) else x - u
    if (!is.finite(gap) || gap * q <= 0) stop("Bellman branch crossed u=x")
    h <- latest_max_kl(x, t)
    g2 <- theta * q / gap
    common <- kap + theta * gap * q - g2 / 2 * (u * (1 - u) + gap^2)
    dA <- -(common + theta * c(-.5, .5) * h) / gap
    list(c(theta * q * x * (1 - x) / gap, dA * x * (1 - x)))
  }
  solve_side <- function(sg) {
    q0 <- sg * 1e-4
    x0 <- plogis(qlogis(t) + q0)
    dx <- x0 - t
    u0 <- t + aa * dx + b2 * dx^2
    init <- c(qlogis(u0) - qlogis(t), rep(dlogA0 * dx, 2))
    qs <- c(1e-4, seq(.01, q_max, by = .01))
    side_rhs <- function(z, state, parms) list(sg * rhs(sg * z, state, parms)[[1]])
    ans <- deSolve::ode(init, qs, side_rhs, parms = NULL, method = "ode45",
      rtol = ode_tol, atol = ode_tol, maxsteps = 200000)
    x <- plogis(qlogis(t) + sg * ans[, 1])
    u <- plogis(qlogis(t) + ans[, 2])
    g <- (1 - theta) * latest_max_kl(x, t) -
      (x * log(x / u) + (1 - x) * log((1 - x) / (1 - u)))
    cbind(x, g, logAm = ans[, 3], logAp = ans[, 4])
  }
  vals <- rbind(solve_side(-1), c(t, 0, 0, 0), solve_side(1))
  vals <- vals[order(vals[, 1]), , drop = FALSE]
  vals <- rbind(c(0, vals[1, -1]), vals, c(1, vals[nrow(vals), -1]))
  list(g = splinefun(vals[, 1], vals[, 2], method = "natural"),
    logAm = splinefun(vals[, 1], vals[, 3], method = "natural"),
    logAp = splinefun(vals[, 1], vals[, 4], method = "natural"),
    central_a = aa, kappa = kap, theta = theta, t = t,
    ode_tol = ode_tol, q_max = q_max)
}

latest_max_compile <- function() invisible(TRUE)


latest_max_logsumexp <- function(x) { mm <- max(x); mm + log(sum(exp(x-mm))) }

latest_max_psi <- function(t, theta, checkpoints = c(128L,256L,512L,1024L),
                           ode_tol = 1e-10, q_max = 28) {
  latest_max_compile()
  bb <- latest_max_bellman(t, theta, ode_tol, q_max)
  dp <- hdzacp_max_endpoint_dp_cpp(t, theta, as.integer(checkpoints))
  calc <- function(dd) {
    x <- (0:dd$m) / dd$m
    wt <- dd$m * bb$g(x) - bb$kappa * log(dd$m)
    km <- latest_max_logsumexp(dd$minus + wt + bb$logAm(x))
    kp <- latest_max_logsumexp(dd$plus + wt + bb$logAp(x))
    c(m = dd$m, logKm = km, logKp = kp,
      logPsi = .25 * log1p(-4 * theta) + km + kp)
  }
  tab <- as.data.frame(do.call(rbind, lapply(dp, calc)))
  list(psi = exp(tail(tab$logPsi, 1)), logPsi = tail(tab$logPsi, 1),
    convergence = tab, bellman = bb)
}


# Independent deterministic check of m(t), v0(t), 2026-09-09.
# The covariance series is evaluated as two Markov additive variances.
latest_max_moment_series <- function(t, checkpoints = c(256L,512L,1024L,2048L,4096L)) {
  stopifnot(t>0,t<1, all(checkpoints>=1))
  M <- max(checkpoints)
  prob <- 1
  first <- second <- list(0,0)
  shifts <- c(-.5,.5)
  mm <- -digamma(1)+2*log(2)-1
  hh <- h2 <- 0
  result <- list()
  for (ell in seq_len(M)) {
    probnew <- (1-t)*c(prob,0)+t*c(0,prob)
    x <- (0:ell)/ell
    kl <- numeric(ell+1L)
    inside <- x>0 & x<1
    kl[inside] <- x[inside]*(log(x[inside])-log(t)) +
      (1-x[inside])*(log1p(-x[inside])-log1p(-t))
    kl[1] <- -log1p(-t)
    kl[ell+1L] <- -log(t)
    gm <- sum(probnew * ell * kl)
    mm <- mm + (1/(ell-.5)+1/(ell+.5))*(gm-.5)
    for (jj in 1:2) {
      ff <- (1-t)*c(first[[jj]],0)+t*c(0,first[[jj]])
      ss <- (1-t)*c(second[[jj]],0)+t*c(0,second[[jj]])
      increment <- ell/(ell+shifts[jj])*kl
      first[[jj]] <- ff+increment*probnew
      second[[jj]] <- ss+2*increment*ff+increment^2*probnew
    }
    prob <- probnew
    hh <- hh+1/ell
    h2 <- h2+1/ell^2
    if (ell %in% checkpoints) {
      means <- vapply(first,sum,numeric(1))
      variances <- vapply(second,sum,numeric(1))-means^2
      vv <- -2*digamma(1)+1-pi^2/6+sum(variances)-2*hh+h2
      aa <- (1-t*(1-t))/(12*t*(1-t))
      result[[length(result)+1L]] <- data.frame(t=t,M=ell,m_raw=mm,
        m_tail1=mm+2*aa/(ell+.5),v0_raw=vv,
        probability=sum(prob),var_minus=variances[1],var_plus=variances[2])
    }
  }
  ans <- do.call(rbind,result)
  ans$v0_extrapolated <- NA_real_
  if (nrow(ans)>=3L) for (i in 3:nrow(ans)) {
    if (ans$M[i]==2*ans$M[i-1] && ans$M[i-1]==2*ans$M[i-2]) {
      ans$v0_extrapolated[i] <- ans$v0_raw[i-2]-4*ans$v0_raw[i-1]+4*ans$v0_raw[i]
    }
  }
  ans
}



latest_max_simpson <- function(x, y) {
  nn <- length(x)
  stopifnot(nn >= 3L, nn %% 2L == 1L,
    max(abs(diff(x) - (x[nn] - x[1]) / (nn - 1L))) < 1e-10)
  (x[nn] - x[1]) / (3 * (nn - 1L)) *
    sum(y * c(1, rep(c(4, 2), (nn - 3L) / 2L), 4, 1))
}

latest_max_constants <- function(n, p, eta = .1,
    cache_dir = file.path("output", "latest_max_calibration_cache"),
    t_nodes = 17L, endpoint_M = 1024L, moment_M = 4096L,
    ode_tol = 1e-10, verbose = TRUE) {
  stopifnot(n > 2, p > 1, eta > 0, eta < .5,
    t_nodes >= 5L, (t_nodes - 1L) %% 4L == 0,
    endpoint_M >= 128L, moment_M >= 256L)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache <- function(key, calc) {
    fn <- file.path(cache_dir, paste0(key, ".rds"))
    if (file.exists(fn)) return(readRDS(fn))
    obj <- calc()
    saveRDS(obj, fn)
    obj
  }
  hardy <- cache(sprintf("hardy_v2_eta%.8f", eta),
    function() latest_max_hardy_b(eta, q_max = 24, tol = 1e-11))
  branch <- latest_max_branch(n, p, eta, hardy)
  if (branch$branch != "subcritical") stop(
    "Current deterministic implementation supports only the manuscript subcritical branch; required branch is ",
    branch$branch, ". No surrogate critical or boundary constants are substituted.")
  N <- log(n); gam <- log(p) / N
  q <- (sqrt(gam) + sqrt(1 + gam))^2
  theta <- (1 - q^(-2)) / 4
  ts <- seq(eta, .5, length.out = t_nodes)
  rows <- lapply(seq_along(ts), function(i) {
    tt <- ts[i]
    mom <- cache(sprintf("mom_v1_t%.10f_M%d", tt, moment_M),
      function() latest_max_moment_series(tt,
        as.integer(moment_M / c(16,8,4,2,1))))
    ps <- cache(sprintf("psi_v2_t%.10f_th%.14f_M%d_tol%.1e", tt, theta, endpoint_M, ode_tol),
      function() {
        if (verbose) message(sprintf("MAX Bellman: theta=%.8f, t=%.4f (%d/%d)", theta,tt,i,length(ts)))
        zz <- latest_max_psi(tt, theta,
          as.integer(endpoint_M / c(8,4,2,1)), ode_tol = ode_tol)
        zz$bellman <- NULL
        zz
      })
    last <- nrow(mom); lastp <- nrow(ps$convergence)
    data.frame(t = tt, m = mom$m_tail1[last], v0 = mom$v0_extrapolated[last],
      m_previous = mom$m_tail1[last-1L], v0_previous = mom$v0_extrapolated[last-1L],
      logPsi = ps$logPsi, logPsi_previous = ps$convergence$logPsi[lastp-1L],
      moment_probability = mom$probability[last])
  })
  tab <- do.call(rbind, rows)
  d0 <- tab$m + (q - 1) * tab$v0 / 4
  integrand <- exp(tab$logPsi - theta * d0) / (tab$t * (1 - tab$t))
  Theta <- 2 * latest_max_simpson(tab$t, integrand)
  co <- seq(1L, nrow(tab), by = 2L)
  Theta_coarse_grid <- 2 * latest_max_simpson(tab$t[co], integrand[co])
  d0_previous <- tab$m_previous + (q - 1) * tab$v0_previous / 4
  Theta_previous_series <- 2 * latest_max_simpson(tab$t,
    exp(tab$logPsi_previous - theta * d0_previous) / (tab$t * (1 - tab$t)))
  lambda2 <- 2 * q^3
  a <- 1 / (theta * sqrt(2*N))
  Cn <- (q - 1) * sqrt(N) * Theta / sqrt(2*pi*lambda2)
  center <- (q - 1) * sqrt(N/2) + a * log(Cn)
  list(a = a, b = center, branch = "subcritical", n = n, p = p, eta = eta,
    theta = theta, q = q, Theta = Theta, Cn = Cn, b_eta = hardy$b,
    q_c = branch$q_c, gamma_c = branch$gamma_c, gamma_n = gam,
    cutoff_05 = center + a * (-log(-log(.95))),
    diagnostics = list(t_nodes_half_interval = t_nodes,
      endpoint_M = endpoint_M, moment_M = moment_M, ode_tol = ode_tol,
      max_logPsi_truncation_change = max(abs(tab$logPsi-tab$logPsi_previous)),
      max_m_series_change = max(abs(tab$m-tab$m_previous)),
      max_v0_series_change = max(abs(tab$v0-tab$v0_previous)),
      center_change_fine_minus_coarse_grid = a * log(Theta / Theta_coarse_grid),
      center_change_fine_minus_previous_series = a * log(Theta / Theta_previous_series),
      probability_mass_error = max(abs(tab$moment_probability-1))),
    constants_grid = tab,
    theory_source = "draft/hdcp_terminal_anchor_corrected_20260909/parts/body.tex equations main-subcritical-normalization and main-endpoint-coefficients")
}

