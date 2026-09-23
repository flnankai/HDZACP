enumerated_rank_statistic <- function(n, k, selected_ranks) {
  entropy <- function(u) {
    if (u <= 0 || u >= 1) return(0)
    -(u * log(u) + (1 - u) * log1p(-u))
  }

  cumulative <- cumsum(seq_len(n) %in% selected_ranks)
  terms <- vapply(seq_len(n), function(ell) {
    count <- cumulative[[ell]]
    numerator <- k * entropy(count / k) +
      (n - k) * entropy((ell - count) / (n - k))
    numerator / ((ell - 0.5) * (n - ell + 0.5))
  }, numeric(1L))
  sum(terms)
}

enumerated_null_moments <- function(n, k) {
  subsets <- utils::combn(n, k)
  values <- apply(subsets, 2L, function(selected) {
    enumerated_rank_statistic(n, k, selected)
  })
  center <- mean(values)
  c(mu = center, sd = sqrt(mean((values - center)^2)))
}

test_that("exact null moments agree with complete finite-sample enumeration", {
  for (n in c(6L, 8L)) {
    moments <- hdzacp_null_moments(n, eta = 0.2)

    expect_s3_class(moments, "data.frame")
    expect_named(moments, c("k", "mu", "sd"))
    expect_true(all(is.finite(moments$mu)))
    expect_true(all(is.finite(moments$sd) & moments$sd > 0))

    for (i in seq_len(nrow(moments))) {
      expected <- enumerated_null_moments(n, moments$k[[i]])
      expect_equal(moments$mu[[i]], unname(expected[["mu"]]), tolerance = 1e-12)
      expect_equal(moments$sd[[i]], unname(expected[["sd"]]), tolerance = 1e-12)
    }

    expect_equal(moments$mu, rev(moments$mu), tolerance = 0)
    expect_equal(moments$sd, rev(moments$sd), tolerance = 0)

    mu_full <- attr(moments, "mu_full")
    sd_full <- attr(moments, "sd_full")
    expect_length(mu_full, n + 1L)
    expect_length(sd_full, n + 1L)
    expect_equal(mu_full[moments$k + 1L], moments$mu, tolerance = 0)
    expect_equal(sd_full[moments$k + 1L], moments$sd, tolerance = 0)
  }
})
