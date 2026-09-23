test_that("scan and permutation test return stable public structures", {
  x <- simulate_hdzacp(
    n = 16, p = 5, change_points = 8, s = 2,
    dependence = "AR1", rho = 0.35, alternative = "student",
    df = 3, seed = 104
  )

  scan <- hdzacp_scan(x, eta = 0.2)
  expect_s3_class(scan, "hdzacp_scan")
  expect_true(all(c("statistics", "estimates", "estimate_fraction",
                    "diagnostics", "n", "p", "eta", "m") %in% names(scan)))
  expect_named(scan$statistics, c("ZAS", "ZAM"))
  expect_named(scan$estimates, c("ZAS", "ZAM"))
  expect_true(all(is.finite(scan$statistics)))
  expect_true(all(scan$estimates >= scan$m & scan$estimates <= scan$n - scan$m))

  fit_one <- hdzacp_test(
    x, eta = 0.2, permutations = 19, seed = 817,
    keep_permutations = TRUE
  )
  fit_two <- hdzacp_test(
    x, eta = 0.2, permutations = 19, seed = 817,
    keep_permutations = TRUE
  )

  expect_s3_class(fit_one, "hdzacp_test")
  expect_equal(fit_one$results, fit_two$results, tolerance = 0)
  expect_equal(fit_one$permutation_statistics,
               fit_two$permutation_statistics, tolerance = 0)
  expect_identical(fit_one$results$method, c("ZAC", "ZAM", "ZAS"))
  expect_true(all(c("statistic", "p.value", "reject", "estimate",
                    "estimate_fraction", "calibration") %in%
                  names(fit_one$results)))
  expect_true(all(fit_one$results$p.value >= 1 / 20 &
                  fit_one$results$p.value <= 1))
  expect_true(all(fit_one$results$estimate >= fit_one$settings$m &
                  fit_one$results$estimate <= fit_one$settings$n - fit_one$settings$m))
  expect_false(is.null(fit_one$permutation_statistics))
})

test_that("method-specific test helpers select exactly one method", {
  x <- simulate_hdzacp(n = 12, p = 3, dependence = "independent",
                       alternative = "normal", seed = 22)
  helpers <- list(ZAC = zac_test, ZAM = zam_test, ZAS = zas_test)

  for (method in names(helpers)) {
    fit <- helpers[[method]](x, eta = 0.2, permutations = 7, seed = 90)
    expect_s3_class(fit, "hdzacp_test")
    expect_identical(fit$results$method, method)
  }
})

test_that("ties are rejected by default and randomized reproducibly", {
  tied <- cbind(
    rep(1:4, each = 3),
    rep(c(0, 1), each = 6),
    rep(1:3, times = 4)
  )

  expect_error(hdzacp_scan(tied, eta = 0.2), "Ties were found")
  expect_error(
    hdzacp_test(tied, calibration = "asymptotic", ties = "random",
                tie_seed = 44),
    "only with permutation"
  )

  fit_one <- hdzacp_test(
    tied, eta = 0.2, permutations = 11, seed = 73,
    ties = "random", tie_seed = 901
  )
  fit_two <- hdzacp_test(
    tied, eta = 0.2, permutations = 11, seed = 73,
    ties = "random", tie_seed = 901
  )
  expect_equal(fit_one$results, fit_two$results, tolerance = 0)
  expect_identical(fit_one$settings$ties, "random")
  expect_identical(fit_one$settings$tie_seed, 901)
})

