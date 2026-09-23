test_that("WBS honors supplied intervals and returns valid change-point objects", {
  x <- simulate_hdzacp(
    n = 30, p = 4, change_points = 15, s = 4,
    dependence = "independent", alternative = "student", df = 2.05,
    seed = 715
  )
  intervals <- matrix(c(0L, 30L), nrow = 1L,
                      dimnames = list(NULL, c("start", "end")))

  fit <- hdzacp_wbs(
    x, methods = c("ZAC", "ZAM", "ZAS"), eta = 0.2,
    alpha = 0.5, intervals = intervals, min_interval = 20,
    permutations = 31, stopping = "unadjusted", seed = 51
  )

  expect_s3_class(fit, "hdzacp_wbs")
  expect_identical(fit$intervals, intervals)
  expect_named(fit$changepoints, c("ZAC", "ZAM", "ZAS"))
  expect_named(fit$segments, c("ZAC", "ZAM", "ZAS"))
  expect_equal(nrow(fit$evidence), 1L)
  expect_true(all(c("start", "end", "p_sum", "p_max", "p_cauchy",
                    "cp_sum", "cp_max", "cp_cauchy") %in%
                  names(fit$evidence)))
  expect_true(all(vapply(fit$changepoints, function(cp) {
    is.integer(cp) && all(cp >= 1L & cp < nrow(x))
  }, logical(1L))))
  expect_true(any(lengths(fit$changepoints) > 0L))

  for (method in names(fit$segments)) {
    segments <- fit$segments[[method]]
    expect_s3_class(segments, "data.frame")
    expect_equal(segments[1L, "start"], 1L)
    expect_equal(segments[nrow(segments), "end"], nrow(x))
  }
})

test_that("simulation has the requested design, dimensions, and seed stability", {
  designs <- c("independent", "AR1", "CS")
  for (dependence in designs) {
    x_one <- simulate_hdzacp(
      n = 24, p = 6, change_points = c(7, 17), s = 3,
      dependence = dependence, rho = if (dependence == "independent") 0 else 0.3,
      alternative = "gamma", gamma_shape = 3,
      pattern = "epidemic", seed = 810
    )
    x_two <- simulate_hdzacp(
      n = 24, p = 6, change_points = c(7, 17), s = 3,
      dependence = dependence, rho = if (dependence == "independent") 0 else 0.3,
      alternative = "gamma", gamma_shape = 3,
      pattern = "epidemic", seed = 810
    )

    expect_identical(dim(x_one), c(24L, 6L))
    expect_true(all(is.finite(x_one)))
    expect_equal(unclass(x_one), unclass(x_two), tolerance = 0)
    design <- attr(x_one, "hdzacp_design")
    expect_identical(design$dependence, dependence)
    expect_identical(design$change_points, c(7L, 17L))
    expect_identical(design$s, 3L)
    expect_identical(design$pattern, "epidemic")
  }
})

test_that("segmentation metrics attain their exact-match optimum", {
  exact <- hdzacp_segmentation_metrics(c(3, 7), c(3, 7), n = 10)
  expect_equal(unname(exact[c("n_true", "n_estimated")]), c(2, 2))
  expect_equal(unname(exact[["count_error"]]), 0)
  expect_equal(unname(exact[["exact_count"]]), 1)
  expect_equal(unname(exact[["scaled_hausdorff"]]), 0)
  expect_equal(unname(exact[["ARI"]]), 1)

  shifted <- hdzacp_segmentation_metrics(c(3, 7), c(4, 7), n = 10)
  expect_equal(unname(shifted[["scaled_hausdorff"]]), 0.25)
  expect_equal(unname(shifted[["true_to_estimated"]]), 0.25)
  expect_equal(unname(shifted[["estimated_to_true"]]), 0.25)
})

