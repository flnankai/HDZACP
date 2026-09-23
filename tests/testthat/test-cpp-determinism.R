test_that("compiled analytic null simulation is invariant to thread count", {
  simulate_cpp <- getFromNamespace("hdzacp_size_simulate_cpp", "HDZACP")
  one <- simulate_cpp(
    n = 16L, p = 5L, eta = 0.2, n_rep = 8L, seed = 931L,
    n_threads = 1L, design = "ar1", rho = 0.3, block_size = 3L
  )
  two <- simulate_cpp(
    n = 16L, p = 5L, eta = 0.2, n_rep = 8L, seed = 931L,
    n_threads = 2L, design = "ar1", rho = 0.3, block_size = 3L
  )

  invariant <- setdiff(names(one), "threads_used")
  expect_equal(one[invariant], two[invariant], tolerance = 0)
  expect_identical(one$threads_used, 1L)
  expect_true(two$threads_used %in% c(1L, 2L))
})

test_that("OU reference paths are invariant to thread count", {
  one <- hdzacp_ou_reference(
    eta = 0.2, paths = 12L, grid = 64L, threads = 1L, seed = 55L
  )
  two <- hdzacp_ou_reference(
    eta = 0.2, paths = 12L, grid = 64L, threads = 2L, seed = 55L
  )

  expect_s3_class(one, "hdzacp_ou_reference")
  expect_equal(one$coarse, two$coarse, tolerance = 0)
  expect_equal(one$fine, two$fine, tolerance = 0)
  expect_true(all(one$fine >= one$coarse))
})
