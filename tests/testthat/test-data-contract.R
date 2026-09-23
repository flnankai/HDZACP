test_that("empirical data are referenced by official HTTPS links, not bundled", {
  urls <- hdzacp_data_urls()

  expect_s3_class(urls, "data.frame")
  expect_named(urls, c("data", "landing_page", "download_url"))
  expect_equal(nrow(urls), 2L)
  expect_true(all(nzchar(urls$data)))
  expect_true(all(grepl("^https://", urls$landing_page)))
  expect_true(all(grepl("^https://", urls$download_url)))
  expect_true(any(grepl("GSE30272", urls$download_url, fixed = TRUE)))
  expect_true(any(grepl("archive.ics.uci.edu", urls$download_url, fixed = TRUE)))

  expect_identical(system.file("data", package = "HDZACP"), "")
  expect_identical(system.file("extdata", package = "HDZACP"), "")
})
