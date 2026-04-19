test_that("bailout flags low LZ, chance accuracy, and both", {
  bo <- inough:::bailout(
    data.frame(
      id            = c("A", "B", "C"),
      lz            = c(0.2, 0.8, 0.1),
      n_trials      = c(100, 100, 100),
      mean_accuracy = c(0.7, 0.52, 0.48),
      stringsAsFactors = FALSE
    ),
    lz_threshold = 0.3
  )

  expect_true("A" %in% bo$id)  # low LZ
  expect_true("B" %in% bo$id)  # chance accuracy
  expect_true("C" %in% bo$id)  # both
  expect_equal(bo$reason[bo$id == "A"], "lz")
  expect_equal(bo$reason[bo$id == "C"], "both")
})
