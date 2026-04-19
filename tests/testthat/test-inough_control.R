test_that("inough_control() returns an object of the expected class", {
  ctrl <- inough_control()
  expect_s3_class(ctrl, "inough_control")
})

test_that("default values match the documented defaults", {
  ctrl <- inough_control()
  expect_equal(ctrl$lz_threshold, 0.2)
  expect_equal(ctrl$window_size,  3L)
  expect_equal(ctrl$sd_threshold, 2)
  expect_identical(ctrl$window_weight, "uniform")
  expect_equal(ctrl$min_chunk, 6L)
  expect_identical(ctrl$comparison, "clean")
})

test_that("screening_threshold is computed from sd_threshold and window weights", {
  ctrl <- inough_control()
  expect_equal(ctrl$screening_threshold, ctrl$sd_threshold * ctrl$window_sd)
})

test_that("triangular weights yield a larger window SD than uniform", {
  # Kish design effect: concentrating weight on the centre increases the
  # variance of the weighted mean per unit weight, so SD is larger than the
  # equal-weighted case.
  uni <- inough_control(window_weight = "uniform")
  tri <- inough_control(window_weight = "triangular")
  expect_gt(tri$window_sd, uni$window_sd)
})

test_that("invalid arguments are rejected", {
  expect_error(inough_control(lz_threshold = -0.1))
  expect_error(inough_control(window_size  = 0))
  expect_error(inough_control(sd_threshold = 0))
  expect_error(inough_control(window_weight = "bogus"))
  expect_error(inough_control(comparison    = "bogus"))
})
