test_that("rolling_wmean handles edges and centre with uniform weights", {
  out <- inough:::rolling_wmean(c(1, 2, 3, 4, 5), weights = rep(1, 5))

  expect_length(out, 5)
  expect_equal(out[1], mean(c(1, 2, 3)))     # left edge: window cropped to 3
  expect_equal(out[3], mean(c(1, 2, 3, 4, 5)))  # centre: full window
  expect_equal(out[5], mean(c(3, 4, 5)))     # right edge: window cropped to 3
})

test_that("rolling_wmean respects unequal weights", {
  out <- inough:::rolling_wmean(c(0, 0, 1, 0, 0), weights = c(1, 2, 1))
  expect_equal(out[3], (0 + 2 * 1 + 0) / 4)  # centre weighted twice
})

test_that("robust_zscore returns values when MAD > 0", {
  zs <- inough:::robust_zscore(c(1, 2, 3, 4, 5))
  expect_false(is.null(zs))
})

test_that("robust_zscore returns NULL when MAD is zero", {
  expect_null(inough:::robust_zscore(c(0, 0, 0, 0, 0)))
  expect_null(inough:::robust_zscore(c(0, 0, 0, 0, 10)))
})
