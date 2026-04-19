test_that("constant sequences have low complexity", {
  expect_lt(lz_complexity(rep(0, 100)), 0.3)
})

test_that("simple alternation has low complexity", {
  expect_lt(lz_complexity(rep(c(0, 1), 50)), 0.4)
})

test_that("random binary sequences have high complexity", {
  set.seed(42)
  expect_gt(lz_complexity(sample(0:1, 200, TRUE)), 0.7)
})

test_that("single-element sequences return 1.0", {
  expect_equal(lz_complexity(c(0)), 1.0)
  expect_equal(lz_complexity(c(1)), 1.0)
})

test_that("invalid inputs raise an error", {
  expect_error(lz_complexity(c()))
  expect_error(lz_complexity(c(0, 2)))
})
