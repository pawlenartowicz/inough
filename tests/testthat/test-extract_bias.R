test_that("extract_bias encodes lag-1 transitions correctly", {
  b <- inough:::extract_bias(data.frame(
    id = c(rep("A", 5), rep("B", 3)),
    response = c(0, 0, 1, 0, 1,  0, 0, 0)
  ))

  expect_equal(b$trial$resp_lag1[1], 0)   # first trial of A
  expect_equal(b$trial$resp_lag1[2], 1)   # 0 -> 0 (same)
  expect_equal(b$trial$resp_lag1[3], -1)  # 0 -> 1 (different)
  expect_equal(b$trial$resp_lag1[6], 0)   # first trial of B
})

test_that("extract_bias returns one LZ value per participant", {
  b <- inough:::extract_bias(data.frame(
    id = c(rep("A", 5), rep("B", 3)),
    response = c(0, 0, 1, 0, 1,  0, 0, 0)
  ))

  expect_length(b$participant$lz, 2)
})
