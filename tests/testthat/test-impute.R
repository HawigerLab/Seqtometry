test_that("PCA (.calc_pca and .inv_pca) works", {
  m <- matrix(ncol = 3, data = c(
    -1, -2, -3, 1, 2, 3,
    -1, -1, -2, 1, 1, 2,
    4, 5, 6, -4, -5, -6))
  p <- prcomp(m, scale. = TRUE)
  r <- .calc_pca(m, 2, TRUE)
  expect_equal(norm(r$v, "F"), norm(p$rotation[, -3], "F"))

  i <- .invert_pca(p$x, p$rotation, p$center, p$scale, FALSE) |> t()
  expect_equal(i, m)
})

test_that(".procrustes works", {
  m1 <- matrix(c(1, 1, 1, 2, 3, 2, 1, 1), nrow = 2)
  m2 <- matrix(c(4, 4, 4, 2, -2, -4, -6, -6), nrow = 2)
  expect_equal(.procrustes(m1, m2), 0)
})

test_that("impute works", {
  m <- matrix(ncol = 14, c(
    1.1, 2, 3, 0, 0, 0,
    1, 2.1, 3, 0, 0, 0,
    1, 2, 3.1, 0, 0, 0,
    1.1, 2.1, 3, 0, 0, 0,
    1.1, 2, 3.1, 0, 0, 0,
    1, 2.1, 3.1, 0, 0, 0,
    0, 0, 0, 1.1, 2, 3,
    0, 0, 0, 1, 2.1, 3,
    0, 0, 0, 1, 2, 3.1,
    0, 0, 0, 1.1, 2.1, 3,
    0, 0, 0, 1.1, 2, 3.1,
    0, 0, 0, 1, 2.1, 3.1,
    1, 0, 3, 0, 0, 0,
    0, 0, 0, 1, 0, 3))
  i <- impute(m, npc = 2L, knn = 6L, ka = 2L)
  expect_gt(i[2, 13], 1.5)
  expect_equal(sum(i[4:6, 13]), 0)
  expect_gt(i[5, 14], 1.5)
  expect_equal(sum(i[1:3, 14]), 0)
})