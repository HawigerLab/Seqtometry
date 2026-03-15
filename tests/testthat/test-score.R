test_that(".gene_indices works", {
  m <- matrix(nrow = 5, ncol = 1, dimnames = list(letters[1:5], "Cell"))
  g <- list(Sig1 = c("a", "c"), Sig2 = c("b", "z"))
  expect_equal(
    .gene_indices(m, g),
    list(Sig1 = c(0L, 2L), Sig2 = c(1L)))
})

test_that(".minmax_scale works", {
  expect_equal(
    .minmax_scale(c(-2, 0, 1, -1, 2)),
    c(0, 0.5, 0.75, 0.25, 1))
})

test_that("score works", {
  m <- matrix(ncol = 3, dimnames = list(letters[1:5], 1:3), c(
    1, 2, 0, 4, 5,
    8, 7, 9, 0, 2,
    3, 4, 2, 0, 1))
  g <- list(Sig = c("a", "c"))
  expect_equal(score(m, g)[["Sig"]], c(0, 1, 0.5))
})
