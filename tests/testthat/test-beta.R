test_that("compute.alpha works", {
  # alpha for identity matrix
  R <- diag(10)
  expect_equal(compute.alpha(R), 0)

  # more interesting case
  R <- matrix(c(1, .3, .3,
                 .3, 1, .2,
                 .3, .2, 1), 3, 3)
  alpha <- compute.alpha(R)
  expect_true(is.numeric(alpha))
  expect_true(alpha >= 0)
})

test_that("compute.beta works", {
  R <- matrix(c(1, .4, .2,
                .4, 1, .3,
                .2, .3, 1), 3, 3)
  rownames(R) <- colnames(R) <- as.character(1:3)

  # split sets
  beta <- compute.beta(R, "1", "3")
  expect_true(beta >= 0 && beta <= 1)

  # larger groups
  beta <- compute.beta(R, c("1","3"), "2")
  expect_true(beta >= 0 && beta <= 1)
})


test_that("hc.beta fails on invalid inputs", {
  R <- matrix(rnorm(9), 3, 3)

  expect_error(hc.beta(1:3))
  expect_error(hc.beta(R, is.corr = TRUE))
  expect_error(hc.beta(R, is.corr = FALSE, n.splits = - 2))
  expect_error(hc.beta(R, is.corr = FALSE, splits = "HCSVD"))
})


test_that("is.corr.matrix works", {
  # Valid cases
  expect_true(is.corr.matrix(diag(3)))
  expect_true(is.corr.matrix(cor(mtcars[, 1:3])))

  # Not a matrix
  expect_false(is.corr.matrix(1:4))

  # Asymmetric matrix
  A <- matrix(c(1, 0.2, 0.3,
                0.1, 1, 0.4,
                0.3, 0.4, 1), 3, 3)
  expect_false(is.corr.matrix(A))

  # Wrong diagonal
  B <- diag(c(1, 2, 1))
  expect_false(is.corr.matrix(B))

  # Negative eigenvalue
  C <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_false(is.corr.matrix(C))
})

