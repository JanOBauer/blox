test_that("is.sim.matrix works", {
  # Valid cases
  expect_true(is.sim.matrix(diag(3)))
  expect_true(is.sim.matrix(cor(mtcars[, 1:3])))

  # Not a matrix
  expect_error(is.sim.matrix(1:4))

  # Asymmetric matrix
  A <- matrix(c(1, 0.2, 0.3,
                0.1, 1, 0.4,
                0.3, 0.4, 1), 3, 3)
  expect_error(is.sim.matrix(A))

  #Value one on off-diagonal
  A <- matrix(c(1, 0, 1,
                0, 1, 0,
                1, 0, 1), 3, 3)
  expect_warning(is.sim.matrix(A))

  # Wrong diagonal
  B <- diag(c(1, 2, 1))
  expect_error(is.sim.matrix(B))

  # Negative eigenvalue
  C <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(is.sim.matrix(C))
})



test_that("calc.distance works", {
  set.seed(123)
  S <- cor(matrix(rnorm(100), 10, 10))
  cluster1 <- 1:3
  cluster2 <- 4:10
  cluster <- list(as.character(cluster1), as.character(cluster2))

  d_avg <- calc.distance(S, cluster, as.character(1:10), "average")
  expect_equal(d_avg, 1 - mean(abs(S[cluster1, cluster2])))

  d_single <- calc.distance(S, cluster, as.character(1:10), "single")
  expect_equal(d_single, 1 - max(abs(S[cluster1, cluster2])))

  d_rv <- calc.distance(S, cluster, as.character(1:10), "RV")
  expect_true(d_rv >= 0 && d_rv <= 1)

  d_ncut <- calc.distance(S, cluster, as.character(1:10), "ncut")
  expect_true(is.numeric(d_ncut))

  d_rcut <- calc.distance(S, cluster, as.character(1:10), "rcut")
  expect_true(is.numeric(d_rcut))
})



test_that("hcsvd returns hclust object", {
  set.seed(123)
  S <- cor(matrix(rnorm(200), 20, 10))

  hc <- hcsvd(S)

  expect_s3_class(hc$hclust, "hclust")
  expect_true(all(hc$dist.mat >= 0))
  expect_equal(length(hc$q.p), ncol(S) - 1)
})



test_that("hcsvd fails on invalid inputs", {
  S <- matrix(rnorm(9), 3, 3)

  expect_error(hcsvd(1:3))
  expect_error(hcsvd(S))

  S <- cor(matrix(rnorm(200), 20, 10))
  expect_error(hcsvd(S, q = "one"))
  expect_error(hcsvd(S, h.power = -1))
  expect_error(hcsvd(S, linkage = "avg"))
})



test_that("bd.approx produces block diagonal structure", {
  set.seed(123)
  S <- cor(matrix(rnorm(200), 20, 10))

  bd <- bd.approx(S, balance = 0.5, linkage = "average")

  # Check number of split elements
  expect_equal(length(bd$split$B1) + length(bd$split$B2), ncol(S))

  # Check balance structure
  expect_equal(sum(bd$BD[1, ] == 0), ncol(S)/2)
})


test_that("bd.approx handles errors", {
  S <- cor(matrix(rnorm(100), 10, 10))

  expect_error(bd.approx(1:5))
  expect_error(bd.approx(S, balance = -0.1))
})


