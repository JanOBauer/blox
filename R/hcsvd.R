## usethis namespace: start
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib blox, .registration = TRUE
## usethis namespace: end


#' @importFrom matrixStats mean2
calc.distance <- function(S, idx.cluster, feature.names, linkage) {
  cluster1 <- feature.names %in% idx.cluster[[1]]
  cluster2 <- feature.names %in% idx.cluster[[2]]

  if (linkage == "average") {
    return(1 - mean2(abs(S[cluster1, cluster2])))
  }

  if (linkage == "ncut") {
    S.abs <- abs(S)
    cut <- sum(S.abs[cluster1, cluster2])
    ass1 <- sum(S.abs[cluster1, ])
    ass2 <- sum(S.abs[cluster2, ])

    return(- cut/ass1 - cut/ass2)
  }

  if (linkage == "rcut") {
    S.abs <- abs(S)
    cut <- sum(S.abs[cluster1, cluster2])

    return(- cut/length(idx.cluster[[1]]) - cut/length(idx.cluster[[2]]) )
  }

  if (linkage == "RV") {
    S1  <- S[cluster1, cluster1]
    S2  <- S[cluster2, cluster2]
    S12 <- S[cluster1, cluster2]

    return(1 - sum(S12^2) / sqrt(sum(S1^2) * sum(S2^2)))
  }

  if (linkage == "single") {
    return(1 - max(abs(S[cluster1, cluster2])))
  }
}



hcsvd.ht <- function(S,
                     S.h,
                     q,
                     linkage,
                     labels,
                     max.iter
) {

  p <- ncol(S)
  if (p == 2) {
    idx.cluster <- list(labels[1], labels[2])
    cluster.distance <- calc.distance(S, idx.cluster, labels, linkage)

    return(list(cluster1 = idx.cluster[1], cluster2 = idx.cluster[2], cluster.distance = cluster.distance, q.p = NA))
  }


  SVD <- eigen(S.h, symmetric = TRUE)

  if (is.numeric(q)) {
    if (q == 1){
      q <- p - 1
    } else {
      q <- min(round(q * p), p - 1)
    }
  }
  if (q == "Kaiser") {
    q <- min(sum(SVD$values >= 1) + 2, p - 1)
  }

  SVD$vectors <- SVD$vectors[, 1:q]

  dof.grid <- 1 : (p - 1)
  distance <- -Inf

  for (dof in dof.grid) {
    V <- calc_sparse_v_cpp(X = S.h, v = SVD$vectors, q = q, n = dof, maxit = max.iter)

    for (i in seq_len(q)) {
      v <- V[, i]
      idx.cluster <- list(labels[v != 0], labels[v == 0])
      distance.ht <- calc.distance(S, idx.cluster, labels, linkage)


      if (distance.ht > distance) {
        distance <- distance.ht
        cluster1 <- idx.cluster[[1]]
        cluster2 <- idx.cluster[[2]]

      }
    }
  }


  cluster.distance <- distance

  return(list(cluster1 = cluster1, cluster2 = cluster2, cluster.distance = cluster.distance, q.p = q / p))
}



is.sim.matrix <- function(S) {

  if (!is.matrix(S)) {
    stop("The provided similarity matrix is not a matrix.")
  }

  if (!all(diag(S) == 1)) {
    stop("The provided similarity matrix does not have ones on the diagonal.")
  }

  if (any(abs(S) > 1)) {
    stop("The provided similarity matrix is not scaled.")
  }

  if (nrow(S) != ncol(S)) {
    stop("The provided similarity matrix is not quadratic.")
  }

  if (!isSymmetric(S)) {
    stop("The provided similarity matrix is not symmetric")
  }

  if (any(abs(S[upper.tri(S)]) == 1)) {
    warning("The provided similarity matrix contains off-diagonal entries equal to 1 or -1.")
  }

  eigenvalues <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if(any(eigenvalues < 0)) {
    stop("The provided similarity matrix is not positive semi-definite.")
  }

  return(TRUE)
}



#' @title Hierarchical Clustering Using Singular Vectors (HC-SVD).
#'
#' @description Performs HC-SVD to reveal the hierarchical structure as described in Bauer (202Xa). This divisive approach iteratively splits each cluster into two subclusters.
#' Candidate splits are determined by the first sparse eigenvectors (sparse approximations of the first eigenvectors, i.e., vectors with many zero entries) of the similarity matrix.
#' The selected split is the one that yields the best block-diagonal approximation of the similarity matrix according to a specified linkage function. The procedure continues until each object is assigned to its own cluster.
#'
#' @param S A scaled \eqn{p}x\eqn{p} similarity matrix. For example, this may be a correlation matrix.
#'
#' @param linkage The linkage function to be used. This should be one of \code{"average"}, \code{"single"}, or
#' \code{"RV"} (for RV-coefficient). Note that the RV-coefficient might not yield an ultrametric distance.
#'
#' @param q Number of sparse eigenvectors to be used. This should be either a numeric value between zero and one to indicate percentages, or \code{"Kaiser"} for as many sparse eigenvectors as
#' there are eigenvalues larger or equal to one. For a numerical value between zero and one, the number of sparse eigenvectors is determined as the corresponding share of the total number of eigenvectors.
#' E.g., \code{q = 1} (100%) uses all sparse eigenvectors and \code{q = 0.5} (50%) will use half of all sparse eigenvectors. For \code{q = 1}, identification is best (see Bauer (202Xa) for details).
#'
#' @param h.power \code{h}-th Hadamard power of \code{S}. This should be a positive integer and increases robustness of the method, as described in Bauer (202Xa).
#'
#' @param max.iter How many iterations should be performed for computing the sparse eigenvectors.
#' Default is \code{500}.
#'
#' @param verbose Print out progress as \eqn{p-1} iterations for divisive hierarchical clustering are performed.
#' Default is \code{TRUE}.
#'
#' @details
#' The sparse loadings are computed using the method proposed by Shen & Huang (2008). The corresponding implementation is written in \code{Rcpp}/\code{RcppArmadillo}
#' for computational efficiency and is based on the \code{R} implementation by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} (\pkg{irlba}).
#' However, the implementation has been adapted to better align with the scope of the \pkg{bdsvd} package which is the base for the \pkg{blox} package.
#'
#' Supplementary details are in \code{\link[blox]{hc.beta}} and in Bauer (202Xb).
#'
#' @return
#' A list with four components:
#' \item{hclust}{
#'  The clustering structure identified by HC-SVD as an object of type \code{hclust}.
#' }
#' \item{dist.matrix}{
#'  The ultrametric distance matrix (cophenetic matrix) of the HC-SVD structure as an object of class \code{dist}.
#' }
#' \item{u.sim}{
#'  The ultrametric similarity matrix of \eqn{S} obtained by HC-SVD as an object of class \code{matrix}. The ultrametric similarity matrix
#'  is calculated as \code{1-dist.matrix}.
#' }
#' \item{q.p}{
#'  A vector of length \eqn{p-1} containing the ratio \eqn{q_i/p_i} of the \eqn{q_i} sparse eigenvectors used relative to all sparse
#'  eigenvectors \eqn{q_i} for the split of each cluster. The ratio is set to \code{NA} if the cluster contains only two variables as the search
#'  for sparse eigenvectors that reflect this obvious split is not required in this case.
#' }
#'
#' @seealso \code{\link[bdsvd]{bdsvd}} \{\link[bdsvd]{bdsvd}\}
#'
#' @references \cite{Bauer, J.O. (202Xa). Divisive hierarchical clustering using block diagonal matrix approximations. Working paper.}
#' @references \cite{Bauer, J.O. (202Xb). Revelle's beta: The wait is over - we can compute it!. Working paper.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #We give one example for variable clustering directly on a correlation matrix,
#' #and we replicate the USArrest example in Bauer (202Xa) for observation clustering.
#' #More elaborate code alongside a different example for variable clustering can be
#' #found in the corresponding supplementary material of that manuscripts.
#'
#' \donttest{
#' ### VARIABLE CLUSTERING
#'
#' #Load the correlation matrix Bechtoldt from the psych
#' #package (see ?Bechtoldt for more information).
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   data("Bechtoldt", package = "psych")
#' }
#'
#' #Compute HC-SVD (with average linkage).
#' hcsvd.obj <- hcsvd(Bechtoldt)
#'
#' #The object of type hclust with corresponding dendrogram can be obtained
#' #directly from hcsvd(...):
#' hc.div <- hcsvd.obj$hclust
#' plot(hc.div, ylab = "")
#'
#' #The dendrogram can also be obtained from the ultrametric distance matrix:
#' plot(hclust(hcsvd.obj$dist.matrix), main = "HC-SVD", sub = "", xlab = "")
#'
#'
#' ### OBSERVATION CLUSTERING
#'
#' #Correct for the known transcription error
#' data("USArrests")
#' USArrests["Maryland", "UrbanPop"] <- 76.6
#'
#' #The distance matrix is scaled (divided by max(D)) to later allow a
#' #transformation to a matrix S that fulfills the properties of a similarity
#' #matrix.
#' D <- as.matrix(dist(USArrests))
#' D <- D / max(D)
#' S <- 1 - D
#'
#' #Compute HC-SVD (with average linkage).
#' hcsvd.obj <- hcsvd(S)
#'
#' #The object of type hclust with corresponding dendrogram can be obtained
#' #directly from hcsvd(...):
#' hc.div <- hcsvd.obj$hclust
#' plot(hc.div, ylab = "")
#'
#' #The dendrogram can also be obtained from the ultrametric distance matrix:
#' plot(hclust(hcsvd.obj$dist.matrix), main = "HC-SVD", sub = "", xlab = "")
#' }
#'
#'
#' @importFrom stats cov
#' @importFrom stats as.dist
#'
#' @export
hcsvd <- function(S,
                  linkage = "average",
                  q = 1,
                  h.power = 2,
                  max.iter,
                  verbose = TRUE
) {

  if (missing(S)) {
    stop("Provide S.")
  }

  if (!is.sim.matrix(S)) {
    stop("S must be a similarity matrix.")
  }

  if (!((q == "Kaiser") | (is.numeric(q) & (q > 0) & (q <= 1)))) {
    stop(sprintf("%s is an invalid argument for q.", q))
  }

  if(h.power <= 0){
    stop("h.power must be positive.")
  }
  S.h <- S^h.power

  if(!(is.numeric(h.power) && h.power == floor(h.power))){
    warning("h.power should be an integer.")
  }

  LINKAGE <- c("average", "RV", "single")
  if (!(linkage %in% LINKAGE)) {
    stop(sprintf("%s is an invalid linkage function.", linkage))
  }

  if (missing(max.iter)) {
    max.iter <- 500
  }

  p <- ncol(S)
  if (length(colnames(S)) == 0 || length(rownames(S)) == 0) {
    colnames(S) <- rownames(S) <- seq_len(p)
  }

  if (length(unique(colnames(S))) != p) {
    stop("Variable names are not unique.")
  }

  q.p <- c()

  labels <- colnames(S)
  dist.matrix <- matrix(0, p, p, dimnames = list(labels, labels))

  merge <- matrix(0, p - 1, 2)
  height <- vector(length = p - 1)
  order <- labels

  sub.matrices <- list(colnames(S))
  cluster.count <- p - 2
  for (iter in 1:(p - 1)) {
    current.labels <- labels %in% sub.matrices[[1]]
    hcsvd.ht <- hcsvd.ht(S = S[current.labels, current.labels],
                         S.h = S.h[current.labels, current.labels],
                         q = q,
                         linkage = linkage,
                         labels = labels[current.labels],
                         max.iter = max.iter)

    q.p <- c(q.p, hcsvd.ht$q.p)

    cluster.rows <- labels %in% hcsvd.ht$cluster1
    cluster.cols <- labels %in% hcsvd.ht$cluster2
    dist.matrix[cluster.rows, cluster.cols] <- hcsvd.ht$cluster.distance
    dist.matrix[cluster.cols, cluster.rows] <- hcsvd.ht$cluster.distance

    sub.matrices <- sub.matrices[-1]

    height[p - iter] <- hcsvd.ht$cluster.distance
    for (i in 1:2) {
      cluster <- hcsvd.ht[[paste0("cluster", i)]]
      if (length(cluster) != 1) {
        merge[p - iter, i] <- cluster.count
        cluster.count <- cluster.count - 1
        sub.matrices <- append(sub.matrices, list(cluster))
      } else {
        merge[p - iter, i] <- -which(labels == cluster)
      }
    }

    order.cluster1 <- order %in% hcsvd.ht$cluster1
    order.cluster2 <- order %in% hcsvd.ht$cluster2
    order[order.cluster1 | order.cluster2] <- c(order[order.cluster2], order[order.cluster1])

    if(verbose){
      message(sprintf("\rSplit %d out of %d (%.2f%%)           ", iter, p - 1, iter / (p - 1) * 100))
    }
  }

  ordered.height <- order(height)

  merge <- merge[ordered.height, ]
  height <- height[ordered.height]

  not.changed <- matrix(TRUE, p - 1, 2)
  for (i in seq_len(p - 1)) {
    change.idx <- which(merge == i)
    merge[merge == i & not.changed] <- which(ordered.height == i)
    not.changed[change.idx] <- FALSE
  }

  hclust <- list(merge = merge, height = height, order = match(order, labels), labels = labels, method = linkage)
  class(hclust) <- "hclust"

  u.sim <- 1 - dist.matrix
  dist.matrix <- stats::as.dist(dist.matrix)
  attr(dist.matrix, "Size") <- p

  if(verbose){
    message("\r======== FINISHED ========                    ")
    message("\n")
  }

  result <- list(hclust = hclust, dist.matrix = dist.matrix, u.sim = u.sim, q.p = q.p, linkage = linkage)
  class(result) <- "hcsvd"
  return(result)
}



#' @title Two Block Diagonal Matrix Approximation.
#'
#' @description Finds the best two block diagonal matrix approximation of a similarity matrix according to some distance (linkage) function as described in Bauer (202Xa).
#' Candidate splits are determined by the first sparse eigenvectors (sparse approximations of the first eigenvectors, i.e., vectors with many zero entries) of the similarity matrix.
#'
#' @param S A scaled \eqn{p}x\eqn{p} similarity matrix. For example, this may be a correlation matrix.
#'
#' @param linkage The linkage function to be used. This should be one of \code{"average"}, \code{"ncut"}, \code{"rcut"}, \code{"RV"} (for RV-coefficient), or \code{"single"}.
#'
#' @param q Number of sparse eigenvectors to be used. This should be either a numeric value between zero and one to indicate percentages, or \code{"Kaiser"} for as many sparse eigenvectors as
#' there are eigenvalues larger or equal to one. For a numerical value between zero and one, the number of sparse eigenvectors is determined as the corresponding share of the total number of eigenvectors.
#' E.g., \code{q = 1} (100%) uses all sparse eigenvectors and \code{q = 0.5} (50%) will use half of all sparse eigenvectors. For \code{q = 1}, identification is best (see Bauer (202Xa) for details).
#'
#' @param h.power \code{h}-th Hadamard power of \code{S}. This should be a positive integer and increases robustness of the method, as described in Bauer (202Xa).
#'
#' @param balance Minimum proportion of the smaller block when splitting into two. Must be a numeric value in \eqn{(0, 0.5]}. For example, \code{balance = 0.5} enforces an exact 50:50 split,
#' while \code{balance = 0.2} allows splits as unbalanced as 20:80, with more balanced splits such as 30:70 or 40:60 also permitted. If an exact split is not possible (e.g., \code{balance = 0.5} when
#' \eqn{p = 9}), the closest integer partition is used (e.g., 4 and 5 per block).
#'
#' @param max.iter How many iterations should be performed for computing the sparse eigenvectors.
#' Default is \code{500}.
#'
#' @details
#' The sparse eigenvectors are computed using the method of Shen and Huang (2008). The method is implemented by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} (\code{irlba}).
#' Here, we use a \code{Rcpp}/\code{RcppArmadillo} implementation based on \code{\link[irlba]{ssvd}} with slight modifications to suit our method and for faster performance.
#'
#' @return
#' A list with four components:
#' \item{B}{
#'  The best two block diagonal matrix approximation.
#' }
#' \item{BD}{
#'  The best two block diagonal matrix approximation permuted to a block diagonal shape: \eqn{BD = P B t(P)}.
#' }
#' \item{P}{
#'  The permutation matrix \eqn{P}: \eqn{BD = P B t(P)}.
#' }
#' \item{clustering}{
#'  The clustering vector as an integer vector of length \eqn{p}, which gives for each component the number \code{1} or \code{2} of the cluster/split to which it belongs.
#' }
#' \item{split}{
#'  A \code{list} containing the two splits.
#' }
#' \item{distance}{
#'  The approximation error (distance) according to the selected linkage function.
#' }
#'
#' @references \cite{Bauer, J.O. (202Xa). Divisive hierarchical clustering using block diagonal matrix approximations. Working paper.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #We give a trivial example for a block diagonal matrix perturbed by
#' #noise, for adapting clustering objectives of spectral clustering,
#' #and for balanced clustering.
#'
#' \donttest{
#' ### TOY EXAMPLE
#'
#' A <- matrix(c(2,1,1,3), 2, 2)  # 2x2 block
#' B <- matrix(c(5,4,4,6), 2, 2)  # 2x2 block
#'
#' # Create a 5x5 zero matrix and insert blocks at right positions.
#' M <- matrix(0, 4, 4)
#' M[1:2, 1:2] <- A
#' M[3:4, 3:4] <- B
#'
#' M.tilde <- M + matrix(rnorm(4^2, 0, 0.2), 4, 4)
#'
#' #Construct a similaritiy matrix with same block structure
#' S <- cov2cor(t(M.tilde) %*% M.tilde)
#' bd <- bd.approx(S)
#'
#' #Block diagonal approximation:
#' bd$B
#'
#' #We can also permute the block diagonal shape:
#' S2 <- S[c(1, 3, 2, 4), c(1, 3, 2, 4)]
#' bd2 <- bd.approx(S2)
#'
#' #bd2$B gives us again the block diagonal approximation
#' bd2$B
#'
#' #And bd2$BD gives us the block diagonal approximation permuted to
#' #block diagonal shape
#' bd2$BD
#'
#'
#' ### ADAPTING CLUSTERING OBJECTIVES
#'
#' #We will use the USArrests example (see ?hcsvd).
#' data("USArrests")
#' USArrests["Maryland", "UrbanPop"] <- 76.6
#' D <- as.matrix(dist(USArrests))
#' S <- 1 - D / max(D)
#'
#' #We compute k = 2 clusters adapting the objective of spectral clustering
#' #with the ratio cut.
#' bd.approx(S, linkage = "rcut")
#'
#'
#' ### BALANCED CLUSTERING
#'
#' #We can also enforce balanced clustering, such as two clusters of equal
#' #size (50:50). We will do this for the USArrests example from above.
#' bd.approx(S, linkage = "rcut", balance = 0.5)
#' }
#'
#'
#' @importFrom stats cov
#' @importFrom stats as.dist
#'
#' @export
bd.approx <- function(S,
                  linkage = "average",
                  q = 1,
                  h.power = 2,
                  balance,
                  max.iter
) {

  if (!is.sim.matrix(S)) {
    stop("S must be a similarity matrix.")
  }

  if (!((q == "Kaiser") | (is.numeric(q) & (q > 0) & (q <= 1)))) {
    stop(sprintf("%s is an invalid argument for q.", q))
  }

  if(h.power <= 0){
    stop("h.power must be positive.")
  }
  S.h <- S^h.power

  if(!(is.numeric(h.power) && h.power == floor(h.power))){
    warning("h.power should be an integer.")
  }

  LINKAGE <- c("average", "ncut", "rcut", "RV", "single")
  if (!(linkage %in% LINKAGE)) {
    stop(sprintf("%s is an invalid linkage function.", linkage))
  }

  if (missing(max.iter)) {
    max.iter <- 500
  }

  p <- ncol(S)
  if (length(colnames(S)) == 0 | length(rownames(S)) == 0) {
    colnames(S) <- rownames(S) <- seq_len(p)
  }

  if (length(unique(colnames(S))) != p) {
    stop("Variable names are not unique.")
  }

  p <- ncol(S)
  if (p == 2) {
    idx.cluster <- list(labels[1], labels[2])
    cluster.distance <- calc.distance(S, idx.cluster, labels, linkage)

    return(list(cluster1 = idx.cluster[1], cluster2 = idx.cluster[2], cluster.distance = cluster.distance, q.p = NA))
  }

  SVD <- eigen(S.h, symmetric = TRUE)

  if (is.numeric(q)) {
    if (q == 1){
      q <- p - 1
    } else {
      q <- min(round(q * p), p - 1)
    }
  }
  if (q == "Kaiser") {
    q <- min(sum(SVD$values >= 1) + 2, p - 1)
  }

  SVD$vectors <- SVD$vectors[, 1:q]
  labels <- colnames(S)

  dof.grid <- 1 : (p - 1)
  if (!missing(balance)) {
    if (balance < 0 | balance > 0.5) {
      stop(sprintf("%s is an invalid value for balance. Must be between 0 and 0.5.", balance))
    }
    min.grid <- max(floor(p * balance), 1)
    dof.grid <- min.grid : (p - min.grid)
  }

  distance <- -Inf
  i <- 1

  for (dof in dof.grid) {
    V <- calc_sparse_v_cpp(X = S.h, v = SVD$vectors, q = q, n = dof, maxit = max.iter)

    for (i in seq_len(q)) {
      v <- V[, i]
      idx.cluster <- list(labels[v != 0], labels[v == 0])
      distance.ht <- calc.distance(S, idx.cluster, labels, linkage)


      if (distance.ht > distance) {
        distance <- distance.ht
        cluster1 <- idx.cluster[[1]]
        cluster2 <- idx.cluster[[2]]

      }
    }
  }

  if (linkage %in% c("ncut", "rcut")) {
    distance <- -distance
  }

  B1 <- labels %in% cluster1
  B2 <- labels %in% cluster2

  B <- S
  B[B1, B2]  <- 0
  B[B2, B1]  <- 0

  BD <- matrix(0, p, p)
  colnames(BD) <- rownames(BD) <- colnames(B)
  sum.B1 <- sum(B1)
  sum.B2 <- sum(B2)
  BD[seq_len(sum.B1), seq_len(sum.B1)] <- B[B1, B1]
  BD[sum.B1 + seq_len(sum.B2), sum.B1 + seq_len(sum.B2)] <- B[B2, B2]

  order.idx <- c(which(B1), which(B2))
  P <- diag(p)[order.idx, ]

  clustering <- rep(1L, p)
  clustering[B2] <- 2L

  split <- list(B1 = cluster1, B2 = cluster2)
  class(split) <- "split"

  result <- list(BD = BD, B = B, P = P, clustering = clustering, split = split, distance = distance)
  class(result) <- "bdapprox"
  return(result)
}





#' @title Ultrametric Distance Property
#'
#' @description This function checks the ultrametric property of a distance matrix obtained by HC-SVD.
#'
#' @param hcsvd.obj An object of type hcsvd(...)
#'
#' @return
#' Returns TRUE if the ultrametric property is fulfilled. Otherwise FALSE.
#'
#' @seealso \code{\link{hcsvd}}
#'
#' @references \cite{Bauer, J.O. (202Xa). Divisive hierarchical clustering using block diagonal matrix approximations. Working paper.}
#'
#' @examples
#' \donttest{
#' #Load the correlation matrix Bechtoldt from the psych
#' #package (see ?Bechtoldt for more information).
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   data("Bechtoldt", package = "psych")
#' }
#'
#' #Compute HC-SVD (with average linkage).
#' hcsvd.obj <- hcsvd(Bechtoldt)
#'
#' #Check the ultrametric property
#' is.ultrametric(hcsvd.obj)
#' }
#'
#' @export
is.ultrametric <- function(hcsvd.obj
) {

  if (!inherits(hcsvd.obj$dist.matrix, "dist")) {
    stop("Please provide an hclust(...) object.")
  }

  dist.mat <- as.matrix(hcsvd.obj$dist.matrix)
  n <- nrow(dist.mat)

  for (i in 1:(n-2)) {
    for (j in (i+1):(n-1)) {
      for (k in (j+1):n) {
        dij <- dist.mat[i,j]
        dik <- dist.mat[i,k]
        djk <- dist.mat[j,k]

        if (dij > max(dik, djk) + 1e-10) {
          return(FALSE)
        }
        if (dik > max(dij, djk) + 1e-10) {
          return(FALSE)
        }
        if (djk > max(dij, dik) + 1e-10) {
          return(FALSE)
        }
      }
    }
  }

  return(TRUE)
}









#' @export
kmeans.svd <- function(S,
                       D,
                      linkage = "average",
                      q = 1,
                      h.power = 2,
                      balance,
                      max.iter
) {

  if (!is.sim.matrix(S)) {
    stop("S must be a similarity matrix.")
  }

  if (!((q == "Kaiser") | (is.numeric(q) & (q > 0) & (q <= 1)))) {
    stop(sprintf("%s is an invalid argument for q.", q))
  }

  if(h.power <= 0){
    stop("h.power must be positive.")
  }
  S.h <- S^h.power

  if(!(is.numeric(h.power) && h.power == floor(h.power))){
    warning("h.power should be an integer.")
  }

  LINKAGE <- c("average", "ncut", "rcut", "RV", "single")
  if (!(linkage %in% LINKAGE)) {
    stop(sprintf("%s is an invalid linkage function.", linkage))
  }

  if (missing(max.iter)) {
    max.iter <- 500
  }

  p <- ncol(S)
  if (length(colnames(S)) == 0 | length(rownames(S)) == 0) {
    colnames(S) <- rownames(S) <- seq_len(p)
  }

  if (length(unique(colnames(S))) != p) {
    stop("Variable names are not unique.")
  }

  p <- ncol(S)
  #if (p == 2) {
  #  idx.cluster <- list(labels[1], labels[2])
  #  cluster.distance <- calc.distance(S, idx.cluster, labels, linkage)
  #
  #  return(list(cluster1 = idx.cluster[1], cluster2 = idx.cluster[2], cluster.distance = cluster.distance, q.p = NA))
  #}

  SVD <- eigen(S.h, symmetric = TRUE)

  if (is.numeric(q)) {
    if (q == 1){
      q <- p - 1
    } else {
      q <- min(round(q * p), p - 1)
    }
  }
  if (q == "Kaiser") {
    q <- min(sum(SVD$values >= 1) + 2, p - 1)
  }

  SVD$vectors <- SVD$vectors[, 1:q]
  labels <- colnames(S)

  dof.grid <- 1 : (p - 1)
  if (!missing(balance)) {
    if (balance < 0 | balance > 0.5) {
      stop(sprintf("%s is an invalid value for balance. Must be between 0 and 0.5.", balance))
    }
    min.grid <- max(floor(p * balance), 1)
    dof.grid <- min.grid : (p - min.grid)
  }

  distance <- -Inf
  i <- 1

  for (dof in dof.grid) {
    V <- calc_sparse_v_cpp(X = S.h, v = SVD$vectors, q = q, n = dof, maxit = max.iter)

    for (i in seq_len(q)) {
      v <- V[, i]
      idx.cluster <- list(labels[v != 0], labels[v == 0])
      #distance.ht <- calc.distance(S, idx.cluster, labels, linkage)


      n.c1 <- length(idx.cluster[[1]])
      n.c2 <- length(idx.cluster[[2]])
      distance.ht <- 1/(2*p) * sum(D) -
        sum(D[idx.cluster[[1]], idx.cluster[[1]] ]^2) / (2*n.c1) -
        sum(D[idx.cluster[[2]], idx.cluster[[2]] ]^2) / (2*n.c2)



      if (distance.ht > distance) {
        distance <- distance.ht
        cluster1 <- idx.cluster[[1]]
        cluster2 <- idx.cluster[[2]]

      }
    }
  }

  if (linkage %in% c("ncut", "rcut")) {
    distance <- -distance
  }

  B1 <- labels %in% cluster1
  B2 <- labels %in% cluster2

  B <- S
  B[B1, B2]  <- 0
  B[B2, B1]  <- 0

  BD <- matrix(0, p, p)
  colnames(BD) <- rownames(BD) <- colnames(B)
  sum.B1 <- sum(B1)
  sum.B2 <- sum(B2)
  BD[seq_len(sum.B1), seq_len(sum.B1)] <- B[B1, B1]
  BD[sum.B1 + seq_len(sum.B2), sum.B1 + seq_len(sum.B2)] <- B[B2, B2]

  order.idx <- c(which(B1), which(B2))
  P <- diag(p)[order.idx, ]

  clustering <- rep(1L, p)
  clustering[B2] <- 2L



  ############
  n.c1 <- length(cluster1)
  n.c2 <- length(cluster2)
  d1 <- sum(D[cluster1, cluster1 ]^2) / (2*n.c1)
  d2 <- sum(D[cluster2, cluster2 ]^2) / (2*n.c2)
  print(paste0("d1 = ", d1))
  print(paste0("d2 = ", d2 ))
  distance <-  d1 + d2
  #############


  split <- list(B1 = cluster1, B2 = cluster2)
  class(split) <- "split"

  result <- list(BD = BD, B = B, P = P, clustering = clustering, split = split, distance = distance)
  class(result) <- "bdapprox"
  return(result)
}












#' @exportS3Method
print.hcsvd <- function(x, ...) {
  cat("hcsvd result\n")
  cat("Linkage:", x$linkage, "\n\n")
  cat("Available components:\n")
  cat(paste0("$", names(x)[-5], collapse = "\n"), "\n")
  invisible(x)
}





#' @exportS3Method
print.bdapprox <- function(x, ...) {
  cat("bd.approx result\n")
  cat("Cluster/Block sizes:", table(x$clustering), "\n")
  cat("Clustering:\n")
  cat(x$clustering, "\n")
  cat("Distance:", x$distance, "\n\n")
  cat("Available components:\n")
  cat(paste0("$", names(x), collapse = "\n"), "\n")
  invisible(x)
}

