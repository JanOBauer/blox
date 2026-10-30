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



exhaustive.split <- function(S,
                             linkage,
                             labels,
                             min.block.size = 1L
) {

  exhaustive <- exhaustive_bd_cpp(
    S = S,
    linkage = linkage,
    min_block_size = as.integer(min.block.size)
  )

  clustering <- as.integer(exhaustive$cluster)

  cluster1 <- labels[clustering == 1]
  cluster2 <- labels[clustering == 2]

  cluster.distance <- exhaustive$objective
  if (linkage %in% c("average", "RV", "single")) {
    cluster.distance <- 1 - cluster.distance
  }

  return(list(
    cluster1 = cluster1,
    cluster2 = cluster2,
    cluster.distance = cluster.distance,
    n.candidates = exhaustive$n.candidates
  ))
}



hcsvd.ht <- function(S,
                     S.h,
                     q,
                     linkage,
                     labels,
                     max.iter,
                     search,
                     exhaustive.max.p
) {

  p <- ncol(S)
  if (p == 2) {
    idx.cluster <- list(labels[1], labels[2])
    cluster.distance <- calc.distance(S, idx.cluster, labels, linkage)

    return(list(cluster1 = idx.cluster[1], cluster2 = idx.cluster[2],
                cluster.distance = cluster.distance, q.p = NA,
                search = "unique", n.candidates = 1))
  }


  search.used <- search
  if (search == "auto") {
    if (p <= exhaustive.max.p) {
      search.used <- "exhaustive"
    } else {
      search.used <- "hcsvd"
    }
  }


  if (search.used == "exhaustive") {
    exact <- exhaustive.split(
      S = S,
      linkage = linkage,
      labels = labels,
      min.block.size = 1L
    )

    return(list(cluster1 = exact$cluster1, cluster2 = exact$cluster2,
                cluster.distance = exact$cluster.distance, q.p = NA,
                search = "exhaustive", n.candidates = exact$n.candidates))
  }


  SVD <- eigen(S.h, symmetric = TRUE)

  if (is.numeric(q)) {
    if (q == 1){
      q <- p - 1
    } else {
      q <- max(1, min(round(q * p), p - 1))
    }
  }
  if (q == "Kaiser") {
    q <- min(sum(SVD$values >= 1) + 2, p - 1)
  }

  SVD$vectors <- SVD$vectors[, 1:q, drop = FALSE]

  hcsvd.ht.cpp <- hcsvd_ht_cpp(
    S = S,
    S_h = S.h,
    V = SVD$vectors,
    q = q,
    linkage = linkage,
    maxit = max.iter
  )

  cluster1 <- labels[hcsvd.ht.cpp$cluster == 1]
  cluster2 <- labels[hcsvd.ht.cpp$cluster == 2]
  cluster.distance <- hcsvd.ht.cpp$distance

  return(list(cluster1 = cluster1, cluster2 = cluster2,
              cluster.distance = cluster.distance, q.p = q / p,
              search = "hcsvd", n.candidates = q * (p - 1)))
}



is.sim.matrix <- function(S) {

  if (!is.matrix(S)) {
    stop("The provided similarity matrix is not a matrix.")
  }

  tol <- sqrt(.Machine$double.eps)

  if (any(abs(diag(S) - 1) > tol)) {
    stop("The provided similarity matrix does not have ones on the diagonal.")
  }

  if (any(abs(S) > 1 + tol)) {
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
  if(any(eigenvalues < -tol)) {
    stop("The provided similarity matrix is not positive semi-definite.")
  }

  return(TRUE)
}



#' @title Hierarchical Clustering Using Singular Vectors (HC-SVD).
#'
#' @description Performs HC-SVD to reveal the hierarchical structure as described in Bauer (202X). This divisive approach iteratively splits each cluster into two subclusters.
#' By default, candidate splits are determined by sparse approximations of the first eigenvectors of the similarity matrix. For small clusters, an exhaustive search over all unique two-way partitions can alternatively be used.
#' The selected split is the one that yields the best block-diagonal approximation of the similarity matrix according to a specified linkage function. The procedure continues until each object is assigned to its own cluster.
#'
#' @param S A scaled \eqn{p}x\eqn{p} similarity matrix. For example, this may be a correlation matrix.
#'
#' @param linkage The linkage function to be used. This should be one of \code{"average"}, \code{"single"}, or
#' \code{"RV"} (for RV-coefficient). Note that the RV-coefficient might not yield an ultrametric distance.
#'
#' @param q Number of sparse eigenvectors to be used whenever the HC-SVD search is applied. This should be either a numeric value between zero and one to indicate percentages, or \code{"Kaiser"}.
#' For a numerical value between zero and one, the number of sparse eigenvectors is determined as the corresponding share of the total number of eigenvectors, capped at \eqn{p-1}. E.g., \code{q = 1} uses all \eqn{p-1} available sparse eigenvectors.
#' The Kaiser rule uses the number of eigenvalues of the Hadamard-transformed similarity matrix that are at least one, plus two additional eigenvectors, capped at \eqn{p-1}.
#'
#' @param h.power \code{h}-th Hadamard power of \code{S} used whenever sparse-eigenvector candidate generation is applied. This should be a positive integer. Larger powers shrink similarities with magnitude smaller than one towards zero and can increase the contrast between within- and between-block similarities, although their effect depends on the similarity structure and clustering objective.
#'
#' @param max.iter How many iterations should be performed for computing the sparse eigenvectors when the HC-SVD search is used.
#' Default is \code{500}.
#'
#' @param verbose Print out progress as \eqn{p-1} iterations for divisive hierarchical clustering are performed.
#' Default is \code{TRUE}.
#'
#' @param search Search strategy. \code{"hcsvd"} uses sparse-eigenvector candidate generation, \code{"exhaustive"} evaluates all unique two-way partitions for the selected linkage function only, and \code{"auto"} uses exhaustive search whenever the current cluster contains at most \code{exhaustive.max.p} objects and HC-SVD otherwise.
#'
#' @param exhaustive.max.p Largest cluster size for which exhaustive search is used when \code{search = "auto"}. Default is \code{20}. The exhaustive C++ implementation is restricted to at most 30 objects.
#'
#' @details
#' The sparse loadings are computed using the method proposed by Shen & Huang (2008). The corresponding implementation is written in \code{Rcpp}/\code{RcppArmadillo}
#' for computational efficiency and is based on the \code{R} implementation by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} \{\pkg{irlba}\}.
#' However, the implementation has been adapted to better align with the scope of the \pkg{bdsvd} package which is the base for the \pkg{blox} package.
#'
#' Supplementary details are in \code{\link[blox]{hc.beta}} and in Bauer (2026).
#'
#' @return
#' A list containing:
#' \item{hclust}{
#'  The clustering structure identified by HC-SVD as an object of type \code{hclust}. For \code{linkage = "RV"}, dendrogram heights are
#'  monotonically adjusted if required because the RV-coefficient need not yield an ultrametric distance.
#' }
#' \item{dist.matrix}{
#'  The distance matrix of the HC-SVD structure as an object of class \code{dist}. This matrix need not be ultrametric for \code{linkage = "RV"}.
#' }
#' \item{u.sim}{
#'  The corresponding similarity matrix of \eqn{S}, calculated as \code{1-dist.matrix}.
#' }
#' \item{q.p}{
#'  A vector of length \eqn{p-1} containing the ratio \eqn{q_i/p_i} of the number of sparse eigenvectors used relative to the total number of eigenvectors for the split of each cluster.
#'  Since HC-SVD uses at most \eqn{p_i-1} sparse eigenvectors, this ratio is always smaller than one. The ratio is set to \code{NA} when exhaustive search is used or when
#'  the cluster contains only two objects, since no sparse eigenvectors are computed in these cases.
#' }
#' \item{search}{
#'  A vector indicating whether \code{"hcsvd"}, \code{"exhaustive"}, or the unique two-object split was used at each recursive step.
#' }
#'
#' @seealso \code{\link[bdsvd]{bdsvd}} \{\pkg{bdsvd}\}
#'
#' @references \cite{Bauer, J.O. (2026). Revelle’s beta: The wait is over—computation becomes possible. Psychometrika.}
#' @references \cite{Bauer, J.O. (202X). Divisive hierarchical clustering using block diagonal matrix approximations. Working paper.}
#' @references \cite{Shen, H. and Huang, J.Z. (2008). Sparse principal component analysis via regularized low rank matrix approximation, J. Multivar. Anal. 99, 1015–1034.}
#'
#' @examples
#' #We give one example for variable clustering directly on a correlation matrix,
#' #and we replicate the USArrest example in Bauer (202X) for observation clustering.
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
                  verbose = TRUE,
                  search = "hcsvd",
                  exhaustive.max.p = 20
) {

  if (missing(S)) {
    stop("Provide S.")
  }

  if (!is.sim.matrix(S)) {
    stop("S must be a similarity matrix.")
  }


  if (!is.numeric(exhaustive.max.p) || length(exhaustive.max.p) != 1L ||
      !is.finite(exhaustive.max.p) || exhaustive.max.p < 2 ||
      exhaustive.max.p > 30 || exhaustive.max.p != floor(exhaustive.max.p)) {
    stop("exhaustive.max.p must be an integer between 2 and 30.")
  }

  p <- ncol(S)


  SEARCH <- c("hcsvd", "exhaustive", "auto")
  if (!(search %in% SEARCH)) {
    stop(sprintf("%s is an invalid search parameter", search))
  }

  uses.hcsvd <- (search == "hcsvd" || (search == "auto" && p > exhaustive.max.p))
  if (uses.hcsvd) {
    if (!((q == "Kaiser") | (is.numeric(q) & (q > 0) & (q <= 1)))) {
      stop(sprintf("%s is an invalid argument for q.", q))
    }

    if (length(h.power) != 1L || !is.numeric(h.power) || is.na(h.power) || !is.finite(h.power) ||
        h.power < 1 || h.power != floor(h.power)) {
      stop("h.power should be a positive integer.")
    }
  }

  if (uses.hcsvd) {
    S.h <- S^h.power
  } else {
    S.h <- S
  }


  LINKAGE <- c("average", "RV", "single")
  if (!(linkage %in% LINKAGE)) {
    stop(sprintf("%s is an invalid linkage function.", linkage))
  }

  if (missing(max.iter)) {
    max.iter <- 500
  }

  if (search == "exhaustive" && p > 30) {
    stop("Exhaustive search is available only for p <= 30.")
  }

  if (length(colnames(S)) == 0 || length(rownames(S)) == 0) {
    colnames(S) <- rownames(S) <- seq_len(p)
  }

  if (length(unique(colnames(S))) != p) {
    stop("Variable names are not unique.")
  }

  q.p <- c()
  search.p <- c()

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
                         max.iter = max.iter,
                         search = search,
                         exhaustive.max.p = exhaustive.max.p)

    q.p <- c(q.p, hcsvd.ht$q.p)
    search.p <- c(search.p, hcsvd.ht$search)

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

  if (linkage == "RV") {
    # RV distances need not be ultrametric. Adjust only the dendrogram heights
    # when a parent would otherwise lie below one of its children.
    for (i in seq_len(p - 1)) {
      previous.clusters <- merge[i, merge[i, ] > 0]

      if (length(previous.clusters) > 0) {
        height[i] <- max(height[i], height[previous.clusters])
      }
    }
  }

  ordered.height <- order(height, seq_along(height))

  merge <- merge[ordered.height, , drop = FALSE]
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

  result <- list(hclust = hclust, dist.matrix = dist.matrix, u.sim = u.sim,
                 q.p = q.p, linkage = linkage, search = search.p)
  class(result) <- "hcsvd"
  return(result)
}



#' @title Two Block Diagonal Matrix Approximation.
#'
#' @description Finds the best two block diagonal matrix approximation of a similarity matrix according to some distance (linkage) function as described in Bauer (202X).
#' By default, candidate splits are determined by sparse approximations of the first eigenvectors of the similarity matrix. Alternatively, all unique two-way partitions can be evaluated by exhaustive search.
#'
#' @param S A scaled \eqn{p}x\eqn{p} similarity matrix. For example, this may be a correlation matrix.
#'
#' @param linkage The linkage function to be used. This should be one of \code{"average"}, \code{"ncut"}, \code{"rcut"}, \code{"RV"} (for RV-coefficient), or \code{"single"}.
#'
#' @param q Number of sparse eigenvectors to be used whenever the HC-SVD search is applied. This should be either a numeric value between zero and one to indicate percentages, or \code{"Kaiser"}.
#' For a numerical value between zero and one, the number of sparse eigenvectors is determined as the corresponding share of the total number of eigenvectors, capped at \eqn{p-1}. E.g., \code{q = 1} uses all \eqn{p-1} available sparse eigenvectors.
#' The Kaiser rule uses the number of eigenvalues of the Hadamard-transformed similarity matrix that are at least one, plus two additional eigenvectors, capped at \eqn{p-1}.
#'
#' @param h.power \code{h}-th Hadamard power of \code{S} used whenever sparse-eigenvector candidate generation is applied. This should be a positive integer. Larger powers shrink similarities with magnitude smaller than one towards zero and can increase the contrast between within- and between-block similarities, although their effect depends on the similarity structure and clustering objective.
#'
#' @param balance Minimum proportion of the smaller block when splitting into two. Must be a numeric value in \eqn{(0, 0.5]}. For example, \code{balance = 0.5} enforces an exact 50:50 split,
#' while \code{balance = 0.2} allows splits as unbalanced as 20:80, with more balanced splits such as 30:70 or 40:60 also permitted. If an exact split is not possible (e.g., \code{balance = 0.5} when
#' \eqn{p = 9}), the closest integer partition is used (e.g., 4 and 5 per block).
#'
#' @param max.iter How many iterations should be performed for computing the sparse eigenvectors when the HC-SVD search is used.
#' Default is \code{500}.
#'
#' @param search Search strategy. \code{"hcsvd"} uses sparse-eigenvector candidate generation, \code{"exhaustive"} evaluates all unique two-way partitions for the selected linkage function only, and \code{"auto"} uses exhaustive search when \eqn{p} is at most \code{exhaustive.max.p} and HC-SVD otherwise.
#'
#' @param exhaustive.max.p Largest problem size for which exhaustive search is used when \code{search = "auto"}. Default is \code{20}. The exhaustive C++ implementation is restricted to at most 30 objects.
#'
#' @details
#' The sparse eigenvectors are computed using the method of Shen and Huang (2008). The method is implemented by Baglama, Reichel, and Lewis in \code{\link[irlba]{ssvd}} (\code{irlba}).
#' Here, we use a \code{Rcpp}/\code{RcppArmadillo} implementation based on \code{\link[irlba]{ssvd}} with slight modifications to suit our method and for faster performance.
#'
#' @return
#' A list containing:
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
#' \item{search}{
#'  The search strategy used for the split.
#' }
#' \item{n.candidates}{
#'  The number of candidate evaluations performed.
#' }
#'
#' @references \cite{Bauer, J.O. (202X). Divisive hierarchical clustering using block diagonal matrix approximations. Working paper.}
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
                      max.iter,
                      search = "hcsvd",
                      exhaustive.max.p = 20
) {

  if (!is.sim.matrix(S)) {
    stop("S must be a similarity matrix.")
  }

  SEARCH <- c("hcsvd", "exhaustive", "auto")
  if (!(search %in% SEARCH)) {
    stop(sprintf("%s is an invalid search parameter", search))
  }

  if (!is.numeric(exhaustive.max.p) || length(exhaustive.max.p) != 1L ||
      !is.finite(exhaustive.max.p) || exhaustive.max.p < 2 ||
      exhaustive.max.p > 30 || exhaustive.max.p != floor(exhaustive.max.p)) {
    stop("exhaustive.max.p must be an integer between 2 and 30.")
  }

  p <- ncol(S)
  search.used <- search
  if (search == "auto") {
    if (p <= exhaustive.max.p) {
      search.used <- "exhaustive"
    } else {
      search.used <- "hcsvd"
    }
  }

  if (search.used == "exhaustive" && p > 30) {
    stop("Exhaustive search is available only for p <= 30.")
  }

  if (search.used == "hcsvd") {
    if (!((q == "Kaiser") | (is.numeric(q) & (q > 0) & (q <= 1)))) {
      stop(sprintf("%s is an invalid argument for q.", q))
    }

    if (length(h.power) != 1L || !is.numeric(h.power) || is.na(h.power) || !is.finite(h.power) ||
        h.power < 1 || h.power != floor(h.power)) {
      stop("h.power should be a positive integer.")
    }
  }

  if (search.used == "hcsvd") {
    S.h <- S^h.power
  } else {
    S.h <- S
  }


  LINKAGE <- c("average", "ncut", "rcut", "RV", "single")
  if (!(linkage %in% LINKAGE)) {
    stop(sprintf("%s is an invalid linkage function.", linkage))
  }

  if (missing(max.iter)) {
    max.iter <- 500
  }

  if (length(colnames(S)) == 0 | length(rownames(S)) == 0) {
    colnames(S) <- rownames(S) <- seq_len(p)
  }

  if (length(unique(colnames(S))) != p) {
    stop("Variable names are not unique.")
  }

  labels <- colnames(S)

  min.grid <- 1L
  if (!missing(balance)) {
    if (balance < 0 | balance > 0.5) {
      stop(sprintf("%s is an invalid value for balance. Must be between 0 and 0.5.", balance))
    }
    min.grid <- max(floor(p * balance), 1)
  }

  if (p == 2) {
    cluster1 <- labels[1]
    cluster2 <- labels[2]
    idx.cluster <- list(cluster1, cluster2)
    distance <- calc.distance(S, idx.cluster, labels, linkage)

    if (linkage %in% c("ncut", "rcut")) {
      distance <- -distance
    }

    B <- S
    B[1, 2] <- 0
    B[2, 1] <- 0

    BD <- B
    P <- diag(p)
    rownames(P) <- colnames(P) <- labels

    clustering <- c(1L, 2L)
    names(clustering) <- labels

    split <- list(B1 = cluster1, B2 = cluster2)
    class(split) <- "split"

    result <- list(BD = BD, B = B, P = P, clustering = clustering,
                   split = split, distance = distance,
                   search = "unique", n.candidates = 1)
    class(result) <- "bdapprox"
    return(result)
  }


  if (search.used == "exhaustive") {

    exact <- exhaustive.split(
      S = S,
      linkage = linkage,
      labels = labels,
      min.block.size = min.grid
    )

    cluster1 <- exact$cluster1
    cluster2 <- exact$cluster2
    distance <- exact$cluster.distance
    n.candidates <- exact$n.candidates

  } else {

    SVD <- eigen(S.h, symmetric = TRUE)

    if (is.numeric(q)) {
      if (q == 1){
        q <- p - 1
      } else {
        q <- max(1, min(round(q * p), p - 1))
      }
    }
    if (q == "Kaiser") {
      q <- min(sum(SVD$values >= 1) + 2, p - 1)
    }

    SVD$vectors <- SVD$vectors[, 1:q, drop = FALSE]

    dof.grid <- min.grid : (p - min.grid)

    distance <- -Inf
    i <- 1

    for (dof in dof.grid) {
      V <- calc_sparse_v_cpp(X = S.h, v = SVD$vectors, q = q, m = dof, maxit = max.iter)

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

    n.candidates <- length(dof.grid) * q
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
  names(clustering) <- labels

  split <- list(B1 = cluster1, B2 = cluster2)
  class(split) <- "split"

  result <- list(BD = BD, B = B, P = P, clustering = clustering,
                 split = split, distance = distance,
                 search = search.used, n.candidates = n.candidates)
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
#' @references \cite{Bauer, J.O. (202X). Divisive hierarchical clustering using block diagonal matrix approximations. Working paper.}
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

  if (n <= 2) {
    return(TRUE)
  }

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

