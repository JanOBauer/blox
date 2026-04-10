compute.alpha <- function(R
) {

  p <- ncol(R)
  alpha <- p/(p - 1) * (1  - p/sum(abs(R)))

  return(alpha)
}



compute.beta <- function(R,
                         one,
                         two
) {

  R12 <- R[colnames(R) %in% one, colnames(R) %in% two]
  beta <- mean(abs(R12)) / mean(abs(R))

  return(beta)
}



get.vars <- function(merge,
                     node
) {
  if (node < 0) {
    return(-node)
  } else {
    return(c(get.vars(merge, merge[node, 1]), get.vars(merge, merge[node, 2])))
  }
}



is.corr.matrix <- function(R) {

  if (!is.matrix(R)) {
    return(FALSE)
  }

  if (any(abs(R) > 1)) {
    return(FALSE)
  }

  if (nrow(R) != ncol(R)) {
    return(FALSE)
  }

  if (!isSymmetric(R)) {
    return(FALSE)
  }

  if (!all(diag(R) == 1)) {
    return(FALSE)
  }

  if (any(abs(R[upper.tri(R)]) == 1)) {
    warning("Correlation matrix contains off-diagonal entries equal to 1 or -1.")
  }

  eigenvalues <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if(any(eigenvalues < 0)) {
    return(FALSE)
  }

  return(TRUE)
}



#' @title Compute Revelle's Beta for all worst split-halves using HC-SVD.
#'
#' @description Performs HC-SVD to reveal the hierarchical variable structure using average linkage as described in Bauer (202Xa). For a data matrix comprising \eqn{p} items,
#' this means that \eqn{p-1} splits are identified. The obtained structure aligns with the structure according to the worst split-half reliability and is thus used to compute
#' a hierarchy of all Revelle's beta as described in Bauer (202Xb).
#'
#' @param R A correlation matrix of dimension \eqn{p}x\eqn{p} or a data matrix of dimension \eqn{n}x\eqn{p} can be provided. If a data matrix is supplied, it must be indicated by setting
#' \code{is.corr = FALSE}, and the correlation matrix will then be calculated as \code{cor(X)}.
#'
#' @param splits An object containing the splits identified by HC-SVD. This can either be the result of \code{\link{hcsvd}} (for all splits)
#' or \code{\link{bd.approx}} (for a single split). If omitted, \code{hc.beta} will internally call \code{hcsvd(R)} and
#' compute Revelle's beta for all \eqn{p - 1} splits.
#'
#' @param n.splits Number of splits for which Revelle's beta is computed. If \code{splits} is from \code{\link{hcsvd}}, the default is all \eqn{p-1}
#'   splits. If \code{splits} is from \code{\link{bd.approx}}, only a single split is available and \code{n.splits} is set to 1.
#'
#' @param is.corr Is the supplied object a correlation matrix. Default is \code{TRUE} and this parameter must be set to \code{FALSE} if
#' a data matrix instead of a correlation matrix is supplied.
#'
#' @param verbose Print out progress as \eqn{p-1} iterations for divisive hierarchical clustering are performed.
#' Default is \code{TRUE}.
#'
#' @details
#' Supplementary details are in Bauer (202Xb).
#'
#' @return
#' A list with \code{n.splits} components. Each split is a list of four components:
#' \item{split}{
#'  The split number.
#' }
#' \item{beta}{
#'  Revelle's beta for this split.
#' }
#' \item{A}{
#'  One of the two sub-scales that has been split.
#' }
#' \item{B}{
#'  One of the two sub-scales that has been split.
#' }
#' \item{beta.alpha}{
#'  Computes the ratio of Revelle's beta and Cronbach's alpha.
#' }
#'
#' @seealso \code{\link{bd.approx}} \code{\link{hcsvd}}
#'
#' @references \cite{Bauer, J.O. (202Xa). Divisive hierarchical clustering using block diagonal matrix approximations. Working paper.}
#' @references \cite{Bauer, J.O. (202Xb). Revelle's beta: The wait is over - we can compute it! Working paper.}
#'
#' @examples
#' #We compute the worst split-half reliabilities on a correlation matrix.
#'
#' \donttest{
#' #Load the correlation matrix Bechtoldt from the psych
#' #package (see ?Bechtoldt for more information).
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   data("Bechtoldt", package = "psych")
#' }
#' R <- Bechtoldt
#'
#'
#' ### RUN HC-SVD FOR HIERARCHICAL VARIABLE CLUSTERING
#'
#' #Compute HC-SVD (with average linkage).
#' hcsvd.obj <- hcsvd(R)
#'
#' #The object of type hclust with corresponding dendrogram can be obtained
#' #directly from hcsvd(...):
#' hc.div <- hcsvd.obj$hclust
#' plot(hc.div, axes = FALSE, ylab = "", main = "Revelle's Beta Splits")
#'
#'
#' ### COMPUTE REVELLE'S BETA FOR ALL IDENTIFIED SPLITS
#'
#' #Compute Revelle's beta
#' betas <- hc.beta(R = R)
#'
#' #Alternatively, you can submit the object obtained from hcsvd(). Thus,
#' #the hiearchy needs not to be computed again using hcsvd().
#' betas <- hc.beta(R = R, splits = hcsvd.obj)
#'
#' #Visualize the splits, e.g., as
#' splits <- sapply(betas, `[[`, "split")
#' beta.values <- sapply(betas, `[[`, "beta")
#'
#' plot(splits, beta.values,
#'   type = "b",
#'   xlab = "Split",
#'   ylab = "Revelle's Beta",
#'   main = "Revelle's Beta Across Splits",
#'   pch = 19)
#'
#' #Visualize the ratio of Revelle's beta and Cronbach's alpha
#' beta.alpha <- sapply(betas, `[[`, "beta.alpha")
#' plot(splits, beta.alpha,
#'   type = "b",
#'   xlab = "Split",
#'   ylab = "Beta/Alpha",
#'   main = "Ratio of Beta and Alpha Across Splits",
#'   pch = 19)
#'
#'
#' ### COMPUTE REVELLE'S BETA FOR THE FIRST IDENTIFIED SPLIT
#'
#' #The first split can be identified using bd.approx()
#' #This is computationally faster, as only the first split
#' #is identified
#' hc.beta(R = R, splits = bd.approx(R))
#'
#' }
#'
#'
#' @importFrom stats cov
#' @importFrom stats cutree
#'
#' @export
hc.beta <- function(R,
                    splits = NULL,
                    n.splits = NULL,
                    is.corr = TRUE,
                    verbose = TRUE
) {

  if (missing(R)) {
    stop("Provide R.")
  }


  if (!is.null(splits)) {
    if (!inherits(splits$split, "split") && !inherits(splits$hclust, "hclust")) {
      stop("'splits' must either be omitted, or a result of hcsvd() or bd.approx().")
    }
  }

  if (is.corr && !is.corr.matrix(R)) {
    stop("R must be a correlation matrix. Set 'is.corr = FALSE' if you want to supply a data matrix.")
  }

  if (!is.corr) {
    X <- R
    if (anyNA(X)) {
      stop("X contains missing value indicator (NA).")
    }
    R <- stats::cor(X)
  }


  if (length(colnames(R)) == 0 || length(rownames(R)) == 0) {
    colnames(R) <- rownames(R) <- seq_len(ncol(R))
  }


  if (!is.null(n.splits) && n.splits < 0) {
    stop("n.splits cannot be negative.")
  }

  if (is.null(splits)) {
    cat("Full hierarchy for beta is being computed using hcsvd(...)\n")
    splits <- hcsvd(R, verbose = verbose)
  }


  if (inherits(splits$split, "split")) {
    if (!setequal(c(splits$split$B1, splits$split$B2), colnames(R))) {
      stop(sprintf("Provide the correlation matrix R used for bd.approx(...)."))
    }

    if (!isTRUE(n.splits == 1 || is.null(n.splits))) {
      warning("bd.approx object was supplied: only one split can be computed; using n.splits = 1")
    }
    n.splits <- 1

    beta <- compute.beta(R, splits$split$B1, splits$split$B2)
    alpha <- compute.alpha(R)

    result <- list(betas = beta, A = splits$split$B1, B = splits$split$B2, beta.alpha = beta/alpha,
                   clustering = splits$clustering)
    class(result) <- "bdapproxbeta"
    return(result)
  }


  if (inherits(splits$hclust, "hclust")) {
    if (!setequal(names(cutree(splits$hclust, k = 1)), colnames(R))) {
      stop(sprintf("Provide the correlation matrix R used for hcsvd(...)."))
    }

    if (is.null(n.splits)) {
      n.splits <- length(splits$hclust$height)
    }

    if (n.splits > length(splits$hclust$height)) {
      stop(sprintf("n.splits must be smaller or equal to p - 1."))
    }

    betas <- lapply(seq_len(n.splits), function(i) {
      prev.cluster <- splits$hclust$merge[length(splits$hclust$height) - i + 1, ]

      A <- splits$hclust$labels[get.vars(splits$hclust$merge, prev.cluster[1])]
      B <- splits$hclust$labels[get.vars(splits$hclust$merge, prev.cluster[2])]

      vars <- colnames(R) %in% c(A, B)
      R.sub <- R[vars, vars]

      beta <- compute.beta(R.sub, A, B)
      alpha <- compute.alpha(R.sub)

      list(split = i, beta = beta, A = A, B = B, beta.alpha = beta/alpha)
    })


    class(betas) <- "hcsvdbeta"
    return(betas)

  }

}





#' @exportS3Method
print.bdapproxbeta <- function(x, ...) {
  cat("Revelle's beta result\n")
  cat("Cluster/Block sizes:", table(x$clustering), "\n")
  cat("Clustering:\n")
  cat(x$clustering, "\n")
  cat("Revelle's beta:", x$betas, "\n\n")
  cat("Available components:\n")
  cat(paste0("$", names(x)[-5], collapse = "\n"), "\n")
  invisible(x)
}


#' @exportS3Method
print.hcsvdbeta <- function(x, ...) {
  cat("\rRevelle's beta has been computed for each split")
  invisible(x)
}

