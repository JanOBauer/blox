// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// Soft-thresholding function
// [[Rcpp::export]]
arma::mat soft_cpp(const arma::mat& X, const arma::mat& u, int q, int n) {
  arma::mat y = X.t() * u;
  arma::mat a = abs(y);
  arma::mat z = sort(a, "ascend", 0);

  arma::rowvec lambda(q);
  for (int j = 0; j < q; j++) {
    lambda(j) = z(n - 1, j);
  }

  a.each_row() -= lambda;
  a.transform([](double val) { return std::max(val, 0.0); });

  return sign(y) % a;
}



// Calculate Sparse Loadings
// [[Rcpp::export]]
arma::mat calc_sparse_v_cpp(const arma::mat& X, arma::mat v, double q, int n, int maxit = 500, double tol = 0.001, bool robust = false) {

  if (n == 0) {
    return v;
  }

  int iter = 0;
  double delta = std::numeric_limits<double>::infinity();
  arma::mat v_check;

  while (delta > tol && iter < maxit) {
    arma::mat v_tilde = v;
    v_check = soft_cpp(X, v_tilde, q, n);
    arma::mat xsv = X * v_check;

    arma::mat R;
    arma::qr_econ(v, R, xsv);

    for (int j = 0; j < q; j++) {
      arma::uvec nonzero_idx = find(xsv.col(j) != 0, 1);
      if (!nonzero_idx.is_empty()) {
        v.col(j) *= (sign(xsv(nonzero_idx(0), j)) * sign(v(nonzero_idx(0), j)));

      }
    }

    arma::vec diag_vals = abs(diagvec(v_tilde.t() * v));
    delta = max(1 - diag_vals);

    ++iter;
  }

  if (iter >= maxit) {
    Rcpp::Rcerr << "Warning: Maximum number of iterations reached before convergence: solution may not be optimal. Consider increasing 'maxit'.\n";
  }

  return v_check * arma::diagmat(1.0 / arma::sqrt(arma::sum(arma::square(v_check), 0)));
}
