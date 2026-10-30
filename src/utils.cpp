// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <utility>

using namespace arma;
using namespace Rcpp;
using namespace std;

// Soft-thresholding function
// [[Rcpp::export]]
arma::mat soft_cpp(const arma::mat& X,
                   const arma::mat& u,
                   int q,
                   int m) {

  arma::mat y = X.t() * u;
  int p = y.n_rows;

  arma::mat out(p, q, arma::fill::zeros);
  std::vector<double> tmp(p);

  for (int j = 0; j < q; ++j) {

    for (int i = 0; i < p; ++i) {
      tmp[i] = std::abs(y(i, j));
    }

    std::nth_element(tmp.begin(), tmp.begin() + m - 1, tmp.end());
    double threshold = tmp[m - 1];

    for (int i = 0; i < p; ++i) {
      double val = std::abs(y(i, j)) - threshold;

      if (val > 0.0) {
        out(i, j) = (y(i, j) >= 0.0) ? val : -val;
      }
    }
  }

  return out;
}



// Calculate Sparse Loadings
// [[Rcpp::export]]
arma::mat calc_sparse_v_cpp(const arma::mat& X, arma::mat v, int q, int m, int maxit = 500, double tol = 0.001, bool robust = false) {

  if (m == 0) {
    return v;
  }

  int iter = 0;
  double delta = std::numeric_limits<double>::infinity();
  arma::mat v_check;

  while (delta > tol && iter < maxit) {
    arma::mat v_tilde = v;
    v_check = soft_cpp(X, v_tilde, q, m);
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

  arma::rowvec col_norms = arma::sqrt(arma::sum(arma::square(v_check), 0));

  v_check.each_row() /= col_norms;

  return v_check;
}




// Calculate one HC-SVD split
// [[Rcpp::export]]
Rcpp::List hcsvd_ht_cpp(const arma::mat& S,
                        const arma::mat& S_h,
                        arma::mat V,
                        int q,
                        std::string linkage,
                        int maxit = 500,
                        double tol = 0.001) {

  int p = S.n_cols;

  arma::mat S_abs;
  arma::mat S_sq;
  arma::vec S_sq_row_sum;

  if (linkage == "average" || linkage == "single") {
    S_abs = arma::abs(S);
  }

  if (linkage == "RV") {
    S_sq = arma::square(S);
    S_sq_row_sum = arma::sum(S_sq, 1);
  }

  double distance = -std::numeric_limits<double>::infinity();
  arma::ivec cluster(p, arma::fill::value(2));

  for (int dof = 1; dof <= p - 1; ++dof) {

    arma::mat V_sparse = calc_sparse_v_cpp(
      S_h,
      V,
      q,
      dof,
      maxit,
      tol,
      false
    );

    for (int i = 0; i < q; ++i) {

      arma::vec v = V_sparse.col(i);
      arma::uvec idx_cluster1 = arma::find(v != 0);
      arma::uvec idx_cluster2 = arma::find(v == 0);

      double distance_ht = NA_REAL;

      if (linkage == "average") {

        double cut = 0.0;

        for (arma::uword j = 0; j < idx_cluster1.n_elem; ++j) {
          for (arma::uword k = 0; k < idx_cluster2.n_elem; ++k) {
            cut += S_abs(idx_cluster1[j], idx_cluster2[k]);
          }
        }

        distance_ht = 1.0 - cut /
          (
              static_cast<double>(idx_cluster1.n_elem) *
                static_cast<double>(idx_cluster2.n_elem)
          );

      } else if (linkage == "RV") {

        double S12 = 0.0;

        for (arma::uword j = 0; j < idx_cluster1.n_elem; ++j) {
          for (arma::uword k = 0; k < idx_cluster2.n_elem; ++k) {
            S12 += S_sq(idx_cluster1[j], idx_cluster2[k]);
          }
        }

        double S1 = 0.0;
        double S2 = 0.0;

        for (arma::uword j = 0; j < idx_cluster1.n_elem; ++j) {
          S1 += S_sq_row_sum[idx_cluster1[j]];
        }

        for (arma::uword j = 0; j < idx_cluster2.n_elem; ++j) {
          S2 += S_sq_row_sum[idx_cluster2[j]];
        }

        S1 -= S12;
        S2 -= S12;

        distance_ht = 1.0 - S12 / std::sqrt(S1 * S2);

      } else if (linkage == "single") {

        double S12 = -std::numeric_limits<double>::infinity();

        for (arma::uword j = 0; j < idx_cluster1.n_elem; ++j) {
          for (arma::uword k = 0; k < idx_cluster2.n_elem; ++k) {
            double val = S_abs(idx_cluster1[j], idx_cluster2[k]);

            if (val > S12) {
              S12 = val;
            }
          }
        }

        distance_ht = 1.0 - S12;

      } else {
        Rcpp::stop("Invalid linkage function.");
      }


      if (distance_ht > distance) {
        distance = distance_ht;
        cluster.fill(2);

        for (arma::uword j = 0; j < idx_cluster1.n_elem; ++j) {
          cluster[idx_cluster1[j]] = 1;
        }
      }
    }
  }


  return Rcpp::List::create(
    Rcpp::Named("cluster") = cluster,
    Rcpp::Named("distance") = distance
  );
}




// Exhaustive search over all unique two-way partitions for one objective.
//
// Object 1 is fixed in the first block, so complementary partitions are
// evaluated only once. Gray-code enumeration changes one membership at a time.
// Only the quantities required for the requested linkage function are computed.
// [[Rcpp::export]]
Rcpp::List exhaustive_bd_cpp(const arma::mat& S,
                             std::string linkage,
                             int min_block_size = 1) {

  int p = S.n_cols;

  if (S.n_rows != static_cast<arma::uword>(p)) {
    Rcpp::stop("S must be a square matrix.");
  }

  if (p < 2) {
    Rcpp::stop("S must contain at least two objects.");
  }

  if (p > 30) {
    Rcpp::stop("Exhaustive search is available only for p <= 30.");
  }

  if (min_block_size < 1 || min_block_size > p / 2) {
    Rcpp::stop("min_block_size must be between 1 and floor(p / 2).");
  }

  if (linkage != "average" && linkage != "RV" &&
      linkage != "ncut" && linkage != "rcut" &&
      linkage != "single") {
    Rcpp::stop("Invalid linkage function.");
  }


  bool use_cut = linkage == "average" ||
    linkage == "ncut" ||
    linkage == "rcut";

  bool use_ncut = linkage == "ncut";
  bool use_rv = linkage == "RV";
  bool use_single = linkage == "single";


  arma::mat W;
  arma::mat Q;

  arma::vec offdiag;
  arma::vec toA;

  arma::vec degree;
  double total_degree = 0.0;

  arma::vec offdiag_q;
  arma::vec toAq;


  if (use_cut || use_single) {
    W = arma::abs(S);
  }

  if (use_cut) {
    offdiag = arma::sum(W, 1) - W.diag();
    toA = arma::vec(p, arma::fill::zeros);
  }

  if (use_ncut) {
    degree = arma::sum(W, 1);
    total_degree = arma::accu(degree);
  }

  if (use_rv) {
    Q = arma::square(S);
    offdiag_q = arma::sum(Q, 1) - Q.diag();
    toAq = arma::vec(p, arma::fill::zeros);
  }


  // Initial partition: block 1 = {1}, block 2 = {2, ..., p}.
  std::vector<bool> inA(p, false);
  inA[0] = true;

  int nA = 1;
  int nB = p - 1;

  double cut = 0.0;
  double volA = 0.0;
  double volB = 0.0;

  double crossq = 0.0;
  double withinAq = 0.0;
  double withinBq = 0.0;


  if (use_cut) {
    for (int i = 1; i < p; ++i) {
      toA[i] = W(i, 0);
      cut += W(0, i);
    }
  }

  if (use_ncut) {
    volA = degree[0];
    volB = total_degree - volA;
  }

  if (use_rv) {
    withinAq = Q(0, 0);

    for (int i = 1; i < p; ++i) {
      toAq[i] = Q(i, 0);
      crossq += Q(0, i);

      for (int j = 1; j < p; ++j) {
        withinBq += Q(i, j);
      }
    }
  }


  struct Edge {
    double weight;
    int i;
    int j;
  };

  std::vector<Edge> edges;
  std::vector<int> edge_rank;
  std::vector<bool> crossing;
  int max_cross_rank = 0;

  if (use_single) {
    edges.reserve(static_cast<std::size_t>(p * (p - 1) / 2));

    for (int i = 0; i < p - 1; ++i) {
      for (int j = i + 1; j < p; ++j) {
        Edge edge;
        edge.weight = W(i, j);
        edge.i = i;
        edge.j = j;
        edges.push_back(edge);
      }
    }

    std::sort(
      edges.begin(),
      edges.end(),
      [](const Edge& a, const Edge& b) {
        return a.weight > b.weight;
      }
    );

    edge_rank.assign(static_cast<std::size_t>(p * p), -1);
    crossing.assign(edges.size(), false);

    for (std::size_t r = 0; r < edges.size(); ++r) {
      int i = edges[r].i;
      int j = edges[r].j;

      edge_rank[static_cast<std::size_t>(i * p + j)] = static_cast<int>(r);
      edge_rank[static_cast<std::size_t>(j * p + i)] = static_cast<int>(r);

      if (inA[i] != inA[j]) {
        crossing[r] = true;
      }
    }

    max_cross_rank = static_cast<int>(edges.size());
    for (std::size_t r = 0; r < crossing.size(); ++r) {
      if (crossing[r]) {
        max_cross_rank = static_cast<int>(r);
        break;
      }
    }
  }


  double best_objective = std::numeric_limits<double>::infinity();
  arma::ivec best_cluster(p, arma::fill::zeros);
  double n_candidates = 0.0;


  auto save_cluster = [&]() {
    for (int i = 0; i < p; ++i) {
      best_cluster[i] = inA[i] ? 1 : 2;
    }
  };


  auto evaluate = [&]() {

    if (nA < min_block_size || nB < min_block_size) {
      return;
    }

    ++n_candidates;

    double objective = std::numeric_limits<double>::infinity();

    if (linkage == "average") {
      objective = cut /
        (static_cast<double>(nA) * static_cast<double>(nB));
    }

    if (linkage == "ncut") {
      objective = cut / volA + cut / volB;
    }

    if (linkage == "rcut") {
      objective = cut / static_cast<double>(nA) +
        cut / static_cast<double>(nB);
    }

    if (linkage == "RV") {
      objective = crossq / std::sqrt(withinAq * withinBq);
    }

    if (linkage == "single") {
      objective = edges[static_cast<std::size_t>(max_cross_rank)].weight;
    }

    if (objective < best_objective) {
      best_objective = objective;
      save_cluster();
    }
  };


  evaluate();

  std::uint64_t n_states = (1ULL << (p - 1));
  std::uint64_t previous_gray = 0ULL;


  for (std::uint64_t state = 1ULL; state < n_states; ++state) {

    std::uint64_t gray = state ^ (state >> 1ULL);
    std::uint64_t changed = gray ^ previous_gray;
    int bit = __builtin_ctzll(changed);
    int v = bit + 1;

    bool move_to_A = ((gray >> bit) & 1ULL) != 0ULL;
    bool old_in_A = inA[v];


    if (use_single) {

      int old_max_cross_rank = max_cross_rank;
      int added_cross_rank = static_cast<int>(edges.size());
      bool removed_old_max = false;

      for (int u = 0; u < p; ++u) {

        if (u == v) {
          continue;
        }

        bool old_cross = inA[u] != old_in_A;
        bool new_cross = inA[u] != move_to_A;

        if (old_cross != new_cross) {

          int rank = edge_rank[static_cast<std::size_t>(v * p + u)];
          crossing[static_cast<std::size_t>(rank)] = new_cross;

          if (new_cross && rank < added_cross_rank) {
            added_cross_rank = rank;
          }

          if (!new_cross && rank == old_max_cross_rank) {
            removed_old_max = true;
          }
        }
      }

      if (added_cross_rank < old_max_cross_rank) {
        max_cross_rank = added_cross_rank;
      } else if (removed_old_max) {

        max_cross_rank = old_max_cross_rank + 1;

        while (max_cross_rank < static_cast<int>(crossing.size()) &&
               !crossing[static_cast<std::size_t>(max_cross_rank)]) {
               ++max_cross_rank;
        }
      }
    }


    if (move_to_A) {

      if (use_cut) {

        double sumA = toA[v];
        double sumB = offdiag[v] - sumA;

        cut += sumB - sumA;

        for (int u = 0; u < p; ++u) {
          if (u != v) {
            toA[u] += W(u, v);
          }
        }
      }

      if (use_ncut) {
        volA += degree[v];
        volB -= degree[v];
      }

      if (use_rv) {

        double sumAq = toAq[v];
        double sumBq = offdiag_q[v] - sumAq;

        crossq += sumBq - sumAq;
        withinAq += Q(v, v) + 2.0 * sumAq;
        withinBq -= Q(v, v) + 2.0 * sumBq;

        for (int u = 0; u < p; ++u) {
          if (u != v) {
            toAq[u] += Q(u, v);
          }
        }
      }

      ++nA;
      --nB;
      inA[v] = true;

    } else {

      if (use_cut) {

        double sumA = toA[v];
        double sumB = offdiag[v] - sumA;

        cut += sumA - sumB;

        for (int u = 0; u < p; ++u) {
          if (u != v) {
            toA[u] -= W(u, v);
          }
        }
      }

      if (use_ncut) {
        volA -= degree[v];
        volB += degree[v];
      }

      if (use_rv) {

        double sumAq = toAq[v];
        double sumBq = offdiag_q[v] - sumAq;

        crossq += sumAq - sumBq;
        withinAq -= Q(v, v) + 2.0 * sumAq;
        withinBq += Q(v, v) + 2.0 * sumBq;

        for (int u = 0; u < p; ++u) {
          if (u != v) {
            toAq[u] -= Q(u, v);
          }
        }
      }

      --nA;
      ++nB;
      inA[v] = false;
    }


    evaluate();
    previous_gray = gray;
  }


  return Rcpp::List::create(
    Rcpp::Named("cluster") = best_cluster,
    Rcpp::Named("objective") = best_objective,
    Rcpp::Named("n.candidates") = n_candidates
  );
}
