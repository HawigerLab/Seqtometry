// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

//' Helper function for performing a weighted Kolmogorov-Smirnov-like procedure.
//' @param gex numeric: normalized gene expression values
//' @param gss list: indices of genes in each gene set
//' @param mus numeric: means of all genes
//' @param sds numeric: standard deviations of all genes
//' @returns Modified Kuiper statistic (sum of minimal and maximal deviations during running sum)
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector wks(
    const arma::vec& gex, const Rcpp::List& gss, 
    const arma::vec& mus, const arma::vec& sds) {
  // Sorting permutation of Z score transformed values
  arma::uvec rank_idx = arma::sort_index((gex - mus) / sds);

  // Compute centered ranks of genes
  int center = gex.size() / 2;
  int start_index = gex.size() % 2 == 0 ? center + 1 : ++center; // Used later in running sum procedure
  std::vector<int> gene_ranks(gex.size());
  for (int x = 0; x < rank_idx.size(); ++x)
    gene_ranks[rank_idx[x]] = x + 1 - center;
  
  Rcpp::NumericVector ret(gss.size()); // Result vector of scores to be returned
  // Calculate score for each signature and place it in the result vector
  for (int x = 0; x < gss.size(); ++x) {
    const Rcpp::IntegerVector& sig_idx = gss[x];
    if (sig_idx.size() == 0) { // Empty signature
      ret[x] = -1;
      continue;
    }

    // Number of genes outside of the signature
    double non_sig_size = gene_ranks.size() - sig_idx.size(); 
    if (non_sig_size == 0) { // Signature contains all genes in the cell
      ret[x] = 1;
      continue;
    }
    double dec_norm = 1 / non_sig_size; // Normalizing factor for decrements in the running sum procedure
    
    std::vector<int> sig_ranks;         // Centered ranks of signature genes
    sig_ranks.reserve(sig_idx.size());
    double inc_norm = 0;                // Normalizing factor for increments in the running sum procedure
    for (const int& i : sig_idx) {
      int k = gene_ranks[i];
      inc_norm += std::abs(k);
      sig_ranks.push_back(k);
    }

    // Sort ranks to place them in descending order for the strided running sum
    std::sort(sig_ranks.begin(), sig_ranks.end(), std::greater<>());
    int stride;          // Number of genes in between two successive members of the gene set
    int i = start_index; // Current index in strided running sum
    double mn = 0;       // Minimum value of running sum
    double mx = 0;       // Maximum value of running sum
    double rs = 0;       // Current value of running sum
    for (const int& k : sig_ranks) {
      // Update stride and index for the current gene
      stride = i - k - 1;
      i = k;

      // Update running sum, keeping track of the associated min and max
      rs -= stride * dec_norm;
      if (rs < mn) mn = rs;
      rs += std::abs(k) / inc_norm;
      if (rs > mx) mx = rs;
    }
    // Last decrement (after last member of signature)
    stride = i + center - 1;
    rs -= stride * dec_norm;
    if (rs < mn) mn = rs;

    // Modified Kuiper statistic instead of classical Kolmogorov-Smirnov statistic
    ret[x] = mn + mx;
  }

  return ret;
}
