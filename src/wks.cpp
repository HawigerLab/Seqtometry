#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Helper function for performing a weighted Kolmogorov-Smirnov like procedure.
//' @param gene_ranks numeric: centered gene expression ranks
//' @param sig_idx numeric: indices of genes belonging to the signature (gene set) being scored
//' @param aux numeric: two element vector containing 
//'   0. Centering value that was used for ranks
//'   1. Start index (after centering) for the running sum procedure
//' @returns Modified Kuiper statistic (sum of minimal and maximal deviations during running sum)
// [[Rcpp::export]]
double wks(const NumericVector& gene_ranks, const NumericVector& sig_idx, const NumericVector& aux) {
  if (sig_idx.size() == 0) return -1;  // Early return due to empty signature

  double non_sig_size = gene_ranks.size() - sig_idx.size(); // Number of genes outside of the signature
  if (non_sig_size == 0) return 1;     // Early return due to signature containing all genes in a cell
  double dec_norm = 1 / non_sig_size;  // Normalizing factor for decrements in the running sum procedure

  std::vector<double> sig_ranks;       // Ranks of signature genes
  sig_ranks.reserve(sig_idx.size());
  double inc_norm = 0;                 // Normalizing factor for increments in the running sum procedure
  for (const auto& i : sig_idx) {
    auto k = gene_ranks[i];
    inc_norm += std::abs(k);
    sig_ranks.push_back(k);
  }
  
  // Sort ranks to place them in descending order for the strided running sum
  std::sort(sig_ranks.begin(), sig_ranks.end(), std::greater<>());
  double stride;     // Number of genes in between two successive members of the gene set
  double i = aux[1]; // Current index in strided running sum
  double mn = 0;     // Minimum value of running sum
  double mx = 0;     // Maximum value of running sum
  double rs = 0;     // Current value of running sum
  for (const auto& k : sig_ranks) {
    stride = i - k - 1;
    i = k;
    rs -= stride * dec_norm;
    if (rs < mn) mn = rs;
    rs += std::abs(k) / inc_norm;
    if (rs > mx) mx = rs;
  }
  // Last decrement (after last member of signature)
  stride = i + aux[0] - 1;
  rs -= stride * dec_norm;
  if (rs < mn) mn = rs;
  
  // Modified Kuiper statistic instead of classical Kolmogorov-Smirnov statistic
  return mn + mx;
}
