// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rcpp, RcppParallel)]]
#include <Rcpp.h>
using namespace Rcpp;

// Build a column-major bitset matrix for presence:
// - hits_by_peptide: List of IntegerVector (1-based subject indices with presence), length = m peptides
// - N: number of subjects (rows)
// Returns: list(data=<raw bytes>, m=<int>, n_words=<int>)
//   data layout: columns are contiguous; each column has n_words uint64 words (little-endian in raw)

// [[Rcpp::export]]
List build_bitset_unpaired(const List& hits_by_peptide, const int N) {
  const int m = hits_by_peptide.size();
  if (N <= 0 || m <= 0) stop("Invalid dimensions");

  const int n_words = (N + 63) / 64;
  const size_t total_words = static_cast<size_t>(m) * static_cast<size_t>(n_words);

  std::vector<uint64_t> words(total_words, 0ULL);

  for (int j = 0; j < m; ++j) {
    SEXP col = hits_by_peptide[j];
    if (TYPEOF(col) == INTSXP) {
      IntegerVector idx(col);
      uint64_t* col_ptr = &words[ static_cast<size_t>(j) * static_cast<size_t>(n_words) ];
      for (int k = 0; k < idx.size(); ++k) {
        int s = idx[k];
        if (s <= 0) continue;              // skip invalid
        int z = s - 1;                     // 0-based
        int w = z >> 6;                    // /64
        int b = z & 63;                    // %64
        col_ptr[w] |= (1ULL << b);
      }
    } else if (!Rf_isNull(col)) {
      stop("hits_by_peptide[%d] must be integer or NULL", j + 1);
    }
  }

  // Pack into a RawVector (bytes)
  RawVector out(static_cast<R_xlen_t>(total_words * 8ULL));
  uint8_t* p = reinterpret_cast<uint8_t*>(RAW(out));
  std::memcpy(p, words.data(), total_words * sizeof(uint64_t));

  return List::create(
    _["data"] = out,
    _["m"] = m,
    _["n_words"] = n_words
  );
}
