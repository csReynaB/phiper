// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <random>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// --------- helpers ------------------------------------------------------------
inline uint32_t pc64(uint64_t x) { return (uint32_t)__builtin_popcountll(x); }

// Column-major bitset: m columns × n_words 64-bit words per column
struct BitsetMat {
  const uint64_t* data;
  int m;
  int n_words;
  inline const uint64_t* col(int j) const { return data + (size_t)j * (size_t)n_words; }
};

// --------- build utilities ----------------------------------------------------

// Build bitset (unpaired): hits_by_peptide is List<IntegerVector> of 1-based subject ids (size N)
static std::vector<uint64_t> build_bitset_unpaired(const List& hits_by_peptide, int N, int &n_words_out) {
  n_words_out = (N + 63) / 64;
  const int m = hits_by_peptide.size();
  std::vector<uint64_t> X((size_t)m * (size_t)n_words_out, (uint64_t)0ULL);

  for (int j = 0; j < m; ++j) {
    IntegerVector hits = hits_by_peptide[j];
    uint64_t* col = X.data() + (size_t)j * (size_t)n_words_out;
    for (int k = 0; k < hits.size(); ++k) {
      int idx1 = hits[k]; if (idx1 <= 0) continue;
      int idx0 = idx1 - 1;
      int w = idx0 >> 6, b = idx0 & 63;
      col[w] |= (uint64_t(1) << b);
    }
  }
  return X;
}

// Build bitset (paired): two lists, both indexed by the SAME 1..P subject order
static std::vector<uint64_t> build_bitset_paired(const List& hits_by_peptide, int P, int &n_words_out) {
  return build_bitset_unpaired(hits_by_peptide, P, n_words_out);
}

static std::vector<uint64_t> build_mask_from_indices(const IntegerVector& idx_1based, int /*N*/, int n_words) {
  std::vector<uint64_t> M((size_t)n_words, (uint64_t)0ULL);
  for (int k = 0; k < idx_1based.size(); ++k) {
    int idx1 = idx_1based[k]; if (idx1 <= 0) continue;
    int idx0 = idx1 - 1;
    int w = idx0 >> 6, b = idx0 & 63;
    M[(size_t)w] |= (uint64_t(1) << b);
  }
  return M;
}

// --------- Parallel reducers: compute T given a mask --------------------------

struct TReduceUnpaired : public RcppParallel::Worker {
  BitsetMat X;
  const uint64_t* maskA;
  const uint32_t* totals; // precomputed totals per column (pc(X))
  int m, n_words, nA, nB;
  double eps, den_mult, min_prev, winsor_z;
  std::string stat_mode, weight_mode;

  // outputs (reduced)
  double num, den;
  int used;

  // observed-only summaries
  bool obs;
  double sum_delta, sum_delta_pos, sum_delta_w, sum_delta_pos_w, sum_w, sum_w2;

  TReduceUnpaired(const BitsetMat& X_, const uint64_t* maskA_,
                  const uint32_t* totals_,
                  int m_, int n_words_, int nA_, int nB_,
                  double eps_, double den_mult_, double min_prev_, double winsor_z_,
                  const std::string& stat_mode_, const std::string& weight_mode_,
                  bool obs_pass)
    : X(X_), maskA(maskA_), totals(totals_),
      m(m_), n_words(n_words_), nA(nA_), nB(nB_),
      eps(eps_), den_mult(den_mult_), min_prev(min_prev_), winsor_z(winsor_z_),
      stat_mode(stat_mode_), weight_mode(weight_mode_),
      num(0.0), den(0.0), used(0),
      obs(obs_pass), sum_delta(0.0), sum_delta_pos(0.0),
      sum_delta_w(0.0), sum_delta_pos_w(0.0), sum_w(0.0), sum_w2(0.0) {}

  TReduceUnpaired(TReduceUnpaired& other, Split)
    : X(other.X), maskA(other.maskA), totals(other.totals),
      m(other.m), n_words(other.n_words), nA(other.nA), nB(other.nB),
      eps(other.eps), den_mult(other.den_mult), min_prev(other.min_prev), winsor_z(other.winsor_z),
      stat_mode(other.stat_mode), weight_mode(other.weight_mode),
      num(0.0), den(0.0), used(0),
      obs(other.obs), sum_delta(0.0), sum_delta_pos(0.0),
      sum_delta_w(0.0), sum_delta_pos_w(0.0), sum_w(0.0), sum_w2(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      const uint64_t* col = X.col((int)j);
      uint32_t x1 = 0u;
      for (int w = 0; w < n_words; ++w) x1 += pc64(col[w] & maskA[w]);
      uint32_t x2 = totals[j] - x1;

      double p1 = (x1 + eps) / (nA + den_mult * eps);
      double p2 = (x2 + eps) / (nB + den_mult * eps);
      if ((p1 >= min_prev) || (p2 >= min_prev)) {
        double z;
        if (stat_mode == "asin") {
          double se = std::sqrt(1.0/(4.0*std::max(1, nA)) + 1.0/(4.0*std::max(1, nB)));
          z = (std::asin(std::sqrt(p2)) - std::asin(std::sqrt(p1))) / se;
        } else {
          double se = std::sqrt(p1*(1.0-p1)/std::max(1,nA) + p2*(1.0-p2)/std::max(1,nB));
          if (se < 1e-12) se = 1e-12;
          z = (p2 - p1) / se;
        }
        if (winsor_z < 1e9) {
          if (z >  winsor_z) z =  winsor_z;
          if (z < -winsor_z) z = -winsor_z;
        }

        double wght = 1.0;
        if (weight_mode == "se_invvar") {
          if (stat_mode == "asin") {
            double se = std::sqrt(1.0/(4.0*std::max(1,nA)) + 1.0/(4.0*std::max(1,nB)));
            wght = 1.0/std::max(se, 1e-6);
          } else {
            double se = std::sqrt(p1*(1.0-p1)/std::max(1,nA) + p2*(1.0-p2)/std::max(1,nB));
            wght = 1.0/std::max(se, 1e-6);
          }
        } else if (weight_mode == "n_eff_sqrt") {
          wght = std::sqrt(std::max((double)nA * p1 + (double)nB * p2, 1e-12));
        }

        num += wght * z;
        den += wght * wght;
        ++used;

        if (obs) {
          double d = (p2 - p1);
          sum_delta       += d;
          sum_delta_pos   += (d > 0.0) ? 1.0 : 0.0;
          sum_delta_w     += wght * d;
          sum_delta_pos_w += wght * ((d > 0.0) ? 1.0 : 0.0);
          sum_w           += wght;
          sum_w2          += wght * wght;
        }
      }
    }
  }

  void join(const TReduceUnpaired& rhs) {
    num += rhs.num; den += rhs.den; used += rhs.used;
    if (obs) {
      sum_delta       += rhs.sum_delta;
      sum_delta_pos   += rhs.sum_delta_pos;
      sum_delta_w     += rhs.sum_delta_w;
      sum_delta_pos_w += rhs.sum_delta_pos_w;
      sum_w           += rhs.sum_w;
      sum_w2          += rhs.sum_w2;
    }
  }
};

struct TReducePaired : public RcppParallel::Worker {
  BitsetMat X1, X2;             // per-peptide bitsets for group1 and group2 (P subjects)
  const uint64_t* flipMask;     // bits=1 -> flip subject (swap g1<->g2)
  const uint32_t* t1;           // totals per peptide in g1 (pc(X1))
  const uint32_t* t2;           // totals per peptide in g2 (pc(X2))
  int m, n_words, P;
  double eps, den_mult, min_prev, winsor_z;
  std::string stat_mode, weight_mode;

  double num, den;
  int used;

  bool obs;
  double sum_delta, sum_delta_pos, sum_delta_w, sum_delta_pos_w, sum_w, sum_w2;

  TReducePaired(const BitsetMat& X1_, const BitsetMat& X2_,
                const uint64_t* F_,
                const uint32_t* t1_, const uint32_t* t2_,
                int m_, int n_words_, int P_,
                double eps_, double den_mult_, double min_prev_, double winsor_z_,
                const std::string& stat_mode_, const std::string& weight_mode_,
                bool obs_pass)
    : X1(X1_), X2(X2_), flipMask(F_), t1(t1_), t2(t2_),
      m(m_), n_words(n_words_), P(P_),
      eps(eps_), den_mult(den_mult_), min_prev(min_prev_), winsor_z(winsor_z_),
      stat_mode(stat_mode_), weight_mode(weight_mode_),
      num(0.0), den(0.0), used(0),
      obs(obs_pass), sum_delta(0.0), sum_delta_pos(0.0),
      sum_delta_w(0.0), sum_delta_pos_w(0.0), sum_w(0.0), sum_w2(0.0) {}

  TReducePaired(TReducePaired& other, Split)
    : X1(other.X1), X2(other.X2), flipMask(other.flipMask),
      t1(other.t1), t2(other.t2),
      m(other.m), n_words(other.n_words), P(other.P),
      eps(other.eps), den_mult(other.den_mult), min_prev(other.min_prev), winsor_z(other.winsor_z),
      stat_mode(other.stat_mode), weight_mode(other.weight_mode),
      num(0.0), den(0.0), used(0),
      obs(other.obs), sum_delta(0.0), sum_delta_pos(0.0),
      sum_delta_w(0.0), sum_delta_pos_w(0.0), sum_w(0.0), sum_w2(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      const uint64_t* c1 = X1.col((int)j);
      const uint64_t* c2 = X2.col((int)j);
      // x1 after flips: t1 - pc(X1&F) + pc(X2&F)
      uint32_t pc1F = 0u, pc2F = 0u;
      for (int w = 0; w < n_words; ++w) {
        uint64_t Fw = flipMask[w];
        pc1F += pc64(c1[w] & Fw);
        pc2F += pc64(c2[w] & Fw);
      }
      uint32_t x1 = t1[j] - pc1F + pc2F;
      uint32_t x2 = t2[j] - pc2F + pc1F;

      double p1 = (x1 + eps) / (P + den_mult * eps);
      double p2 = (x2 + eps) / (P + den_mult * eps);
      if ((p1 >= min_prev) || (p2 >= min_prev)) {
        double z;
        if (stat_mode == "asin") {
          double se = std::sqrt(1.0/(4.0*std::max(1,P)) + 1.0/(4.0*std::max(1,P)));
          z = (std::asin(std::sqrt(p2)) - std::asin(std::sqrt(p1))) / se;
        } else {
          double se = std::sqrt(p1*(1.0-p1)/std::max(1,P) + p2*(1.0-p2)/std::max(1,P));
          if (se < 1e-12) se = 1e-12;
          z = (p2 - p1) / se;
        }
        if (winsor_z < 1e9) {
          if (z >  winsor_z) z =  winsor_z;
          if (z < -winsor_z) z = -winsor_z;
        }

        double wght = 1.0;
        if (weight_mode == "se_invvar") {
          if (stat_mode == "asin") {
            double se = std::sqrt(1.0/(4.0*std::max(1,P)) + 1.0/(4.0*std::max(1,P)));
            wght = 1.0/std::max(se, 1e-6);
          } else {
            double se = std::sqrt(p1*(1.0-p1)/std::max(1,P) + p2*(1.0-p2)/std::max(1,P));
            wght = 1.0/std::max(se, 1e-6);
          }
        } else if (weight_mode == "n_eff_sqrt") {
          wght = std::sqrt(std::max((double)P * (p1 + p2), 1e-12)); // both groups size=P
        }

        num += wght * z;
        den += wght * wght;
        ++used;

        if (obs) {
          double d = (p2 - p1);
          sum_delta       += d;
          sum_delta_pos   += (d > 0.0) ? 1.0 : 0.0;
          sum_delta_w     += wght * d;
          sum_delta_pos_w += wght * ((d > 0.0) ? 1.0 : 0.0);
          sum_w           += wght;
          sum_w2          += wght * wght;
        }
      }
    }
  }

  void join(const TReducePaired& rhs) {
    num += rhs.num; den += rhs.den; used += rhs.used;
    if (obs) {
      sum_delta       += rhs.sum_delta;
      sum_delta_pos   += rhs.sum_delta_pos;
      sum_delta_w     += rhs.sum_delta_w;
      sum_delta_pos_w += rhs.sum_delta_pos_w;
      sum_w           += rhs.sum_w;
      sum_w2          += rhs.sum_w2;
    }
  }
};

// Compute totals per column once
static std::vector<uint32_t> col_totals(const BitsetMat& X) {
  std::vector<uint32_t> tt((size_t)X.m, 0u);
  for (int j = 0; j < X.m; ++j) {
    const uint64_t* c = X.col(j);
    uint32_t s = 0u; for (int w = 0; w < X.n_words; ++w) s += pc64(c[w]);
    tt[(size_t)j] = s;
  }
  return tt;
}

// --------- unified entry point ------------------------------------------------

// [[Rcpp::export]]
Rcpp::List perm_bitset_T_parallel(
    bool paired,
    Rcpp::List hits_by_peptide_unpaired,
    int N_unpaired,
    Rcpp::IntegerVector g1_idx_unpaired,
    int nA,
    int nB,
    Rcpp::List hits_g1_paired,
    Rcpp::List hits_g2_paired,
    int P,
    int B,
    double eps,
    double den_mult,
    double min_max_prev,
    double winsor_z,
    std::string stat_mode,
    std::string weight_mode,
    int seed,
    int grain = 512
) {
  std::mt19937 rng(seed);

  double T_obs = 0.0;
  int n_used = 0;
  double m_eff = 0.0, mean_delta = 0.0, frac_delta_pos = 0.0, mean_delta_w = 0.0, frac_delta_pos_w = 0.0;

  int m = 0, n_words = 0;
  int b_hits = 0;

  if (!paired) {
    // ---------- UNPAIRED ------------------------------------------------------
    m = hits_by_peptide_unpaired.size();
    std::vector<uint64_t> Xbuf = build_bitset_unpaired(hits_by_peptide_unpaired, N_unpaired, n_words);
    BitsetMat X{ Xbuf.data(), m, n_words };

    // totals per column
    std::vector<uint32_t> tot = col_totals(X);

    // observed mask (group1 indices)
    std::vector<uint64_t> Mobs = build_mask_from_indices(g1_idx_unpaired, N_unpaired, n_words);

    // observed T + summaries (parallel)
    TReduceUnpaired red_obs(X, Mobs.data(), tot.data(), m, n_words, nA, nB,
                            eps, den_mult, min_max_prev, winsor_z, stat_mode, weight_mode, true);
    parallelReduce(0, (size_t)m, red_obs, (size_t)grain);
    if (red_obs.den > 0.0) T_obs = red_obs.num / std::sqrt(red_obs.den);
    n_used          = red_obs.used;
    if (red_obs.sum_w2 > 0.0) m_eff = (red_obs.sum_w * red_obs.sum_w) / red_obs.sum_w2;
    mean_delta      = (n_used > 0) ? red_obs.sum_delta / (double)n_used : 0.0;
    frac_delta_pos  = (n_used > 0) ? red_obs.sum_delta_pos / (double)n_used : 0.0;
    mean_delta_w    = (red_obs.sum_w > 0.0) ? red_obs.sum_delta_w / red_obs.sum_w : 0.0;
    frac_delta_pos_w= (red_obs.sum_w > 0.0) ? red_obs.sum_delta_pos_w / red_obs.sum_w : 0.0;

    // permute (balanced masks)
    std::vector<int> tmp_idx(N_unpaired);
    std::vector<uint64_t> M((size_t)n_words, (uint64_t)0ULL);

    for (int b = 0; b < B; ++b) {
      // partial shuffle for first nA
      for (int i = 0; i < N_unpaired; ++i) tmp_idx[i] = i;
      for (int i = 0; i < nA; ++i) {
        std::uniform_int_distribution<int> dist(i, N_unpaired - 1);
        int j = dist(rng);
        std::swap(tmp_idx[i], tmp_idx[j]);
      }
      std::fill(M.begin(), M.end(), (uint64_t)0ULL);
      for (int i = 0; i < nA; ++i) {
        int idx0 = tmp_idx[i];
        int w = idx0 >> 6, bb = idx0 & 63;
        M[(size_t)w] |= (uint64_t(1) << bb);
      }

      TReduceUnpaired red(X, M.data(), tot.data(), m, n_words, nA, nB,
                          eps, den_mult, min_max_prev, winsor_z, stat_mode, weight_mode, false);
      parallelReduce(0, (size_t)m, red, (size_t)grain);
      double Tb = (red.den > 0.0) ? (red.num / std::sqrt(red.den)) : 0.0;
      if (std::fabs(Tb) >= std::fabs(T_obs)) ++b_hits;
    }
  } else {
    // ---------- PAIRED --------------------------------------------------------
    m = hits_g1_paired.size();
    std::vector<uint64_t> X1buf = build_bitset_paired(hits_g1_paired, P, n_words);
    std::vector<uint64_t> X2buf = build_bitset_paired(hits_g2_paired, P, n_words);
    BitsetMat X1{ X1buf.data(), m, n_words }, X2{ X2buf.data(), m, n_words };

    std::vector<uint32_t> t1 = col_totals(X1);
    std::vector<uint32_t> t2 = col_totals(X2);

    // observed pass: flipMask = 0 (no flips)
    std::vector<uint64_t> Fobs((size_t)n_words, (uint64_t)0ULL);
    TReducePaired red_obs(X1, X2, Fobs.data(), t1.data(), t2.data(),
                          m, n_words, P,
                          eps, den_mult, min_max_prev, winsor_z, stat_mode, weight_mode, true);
    parallelReduce(0, (size_t)m, red_obs, (size_t)grain);
    if (red_obs.den > 0.0) T_obs = red_obs.num / std::sqrt(red_obs.den);
    n_used          = red_obs.used;
    if (red_obs.sum_w2 > 0.0) m_eff = (red_obs.sum_w * red_obs.sum_w) / red_obs.sum_w2;
    mean_delta      = (n_used > 0) ? red_obs.sum_delta / (double)n_used : 0.0;
    frac_delta_pos  = (n_used > 0) ? red_obs.sum_delta_pos / (double)n_used : 0.0;
    mean_delta_w    = (red_obs.sum_w > 0.0) ? red_obs.sum_delta_w / red_obs.sum_w : 0.0;
    frac_delta_pos_w= (red_obs.sum_w > 0.0) ? red_obs.sum_delta_pos_w / red_obs.sum_w : 0.0;

    // permutations: random independent flips per subject (Bernoulli 0.5)
    std::bernoulli_distribution bern(0.5);
    std::vector<uint64_t> F((size_t)n_words, (uint64_t)0ULL);
    for (int b = 0; b < B; ++b) {
      std::fill(F.begin(), F.end(), (uint64_t)0ULL);
      for (int s = 0; s < P; ++s) {
        if (bern(rng)) {
          int w = s >> 6, bb = s & 63;
          F[(size_t)w] |= (uint64_t(1) << bb);
        }
      }
      TReducePaired red(X1, X2, F.data(), t1.data(), t2.data(),
                        m, n_words, P,
                        eps, den_mult, min_max_prev, winsor_z, stat_mode, weight_mode, false);
      parallelReduce(0, (size_t)m, red, (size_t)grain);
      double Tb = (red.den > 0.0) ? (red.num / std::sqrt(red.den)) : 0.0;
      if (std::fabs(Tb) >= std::fabs(T_obs)) ++b_hits;
    }
  }

  double p_perm = (1.0 + (double)b_hits) / (1.0 + (double)B);

  return Rcpp::List::create(
    _["T_obs"]            = T_obs,
    _["p_perm"]           = p_perm,
    _["b"]                = b_hits,
    _["n_peptides_used"]  = n_used,
    _["m_eff"]            = m_eff,
    _["mean_delta"]       = mean_delta,
    _["frac_delta_pos"]   = frac_delta_pos,
    _["mean_delta_w"]     = mean_delta_w,
    _["frac_delta_pos_w"] = frac_delta_pos_w
  );
}
