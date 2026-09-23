// Memory-efficient computational kernels for HDZACP.
//
// The public R interface performs data validation, tie handling, and method
// labeling.  These kernels operate on strict pooled-rank permutations and use
// the exact hypergeometric null moments from the manuscript.
//
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr double PI = 3.141592653589793238462643383279502884;

inline double entropy_h(const double u) {
  if (u <= 0.0 || u >= 1.0) return 0.0;
  return -(u * std::log(u) + (1.0 - u) * std::log1p(-u));
}

// Store a*h(c/a), rather than all g(k, ell, c).  Since
//
// g(k,ell,c) = {k h(c/k) + (n-k) h((ell-c)/(n-k))}
//              / { (ell-.5)(n-ell+.5) },
//
// this representation is exact and uses O(n^2), rather than O(n^3), memory.
struct EntropyTable {
  int n;
  std::vector<std::size_t> offset;
  std::vector<double> h;
  std::vector<double> inv_d;

  explicit EntropyTable(const int n_)
      : n(n_), offset(n_ + 1), inv_d(n_ + 1, 0.0) {
    std::size_t total = 0;
    for (int a = 0; a <= n; ++a) {
      offset[a] = total;
      total += static_cast<std::size_t>(a + 1);
    }
    h.assign(total, 0.0);
    for (int a = 1; a <= n; ++a) {
      const std::size_t base = offset[a];
      for (int c = 1; c < a; ++c) {
        h[base + static_cast<std::size_t>(c)] =
          static_cast<double>(a) * entropy_h(
            static_cast<double>(c) / static_cast<double>(a));
      }
    }
    for (int ell = 1; ell <= n; ++ell) {
      inv_d[ell] = 1.0 /
        ((static_cast<double>(ell) - 0.5) *
         (static_cast<double>(n - ell) + 0.5));
    }
  }

  inline double g(const int k, const int ell, const int c) const {
    return (h[offset[k] + static_cast<std::size_t>(c)] +
            h[offset[n - k] + static_cast<std::size_t>(ell - c)]) *
      inv_d[ell];
  }

  double bytes() const {
    return static_cast<double>(
      h.size() * sizeof(double) +
      offset.size() * sizeof(std::size_t) +
      inv_d.size() * sizeof(double));
  }
};

void validate_n_mn(const int n, const int mn) {
  if (n < 4 || mn < 1 || 2 * mn >= n) {
    Rcpp::stop("Require n >= 4 and 1 <= mn < n/2.");
  }
}

struct NullSetup {
  int n;
  int mn;
  EntropyTable table;
  std::vector<double> mu;
  std::vector<double> sd;

  NullSetup(const int n_, const int mn_)
      : n(n_), mn(mn_), table(n_),
        mu(n_ + 1, NA_REAL), sd(n_ + 1, NA_REAL) {
    // Exact first and second moments of the hypergeometric rank-count path.
    for (int k = mn; k <= n / 2; ++k) {
      std::vector<long double> probability(k + 1, 0.0L);
      std::vector<long double> first(k + 1, 0.0L);
      std::vector<long double> second(k + 1, 0.0L);
      std::vector<long double> next_probability(k + 1, 0.0L);
      std::vector<long double> next_first(k + 1, 0.0L);
      std::vector<long double> next_second(k + 1, 0.0L);
      probability[0] = 1.0L;

      for (int ell = 0; ell < n; ++ell) {
        std::fill(next_probability.begin(), next_probability.end(), 0.0L);
        std::fill(next_first.begin(), next_first.end(), 0.0L);
        std::fill(next_second.begin(), next_second.end(), 0.0L);
        const int lo = std::max(0, k + ell - n);
        const int hi = std::min(k, ell);
        const long double denominator = static_cast<long double>(n - ell);

        for (int c = lo; c <= hi; ++c) {
          if (probability[c] == 0.0L) continue;
          const long double q0 =
            static_cast<long double>(n - ell - k + c) / denominator;
          const long double q1 =
            static_cast<long double>(k - c) / denominator;

          if (q0 > 0.0L) {
            const long double increment =
              static_cast<long double>(table.g(k, ell + 1, c));
            next_probability[c] += q0 * probability[c];
            next_first[c] += q0 *
              (first[c] + increment * probability[c]);
            next_second[c] += q0 *
              (second[c] + 2.0L * increment * first[c] +
               increment * increment * probability[c]);
          }
          if (q1 > 0.0L) {
            const int c1 = c + 1;
            const long double increment =
              static_cast<long double>(table.g(k, ell + 1, c1));
            next_probability[c1] += q1 * probability[c];
            next_first[c1] += q1 *
              (first[c] + increment * probability[c]);
            next_second[c1] += q1 *
              (second[c] + 2.0L * increment * first[c] +
               increment * increment * probability[c]);
          }
        }

        probability.swap(next_probability);
        first.swap(next_first);
        second.swap(next_second);
      }

      const long double mean = first[k];
      long double variance = second[k] - mean * mean;
      if (variance < 0.0L && variance > -1e-24L) variance = 0.0L;
      if (!(variance > 0.0L)) {
        Rcpp::stop("Non-positive exact null variance at k=%d.", k);
      }
      mu[k] = static_cast<double>(mean);
      sd[k] = std::sqrt(static_cast<double>(variance));
    }

    for (int k = n / 2 + 1; k <= n - mn; ++k) {
      mu[k] = mu[n - k];
      sd[k] = sd[n - k];
    }
  }
};

const NullSetup& cached_setup(const int n, const int mn) {
  // R invokes exported routines serially in an ordinary process.  Separate
  // PSOCK workers load separate DLLs and therefore have separate caches.
  static std::unique_ptr<NullSetup> setup;
  if (!setup || setup->n != n || setup->mn != mn) {
    setup = std::make_unique<NullSetup>(n, mn);
  }
  return *setup;
}

Rcpp::List moments_list(const NullSetup& setup) {
  const int count = setup.n - 2 * setup.mn + 1;
  Rcpp::IntegerVector k(count);
  Rcpp::NumericVector mu(count);
  Rcpp::NumericVector sd(count);
  for (int index = 0; index < count; ++index) {
    const int split = setup.mn + index;
    k[index] = split;
    mu[index] = setup.mu[split];
    sd[index] = setup.sd[split];
  }
  return Rcpp::List::create(
    Rcpp::_ ["k"] = k,
    Rcpp::_ ["mu"] = mu,
    Rcpp::_ ["sd"] = sd,
    // In R, vector element k + 1 stores the moment for split k.
    Rcpp::_ ["mu_full"] = Rcpp::wrap(setup.mu),
    Rcpp::_ ["sd_full"] = Rcpp::wrap(setup.sd),
    Rcpp::_ ["table_bytes"] = setup.table.bytes()
  );
}

void validate_moments(const Rcpp::NumericVector& mu_full,
                      const Rcpp::NumericVector& sd_full,
                      const int n,
                      const int mn) {
  if (mu_full.size() != n + 1 || sd_full.size() != n + 1) {
    Rcpp::stop("mu_full and sd_full must have length n + 1.");
  }
  for (int k = mn; k <= n - mn; ++k) {
    if (!std::isfinite(mu_full[k]) || !std::isfinite(sd_full[k]) ||
        !(sd_full[k] > 0.0)) {
      Rcpp::stop("Moments must be finite and sd positive at every candidate split.");
    }
  }
}

std::vector<int> validate_rank_matrix(
    const Rcpp::IntegerMatrix& rank_matrix) {
  const int n = rank_matrix.nrow();
  const int p = rank_matrix.ncol();
  if (n < 1 || p < 1) Rcpp::stop("rank_matrix must be non-empty.");

  std::vector<int> ranks(static_cast<std::size_t>(n) * p);
  std::vector<unsigned char> seen(n, 0);
  for (int j = 0; j < p; ++j) {
    std::fill(seen.begin(), seen.end(), 0);
    for (int i = 0; i < n; ++i) {
      const int rank = rank_matrix(i, j);
      if (rank < 1 || rank > n || seen[rank - 1]) {
        Rcpp::stop("Each rank_matrix column must be a permutation of 1,...,n.");
      }
      seen[rank - 1] = 1;
      ranks[static_cast<std::size_t>(i) * p + j] = rank;
    }
  }
  return ranks;
}

std::vector<int> validate_order(const Rcpp::IntegerVector& order,
                                const int n) {
  if (order.size() != n) Rcpp::stop("order must have length n.");
  std::vector<int> order0(n);
  std::vector<unsigned char> seen(n, 0);
  for (int i = 0; i < n; ++i) {
    const int value = order[i] - 1;
    if (value < 0 || value >= n || seen[value]) {
      Rcpp::stop("order must be a permutation of 1,...,n.");
    }
    seen[value] = 1;
    order0[i] = value;
  }
  return order0;
}

struct ScanResult {
  double sum;
  double maximum;
  double dense_mid;
  double coordinate_mid_deficit;
  double coordinate_mid_raw;
  int cp_sum;
  int cp_max;
};

ScanResult scan_rank_matrix(const std::vector<int>& ranks,
                            const std::vector<int>& order,
                            const int n,
                            const int p,
                            const int mn,
                            const std::vector<double>& mu,
                            const std::vector<double>& sd,
                            const EntropyTable& table) {
  // Layout is rank x coordinate, so the coordinate loop is contiguous.
  std::vector<unsigned char> included(static_cast<std::size_t>(n) * p, 0);
  std::vector<int> cumulative(p, 0);
  std::vector<double> z(p, 0.0);

  ScanResult result{
    -std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity(),
    NA_REAL, NA_REAL, NA_REAL, mn, mn
  };
  int inserted = 0;
  const double root_p = std::sqrt(static_cast<double>(p));

  for (int k = mn; k <= n - mn; ++k) {
    while (inserted < k) {
      const std::size_t row_base =
        static_cast<std::size_t>(order[inserted]) * p;
      for (int j = 0; j < p; ++j) {
        const int rank0 = ranks[row_base + j] - 1;
        included[static_cast<std::size_t>(rank0) * p + j] = 1;
      }
      ++inserted;
    }

    std::fill(cumulative.begin(), cumulative.end(), 0);
    std::fill(z.begin(), z.end(), 0.0);
    const double* h_left = table.h.data() + table.offset[k];
    const double* h_right = table.h.data() + table.offset[n - k];
    for (int ell = 1; ell <= n; ++ell) {
      const std::size_t rank_base = static_cast<std::size_t>(ell - 1) * p;
      const double inv_d = table.inv_d[ell];
      for (int j = 0; j < p; ++j) {
        const int count =
          (cumulative[j] += included[rank_base + static_cast<std::size_t>(j)]);
        z[j] += (h_left[count] + h_right[ell - count]) * inv_d;
      }
    }

    double deficit_sum = 0.0;
    double deficit_max = -std::numeric_limits<double>::infinity();
    for (int j = 0; j < p; ++j) {
      const double deficit = (mu[k] - z[j]) / sd[k];
      deficit_sum += deficit;
      if (deficit > deficit_max) deficit_max = deficit;
    }
    const double dense = deficit_sum / root_p;
    // Strict comparison and increasing k implement the smallest-maximizer rule.
    if (dense > result.sum) {
      result.sum = dense;
      result.cp_sum = k;
    }
    if (deficit_max > result.maximum) {
      result.maximum = deficit_max;
      result.cp_max = k;
    }
    if (k == n / 2) {
      result.dense_mid = dense;
      result.coordinate_mid_deficit = (mu[k] - z[0]) / sd[k];
      result.coordinate_mid_raw = z[0];
    }
  }
  return result;
}

inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31U);
}

inline double clipped_cot_score(double pvalue) {
  const double epsilon = 8.0 * std::numeric_limits<double>::epsilon();
  pvalue = std::max(epsilon, std::min(1.0 - epsilon, pvalue));
  return std::tan(PI * (0.5 - pvalue));
}

std::vector<int> upper_tail_counts(const std::vector<double>& values) {
  std::vector<double> sorted(values);
  std::sort(sorted.begin(), sorted.end());
  std::vector<int> counts(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    const auto first_equal =
      std::lower_bound(sorted.begin(), sorted.end(), values[i]);
    counts[i] = static_cast<int>(sorted.end() - first_equal);
  }
  return counts;
}

int effective_threads(const int requested) {
  int used = 1;
#ifdef _OPENMP
  #pragma omp parallel num_threads(requested)
  {
    #pragma omp single
    used = omp_get_num_threads();
  }
#else
  (void)requested;
#endif
  return used;
}

// Independent uniform rank permutations have exactly the same null law as
// ranks of iid continuous coordinates.  The Gaussian branches add the stated
// AR(1) or block dependence across coordinates, never across rows.
std::vector<int> make_null_ranks(const int n,
                                 const int p,
                                 const std::string& design,
                                 const double rho,
                                 const int block_size,
                                 std::mt19937_64& rng) {
  std::vector<int> ranks(static_cast<std::size_t>(n) * p);
  if (design == "independent" || rho == 0.0) {
    std::vector<int> permutation(n);
    for (int j = 0; j < p; ++j) {
      std::iota(permutation.begin(), permutation.end(), 1);
      std::shuffle(permutation.begin(), permutation.end(), rng);
      for (int i = 0; i < n; ++i) {
        ranks[static_cast<std::size_t>(i) * p + j] = permutation[i];
      }
    }
    return ranks;
  }

  std::normal_distribution<double> normal(0.0, 1.0);
  std::vector<double> data(static_cast<std::size_t>(n) * p, 0.0);
  if (design == "ar1") {
    const double noise_sd = std::sqrt(1.0 - rho * rho);
    for (int i = 0; i < n; ++i) {
      const std::size_t row = static_cast<std::size_t>(i) * p;
      data[row] = normal(rng);
      for (int j = 1; j < p; ++j) {
        data[row + j] = rho * data[row + j - 1] + noise_sd * normal(rng);
      }
    }
  } else {
    const double common_sd = std::sqrt(rho);
    const double noise_sd = std::sqrt(1.0 - rho);
    for (int i = 0; i < n; ++i) {
      const std::size_t row = static_cast<std::size_t>(i) * p;
      for (int block_start = 0; block_start < p;
           block_start += block_size) {
        const double common = common_sd * normal(rng);
        const int block_end = std::min(p, block_start + block_size);
        for (int j = block_start; j < block_end; ++j) {
          data[row + j] = common + noise_sd * normal(rng);
        }
      }
    }
  }

  std::vector<std::pair<double, int>> values(n);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) {
      values[i] = std::make_pair(data[static_cast<std::size_t>(i) * p + j], i);
    }
    std::sort(values.begin(), values.end());
    for (int ell = 0; ell < n; ++ell) {
      ranks[static_cast<std::size_t>(values[ell].second) * p + j] = ell + 1;
    }
  }
  return ranks;
}

void validate_finite_matrix(const Rcpp::NumericMatrix& matrix,
                            const char* name) {
  for (R_xlen_t index = 0; index < matrix.size(); ++index) {
    if (!std::isfinite(matrix[index])) {
      Rcpp::stop("%s must contain only finite values.", name);
    }
  }
}

std::pair<int, double> hdd_stat_from_distance(
    const Rcpp::NumericMatrix& distance_matrix,
    const std::vector<int>& order) {
  const int n = distance_matrix.nrow();
  int cp_r = 1;
  double max_delta = 0.0;  // HDDchangepoint sets the first column to zero.

  for (int position = 1; position < n; ++position) {
    double delta = 0.0;
    const int left = order[position - 1];
    const int right = order[position];
    for (int i = 0; i < n; ++i) {
      delta += std::abs(distance_matrix(i, right) -
                        distance_matrix(i, left));
    }
    delta /= static_cast<double>(n);
    if (delta > max_delta) {
      max_delta = delta;
      cp_r = position + 1;  // R one-based first observation after the split.
    }
  }

  if (cp_r == n || cp_r <= 2) return std::make_pair(cp_r, 0.0);
  const int n_before = cp_r - 1;
  const int n_after = n - n_before;
  double total = 0.0;
  for (int i = 0; i < n; ++i) {
    double sum_before = 0.0;
    double sumsq_before = 0.0;
    double sum_after = 0.0;
    double sumsq_after = 0.0;
    for (int position = 0; position < n_before; ++position) {
      const double value = distance_matrix(i, order[position]);
      sum_before += value;
      sumsq_before += value * value;
    }
    for (int position = n_before; position < n; ++position) {
      const double value = distance_matrix(i, order[position]);
      sum_after += value;
      sumsq_after += value * value;
    }
    const double mean_before = sum_before / n_before;
    const double mean_after = sum_after / n_after;
    total += sumsq_before / n_before + sumsq_after / n_after -
      2.0 * mean_before * mean_after;
  }
  return std::make_pair(cp_r, total / static_cast<double>(n));
}

}  // namespace


// Exact null mean and standard deviation for every retained split.
// R element k + 1 of mu_full/sd_full corresponds to split k.
// [[Rcpp::export]]
Rcpp::List hdzacp_null_moments_cpp(const int n, const int mn) {
  validate_n_mn(n, mn);
  return moments_list(cached_setup(n, mn));
}


// Scan a supplied strict pooled-rank matrix in a supplied observation order.
// Each rank column and order must be a permutation of 1,...,n.
// [[Rcpp::export]]
Rcpp::NumericVector hdzacp_scan_rank_cpp(
    const Rcpp::IntegerMatrix rank_matrix,
    const Rcpp::IntegerVector order,
    const Rcpp::NumericVector mu_full,
    const Rcpp::NumericVector sd_full,
    const int mn) {
  const int n = rank_matrix.nrow();
  const int p = rank_matrix.ncol();
  validate_n_mn(n, mn);
  if (p < 1) Rcpp::stop("rank_matrix must have at least one column.");
  validate_moments(mu_full, sd_full, n, mn);
  const std::vector<int> ranks = validate_rank_matrix(rank_matrix);
  const std::vector<int> order0 = validate_order(order, n);
  const std::vector<double> mu(mu_full.begin(), mu_full.end());
  const std::vector<double> sd(sd_full.begin(), sd_full.end());
  const ScanResult result = scan_rank_matrix(
    ranks, order0, n, p, mn, mu, sd, cached_setup(n, mn).table);

  return Rcpp::NumericVector::create(
    Rcpp::_ ["S_sum"] = result.sum,
    Rcpp::_ ["S_max"] = result.maximum,
    Rcpp::_ ["cp_sum"] = result.cp_sum,
    Rcpp::_ ["cp_max"] = result.cp_max,
    Rcpp::_ ["dense_mid"] = result.dense_mid,
    Rcpp::_ ["coord_mid_deficit"] = result.coordinate_mid_deficit,
    Rcpp::_ ["coord_mid_raw"] = result.coordinate_mid_raw
  );
}


// Whole-vector permutation calibration for ZAS, ZAM, and the symmetric-orbit
// ZAC procedure.  The observed ordering is orbit element zero.
// [[Rcpp::export]]
Rcpp::List hdzacp_test_rank_cpp(
    const Rcpp::IntegerMatrix rank_matrix,
    const int mn,
    const int n_perm,
    const int seed,
    const Rcpp::NumericVector mu_full,
    const Rcpp::NumericVector sd_full,
    const bool keep_orbit = false) {
  const int n = rank_matrix.nrow();
  const int p = rank_matrix.ncol();
  validate_n_mn(n, mn);
  if (p < 1 || n_perm < 1) {
    Rcpp::stop("Require p >= 1 and n_perm >= 1.");
  }
  validate_moments(mu_full, sd_full, n, mn);
  const std::vector<int> ranks = validate_rank_matrix(rank_matrix);
  const std::vector<double> mu(mu_full.begin(), mu_full.end());
  const std::vector<double> sd(sd_full.begin(), sd_full.end());
  const EntropyTable& table = cached_setup(n, mn).table;

  const int orbit_size = n_perm + 1;
  std::vector<double> orbit_sum(orbit_size);
  std::vector<double> orbit_max(orbit_size);
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::mt19937_64 rng(splitmix64(static_cast<std::uint64_t>(
    static_cast<std::uint32_t>(seed))));
  int observed_cp_sum = mn;
  int observed_cp_max = mn;

  for (int b = 0; b < orbit_size; ++b) {
    if (b > 0) {
      std::iota(order.begin(), order.end(), 0);
      std::shuffle(order.begin(), order.end(), rng);
    }
    const ScanResult statistic =
      scan_rank_matrix(ranks, order, n, p, mn, mu, sd, table);
    orbit_sum[b] = statistic.sum;
    orbit_max[b] = statistic.maximum;
    if (b == 0) {
      observed_cp_sum = statistic.cp_sum;
      observed_cp_max = statistic.cp_max;
    }
  }

  const std::vector<int> count_sum = upper_tail_counts(orbit_sum);
  const std::vector<int> count_max = upper_tail_counts(orbit_max);
  std::vector<double> orbit_p_sum(orbit_size);
  std::vector<double> orbit_p_max(orbit_size);
  std::vector<double> orbit_cauchy(orbit_size);
  for (int b = 0; b < orbit_size; ++b) {
    orbit_p_sum[b] = static_cast<double>(count_sum[b]) / orbit_size;
    orbit_p_max[b] = static_cast<double>(count_max[b]) / orbit_size;
    orbit_cauchy[b] =
      0.5 * clipped_cot_score(orbit_p_sum[b]) +
      0.5 * clipped_cot_score(orbit_p_max[b]);
  }
  const std::vector<int> count_cauchy = upper_tail_counts(orbit_cauchy);

  Rcpp::List answer = Rcpp::List::create(
    Rcpp::_ ["S_sum"] = orbit_sum[0],
    Rcpp::_ ["S_max"] = orbit_max[0],
    Rcpp::_ ["cp_sum"] = observed_cp_sum,
    Rcpp::_ ["cp_max"] = observed_cp_max,
    Rcpp::_ ["T_cauchy"] = orbit_cauchy[0],
    Rcpp::_ ["p_sum"] = orbit_p_sum[0],
    Rcpp::_ ["p_max"] = orbit_p_max[0],
    Rcpp::_ ["p_cauchy"] =
      static_cast<double>(count_cauchy[0]) / orbit_size
  );

  if (keep_orbit) {
    answer["orbit"] = Rcpp::DataFrame::create(
      Rcpp::_ ["ordering"] = Rcpp::seq(0, n_perm),
      Rcpp::_ ["S_sum"] = Rcpp::wrap(orbit_sum),
      Rcpp::_ ["S_max"] = Rcpp::wrap(orbit_max),
      Rcpp::_ ["p_sum"] = Rcpp::wrap(orbit_p_sum),
      Rcpp::_ ["p_max"] = Rcpp::wrap(orbit_p_max),
      Rcpp::_ ["T_cauchy"] = Rcpp::wrap(orbit_cauchy)
    );
  }
  return answer;
}


// Simulate the original rank scans under the analytic null designs used in
// the manuscript's direct-asymptotic-calibration audit.
// [[Rcpp::export]]
Rcpp::List hdzacp_size_simulate_cpp(
    const int n,
    const int p,
    const double eta,
    const int n_rep,
    const int seed,
    const int n_threads,
    const std::string design = "independent",
    const double rho = 0.0,
    const int block_size = 5) {
  if (!std::isfinite(eta) || !(eta > 0.0 && eta < 0.5)) {
    Rcpp::stop("Require 0 < eta < 0.5.");
  }
  const int mn = static_cast<int>(std::ceil(n * eta));
  validate_n_mn(n, mn);
  if (p < 1 || n_rep < 1 || n_threads < 1 || block_size < 1) {
    Rcpp::stop("p, n_rep, n_threads, and block_size must be positive.");
  }
  if (design != "independent" && design != "ar1" && design != "block") {
    Rcpp::stop("design must be one of 'independent', 'ar1', or 'block'.");
  }
  if (!std::isfinite(rho) || std::abs(rho) >= 1.0 ||
      (design == "block" && rho < 0.0)) {
    Rcpp::stop("Require |rho| < 1; block additionally requires rho >= 0.");
  }

  const NullSetup& setup = cached_setup(n, mn);
  Rcpp::NumericVector sums(n_rep);
  Rcpp::NumericVector maxima(n_rep);
  Rcpp::NumericVector dense_mid(n_rep);
  Rcpp::NumericVector coordinate_mid_deficit(n_rep);
  Rcpp::NumericVector coordinate_mid_raw(n_rep);
  Rcpp::IntegerVector cp_sum(n_rep);
  Rcpp::IntegerVector cp_max(n_rep);
  std::vector<int> identity(n);
  std::iota(identity.begin(), identity.end(), 0);

  // No R API is called inside the OpenMP loop.
  double* sum_ptr = sums.begin();
  double* max_ptr = maxima.begin();
  double* dense_ptr = dense_mid.begin();
  double* deficit_ptr = coordinate_mid_deficit.begin();
  double* raw_ptr = coordinate_mid_raw.begin();
  int* cp_sum_ptr = cp_sum.begin();
  int* cp_max_ptr = cp_max.begin();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int replication = 0; replication < n_rep; ++replication) {
    const std::uint64_t replication_seed = splitmix64(
      static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed)) ^
      (0xd1b54a32d192ed03ULL *
       static_cast<std::uint64_t>(replication + 1)));
    std::mt19937_64 rng(replication_seed);
    const std::vector<int> ranks =
      make_null_ranks(n, p, design, rho, block_size, rng);
    const ScanResult result = scan_rank_matrix(
      ranks, identity, n, p, mn, setup.mu, setup.sd, setup.table);
    sum_ptr[replication] = result.sum;
    max_ptr[replication] = result.maximum;
    dense_ptr[replication] = result.dense_mid;
    deficit_ptr[replication] = result.coordinate_mid_deficit;
    raw_ptr[replication] = result.coordinate_mid_raw;
    cp_sum_ptr[replication] = result.cp_sum;
    cp_max_ptr[replication] = result.cp_max;
  }

  return Rcpp::List::create(
    Rcpp::_ ["S_sum"] = sums,
    Rcpp::_ ["S_max"] = maxima,
    Rcpp::_ ["cp_sum"] = cp_sum,
    Rcpp::_ ["cp_max"] = cp_max,
    Rcpp::_ ["dense_mid"] = dense_mid,
    Rcpp::_ ["coord_mid_deficit"] = coordinate_mid_deficit,
    Rcpp::_ ["coord_mid_raw"] = coordinate_mid_raw,
    Rcpp::_ ["k_mid"] = n / 2,
    Rcpp::_ ["mn"] = mn,
    Rcpp::_ ["moments"] = moments_list(setup),
    Rcpp::_ ["threads_used"] = effective_threads(n_threads),
    Rcpp::_ ["openmp_enabled"] =
#ifdef _OPENMP
      true
#else
      false
#endif
  );
}


// Exact-grid simulation of the stationary OU reference process.  The coarse
// grid is the even-index subset of the fine grid on the same paths.
// [[Rcpp::export]]
Rcpp::NumericMatrix hdzacp_ou_max_cpp(
    const int n_paths,
    const int n_grid,
    const double L_eta,
    const int n_threads,
    const int seed) {
  if (n_paths < 1 || n_grid < 2 || n_grid % 2 != 0 ||
      !std::isfinite(L_eta) || !(L_eta > 0.0) || n_threads < 1) {
    Rcpp::stop(
      "Require positive paths/threads/L_eta and a positive even n_grid.");
  }
  const double delta = 2.0 * L_eta / static_cast<double>(n_grid);
  const double autoregression = std::exp(-delta);
  const double innovation_sd = std::sqrt(-std::expm1(-2.0 * delta));
  Rcpp::NumericMatrix maxima(n_paths, 2);
  double* coarse = maxima.begin();
  double* fine = coarse + n_paths;

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int path = 0; path < n_paths; ++path) {
    const std::uint64_t path_seed = splitmix64(
      static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed)) ^
      (0x94d049bb133111ebULL * static_cast<std::uint64_t>(path + 1)));
    std::mt19937_64 rng(path_seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    double value = normal(rng);
    double fine_max = value;
    double coarse_max = value;
    for (int step = 1; step <= n_grid; ++step) {
      value = autoregression * value + innovation_sd * normal(rng);
      if (value > fine_max) fine_max = value;
      if (step % 2 == 0 && value > coarse_max) coarse_max = value;
    }
    coarse[path] = coarse_max;
    fine[path] = fine_max;
  }
  Rcpp::colnames(maxima) =
    Rcpp::CharacterVector::create("coarse", "fine");
  return maxima;
}


// Fast, exactly equivalent implementation of HDDchangepoint::dist1().
// [[Rcpp::export]]
Rcpp::NumericMatrix hdzacp_hdd_distance_matrix_cpp(
    const Rcpp::NumericMatrix data) {
  const int n = data.nrow();
  const int p = data.ncol();
  if (n < 3 || p < 1) {
    Rcpp::stop("HDD distance requires n >= 3 and p >= 1.");
  }
  validate_finite_matrix(data, "data");

  std::vector<double> row_mean(n, 0.0);
  std::vector<double> row_sd_population(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) row_mean[i] += data(i, j);
    row_mean[i] /= static_cast<double>(p);
    for (int j = 0; j < p; ++j) {
      const double centered = data(i, j) - row_mean[i];
      row_sd_population[i] += centered * centered;
    }
    row_sd_population[i] =
      std::sqrt(row_sd_population[i] / static_cast<double>(p));
  }

  Rcpp::NumericMatrix feature_distance(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double difference_mean = row_mean[i] - row_mean[j];
      const double difference_sd =
        row_sd_population[i] - row_sd_population[j];
      const double value = std::sqrt(
        difference_mean * difference_mean + difference_sd * difference_sd);
      feature_distance(i, j) = value;
      feature_distance(j, i) = value;
    }
  }

  Rcpp::NumericMatrix answer(n, n);
  const double scale = 1.0 / static_cast<double>(n - 2);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double total = 0.0;
      for (int z = 0; z < n; ++z) {
        if (z == i || z == j) continue;
        total += std::abs(feature_distance(i, z) - feature_distance(j, z));
      }
      answer(i, j) = scale * total;
      answer(j, i) = answer(i, j);
    }
  }
  return answer;
}


// Fast counterpart of HDDchangepoint::test_statistic() after its distance
// matrix has been computed. order is a one-based whole-observation ordering.
// [[Rcpp::export]]
Rcpp::NumericVector hdzacp_hdd_stat_from_distance_cpp(
    const Rcpp::NumericMatrix distance_matrix,
    const Rcpp::IntegerVector order) {
  const int n = distance_matrix.nrow();
  if (n < 3 || distance_matrix.ncol() != n) {
    Rcpp::stop("distance_matrix must be square with at least three rows.");
  }
  validate_finite_matrix(distance_matrix, "distance_matrix");
  const std::vector<int> order0 = validate_order(order, n);
  const std::pair<int, double> statistic =
    hdd_stat_from_distance(distance_matrix, order0);
  return Rcpp::NumericVector::create(
    Rcpp::_ ["changepoint"] = statistic.first,
    Rcpp::_ ["statistic"] = statistic.second
  );
}


// Bernoulli endpoint dynamic program used in the subcritical analytic MAX
// normalization.  This is compiled with the package rather than at run time.
// [[Rcpp::export]]
Rcpp::List hdzacp_max_endpoint_dp_cpp(
    const double t,
    const double theta,
    const Rcpp::IntegerVector checkpoints) {
  if (!std::isfinite(t) || !(t > 0.0 && t < 1.0)) {
    Rcpp::stop("t must lie strictly between zero and one.");
  }
  if (!std::isfinite(theta) || !(theta > 0.0 && theta < 0.25)) {
    Rcpp::stop("theta must lie strictly between zero and 1/4.");
  }
  if (checkpoints.size() < 1) {
    Rcpp::stop("checkpoints must be a non-empty integer vector.");
  }
  for (int i = 0; i < checkpoints.size(); ++i) {
    if (checkpoints[i] < 1 ||
        (i > 0 && checkpoints[i] <= checkpoints[i - 1])) {
      Rcpp::stop("checkpoints must be positive and strictly increasing.");
    }
  }

  const int maximum = checkpoints[checkpoints.size() - 1];
  int checkpoint_index = 0;
  std::vector<double> log_minus(1, 0.0);
  std::vector<double> log_plus(1, 0.0);
  Rcpp::List output(checkpoints.size());
  const double log_failure = std::log1p(-t);
  const double log_success = std::log(t);

  for (int m = 1; m <= maximum; ++m) {
    std::vector<double> next_minus(m + 1);
    std::vector<double> next_plus(m + 1);
    for (int successes = 0; successes <= m; ++successes) {
      const double x = static_cast<double>(successes) / m;
      double relative_entropy = 0.0;
      if (successes > 0) relative_entropy += x * std::log(x / t);
      if (successes < m) {
        relative_entropy +=
          (1.0 - x) * std::log((1.0 - x) / (1.0 - t));
      }

      for (int component = 0; component < 2; ++component) {
        const std::vector<double>& previous =
          component == 0 ? log_minus : log_plus;
        double log_probability;
        if (successes == 0) {
          log_probability = log_failure + previous[0];
        } else if (successes == m) {
          log_probability = log_success + previous[m - 1];
        } else {
          const double from_failure = log_failure + previous[successes];
          const double from_success = log_success + previous[successes - 1];
          log_probability = std::max(from_failure, from_success) +
            std::log1p(std::exp(-std::abs(from_failure - from_success)));
        }
        const double shift = component == 0 ? -0.5 : 0.5;
        const double value = log_probability +
          theta * static_cast<double>(m) /
          (static_cast<double>(m) + shift) * relative_entropy;
        if (component == 0) {
          next_minus[successes] = value;
        } else {
          next_plus[successes] = value;
        }
      }
    }
    log_minus.swap(next_minus);
    log_plus.swap(next_plus);

    if (checkpoint_index < checkpoints.size() &&
        m == checkpoints[checkpoint_index]) {
      output[checkpoint_index] = Rcpp::List::create(
        Rcpp::_ ["m"] = m,
        Rcpp::_ ["minus"] = Rcpp::wrap(log_minus),
        Rcpp::_ ["plus"] = Rcpp::wrap(log_plus)
      );
      ++checkpoint_index;
    }
  }
  return output;
}

