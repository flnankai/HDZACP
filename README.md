# HDZACP

`HDZACP` implements the procedures in *Nonparametric Change-Point
Detection and Inference for High-Dimensional Distributions*. The package is
designed for an ordered numeric matrix whose rows are observations and whose
columns are coordinates. It tests for and locates changes in coordinatewise
marginal distributions without imposing a parametric marginal model or
requiring marginal moments.

The three proposed procedures use the same standardized coordinatewise rank
scan and differ only in how evidence is aggregated:

- **ZAS** is the sum scan and targets changes spread over many coordinates.
- **ZAM** is the maximum scan and targets changes concentrated in a few
  coordinates.
- **ZAC** combines the ZAS and ZAM permutation evidence and adapts to unknown
  sparsity.

The default calibration is whole-vector permutation. Each permutation moves an
entire row, rather than permuting coordinates separately, and therefore
preserves contemporaneous dependence among coordinates. Its finite-sample
validity requires exchangeability of the ordered observation vectors under the
null. The methods target marginal distribution changes; a change confined to
the copula while every marginal distribution remains fixed can be missed.

## Installation

Install the development version from GitHub with only the core dependencies:

```r
install.packages("remotes")
remotes::install_github(
  "flnankai/HDZACP",
  dependencies = c("Depends", "Imports")
)
```

The comparison methods are optional. The CRAN dependencies can be installed
with

```r
install.packages(c("ade4", "ecp", "gSeg"))
```

and the two packages currently distributed through GitHub can be installed
with

```r
remotes::install_github("zhangxiany-tamu/KDist")
remotes::install_github("rezadrikvandi/HDDchangepoint")
```

These packages remain the work of their respective authors. `HDZACP` provides
documented adapters so that their outputs can be compared on a common input
matrix; it does not present the competing methods as part of the proposed
methodology.

## Basic use

```r
library(HDZACP)

set.seed(42)
n <- 160
p <- 80
x <- matrix(rnorm(n * p), nrow = n, ncol = p)
x[81:n, 1:5] <- matrix(rt((n - 80) * 5, df = 3) / sqrt(3), n - 80, 5)

fit <- hdzacp_test(
  x,
  methods = c("ZAC", "ZAM", "ZAS"),
  calibration = "permutation",
  eta = 0.1,
  permutations = 199,
  seed = 1
)
fit$results
```

`estimate` in the result is the split-after index: an estimate of 80 divides
rows `1:80` from rows `81:n`. The component functions are convenient when only
one proposed procedure is needed:

```r
zac_test(x, permutations = 199, seed = 1)
zam_test(x, permutations = 199, seed = 1)
zas_test(x, permutations = 199, seed = 1)
```

For ZAC, the reported location is inherited from the component with the
smaller permutation p-value; a tie is resolved in favor of ZAS. ZAC's default
permutation p-value is obtained symmetrically over the same observed-plus-
permuted orbit as the two component p-values. It is not formed by simply
inserting two Monte Carlo p-values into an independent theoretical Cauchy
formula.

Multiple changes are estimated by the paper's componentwise and
Cauchy-adaptive wild binary segmentation:

```r
wbs_fit <- hdzacp_wbs(
  x,
  methods = c("ZAC", "ZAM", "ZAS"),
  calibration = "permutation",
  n_intervals = 100,
  min_interval = 30,
  permutations = 99,
  seed = 2
)
wbs_fit$changepoints
```

The same random interval collection is used for the three proposed methods.
For an exact reproduction of the article's multiple-change simulation, use
`n_intervals = 50`, `min_interval = 40`, `permutations = 99`, and
`stopping = "unadjusted"`; the full interval is included in addition to the
50 random intervals.

## Calibration and ties

Permutation calibration is recommended for routine use. With `B`
permutations, the smallest attainable p-value is `1 / (B + 1)`; therefore `B`
must be large enough for the intended significance level, particularly inside
multiple-change procedures. `calibration = "asymptotic"` exposes the analytic
laws studied in the article. The implementation reports unsupported
dimension-growth branches instead of silently replacing them with a different
limit; in particular, the current analytic MAX implementation is restricted to
the subcritical branch covered by its numerical constants.

The rank theory assumes continuous margins. By default, tied values trigger a
clear error. For permutation inference, `ties = "random"` may be used with a
fixed `tie_seed`; the same randomized tie-breaking rule is then used throughout
the orbit. Random tie breaking is not offered as an analytic calibration.

## Comparison methods

`hdzacp_compare()` standardizes the interfaces and returned columns for the
four methods used in the article:

| Label | External implementation used |
|---|---|
| KDist | `KDist::kcpd_single()` and `KDist::kcpd_wbs()` |
| HDD | `HDDchangepoint::test_single_changepoint()` and `HDDchangepoint::multiple_changepoint_detection_wbs()` |
| E-Divisive | `ecp::e.divisive()` |
| gSeg | `ade4::mstree(stats::dist(x))` followed by `gSeg::gseg1()` |

```r
comparison <- hdzacp_compare(x, permutations = 199, seed = 3)
comparison
```

The article reproduction preserves any article-specific calibration choices
and records them in the output. In particular, its accelerated HDD calculation
uses a statistic numerically equivalent to the `HDDchangepoint` distance and
scan statistic together with the article's conditional plus-one permutation
calibration; this is explicitly distinguished from the external package's
default call.

## Reproducing the article

Every major experiment has a callable function. `quick = TRUE` reduces the
number of Monte Carlo replications and permutations for a smoke test; omitting
it, or setting `quick = FALSE`, uses the article settings.

```r
reproduce_hdzacp_size(quick = TRUE, output_dir = "results/size_smoke")
reproduce_hdzacp_power(quick = TRUE, output_dir = "results/power_smoke")
reproduce_hdzacp_location(quick = TRUE, output_dir = "results/location_smoke")
reproduce_hdzacp_multiple(quick = TRUE, output_dir = "results/multiple_smoke")
reproduce_hdzacp_analytic_size(
  quick = TRUE,
  output_dir = "results/analytic_size_smoke"
)

# Run the complete workflow. This is computationally expensive.
reproduce_hdzacp_article(quick = FALSE, output_dir = "results/article")
```

Command-line entry points are installed under `inst/reproduce/`. For example,
after cloning the repository:

```sh
Rscript inst/reproduce/00_run_all.R --quick
Rscript inst/reproduce/00_run_all.R --output-dir=results/article
```

The full testing experiments use 1,000 Monte Carlo replications and 199
whole-vector permutations. The full multiple-change experiment uses 100
replications and 99 local permutations. Runtime can therefore be substantial;
the smoke mode is intended only to verify installation and execution, not to
reproduce numerical tables.

## Empirical data

No empirical data are included in the package or repository.
`hdzacp_data_urls()` returns the official landing and download URLs. Users must
review and comply with the data providers' terms and cite the original data
sources.

### BrainCloud

- GEO accession: [GSE30272](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30272)
- Normalized Series Matrix:
  [GSE30272_series_matrix.txt.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30272/matrix/GSE30272_series_matrix.txt.gz)

The article retains subjects aged at least 20 years, regresses each probe on an
intercept and the two annotated surrogate variables without adjusting for age,
removes incomplete and zero-variance probes, retains the 2,000 largest residual
variances, applies featurewise median/MAD standardization, and orders the 148
subjects by age. The analysis uses `eta = 0.1` and 999 permutations.

```r
urls <- hdzacp_data_urls()
reproduce_hdzacp_braincloud(
  data_file = "path/to/GSE30272_series_matrix.txt.gz",
  output_dir = "results/braincloud"
)
```

### Gas Sensor Array Drift

- UCI landing page:
  [Gas Sensor Array Drift at Different Concentrations](https://archive.ics.uci.edu/dataset/270/gas%2Bsensor%2Barray%2Bdrift%2Bdataset%2Bat%2Bdifferent%2Bconcentrations)
- Official archive:
  [ZIP download](https://archive.ics.uci.edu/static/public/270/gas+sensor+array+drift+dataset+at+different+concentrations.zip)
- DOI: [10.24432/C5MK6M](https://doi.org/10.24432/C5MK6M)

The article retains class-1 ethanol measurements, adjusts every feature for log
concentration through a four-degree-of-freedom natural-spline location and
scale regression, and applies final median/MAD standardization. At most 50
observations are sampled without replacement from each temporal batch, giving
batch sizes `(50, 50, 50, 50, 28, 50, 50, 30, 50, 50)` and `n = 458`,
`p = 128`. Within-batch rows are randomized with a fixed seed while the ten
batches remain in temporal order. The interpretation is descriptive: exact
within-batch acquisition times are unavailable and concentration support is
confounded with batch.

```r
reproduce_hdzacp_gas_sensor(
  data_zip = "path/to/gas_sensor_array_drift.zip",
  output_dir = "results/gas_sensor"
)
```

See `vignette("methods-and-calibration", package = "HDZACP")` and
`vignette("article-reproduction", package = "HDZACP")` for the statistical
interfaces and the complete reproduction map.

## License

`HDZACP` is released under GPL-3 or later. External comparison packages and
empirical data retain their own licenses and terms.

