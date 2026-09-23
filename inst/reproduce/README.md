# Article reproduction scripts

These scripts reproduce the simulations and empirical analyses in
*Nonparametric Change-Point Detection and Inference for High-Dimensional
Distributions*. They call exported `HDZACP` functions; no method implementation
is duplicated here.

The package must be installed before the scripts are run. The seven-method
comparisons additionally require:

```r
install.packages(c("ade4", "ecp", "gSeg", "remotes"))
remotes::install_github("zhangxiany-tamu/KDist")
remotes::install_github("rezadrikvandi/HDDchangepoint")
```

KDist is provided by `zhangxiany-tamu/KDist`; HDD is provided by
`rezadrikvandi/HDDchangepoint`; E-Divisive is called through `ecp`; and gSeg is
called through `gSeg` with an MST from `ade4`. These are optional external
implementations and are not copied into this package.

## Full and smoke modes

Running a script without `--quick` selects the full article defaults. The main
testing experiments then use 1,000 Monte Carlo replications and 199
permutations; the multiple-change study uses 100 replications and 99 local
permutations. These jobs can take many hours.

Use `--quick` only to verify installation, compiled code, optional dependencies,
and output paths:

```sh
Rscript inst/reproduce/00_run_all.R --quick
Rscript inst/reproduce/01_size.R --quick
```

The smoke results are not estimates of the published size, power, or location
performance.

All scripts accept `--output-dir=PATH`. When the called function supports them,
`--cores=N` and `--seed=N` are forwarded. The defaults remain those recorded in
the package's article-reproduction functions.

## Script map

| Script | Exported function |
|---|---|
| `00_run_all.R` | `reproduce_hdzacp_article()` |
| `01_size.R` | `reproduce_hdzacp_size()` |
| `02_power.R` | `reproduce_hdzacp_power()` |
| `03_location.R` | `reproduce_hdzacp_location()` |
| `04_multiple.R` | `reproduce_hdzacp_multiple()` |
| `05_analytic_size.R` | `reproduce_hdzacp_analytic_size()` |
| `06_braincloud.R` | `reproduce_hdzacp_braincloud()` |
| `07_gas_sensor.R` | `reproduce_hdzacp_gas_sensor()` |

## Empirical data

No empirical data are stored in the package or GitHub repository. The official
sources are:

- BrainCloud GSE30272:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30272
- BrainCloud normalized Series Matrix:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30272/matrix/GSE30272_series_matrix.txt.gz
- UCI Gas Sensor Array Drift at Different Concentrations:
  https://archive.ics.uci.edu/dataset/270/gas%2Bsensor%2Barray%2Bdrift%2Bdataset%2Bat%2Bdifferent%2Bconcentrations
- UCI Gas Sensor archive:
  https://archive.ics.uci.edu/static/public/270/gas+sensor+array+drift+dataset+at+different+concentrations.zip

Pass a local BrainCloud Series Matrix with `--data-file=PATH` to
`06_braincloud.R`, and pass a local Gas Sensor ZIP archive with
`--data-zip=PATH` to `07_gas_sensor.R`. If a reproduction function supports
official on-demand download and no path is supplied, it will print the source
URL before downloading to the user-selected output or cache directory. It must
never write raw data into the installed package.

Examples:

```sh
Rscript inst/reproduce/06_braincloud.R \
  --data-file=data/GSE30272_series_matrix.txt.gz \
  --output-dir=results/braincloud

Rscript inst/reproduce/07_gas_sensor.R \
  --data-zip=data/gas_sensor_array_drift.zip \
  --output-dir=results/gas_sensor
```

Each output directory contains the effective configuration and seeds. Preserve
those files with any numerical result used in a paper or supplement.
