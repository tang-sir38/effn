# effn

**Effective sample size and minimum detectable effect for single-cell comparisons.**

`effn` takes a single-cell count matrix plus donor labels and returns, per gene:

1. the **effective donor count** that survives the sequencing-depth penalty;
2. the **minimum detectable effect** (MDE) at a stated power, using exact *t* quantiles at the design degrees of freedom;
3. the **donors still required** to support a target fold change, solved by bisection.

It answers the question a laboratory actually asks at analysis time — *"how many of the genes on this list can I carry into a second study?"* — from the data in hand, rather than from dispersions assumed in a protocol.

This repository accompanies the Application Note:

> Tang C, Zhan L, Zheng Y, Wu F, Tan H, Lin L. *effn: how many donors and how much depth does a single-cell comparison need?* Bioinformatics (submitted).

---

## What is in this repository

```
effn/         the R package source (install with R CMD INSTALL, or R CMD build then install)
analysis/     the numbered R scripts (00-13) that produced every table and figure in the paper
data/         where the public input data go; nothing is redistributed here
```

The package is base R plus `Matrix`. It has no compiled code.

---

## Install

```bash
R CMD INSTALL effn_1.1.0.tar.gz     # or: R CMD INSTALL effn/
```

Requirements: R >= 4.0, `Matrix`. Optional adapters use `SingleCellExperiment` or `Seurat`, neither of which is required for any number in the paper.

---

## Quick start (no download needed)

An 800-cell subsample of the worked dataset ships inside the package. This reproduces the headline medians in seconds:

```r
library(effn)
f <- system.file("extdata", "gse96583_small.rds", package = "effn")
s <- readRDS(f)
fit <- effn(s$counts, s$group, s$donor, n_boot = 200, n_tech = 40)
print(fit)
```

---

## Full pipeline

The scripts expect a working directory that holds the downloaded GEO files and an output folder. Set it once:

```bash
export EFFN_ROOT=/path/to/work      # default: /root/effn
mkdir -p "$EFFN_ROOT/out" "$EFFN_ROOT/GSE96583"
```

`EFFN_OUT` overrides the output folder alone if you want to split inputs from outputs.

**1. Get the data** — GSE96583 (Kang et al., *Nat Biotechnol* 2018), batch 2:

```bash
BASE=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl
cd "$EFFN_ROOT/GSE96583"
wget $BASE/GSE96583_RAW.tar
wget $BASE/GSE96583_batch2.genes.tsv.gz
tar -xf GSE96583_RAW.tar -C .
```

The per-cell demuxlet annotation that carries the donor label (`GSE96583_batch2.total.tsne.df.tsv.gz`) is in the same series and is required by `01_build_gse96583.R`.

**2. Run the numbered scripts in order:**

| Script | What it does | Needs |
|---|---|---|
| `00_definitions.R` | reference implementation of the estimator, for auditing the algebra; no data | — |
| `01_build_gse96583.R` | assembles the 24,679-cell matrix, attaches donor labels, exports the per-cell table | GEO files |
| `02_effn_main.R` | runs `effn()`, then re-runs it with the same seed to prove determinism | stage 01 |
| `03_depth_and_export.R` | the depth grid, within-gene exchange rate, two-scale comparison | stage 02 |
| `04_sim_calibration.R` | the two calibration arms and the per-effect-size rejection grid | — (simulated) |
| `05_subsample_recheck.R` | rebuilds the 800-cell subsample and re-runs `effn()` on it | stage 02 |
| `06_figure1.R` | Figure 1 | stage 02, 03 |
| `07_unpaired_check.R` | the same counts treated as an unpaired design | stage 02 |
| `08_hierarchicell_compare.R` | Hierarchicell 1.0.0 on the same counts, its own defaults | stage 02, Hierarchicell |
| `09_hc_power_curve.R` | the donor grid for that comparison | stage 08 |
| `10_small_dataset.R` | builds the 800-cell subsample that ships with the package; prints its md5 | stage 02 |
| `11_small_reproduce.R` | reproduces the subsample medians from the shipped object | package only |
| `12_small_ci.R` | jackknife against bootstrap intervals on the subsample | stage 10, 11 |
| `13_quantile_gap.R` | how far the normal approximation understates the MDE, by design df | — |

A convenience runner:

```bash
bash analysis/run_all.sh          # stages 01-04, 07, 10; needs the GEO files first
```

**3. Nothing here requires a cluster.** The 24,679-cell run takes minutes on one core; the 800-cell subsample takes seconds.

---

## Environment used in the paper

| Item | Value |
|---|---|
| Package | effn 1.1.0, GPL-3, no compiled code |
| R (24,679-cell run) | 4.6.1, Linux |
| R (package check) | 4.6.0, Windows 11 x64 |
| Matrix | 1.7.6 |
| `R CMD check` | `--ignore-vignettes --no-manual`, Status: OK |
| Test suite | 16 blocks, base R only, all pass |

The vignette requires LaTeX to build; the Sweave source was executed separately and the medians it prints match the paper.

---

## Data

All input data are public and are **not** redistributed here:

- **GSE96583** — Kang HM, Subramaniam M, Targ S, et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nat Biotechnol* 2018;36:89–94. NCBI GEO, open access.

No controlled-access data were used. This is a secondary analysis of de-identified public data; no ethics approval was required.

An 800-cell subsample (`effn/inst/extdata/gse96583_small.rds`) is included so the headline numbers can be reproduced without any download.

---

## Citation

If you use `effn`, please cite the Application Note above. A snapshot of this
repository (version 1.1.0) is archived on Zenodo; the DOI is given in the paper
and in the release notes of this repository.

`CITATION.cff` in this repository carries the same information in machine-readable form.

---

## Licence

GPL-3. See `LICENSE`.
