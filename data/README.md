# Input data

No data are redistributed in this repository. The single dataset used in the
paper is public and is downloaded from NCBI GEO.

## GSE96583 — 8 donors, control and interferon-β

Kang HM, Subramaniam M, Targ S, et al. Multiplexed droplet single-cell
RNA-sequencing using natural genetic variation. *Nat Biotechnol* 2018;36:89–94.

Open access, no controlled-access component.

| File | What it is |
|---|---|
| `GSE96583_RAW.tar` | the per-library count matrices (batch 2 uses `GSM2560248` and `GSM2560249`) |
| `GSE96583_batch2.genes.tsv.gz` | gene identifiers for batch 2 |
| `GSE96583_batch2.total.tsne.df.tsv.gz` | per-cell annotation carrying the demuxlet donor label and the condition |

All three come from:

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/
```

Download and unpack into `$EFFN_ROOT/GSE96583` as described in the top-level
README. `01_build_gse96583.R` expects `GSE96583_RAW.tar` to be extracted into a
subdirectory named `extract/`.

## What is included instead

`effn/inst/extdata/gse96583_small.rds` is an 800-cell subsample (50 cells per
donor and condition) of the processed object. It ships with the package so that
the headline medians in the paper can be reproduced without downloading
anything:

```r
library(effn)
s <- readRDS(system.file("extdata", "gse96583_small.rds", package = "effn"))
print(effn(s$counts, s$group, s$donor, n_boot = 200, n_tech = 40))
```

## A second dataset that was not used

GSE203191 (endometriosis, eleven samples) was checked and set aside: its series
matrix annotates tissue of origin only, with no case or control arm, so it
supports a between-tissue comparison but not the donor-level case–control
contrast this package is built for.
