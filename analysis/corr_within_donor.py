# -*- coding: utf-8 -*-
"""Within-donor correlation implied by the two analyses of GSE96583.

The paired run reports sigma_obs as the SD of the within-donor difference; the
unpaired run reports sigma_obs as the SD of a single donor pseudobulk. If both
conditions have variance s2 and correlate with r, the difference has variance
2*s2*(1 - r), so r = 1 - sd_diff^2 / (2 * sd_single^2).

Writes within_donor_correlation.txt into $EFFN_OUT and is the source for the
ledger row that quotes the median in the manuscript.
"""
import csv
import gzip
import os
import statistics as st

# the numbered scripts write flat into $EFFN_OUT; this helper reads the same folder
# this file is analysis/aux/<name>.py, so three levels up is the repository root
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
GSE = os.environ.get("EFFN_OUT", os.path.join(os.environ.get("EFFN_ROOT", _REPO), "out"))


def sigma_obs(name):
    out = {}
    with gzip.open(os.path.join(GSE, name), "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["status"] == "ok":
                out[r["gene"]] = float(r["sigma_obs"])
    return out


paired = sigma_obs("gse96583_per_gene.tsv.gz")
unpaired = sigma_obs("gse96583_unpaired_per_gene.tsv.gz")
common = [g for g in paired if g in unpaired]
r = sorted(1.0 - paired[g] ** 2 / (2.0 * unpaired[g] ** 2) for g in common)


def q(p):
    i = min(len(r) - 1, max(0, int(round(p * (len(r) - 1)))))
    return r[i]


lines = [
    "genes paired_ok=%d unpaired_ok=%d both=%d" % (len(paired), len(unpaired), len(common)),
    "median_r=%.6f" % st.median(r),
    "q1_r=%.6f" % q(0.25),
    "q3_r=%.6f" % q(0.75),
    "share_below_0.2=%.6f" % (sum(1 for x in r if x < 0.2) / len(r)),
    "median_sd_diff_paired=%.6f" % st.median([paired[g] for g in common]),
    "median_sd_single_unpaired=%.6f" % st.median([unpaired[g] for g in common]),
]
out = os.path.join(GSE, "within_donor_correlation.txt")
open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("wrote", out)
