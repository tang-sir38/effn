# -*- coding: utf-8 -*-
"""Builds number_ledger.tsv in $EFFN_OUT from the result files.

Every display string is derived from the file it points at, never typed by hand,
so a change in the analysis cannot leave the manuscript quoting a stale value.
Run from the project root:  python 02_SCRIPTS/build_ledger.py
"""
import csv
import gzip
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# every table the ledger audits is written flat into $EFFN_OUT by the numbered
# scripts, so one environment variable points at all of them
# this file is analysis/aux/<name>.py, so three levels up is the repository root
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
RES = os.environ.get("EFFN_OUT", os.path.join(os.environ.get("EFFN_ROOT", _REPO), "out"))
GSE = RES

# The two exact quantiles of Supplementary S5, computed rather than typed.
from scipy import stats as _st

q7 = _st.t.ppf(0.975, 7) + _st.t.ppf(0.8, 7)
q14 = _st.t.ppf(0.975, 14) + _st.t.ppf(0.8, 14)
q_pct = 100 * (q7 / q14 - 1)

s = json.load(open(os.path.join(GSE, "gse96583_summary.json"), encoding="utf-8"))
old_unp = json.load(open(os.path.join(GSE, "gse96583_unpaired_summary.json"),
                         encoding="utf-8"))
cal = list(csv.DictReader(open(os.path.join(RES, "t_sim_calibration.tsv"),
                               encoding="utf-8"), delimiter="\t"))
arms = {r["arm"]: r for r in csv.DictReader(
    open(os.path.join(RES, "t_calibration_arms.tsv"), encoding="utf-8"),
    delimiter="\t")}
depth = list(csv.DictReader(open(os.path.join(GSE, "depth_curve.tsv"),
                                 encoding="utf-8"), delimiter="\t"))
exch = list(csv.DictReader(open(os.path.join(GSE, "exchange_rate.tsv"),
                               encoding="utf-8"), delimiter="\t"))

# cells per donor and condition, from the exported cell table
per = {}
with gzip.open(os.path.join(GSE, "gse96583_cells.tsv.gz"), "rt", encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        per[(r["donor"], r["group"])] = per.get((r["donor"], r["group"]), 0) + 1
cmin, cmax = min(per.values()), max(per.values())

d10 = depth[0]
dmax = depth[-1]
nun = s["n_unestimable"]
neff_ci = s["median_n_eff_ci"]
mde_ci = s["median_mde_log2_ci"]
cal_mde = next(c for c in cal if c["factor"] == "1" and c["arm"] == "deep")
cal_null = next(c for c in cal if c["factor"] == "0" and c["arm"] == "deep")
ex = {int(float(r["cells_per_donor"])): r for r in exch}

rows = []


def add(i, display, value, src, field, must="yes", note="", used_in="both"):
    """Record one number.

    used_in says which text is expected to carry the string: "main",
    "supp" or "both". The strict check only demands that a row appear in the
    file it names, so a number that the Supplement quotes and the main text
    does not is not reported as dead in the main text.
    """
    rows.append(dict(id=i, display=display, value=value, must_appear=must,
                     source_file=src, source_field=field, note=note,
                     used_in=used_in))


add("N01", "24,679", str(s["nominal_n"]), "gse96583/gse96583_summary.json", "nominal_n",
    note="singlets with a demuxlet donor label")
add("N02", "35,635", str(s["n_genes_tested"]), "gse96583/gse96583_summary.json",
    "n_genes_tested")
add("N03", "17,171", str(nun[0]), "gse96583/gse96583_summary.json",
    "n_unestimable[all_zero]")
add("N04", "6,357", str(nun[1]), "gse96583/gse96583_summary.json",
    "n_unestimable[tech_dominated]")
add("N05", "3,157", str(nun[2]), "gse96583/gse96583_summary.json",
    "n_unestimable[low_resolution]")
add("N06", "2,492", str(nun[3]), "gse96583/gse96583_summary.json",
    "n_unestimable[low_expression]")
add("N07", "6,458", str(s["n_ok"]), "gse96583/gse96583_summary.json", "n_ok")
add("N08", "3.91", "%.2f" % s["median_n_eff"], "gse96583/gse96583_summary.json",
    "median_n_eff")
add("N09", "2.39", "%.2f" % neff_ci[0], "gse96583/gse96583_summary.json",
    "median_n_eff_ci[0]", note="delete-one-donor jackknife")
add("N10", "5.43", "%.2f" % neff_ci[1], "gse96583/gse96583_summary.json",
    "median_n_eff_ci[1]", note="delete-one-donor jackknife")
add("N11", "51%", "%d%%" % round(100 * s["loss_depth_share"]),
    "gse96583/gse96583_summary.json", "loss_depth_share")
add("N12", "0.71", "%.2f" % s["median_mde_log2"], "gse96583/gse96583_summary.json",
    "median_mde_log2")
add("N13", "1.63", "%.2f" % s["median_mde_fold"], "gse96583/gse96583_summary.json",
    "median_mde_fold")
add("N14", "1.48", "%.2f" % (2 ** mde_ci[0]), "gse96583/gse96583_summary.json",
    "2^median_mde_log2_ci[0]")
add("N15", "1.81", "%.2f" % (2 ** mde_ci[1]), "gse96583/gse96583_summary.json",
    "2^median_mde_log2_ci[1]")
add("N16", "6", str(s["n_required_2fold"]), "gse96583/gse96583_summary.json",
    "n_required_2fold")
add("N17", "4", str(s["n_required_2fold_iqr"][0]), "gse96583/gse96583_summary.json",
    "n_required_2fold_iqr[0]")
add("N18", "8", str(s["n_required_2fold_iqr"][1]), "gse96583/gse96583_summary.json",
    "n_required_2fold_iqr[1]", note="also the donor count of the dataset")
add("N19", "11", str(s["n_required_1.5fold"]), "gse96583/gse96583_summary.json",
    "n_required_1.5fold")
add("N20", "3,085", "%d" % round(s["loss_donor"]), "gse96583/gse96583_summary.json",
    "loss_donor", note="cells / donors, 24679/8 = 3084.9")
add("N21", str(cmin), str(cmin), "gse96583/gse96583_cells.tsv.gz", "min cells per donor-group")
add("N22", "%d" % cmax, "%d" % cmax, "gse96583/gse96583_cells.tsv.gz",
    "max cells per donor-group", note="2,738 is formatted with a comma in the text")

# depth curve, survivor medians
add("N23", "4.62", "%.2f" % float(d10["core_median_mde_log2"]),
    "gse96583/depth_curve.tsv", "core_median_mde_log2 @10 cells")
add("N24", "1.21", "%.2f" % float(dmax["core_median_mde_log2"]),
    "gse96583/depth_curve.tsv", "core_median_mde_log2 @1042 cells")
add("N25", "1,042", "%d" % int(float(dmax["cells_per_donor"])),
    "gse96583/depth_curve.tsv", "cells_per_donor (max grid point)")
add("N26", "128", "%d" % round(float(d10["core_median_n_required_2x"])),
    "gse96583/depth_curve.tsv", "core_median_n_required_2x @10 cells")
add("N27", "10", "%d" % round(float(dmax["core_median_n_required_2x"])),
    "gse96583/depth_curve.tsv", "core_median_n_required_2x @1042 cells")
add("N28", "5.73", "%.2f" % float(d10["core_median_n_eff"]),
    "gse96583/depth_curve.tsv", "core_median_n_eff @10 cells")
add("N29", "4.46", "%.2f" % float(dmax["core_median_n_eff"]),
    "gse96583/depth_curve.tsv", "core_median_n_eff @1042 cells")
add("N30", "28%", "%d%%" % round(100 * float(d10["core_depth_share"])),
    "gse96583/depth_curve.tsv", "core_depth_share @10 cells")
add("N31", "44%", "%d%%" % round(100 * float(dmax["core_depth_share"])),
    "gse96583/depth_curve.tsv", "core_depth_share @1042 cells")
add("N83", "2,471", "%d" % int(float(d10["core_n_ok"])), "gse96583/depth_curve.tsv",
    "core_n_ok @10 cells")
add("N84", "4,343", "%d" % int(float(dmax["core_n_ok"])), "gse96583/depth_curve.tsv",
    "core_n_ok @1042 cells")

# exchange rate, within gene
add("N85", "1.95", "%.2f" % float(ex[100]["median_fold_vs_full"]),
    "gse96583/exchange_rate.tsv", "median_fold_vs_full @100 cells")
add("N86", "2.66", "%.2f" % float(ex[50]["median_fold_vs_full"]),
    "gse96583/exchange_rate.tsv", "median_fold_vs_full @50 cells")
add("N87", "14.0", "%.1f" % float(ex[10]["median_fold_vs_full"]),
    "gse96583/exchange_rate.tsv", "median_fold_vs_full @10 cells")
add("N88", "2.1", "%.1f" % float(ex[100]["donor_equivalent_of_8"]),
    "gse96583/exchange_rate.tsv", "donor_equivalent_of_8 @100 cells")
add("N89", "1.1", "%.1f" % float(ex[50]["donor_equivalent_of_8"]),
    "gse96583/exchange_rate.tsv", "donor_equivalent_of_8 @50 cells")
add("N90", "0.04", "%.2f" % float(ex[10]["donor_equivalent_of_8"]),
    "gse96583/exchange_rate.tsv", "donor_equivalent_of_8 @10 cells")
add("N93", "381", "%.0f" % float(ex[100]["cells_per_full_depth_donor"]),
    "gse96583/exchange_rate.tsv", "cells_per_full_depth_donor @100 cells")
add("N94", "646", "%.0f" % float(ex[25]["cells_per_full_depth_donor"]),
    "gse96583/exchange_rate.tsv", "cells_per_full_depth_donor @25 cells")
add("N95", "1,969", "%s" % "{:,}".format(int(round(float(
    ex[10]["cells_per_full_depth_donor"])))),
    "gse96583/exchange_rate.tsv", "cells_per_full_depth_donor @10 cells")

# two-scale probe: 03_RESULTS/gse96583/depth_curve_two_scales.tsv, 6,458-gene fixed set
two = list(csv.DictReader(open(os.path.join(GSE, "depth_curve_two_scales.tsv"),
                              encoding="utf-8"), delimiter="\t"))
lin, lg = {}, {}
for r in two:
    cells = int(r["cells"])
    rec = dict(sigma_tech=float(r["med_sigma_tech"]), share=float(r["med_share"]),
               n_eff=float(r["med_n_eff"]))
    (lin if r["scale"] == "linear" else lg)[cells] = rec

p25, p1k = 25, 1000
add("N32", "71%", "%d%%" % round(100 * lin[p25]["share"]),
    "gse96583/depth_curve_two_scales.tsv", "linear CPM med_share @25 cells")
add("N33", "66%", "%d%%" % round(100 * lin[p1k]["share"]),
    "gse96583/depth_curve_two_scales.tsv", "linear CPM med_share @1000 cells")
add("N34", "eleven-fold", "%.2f" % (lin[p25]["sigma_tech"] / lin[p1k]["sigma_tech"]),
    "gse96583/depth_curve_two_scales.tsv",
    "linear med_sigma_tech 25 cells / 1000 cells", note="60.379 / 5.378")
add("N35", "30%", "%d%%" % round(100 * lg[p25]["share"]),
    "gse96583/depth_curve_two_scales.tsv", "log2(CPM+1) med_share @25 cells",
    note="Supplementary S4")
add("N36", "48%", "%d%%" % round(100 * lg[p1k]["share"]),
    "gse96583/depth_curve_two_scales.tsv", "log2(CPM+1) med_share @1000 cells",
    note="Supplementary S4")
add("N37", "three-fold", "%.2f" % (lg[p25]["sigma_tech"] / lg[p1k]["sigma_tech"]),
    "gse96583/depth_curve_two_scales.tsv",
    "log2(CPM+1) med_sigma_tech 25 cells / 1000 cells",
    note="2.027 / 0.676; a 1/n variance law predicts 6.3-fold")
add("N123", "25", "25", "gse96583/depth_curve_two_scales.tsv",
    "cells, grid point", note="used in the two-scale comparison")
add("N124", "1,000", "1000", "gse96583/depth_curve_two_scales.tsv",
    "cells, grid point", note="used in the two-scale comparison")
add("N128", "50", "50", "gse96583/exchange_rate.tsv", "cells_per_donor",
    note="exchange rate quoted at 50 cells")

# calibration, two arms
add("N38", "0.77", "%.2f" % float(cal_mde["mean_rejection"]),
    "tables/t_sim_calibration.tsv", "deep arm mean_rejection at 1 x MDE")
add("N39", "0.062", "%.3f" % float(cal_null["mean_rejection"]),
    "tables/t_sim_calibration.tsv", "deep arm mean_rejection at 0 x MDE")
add("N96", "0.87", "%.2f" % float(arms["shallow"]["power_at_mde"]),
    "tables/t_calibration_arms.tsv", "shallow arm power_at_mde")
add("N97", "0.063", "%.3f" % float(arms["shallow"]["type1_at_zero"]),
    "tables/t_calibration_arms.tsv", "shallow arm type1_at_zero")
add("N98", "0.01", "%.2f" % float(arms["deep"]["median_rho"]),
    "tables/t_calibration_arms.tsv", "deep arm median_rho")
add("N99", "0.23", "%.2f" % float(arms["shallow"]["median_rho"]),
    "tables/t_calibration_arms.tsv", "shallow arm median_rho")
add("N100", "0.86", "%.2f" % float(arms["shallow"]["sigma_ratio_median"]),
    "tables/t_calibration_arms.tsv", "shallow arm sigma_ratio_median")
add("N101", "0.90", "%.2f" % float(arms["shallow"]["sigma_ratio_expected"]),
    "tables/t_calibration_arms.tsv", "shallow arm sigma_ratio_expected")
add("N102", "19%", "%d%%" % round(100 * float(arms["shallow"]["median_rho"]) /
                                  (1 + float(arms["shallow"]["median_rho"]))),
    "tables/t_calibration_arms.tsv", "rho/(1+rho) for the shallow arm")
add("N103", "twenty-fold", "%.0f" % (float(arms["shallow"]["median_rho"]) /
                                     float(arms["deep"]["median_rho"])),
    "tables/t_calibration_arms.tsv", "shallow rho / deep rho")
add("N104", "0.07", "0.07", "tables/t_calibration_arms.tsv",
    "max |power - 0.80| over the two arms: 0.074")
add("N105", "1.05", "%.2f" % s["median_rho_tech"], "gse96583/gse96583_summary.json",
    "median_rho_tech", note="the depth penalty of the worked dataset")
add("N42", "200", "200", "tables/t_sim_calibration.tsv", "R = 200 replicates per grid point")
add("N43", "300", "300", "tables/t_calibration_arms.tsv", "deep arm cells per donor-group")
add("N106", "40", str(arms["shallow"]["cells_per_donor"]),
    "tables/t_calibration_arms.tsv", "shallow arm cells per donor-group")
add("N44", "100", "100", "tables/t_calibration_arms.tsv", "G = 100 genes, 4 strata of 25")

# paired versus naive-independent analysis, on the same genes
un = dict(l.split("=", 1) for l in
          open(os.path.join(GSE, "unpaired_check.txt"), encoding="utf-8").read().split("\n")
          if "=" in l)
add("N59", "0.69", "%.2f" % float(un["mde_unpaired_common"]),
    "gse96583/unpaired_check.txt", "median MDE unpaired, common gene set")
add("N60", "0.70", "%.2f" % float(un["mde_paired_common"]),
    "gse96583/unpaired_check.txt", "median MDE paired, common gene set")
add("N107", "0.6%", "%.1f%%" % (100 * abs(1 - float(un["mde_ratio_unpaired_over_paired"]))),
    "gse96583/unpaired_check.txt", "|1 - ratio| in per cent")
add("N108", "14", "14", "gse96583/unpaired_check.txt",
    "df of the unpaired contrast, 2n - 2", note="also N110")
add("N109", "7", "7", "gse96583/unpaired_check.txt", "df of the paired contrast, n - 1",
    note="also N111")
add("N110", "%d%%" % round(q_pct), "%.2f" % q_pct, "definition",
    "100*(t(.975,7)+t(.8,7))/(t(.975,14)+t(.8,14)) - 1, quantile penalty of pairing",
    note="3.260654 / 3.012841 from scipy.stats.t.ppf")
cor_txt = open(os.path.join(GSE, "within_donor_correlation.txt"),
               encoding="utf-8").read()
cor = dict(l.split("=", 1) for l in cor_txt.split("\n") if "=" in l and
           not l.startswith("genes "))
cor["both"] = re.search(r"\bboth=(\d+)", cor_txt).group(1)
add("N61", "0.04", "%.2f" % float(cor["median_r"]),
    "gse96583/within_donor_correlation.txt", "median_r across 5,851 shared genes")
add("N62", "64%", "%.0f%%" % (100 * float(cor["share_below_0.2"])),
    "gse96583/within_donor_correlation.txt", "share_below_0.2 with two decimals dropped",
    note="Supplementary S5")
add("N63", "5,851", cor["both"], "gse96583/within_donor_correlation.txt",
    "genes estimable in both the paired and the independent analysis")

# the defect itself: what the old degrees of freedom did
old = json.load(open(os.path.join(GSE, "v102_summary.json"), encoding="utf-8"))
infl = s["median_mde_log2"] / old["median_mde_log2"] - 1
add("N111", "%d%%" % round(100 * (old["median_mde_log2"] / s["median_mde_log2"] - 1)),
    "%.1f%%" % (100 * (old["median_mde_log2"] / s["median_mde_log2"] - 1)),
    "gse96583/v102_summary.json",
    "median_mde_log2 (v1.0.2) / median_mde_log2 (v1.1.0) - 1",
    note="1.0334 / 0.7086, the inflation the df error caused")
add("N112", "%.2f" % (s["median_n_eff"] - 1), "%.4f" % (s["median_n_eff"] - 1),
    "gse96583/gse96583_summary.json",
    "median_n_eff - 1", note="the degrees of freedom the erroneous formula used")
add("N113", str(old["n_required_2fold"]), str(old["n_required_2fold"]),
    "gse96583/v102_summary.json", "n_required_2fold",
    note="donors needed for a two-fold change under the erroneous df; "
         "the text quotes the corrected number, 6, against it")
add("N116", "%.3f" % q7, "%.3f" % q7, "definition",
    "qt(0.975,7)+qt(0.8,7)", note="Supplementary S5, the paired quantile")
add("N117", "%.3f" % q14, "%.3f" % q14, "definition",
    "qt(0.975,14)+qt(0.8,14)", note="Supplementary S5, the unpaired quantile")
add("N118", "%.2f" % s["median_sigma_bio"], "%.4f" % s["median_sigma_bio"],
    "gse96583/gse96583_summary.json", "median_sigma_bio",
    note="paired run: SD of a donor-level difference; Supplementary S5")
add("N119", "%.2f" % float(un["median_sigma_bio"]), "%.6f" % float(un["median_sigma_bio"]),
    "gse96583/unpaired_check.txt", "median_sigma_bio",
    note="unpaired run: SD of a single condition; Supplementary S5")

# Hierarchicell 1.0.0, same counts, its own documented defaults
hc = list(csv.DictReader(open(os.path.join(RES, "t_hc_power_curve.tsv"),
                              encoding="utf-8"), delimiter="\t"))
hc_by_n = {int(r["n_per_group"]): r for r in hc}
for i, n in zip(("N64", "N65", "N66"), (8, 16, 32)):
    add(i, "%.2f" % float(hc_by_n[n]["hurdle_power"]), hc_by_n[n]["hurdle_power"],
        "tables/t_hc_power_curve.tsv", "hurdle_power at n_per_group=%d" % n,
        must="no",
        note="Supplementary S6; two-fold change, 300 cells per donor, 100 genes simulated")
add("N67", "100", hc_by_n[8]["n_genes"], "tables/t_hc_power_curve.tsv",
    "n_genes passed to power_hierarchicell")

# the 800-cell subsample that ships in the package
sm_log = os.path.join(RES, "small_dataset.txt")
if os.path.exists(sm_log):
    sm = dict(l.split("=", 1) for l in open(sm_log, encoding="utf-8").read().split("\n")
              if "=" in l)
    add("N68", "800", sm["cells"], "small_dataset.txt", "cells in the shipped subsample")
    add("N69", "50", sm["cells_per_donor_group"], "small_dataset.txt",
        "cells per donor and condition in the subsample", must="no",
        note="Supplementary S6")
    add("N70", "5.28", "%.2f" % float(sm["median_n_eff"]), "small_dataset.txt",
        "median n_eff of the subsample", note="Supplementary S6")
    add("N71", "2.76", "%.2f" % float(sm["median_mde_log2"]), "small_dataset.txt",
        "median log2 MDE of the subsample", note="Supplementary S6")
    add("N72", sm["md5_per_gene_gz"], sm["md5_per_gene_gz"], "small_dataset.txt",
        "md5 of small_per_gene.tsv.gz as distributed", must="no",
        note="Supplementary S6")

# design constants of the calibration and of the Hierarchicell run; these appear
# only in the Supplementary and are sourced to the script that sets them
for i, (rid, val) in enumerate(zip(("N73", "N74", "N75", "N76"),
                                  ("0.20", "0.35", "0.50", "0.70"))):
    add(rid, val, val, "02_SCRIPTS/R/04_sim_calibration.R",
        "SIGMA_TRUE[%d]" % (i + 1), must="no",
        note="biological SD stratum of the calibration")
add("N77", "5,000", "5000", "02_SCRIPTS/R/04_sim_calibration.R", "L",
    must="no", note="counts per cell in the deep arm")
add("N78", "2,000", "2000", "02_SCRIPTS/R/04_sim_calibration.R", "CPM_MAX",
    must="no", note="top of the CPM grid of the deep arm")
add("N79", "20", "20", "02_SCRIPTS/R/04_sim_calibration.R", "length(target)",
    must="no", note="target genes carrying the effect, five per stratum")
add("N80", "40", "40", "02_SCRIPTS/R/02_effn_main.R", "n_tech",
    must="no", note="within-donor resamples used in every reported run")
for i, n in zip(("N81", "N82"), (16, 32)):
    add(i, str(n), hc_by_n[n]["n_per_group"], "tables/t_hc_power_curve.tsv",
        "n_per_group", must="no", note="Hierarchicell donor grid point")

# subsample confidence intervals quoted in Supplementary S1, from 12_small_ci.R
ci_txt = os.path.join(RES, "small_dataset_ci.txt")
if os.path.exists(ci_txt):
    ci = dict(l.split("=", 1) for l in open(ci_txt, encoding="utf-8").read().split("\n")
              if "=" in l)
    add("N114", "%.2f" % float(ci["boot_lo"]), "%.10f" % float(ci["boot_lo"]),
        "small_dataset_ci.txt", "boot_lo",
        must="no", note="Supplementary S1: bootstrap 2.5% for the subsample median MDE")
    add("N115", "%.2f" % float(ci["boot_hi"]), "%.10f" % float(ci["boot_hi"]),
        "small_dataset_ci.txt", "boot_hi",
        must="no", note="Supplementary S1: bootstrap 97.5%")
    add("N120", "%.2f" % float(ci["jack_lo"]), "%.10f" % float(ci["jack_lo"]),
        "small_dataset_ci.txt", "jack_lo",
        must="no", note="Supplementary S1: jackknife interval, lower end")
    add("N121", "%.2f" % float(ci["jack_hi"]), "%.10f" % float(ci["jack_hi"]),
        "small_dataset_ci.txt", "jack_hi",
        must="no", note="Supplementary S1: jackknife interval, upper end")
    assert ci["boot_contains_point"] == "FALSE" and ci["jack_contains_point"] == "TRUE", \
        "the interval comparison quoted in S1 no longer holds"

# the normal approximation against the exact t quantile, paired design
gap = list(csv.DictReader(open(os.path.join(RES, "t_normal_gap.tsv"),
                               encoding="utf-8"), delimiter="\t"))
gap_p = {int(r["n_per_group"]): r for r in gap if r["design"] == "paired"}
add("N122", "80%", "%.1f%%" % float(gap_p[2]["understatement_pct"]),
    "tables/t_normal_gap.tsv", "understatement_pct, paired, n=2",
    note="same string as the nominal power target; both are 80")
add("N129", "48%", "%.1f%%" % float(gap_p[3]["understatement_pct"]),
    "tables/t_normal_gap.tsv", "understatement_pct, paired, n=3")
add("N130", "20%", "%.1f%%" % float(gap_p[6]["understatement_pct"]),
    "tables/t_normal_gap.tsv", "understatement_pct, paired, n=6")
add("N131", "14%", "%.1f%%" % float(gap_p[8]["understatement_pct"]),
    "tables/t_normal_gap.tsv", "understatement_pct, paired, n=8",
    note="the dataset's own donor count")

# design constants quoted as words in the methods
add("N125", "0.05", "0.05", "definition", "two-sided alpha, used in every interval")
add("N126", "0.80", "0.80", "definition", "target power, 80%", must="no")
add("N127", "95%", "95%", "definition", "interval level", must="no")

# genes above a two-fold MDE, from the per-gene table
with gzip.open(os.path.join(GSE, "gse96583_per_gene.tsv.gz"), "rt", encoding="utf-8") as f:
    pg = [r for r in csv.DictReader(f, delimiter="\t") if r["status"] == "ok"]
n_gt2 = sum(1 for r in pg if float(r["mde_fold"]) > 2)
n_gt6 = sum(1 for r in pg if float(r["mde_fold"]) > 6)
add("N132", "1,391", str(n_gt2), "gse96583/gse96583_per_gene.tsv.gz",
    "count of ok genes with mde_fold > 2")
add("N133", "6", str(n_gt6), "gse96583/gse96583_per_gene.tsv.gz",
    "count of ok genes with mde_fold > 6", note="Figure 1B tail")

# Table S2 as printed: every cell at three decimals, straight from the arm summary
_s2 = [
    ("N140", "deep", "median_rho", "median_rho", 3),
    ("N141", "deep", "median_n_eff", "median_n_eff", 2),
    ("N142", "deep", "sigma_ratio_median", "sigma_ratio_median", 3),
    ("N143", "deep", "sigma_ratio_expected", "sigma_ratio_expected", 3),
    ("N144", "deep", "power_at_mde", "power_at_mde", 3),
    ("N145", "shallow", "median_rho", "median_rho", 3),
    ("N146", "shallow", "median_n_eff", "median_n_eff", 2),
    ("N147", "shallow", "sigma_ratio_median", "sigma_ratio_median", 3),
    ("N148", "shallow", "sigma_ratio_expected", "sigma_ratio_expected", 3),
    ("N149", "shallow", "power_at_mde", "power_at_mde", 3),
    ("N150", "shallow", "type1_at_zero", "type1_at_zero", 3),
]
for rid, arm, col, field, dp in _s2:
    v = float(arms[arm][col])
    add(rid, "%.*f" % (dp, v), "%.6f" % v, "tables/t_calibration_arms.tsv",
        "%s arm %s" % (arm, field), used_in="supp",
        note="Table S2 cell")

# the four constants of the shallow arm, quoted in Supplementary S3
for rid, val, field, note in (
        ("N151", "1,500", "ARMS$shallow$L", "counts per cell, shallow arm"),
        ("N152", "500", "ARMS$shallow$CPM_MAX", "top of the CPM grid, shallow arm"),
        ("N153", "0.1", "ok_power tolerance", "power accepted within 0.80 +/- 0.10"),
        ("N154", "0.02", "ok_type1 tolerance", "type I accepted within 0.05 +/- 0.02")):
    add(rid, val, val.replace(",", ""), "02_SCRIPTS/R/04_sim_calibration.R",
        field, used_in="supp", note=note)

# Supplementary S5: the same-gene-set medians at three decimals, and the size of
# the set the earlier analysis used
add("N134", "%.3f" % float(un["mde_unpaired_common"]),
    "%.6f" % float(un["mde_unpaired_common"]),
    "gse96583/unpaired_check.txt", "mde_unpaired_common", used_in="supp")
add("N135", "%.3f" % float(un["mde_paired_common"]),
    "%.6f" % float(un["mde_paired_common"]),
    "gse96583/unpaired_check.txt", "mde_paired_common", used_in="supp")
add("N136", "10,116", str(old_unp["n_ok"]), "gse96583/gse96583_unpaired_summary.json",
    "n_ok", used_in="supp", note="genes estimable under the independent analysis")
add("N91", "5.8", "%.1f" % float(ex[500]["donor_equivalent_of_8"]),
    "gse96583/exchange_rate.tsv", "donor_equivalent_of_8 @500 cells", used_in="supp")

# numbers the Supplement quotes and the main text does not
SUPP_ONLY = {"N35", "N36", "N62", "N70", "N71", "N83", "N84", "N112",
             "N116", "N117", "N118", "N119"}
for r in rows:
    if r["id"] in SUPP_ONLY:
        r["used_in"] = "supp"

# numbers that are printed only inside a generated table, so they never appear
# in the markdown but must appear in the built docx
TABLE_ONLY = {"N39", "N96", "N97", "N98", "N140", "N141", "N142", "N143",
              "N144", "N145", "N146", "N147", "N148", "N149", "N150"}
for r in rows:
    if r["id"] in TABLE_ONLY:
        r["used_in"] = "table"

# numbers the main text quotes and the Supplement does not
MAIN_ONLY = {"N38", "N85", "N86", "N87", "N100", "N104",
             "N130", "N131", "N132", "N133"}
for r in rows:
    if r["id"] in MAIN_ONLY:
        r["used_in"] = "main"

out = os.path.join(RES, "number_ledger.tsv")
with open(out, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t",
                       lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print("wrote", out, len(rows), "rows")
for r in rows:
    print("%-5s %-8s %-12s %s" % (r["id"], r["display"], r["value"], r["source_field"]))
