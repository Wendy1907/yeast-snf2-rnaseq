#!/usr/bin/env python3
"""Compare a pipeline run with the expected results stored in the repository.

Usage: check_results.py --expected expected_results --observed observed --out check_report.txt

'observed' is a flat folder holding this run's gene_counts.tsv, deseq2_results.csv,
*.hisat2.log, the five PNG figures and multiqc_report.html. Standard library only.

What counts as a match:
  * gene counts           identical gene and sample layout; summed absolute difference <= 0.1% of all reads
  * alignment rates       within 0.05 percentage points for every sample
  * DESeq2 results        same genes; fold changes within 0.1% (relative); p-values within 0.01 on the log10 scale
  * figures, MultiQC      present and not empty (they are meant to be looked at, not compared byte by byte)
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

FIGURES = ["ma_plot.png", "pca.png", "sample_distance_heatmap.png",
           "top_genes_heatmap.png", "volcano.png", "multiqc_report.html"]


def read_table(path, sep):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter=sep))
    return rows[0], rows[1:]


def check_counts(expected_dir, observed_dir):
    name = "Gene count matrix"
    exp_path = os.path.join(expected_dir, "counts", "gene_counts.tsv")
    obs_path = os.path.join(observed_dir, "gene_counts.tsv")
    if not os.path.exists(obs_path):
        return False, name, "file missing from this run"
    exp_head, exp_rows = read_table(exp_path, "\t")
    obs_head, obs_rows = read_table(obs_path, "\t")
    if exp_head != obs_head:
        return False, name, "sample columns differ: %s vs %s" % (obs_head[1:], exp_head[1:])
    if [r[0] for r in exp_rows] != [r[0] for r in obs_rows]:
        return False, name, "gene list differs (%d vs %d genes)" % (len(obs_rows), len(exp_rows))
    total = sum(int(v) for r in exp_rows for v in r[1:])
    diff = sum(abs(int(a) - int(b)) for r1, r2 in zip(exp_rows, obs_rows) for a, b in zip(r1[1:], r2[1:]))
    shape = "%d genes x %d samples" % (len(exp_rows), len(exp_head) - 1)
    if diff == 0:
        return True, name, "identical (%s)" % shape
    frac = diff / total
    return frac <= 0.001, name, "%s; summed absolute difference %d reads (%.4f%% of %d)" % (shape, diff, 100 * frac, total)


def alignment_rate(path):
    with open(path) as handle:
        match = re.search(r"([0-9.]+)% overall alignment rate", handle.read())
    return float(match.group(1)) if match else None


def check_alignment(expected_dir, observed_dir):
    name = "HISAT2 alignment rates"
    logs = sorted(glob.glob(os.path.join(expected_dir, "alignment_logs", "*.hisat2.log")))
    bad = []
    for log in logs:
        sample = os.path.basename(log)
        obs = os.path.join(observed_dir, sample)
        exp_rate = alignment_rate(log)
        obs_rate = alignment_rate(obs) if os.path.exists(obs) else None
        if obs_rate is None or abs(obs_rate - exp_rate) > 0.05:
            bad.append("%s (expected %s%%, got %s)" % (sample.replace(".hisat2.log", ""), exp_rate, obs_rate))
    if bad:
        return False, name, "; ".join(bad)
    return True, name, "%d/%d samples within 0.05 percentage points" % (len(logs), len(logs))


def to_float(text):
    return None if text in ("NA", "", "NaN") else float(text)


def same_pvalue(a, b):
    if a is None or b is None:
        return a is None and b is None
    if a == 0 or b == 0:
        return a == b
    return abs(math.log10(a) - math.log10(b)) <= 0.01


def check_deseq2(expected_dir, observed_dir):
    name = "DESeq2 results"
    exp_path = os.path.join(expected_dir, "deseq2", "deseq2_results.csv")
    obs_path = os.path.join(observed_dir, "deseq2_results.csv")
    if not os.path.exists(obs_path):
        return False, name, "file missing from this run"
    exp_head, exp_rows = read_table(exp_path, ",")
    obs_head, obs_rows = read_table(obs_path, ",")
    if exp_head != obs_head:
        return False, name, "columns differ"
    exp = {r[0]: dict(zip(exp_head, r)) for r in exp_rows}
    obs = {r[0]: dict(zip(obs_head, r)) for r in obs_rows}
    if set(exp) != set(obs):
        return False, name, "tested gene sets differ (%d vs %d genes)" % (len(obs), len(exp))
    bad = []
    for gene, e in exp.items():
        o = obs[gene]
        for col in ("baseMean", "log2FoldChange", "log2FoldChange_shrunk"):
            a, b = to_float(e[col]), to_float(o[col])
            if (a is None) != (b is None) or (a is not None and not math.isclose(a, b, rel_tol=1e-3, abs_tol=1e-6)):
                bad.append("%s:%s" % (gene, col))
        for col in ("pvalue", "padj"):
            if not same_pvalue(to_float(e[col]), to_float(o[col])):
                bad.append("%s:%s" % (gene, col))
    if bad:
        return False, name, "%d values outside tolerance, e.g. %s" % (len(bad), ", ".join(bad[:5]))
    return True, name, "%d genes compared, all values within tolerance" % len(exp)


def check_files(observed_dir):
    name = "Figures and MultiQC report"
    missing = [f for f in FIGURES
               if not os.path.exists(os.path.join(observed_dir, f)) or os.path.getsize(os.path.join(observed_dir, f)) == 0]
    if missing:
        return False, name, "missing or empty: %s" % ", ".join(missing)
    return True, name, "all %d files present" % len(FIGURES)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected", required=True)
    parser.add_argument("--observed", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    results = []
    for check in (lambda: check_counts(args.expected, args.observed),
                  lambda: check_alignment(args.expected, args.observed),
                  lambda: check_deseq2(args.expected, args.observed),
                  lambda: check_files(args.observed)):
        try:
            results.append(check())
        except Exception as error:  # a broken input should be reported, not crash the pipeline
            results.append((False, "Check could not run", repr(error)))

    lines = ["Reproducibility check against expected_results/"]
    for ok, name, detail in results:
        lines.append("  [%s] %s: %s" % ("PASS" if ok else "FAIL", name, detail))
    lines.append("OVERALL: %s" % ("PASS" if all(r[0] for r in results) else "FAIL"))
    text = "\n".join(lines) + "\n"
    with open(args.out, "w") as handle:
        handle.write(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
