#!/usr/bin/env python3
"""
plot_qc.py

QC plots for the unified peak set as a whole -- as opposed to
qc_peaks.py's per-input-file quality metrics, this looks at reproducibility
across the entire final peak set.

Produces two plots, both derived directly from align_peaks.py's mapping
output (mapping/all_peaks_mapped.tsv), so no re-running of the alignment
itself is needed:

1. peak_support_distribution.png -- a histogram of how many distinct
   input files/assays support each unified peak (i.e. placed at least one
   original peak in that unified window). Shows how reproducible the
   unified peak set is overall: a set dominated by peaks seen in only one
   or two files suggests a lot of singleton/noisy calls, while a set
   shifted toward higher support counts suggests a more reproducible core
   set of peaks.

2. peaks_retained_by_min_support.png -- the number of unified peaks that
   would remain if you required at least K supporting files, for every K
   from 1 up to the total number of input files. This is a direct,
   monotonically non-increasing readout of "how much of the unified peak
   set survives at each reproducibility bar," computed straight from data
   already on hand (not by re-running the alignment on file subsets).

Usage
-----
    python3 plot_qc.py --mapping-tsv results/mapping/all_peaks_mapped.tsv \
        --outdir results/qc_plots
"""
import argparse
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def compute_peak_support(mapping_df):
    """One row per unified peak: unified_id, n_files_supporting -- the
    number of distinct input files that placed at least one original peak
    in that unified window."""
    return (
        mapping_df.groupby("unified_id")["file_name"]
        .nunique()
        .rename("n_files_supporting")
        .reset_index()
    )


def plot_support_distribution(peak_support, n_total_files, out_path):
    counts = (
        peak_support["n_files_supporting"]
        .value_counts()
        .reindex(range(1, n_total_files + 1), fill_value=0)
        .sort_index()
    )
    fig, ax = plt.subplots(figsize=(max(6, n_total_files * 0.15), 5))
    ax.bar(counts.index, counts.values, color="#3b6fa0", edgecolor="none", width=0.9)
    ax.set_xlabel("Number of datasets/assays supporting the peak")
    ax.set_ylabel("Number of unified peaks")
    ax.set_title("Distribution of peak support across datasets")
    if n_total_files <= 30:
        ax.set_xticks(counts.index)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return counts


def plot_retained_by_threshold(peak_support, n_total_files, out_path):
    thresholds = list(range(1, n_total_files + 1))
    support_values = peak_support["n_files_supporting"].to_numpy()
    retained = [int((support_values >= k).sum()) for k in thresholds]

    fig, ax = plt.subplots(figsize=(max(6, n_total_files * 0.15), 5))
    ax.step(thresholds, retained, where="post", color="#a0473b")
    ax.set_xlabel("Minimum number of supporting datasets required")
    ax.set_ylabel("Number of unified peaks retained")
    ax.set_title("Unified peaks retained as the support threshold increases")
    if n_total_files <= 30:
        ax.set_xticks(thresholds)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return pd.DataFrame({"min_datasets_required": thresholds, "n_peaks_retained": retained})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapping-tsv", required=True,
        help="Path to align_peaks.py's mapping/all_peaks_mapped.tsv output",
    )
    parser.add_argument("--outdir", required=True, help="Output directory for plots/tables")
    args = parser.parse_args()

    if not os.path.isfile(args.mapping_tsv):
        sys.exit(f"Mapping file not found: {args.mapping_tsv}")

    df = pd.read_csv(args.mapping_tsv, sep="\t")
    required = {"file_name", "unified_id"}
    if not required.issubset(df.columns):
        sys.exit(f"{args.mapping_tsv} is missing required columns {required}")

    n_total_files = df["file_name"].nunique()
    peak_support = compute_peak_support(df)

    try:
        os.makedirs(args.outdir, exist_ok=True)
        dist_png = os.path.join(args.outdir, "peak_support_distribution.png")
        retained_png = os.path.join(args.outdir, "peaks_retained_by_min_support.png")

        counts = plot_support_distribution(peak_support, n_total_files, dist_png)
        retained_df = plot_retained_by_threshold(peak_support, n_total_files, retained_png)

        counts_path = os.path.join(args.outdir, "peak_support_distribution.tsv")
        counts.rename("n_peaks").rename_axis("n_files_supporting").reset_index().to_csv(
            counts_path, sep="\t", index=False
        )
        retained_path = os.path.join(args.outdir, "peaks_retained_by_min_support.tsv")
        retained_df.to_csv(retained_path, sep="\t", index=False)
    except PermissionError as exc:
        sys.exit(
            f"\nPermissionError writing output: {exc}\n\n"
            f"Likely a Docker + NFS root-squash mismatch: make sure --outdir "
            f"already exists and is owned by you, and run the container as "
            f"your own user, e.g. --user $(id -u):$(id -g)."
        )

    print(f"Total unified peaks: {len(peak_support)}")
    print(f"Total input files:   {n_total_files}")
    print(f"\nSupport distribution plot:  {dist_png}")
    print(f"Support distribution table: {counts_path}")
    print(f"\nRetained-by-threshold plot:  {retained_png}")
    print(f"Retained-by-threshold table: {retained_path}")


if __name__ == "__main__":
    main()
