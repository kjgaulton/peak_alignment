#!/usr/bin/env python3
"""
qc_peaks.py

Quality metrics for each original peak file, based on how its peaks land
in the unified peak set produced by align_peaks.py.

For every original peak, once it's been mapped to a unified peak, we know
which OTHER input files also placed a peak in that same unified window --
that's a direct measure of reproducibility. This script summarizes that
per file:

  - pct_overlap_other_file: percent of a file's peaks whose unified window
    is also occupied by at least one peak from a different file. A file
    scoring low here is contributing a lot of peaks nobody else supports.
  - mean / median / max_other_files: how many *other* files, on average,
    back up a given peak from this file.
  - A full distribution table: for each file, how many of its peaks fall
    into a unified window supported by exactly 0, 1, 2, ... other files.

This does not recompute peak overlaps from scratch -- it reads the
`all_peaks_mapped.tsv` that align_peaks.py already produced (file_name,
chrom, start, end, unified_id), so it stays fast even at tens of millions
of peaks.

Usage
-----
    python3 qc_peaks.py --mapping-tsv results/mapping/all_peaks_mapped.tsv \
        --outdir results/qc --flag-below 50

Output
------
- <outdir>/file_quality_summary.tsv   One row per input file with the
  metrics above.
- <outdir>/file_overlap_distribution.tsv  Long-format: file_name,
  n_other_files, n_peaks, pct_of_file_peaks.
- <outdir>/flagged_files.txt   One file_name per line, for any file whose
  pct_overlap_other_file fell below --flag-below (only written if
  --flag-below is given). Feed straight into align_peaks.py's --exclude:
      python3 align_peaks.py --input peaks/*.narrowPeak \\
          --exclude $(cat results/qc/flagged_files.txt) \\
          --outdir results_v2 ...
"""
import argparse
import os
import sys

import pandas as pd


def compute_metrics(df):
    # Distinct files contributing to each unified peak. Since this is a
    # per-unified-id count (not per-row), it doesn't matter if a single
    # file contributes more than one of its own peaks to the same unified
    # window -- that file is still counted once, so "n_other_files" is the
    # same for every row of that file at that unified id.
    n_files_at_unified = (
        df.groupby("unified_id")["file_name"].nunique().rename("n_files_at_unified")
    )
    df = df.merge(n_files_at_unified, on="unified_id", how="left")
    df["n_other_files"] = df["n_files_at_unified"] - 1
    df["supported_by_other"] = df["n_other_files"] >= 1
    return df


def build_summary(df):
    summary = df.groupby("file_name").agg(
        n_peaks=("unified_id", "size"),
        n_unified_peaks=("unified_id", "nunique"),
        n_supported_by_other=("supported_by_other", "sum"),
        mean_other_files=("n_other_files", "mean"),
        median_other_files=("n_other_files", "median"),
        max_other_files=("n_other_files", "max"),
    ).reset_index()
    summary["pct_overlap_other_file"] = (
        100 * summary["n_supported_by_other"] / summary["n_peaks"]
    )
    summary = summary[
        ["file_name", "n_peaks", "n_unified_peaks", "n_supported_by_other",
         "pct_overlap_other_file", "mean_other_files", "median_other_files",
         "max_other_files"]
    ]
    return summary.sort_values("pct_overlap_other_file").reset_index(drop=True)


def build_distribution(df):
    dist = (
        df.groupby(["file_name", "n_other_files"])
        .size()
        .rename("n_peaks")
        .reset_index()
    )
    dist["pct_of_file_peaks"] = (
        100 * dist["n_peaks"] / dist.groupby("file_name")["n_peaks"].transform("sum")
    )
    return dist.sort_values(["file_name", "n_other_files"]).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapping-tsv", required=True,
        help="Path to align_peaks.py's mapping/all_peaks_mapped.tsv output",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--flag-below", type=float, default=None,
        help=(
            "If set, write flagged_files.txt listing any file whose "
            "pct_overlap_other_file is below this percentage (e.g. 50)."
        ),
    )
    args = parser.parse_args()

    if not os.path.isfile(args.mapping_tsv):
        sys.exit(f"Mapping file not found: {args.mapping_tsv}")

    df = pd.read_csv(args.mapping_tsv, sep="\t")
    required = {"file_name", "unified_id"}
    if not required.issubset(df.columns):
        sys.exit(f"{args.mapping_tsv} is missing required columns {required}")

    df = compute_metrics(df)
    summary = build_summary(df)
    distribution = build_distribution(df)

    os.makedirs(args.outdir, exist_ok=True)
    summary_path = os.path.join(args.outdir, "file_quality_summary.tsv")
    dist_path = os.path.join(args.outdir, "file_overlap_distribution.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    distribution.to_csv(dist_path, sep="\t", index=False)

    print(summary.to_string(index=False))
    print(f"\nSummary:     {summary_path}")
    print(f"Distribution: {dist_path}")

    if args.flag_below is not None:
        flagged = summary.loc[
            summary["pct_overlap_other_file"] < args.flag_below, "file_name"
        ].tolist()
        flagged_path = os.path.join(args.outdir, "flagged_files.txt")
        with open(flagged_path, "w") as fh:
            fh.write("\n".join(flagged) + ("\n" if flagged else ""))
        print(f"\nFlagged {len(flagged)} file(s) below {args.flag_below}% overlap-with-other-file rate:")
        for f in flagged:
            print(f"  {f}")
        print(f"Flagged list: {flagged_path}")


if __name__ == "__main__":
    main()
