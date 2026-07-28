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

Note: align_peaks.py now computes these same metrics and auto-excludes
low-quality files itself (see its --min-pct-overlap-other-file /
--min-median-other-files / --min-max-other-files flags, mirrored here).
Use this script standalone when you want to inspect or re-threshold
metrics from an existing mapping table without re-running the alignment
(e.g. an older result set, or trying out different thresholds).

Usage
-----
    python3 qc_peaks.py --mapping-tsv results/mapping/all_peaks_mapped.tsv \
        --outdir results/qc

Output
------
- <outdir>/file_quality_summary.tsv   One row per input file with the
  metrics above.
- <outdir>/file_overlap_distribution.tsv  Long-format: file_name,
  n_other_files, n_peaks, pct_of_file_peaks.
- <outdir>/flagged_files.txt   One file_name per line, for any file
  failing a threshold (only written if at least one file is flagged).
  Feed straight into align_peaks.py's --exclude:
      python3 align_peaks.py --input peaks/*.narrowPeak \\
          --exclude $(cat results/qc/flagged_files.txt) \\
          --outdir results_v2 ...
"""
import argparse
import os
import sys

import pandas as pd

from qc_metrics import compute_metrics, build_summary, build_distribution, flag_low_quality_files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapping-tsv", required=True,
        help="Path to align_peaks.py's mapping/all_peaks_mapped.tsv output",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--min-pct-overlap-other-file", type=float, default=75.0,
        help="Flag a file if this percentage of its peaks aren't backed by another file (default: 75).",
    )
    parser.add_argument(
        "--min-median-other-files", type=float, default=2.0,
        help="Flag a file if the median number of other files backing its peaks is below this (default: 2).",
    )
    parser.add_argument(
        "--min-max-other-files", type=float, default=None,
        help=(
            "Flag a file if even its best-supported peak is backed by "
            "fewer than this many other files. Default: half the number "
            "of distinct files in the mapping table."
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

    min_max_other_files = (
        args.min_max_other_files
        if args.min_max_other_files is not None
        else df["file_name"].nunique() / 2
    )

    try:
        os.makedirs(args.outdir, exist_ok=True)
        summary_path = os.path.join(args.outdir, "file_quality_summary.tsv")
        dist_path = os.path.join(args.outdir, "file_overlap_distribution.tsv")
        summary.to_csv(summary_path, sep="\t", index=False)
        distribution.to_csv(dist_path, sep="\t", index=False)
    except PermissionError as exc:
        sys.exit(
            f"\nPermissionError writing output: {exc}\n\n"
            f"Likely a Docker + NFS root-squash mismatch: make sure --outdir "
            f"already exists and is owned by you, and run the container as "
            f"your own user, e.g. --user $(id -u):$(id -g)."
        )

    print(summary.to_string(index=False))
    print(f"\nSummary:     {summary_path}")
    print(f"Distribution: {dist_path}")

    flagged = flag_low_quality_files(
        summary, args.min_pct_overlap_other_file, args.min_median_other_files,
        min_max_other_files,
    )
    print(
        f"\nThresholds: pct_overlap_other_file < {args.min_pct_overlap_other_file}, "
        f"median_other_files < {args.min_median_other_files}, "
        f"max_other_files < {min_max_other_files:g}"
    )
    if flagged:
        flagged_path = os.path.join(args.outdir, "flagged_files.txt")
        with open(flagged_path, "w") as fh:
            fh.write("\n".join(flagged) + "\n")
        print(f"Flagged {len(flagged)} file(s):")
        for f in flagged:
            print(f"  {f}")
        print(f"Flagged list: {flagged_path}")
    else:
        print("No files flagged.")


if __name__ == "__main__":
    main()
