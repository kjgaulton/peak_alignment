#!/usr/bin/env python3
"""
align_peaks.py

Build a unified, non-overlapping consensus peak set across multiple
narrowPeak (MACS2) files, and map every original peak in every input
file to the unified peak it was collapsed into.

Algorithm (iterative greedy summit-window selection)
-----------------------------------------------------
1. Pool all peaks from all input files together.
2. Repeat until the pool is empty:
     a. Take the peak in the pool with the strongest signal (signalValue,
        narrowPeak column 7) -- this is the "seed" peak.
     b. Build a fixed-width window (default 300bp) centered on the seed
        peak's summit (chromStart + summit offset, narrowPeak column 10).
        Clip the window to chromosome boundaries.
     c. Find every other peak still in the pool whose ORIGINAL interval
        overlaps the seed peak's ORIGINAL interval. Remove the seed peak
        and all of these overlapping peaks from the pool. Record that each
        of them (seed included) maps to the new unified window.
     d. Continue with whatever peaks remain in the pool (they did not
        overlap the seed peak).
3. The set of fixed-width windows created in step (b) across all
   iterations is the unified peak set.

Because peaks are always processed in descending signal order and a peak
is only ever consumed once, this is equivalent to literally repeating
"among overlapping peaks, keep the strongest, resize it, continue with the
non-overlapping remainder" until no peaks are left.

Output
------
- <outdir>/unified_peaks.bed        Final non-overlapping 300bp consensus peaks.
- <outdir>/mapping/<file>.mapped.bed  Per input file: every original peak with
  the unified peak it was assigned to appended as extra columns.
- <outdir>/mapping/all_peaks_mapped.tsv  Same mapping, all files combined.

Usage
-----
    python3 align_peaks.py --input peaks/*.narrowPeak --outdir results \
        --window 300 --genome hg38
"""
import argparse
import glob
import os
import sys
from dataclasses import dataclass

import pandas as pd
from intervaltree import IntervalTree

from chrom_sizes import GENOME_SIZES

NARROWPEAK_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "signalValue", "pValue", "qValue", "summit_offset",
]


@dataclass
class Peak:
    global_id: int
    file_index: int
    file_name: str
    orig_line_index: int
    chrom: str
    start: int
    end: int
    name: str
    signalValue: float
    summit: int  # absolute summit position


def load_narrowpeak(path, file_index):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    ncols = df.shape[1]
    if ncols < 10:
        raise ValueError(
            f"{path}: expected narrowPeak with >=10 columns, got {ncols}. "
            "If these are plain BED files without a summit/signal column, "
            "this script needs to be adapted for that format."
        )
    df = df.iloc[:, :10]
    df.columns = NARROWPEAK_COLS
    df["file_index"] = file_index
    df["file_name"] = os.path.basename(path)
    df["orig_line_index"] = df.index
    # Absolute summit position. narrowPeak uses -1 when no summit was called;
    # fall back to the interval midpoint in that case.
    has_summit = df["summit_offset"] >= 0
    df["summit"] = df["start"] + df["summit_offset"]
    df.loc[~has_summit, "summit"] = (
        (df.loc[~has_summit, "start"] + df.loc[~has_summit, "end"]) // 2
    )
    return df


def build_peak_pool(input_files):
    frames = []
    for i, path in enumerate(input_files):
        frames.append(load_narrowpeak(path, i))
    combined = pd.concat(frames, ignore_index=True)
    combined["global_id"] = combined.index
    return combined


def build_chrom_trees(peaks_df):
    trees = {}
    for chrom, sub in peaks_df.groupby("chrom"):
        tree = IntervalTree()
        for row in sub.itertuples():
            if row.end > row.start:
                tree[row.start:row.end] = row.global_id
        trees[chrom] = tree
    return trees


def clip_window(chrom, center, half_width, genome_sizes):
    start = center - half_width
    end = center + half_width
    chrom_len = genome_sizes.get(chrom)
    if start < 0:
        start = 0
    if chrom_len is not None and end > chrom_len:
        end = chrom_len
    if chrom_len is not None and start > chrom_len:
        start = chrom_len
    return start, end


def run_iterative_selection(peaks_df, window, genome_sizes):
    """Returns (unified_peaks_df, mapping_series) where mapping_series maps
    global_id -> unified_peak_id."""
    half = window // 2
    trees = build_chrom_trees(peaks_df)
    by_id = peaks_df.set_index("global_id")

    order = peaks_df.sort_values(
        "signalValue", ascending=False
    )["global_id"].tolist()

    consumed = set()
    mapping = {}
    unified_rows = []
    unified_id = 0

    for gid in order:
        if gid in consumed:
            continue
        seed = by_id.loc[gid]
        chrom = seed["chrom"]
        tree = trees[chrom]

        # Find every still-active peak whose ORIGINAL interval overlaps the
        # seed peak's ORIGINAL interval (this always includes the seed
        # itself, since the seed's own interval is still in the tree).
        hits = list(tree.overlap(seed["start"], seed["end"]))
        if not hits:
            # Seed's own interval was zero-length or otherwise missing from
            # the tree; fall back to treating it as a singleton cluster.
            overlapping_ids = [gid]
        else:
            overlapping_ids = [iv.data for iv in hits]

        start, end = clip_window(chrom, int(seed["summit"]), half, genome_sizes)
        uid = f"unified_peak_{unified_id}"
        unified_id += 1
        unified_rows.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "unified_id": uid,
            "seed_signalValue": seed["signalValue"],
            "seed_file": seed["file_name"],
            "seed_peak_name": seed["name"],
            "n_peaks_merged": len(overlapping_ids),
        })

        # Remove exactly these intervals (by identity, not by a fresh range
        # query) so we never collaterally drop peaks that only overlap a
        # cluster member but not the seed itself.
        for iv in hits:
            tree.remove(iv)
        for oid in overlapping_ids:
            consumed.add(oid)
            mapping[oid] = uid

    unified_df = pd.DataFrame(unified_rows)
    return unified_df, mapping


def write_outputs(peaks_df, unified_df, mapping, outdir):
    os.makedirs(outdir, exist_ok=True)
    mapping_dir = os.path.join(outdir, "mapping")
    os.makedirs(mapping_dir, exist_ok=True)

    unified_sorted = unified_df.sort_values(["chrom", "start"]).reset_index(drop=True)
    unified_path = os.path.join(outdir, "unified_peaks.bed")
    unified_sorted.to_csv(unified_path, sep="\t", header=True, index=False)

    peaks_df = peaks_df.copy()
    peaks_df["unified_id"] = peaks_df["global_id"].map(mapping)

    unified_lookup = unified_df.set_index("unified_id")[["chrom", "start", "end"]]
    unified_lookup = unified_lookup.rename(
        columns={"chrom": "unified_chrom", "start": "unified_start", "end": "unified_end"}
    )
    peaks_df = peaks_df.join(unified_lookup, on="unified_id")

    out_cols = [
        "chrom", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "summit_offset", "summit",
        "unified_id", "unified_chrom", "unified_start", "unified_end",
    ]

    combined_path = os.path.join(mapping_dir, "all_peaks_mapped.tsv")
    peaks_df[["file_name"] + out_cols].to_csv(combined_path, sep="\t", index=False)

    for fname, sub in peaks_df.groupby("file_name"):
        out_path = os.path.join(mapping_dir, f"{os.path.splitext(fname)[0]}.mapped.bed")
        sub[out_cols].sort_values(["chrom", "start"]).to_csv(
            out_path, sep="\t", header=True, index=False
        )

    return unified_path, combined_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", nargs="+", required=True,
        help="narrowPeak files or glob patterns (e.g. peaks/*.narrowPeak)",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--window", type=int, default=300,
        help="Fixed window width centered on each retained summit (default: 300)",
    )
    parser.add_argument(
        "--genome", default="hg38", choices=sorted(GENOME_SIZES.keys()),
        help="Genome build, used to clip windows to chromosome ends (default: hg38)",
    )
    args = parser.parse_args()

    input_files = []
    for pattern in args.input:
        matches = sorted(glob.glob(pattern))
        input_files.extend(matches if matches else [pattern])
    input_files = [f for f in input_files if os.path.isfile(f)]
    if not input_files:
        sys.exit("No input files found.")

    print(f"Loading {len(input_files)} peak file(s)...")
    peaks_df = build_peak_pool(input_files)
    print(f"Total peaks pooled: {len(peaks_df)}")

    genome_sizes = GENOME_SIZES[args.genome]
    unified_df, mapping = run_iterative_selection(peaks_df, args.window, genome_sizes)
    print(f"Unified peaks produced: {len(unified_df)}")

    unified_path, combined_path = write_outputs(peaks_df, unified_df, mapping, args.outdir)
    print(f"Unified peak set:      {unified_path}")
    print(f"Combined mapping:      {combined_path}")
    print(f"Per-file mapped beds:  {os.path.join(args.outdir, 'mapping')}/")


if __name__ == "__main__":
    main()
