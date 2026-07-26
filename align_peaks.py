#!/usr/bin/env python3
"""
align_peaks.py

Build a unified, non-overlapping consensus peak set across multiple
narrowPeak (MACS2) files -- optionally alongside coordinate-only BED files
that carry no signal/summit information -- and map every original peak in
every input file to the unified peak it was collapsed into.

Algorithm (iterative greedy summit-window selection, two priority tiers)
--------------------------------------------------------------------------
Peaks come in two tiers:
  - Tier 0 ("scored"): narrowPeak entries, ranked by signalValue.
  - Tier 1 ("coordinate-only"): plain BED entries with no score or summit,
    ranked by peak width (wider first), used only as a fallback signal in
    the absence of anything better.

1. Pool all peaks from all input files (both tiers) together.
2. Process peaks in priority order -- every tier 0 peak (by descending
   signalValue) before any tier 1 peak (by descending width). Repeat until
   the pool is empty:
     a. Take the next peak in that order that hasn't already been consumed
        -- this is the "seed" peak.
     b. Build a fixed-width window (default 300bp) centered on the seed
        peak's summit (narrowPeak: chromStart + summit offset; coordinate-
        only: interval midpoint). Clip the window to chromosome boundaries.
     c. Find every other peak still in the pool whose ORIGINAL interval
        overlaps the seed peak's ORIGINAL interval, regardless of tier.
        Remove the seed peak and all of these overlapping peaks from the
        pool. Record that each of them (seed included) maps to the new
        unified window.
     d. Continue with whatever peaks remain in the pool (they did not
        overlap the seed peak).
3. The set of fixed-width windows created in step (b) across all
   iterations is the unified peak set.

Because tier 0 peaks are always processed before tier 1 peaks, any
coordinate-only peak that overlaps a scored peak gets absorbed into that
scored peak's window during tier 0 processing and never gets a turn to
seed its own window. A coordinate-only peak only seeds a new unified
window in places no scored peak's original interval touches. Within tier
1, wider peaks get priority over narrower ones under the same rule.

Output
------
- <outdir>/unified_peaks.bed        Final non-overlapping 300bp consensus peaks
  (with seed/provenance columns).
- <outdir>/mapping/<file>.mapped.bed  Per input file, one row per original
  peak in its original order: a headerless BED4 of chrom, start, end,
  unified_id -- the coordinates of the unified peak that original peak was
  collapsed into (original peak's own coordinates/score/etc. are dropped).
- <outdir>/mapping/all_peaks_mapped.tsv  Same mapping for all input files
  combined, with a file_name column added.

Usage
-----
    python3 align_peaks.py --input peaks/*.narrowPeak --outdir results \
        --window 300 --genome hg38

    # with a lower-priority tier of coordinate-only BED files (chrom, start,
    # end[, name]) mixed in:
    python3 align_peaks.py --input peaks/*.narrowPeak \
        --coord-input extra_peaks/*.bed --outdir results \
        --window 300 --genome hg38
"""
import argparse
import glob
import os
import sys

import pandas as pd
from intervaltree import IntervalTree

from chrom_sizes import GENOME_SIZES

NARROWPEAK_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "signalValue", "pValue", "qValue", "summit_offset",
]

# Common schema every peak (scored or coordinate-only) is normalized to.
COMMON_COLS = [
    "chrom", "start", "end", "name", "score", "strand",
    "signalValue", "pValue", "qValue", "summit_offset",
    "summit", "tier", "width", "rank_value",
]

TIER_SCORED = 0
TIER_COORD_ONLY = 1


def load_narrowpeak(path, file_index):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    ncols = df.shape[1]
    if ncols < 10:
        raise ValueError(
            f"{path}: expected narrowPeak with >=10 columns, got {ncols}. "
            "Use --coord-input instead if these are coordinate-only BED files."
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
    df["tier"] = TIER_SCORED
    df["width"] = df["end"] - df["start"]
    df["rank_value"] = df["signalValue"]
    return df[COMMON_COLS + ["file_index", "file_name", "orig_line_index"]]


def load_coord_only(path, file_index):
    """Plain BED with just chrom/start/end (a 4th name column is used if
    present). No score and no summit -- ranked within its own lower-priority
    tier by peak width, and centered on the interval midpoint."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    ncols = df.shape[1]
    if ncols < 3:
        raise ValueError(f"{path}: expected at least 3 columns (chrom, start, end).")
    df = df.iloc[:, :4] if ncols >= 4 else df.iloc[:, :3]
    if ncols >= 4:
        df.columns = ["chrom", "start", "end", "name"]
    else:
        df.columns = ["chrom", "start", "end"]
        df["name"] = [f"coord_peak_{i}" for i in df.index]

    df["file_index"] = file_index
    df["file_name"] = os.path.basename(path)
    df["orig_line_index"] = df.index
    df["score"] = float("nan")
    df["strand"] = "."
    df["signalValue"] = float("nan")
    df["pValue"] = float("nan")
    df["qValue"] = float("nan")
    df["summit_offset"] = float("nan")
    df["summit"] = (df["start"] + df["end"]) // 2
    df["tier"] = TIER_COORD_ONLY
    df["width"] = df["end"] - df["start"]
    df["rank_value"] = df["width"]
    return df[COMMON_COLS + ["file_index", "file_name", "orig_line_index"]]


def build_peak_pool(scored_files, coord_only_files):
    frames = []
    for i, path in enumerate(scored_files):
        frames.append(load_narrowpeak(path, i))
    offset = len(scored_files)
    for j, path in enumerate(coord_only_files):
        frames.append(load_coord_only(path, offset + j))
    if not frames:
        raise ValueError("No input peak files provided.")
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

    # Tier 0 (scored) entirely before tier 1 (coordinate-only); descending
    # rank_value within each tier (signalValue for tier 0, width for tier 1).
    order = peaks_df.sort_values(
        ["tier", "rank_value"], ascending=[True, False]
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
            "seed_tier": int(seed["tier"]),
            "seed_signalValue": seed["signalValue"],
            "seed_width": int(seed["width"]),
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

    # Per-original-peak mapping output: just the unified peak's own
    # coordinates plus its ID (BED4: chrom, start, end, name=unified_id).
    # One row per original input peak, in its original file order, with the
    # original peak's coordinates replaced by the unified window it maps to.
    bed4 = peaks_df[
        ["file_name", "orig_line_index", "unified_chrom", "unified_start", "unified_end", "unified_id"]
    ].rename(columns={
        "unified_chrom": "chrom", "unified_start": "start", "unified_end": "end",
    })

    combined_path = os.path.join(mapping_dir, "all_peaks_mapped.tsv")
    bed4.sort_values(["file_name", "orig_line_index"])[
        ["file_name", "chrom", "start", "end", "unified_id"]
    ].to_csv(combined_path, sep="\t", header=True, index=False)

    for fname, sub in bed4.groupby("file_name"):
        out_path = os.path.join(mapping_dir, f"{os.path.splitext(fname)[0]}.mapped.bed")
        sub.sort_values("orig_line_index")[["chrom", "start", "end", "unified_id"]].to_csv(
            out_path, sep="\t", header=False, index=False
        )

    return unified_path, combined_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", nargs="+", required=True,
        help="narrowPeak files or glob patterns (e.g. peaks/*.narrowPeak)",
    )
    parser.add_argument(
        "--coord-input", nargs="+", default=[],
        help=(
            "Optional coordinate-only BED files or glob patterns "
            "(chrom, start, end[, name] -- no score/summit). Treated as a "
            "strictly lower-priority tier: only ever seeds a unified window "
            "in regions no --input peak occupies, ranked by width within "
            "the tier."
        ),
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

    def resolve(patterns):
        files = []
        for pattern in patterns:
            matches = sorted(glob.glob(pattern))
            files.extend(matches if matches else [pattern])
        return [f for f in files if os.path.isfile(f)]

    input_files = resolve(args.input)
    if not input_files:
        sys.exit("No input files found for --input.")

    coord_files = resolve(args.coord_input) if args.coord_input else []
    if args.coord_input and not coord_files:
        sys.exit("No input files found for --coord-input.")

    print(f"Loading {len(input_files)} scored peak file(s)...")
    if coord_files:
        print(f"Loading {len(coord_files)} coordinate-only peak file(s) (lower-priority tier)...")
    peaks_df = build_peak_pool(input_files, coord_files)
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
