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

1. Pool all peaks from all input files (both tiers) together. For every
   peak, precompute the fixed-width window (default 300bp) it WOULD get if
   chosen as a seed: centered on its own summit (narrowPeak: chromStart +
   summit offset; coordinate-only: interval midpoint), clipped to
   chromosome boundaries.
2. Process peaks in priority order -- every tier 0 peak (by descending
   signalValue) before any tier 1 peak (by descending width). Repeat until
   the pool is empty:
     a. Take the next peak in that order that hasn't already been consumed
        -- this is the "seed" peak. Its precomputed window is the unified
        window.
     b. Find every other peak still in the pool whose OWN precomputed
        window overlaps the seed's window, regardless of tier. Remove the
        seed and all of these from the pool. Record that each of them
        (seed included) maps to the new unified window.
     c. Continue with whatever peaks remain in the pool (their windows did
        not overlap the seed's window).
3. The set of fixed-width windows chosen in step (a) across all
   iterations is the unified peak set.

Clustering is judged on each peak's own fixed-width summit window, not its
raw chromStart/chromEnd -- this is what actually guarantees the final
unified set is non-overlapping, and it's also how multi-summit peaks are
handled: if a peak caller (e.g. MACS2 --call-summits) reports several
summits sharing one broad raw span, each summit is just another peak in
the pool, and two of them only collapse into the same unified peak if
their own 300bp windows would actually overlap -- summits far enough apart
within the same broad call can each seed their own separate unified peak,
while summits close enough together still collapse to whichever has the
strongest signal, exactly like any other overlap.

Because tier 0 peaks are always processed before tier 1 peaks, any
coordinate-only peak whose window overlaps a scored peak's window gets
absorbed into that scored peak's window during tier 0 processing and
never gets a turn to seed its own. A coordinate-only peak only seeds a new
unified window in places no scored peak's window touches. Within tier 1,
wider peaks get priority over narrower ones under the same rule.

Quality control and auto-exclusion
-----------------------------------
After the first alignment pass, every original peak is mapped to a unified
window, which tells us which OTHER input files also placed a peak in that
same window -- a direct measure of reproducibility. By default this script
scores every input file on that basis and re-runs the alignment excluding
any file that falls below either threshold:

  - --min-pct-overlap-other-file (default 75): percent of a file's peaks
    whose unified window is also occupied by a peak from another file.
  - --min-median-other-files (default 2): median number of other files
    backing up a given peak from this file.

Use --qc-report-only to compute and write these metrics without excluding
anything, or --no-qc to skip QC entirely.

Output
------
- <outdir>/unified_peaks.bed        Final non-overlapping 300bp consensus peaks
  (with seed/provenance columns). Reflects the QC-filtered file set unless
  --no-qc or --qc-report-only is given.
- <outdir>/mapping/<file>.mapped.bed  Per input file, one row per original
  peak in its original order: a headerless BED4 of chrom, start, end,
  unified_id -- the coordinates of the unified peak that original peak was
  collapsed into (original peak's own coordinates/score/etc. are dropped).
- <outdir>/mapping/all_peaks_mapped.tsv  Same mapping for all input files
  combined, with a file_name column added.
- <outdir>/qc/file_quality_summary.tsv, <outdir>/qc/file_overlap_distribution.tsv,
  <outdir>/qc/excluded_files.txt  QC metrics from the initial (pre-exclusion)
  pass, and the list of files that were auto-excluded (unless --no-qc).

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
import fnmatch
import glob
import os
import sys

import pandas as pd
from intervaltree import IntervalTree

from chrom_sizes import GENOME_SIZES
from qc_metrics import compute_metrics, build_summary, build_distribution, flag_low_quality_files

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


def compute_windows(peaks_df, window, genome_sizes):
    """Adds window_start/window_end: the fixed-width window each peak would
    get if IT were chosen as a seed, centered on its own summit and clipped
    to chromosome boundaries. Clustering is decided from these windows, not
    from the raw chromStart/chromEnd -- see run_iterative_selection."""
    half = window // 2
    peaks_df = peaks_df.copy()

    starts = (peaks_df["summit"] - half).clip(lower=0)
    ends = peaks_df["summit"] + half

    chrom_len = peaks_df["chrom"].map(genome_sizes)
    unknown = chrom_len.isna()
    ends = ends.where(unknown | (ends <= chrom_len), chrom_len)
    starts = starts.where(unknown | (starts <= chrom_len), chrom_len)

    peaks_df["window_start"] = starts.astype("int64")
    peaks_df["window_end"] = ends.astype("int64")
    return peaks_df


def build_chrom_trees(peaks_df):
    trees = {}
    for chrom, sub in peaks_df.groupby("chrom"):
        tree = IntervalTree()
        for row in sub.itertuples():
            if row.window_end > row.window_start:
                tree[row.window_start:row.window_end] = row.global_id
        trees[chrom] = tree
    return trees


def run_iterative_selection(peaks_df, window, genome_sizes):
    """Returns (unified_peaks_df, mapping_series) where mapping_series maps
    global_id -> unified_peak_id.

    Clustering ("do these peaks overlap") is judged on each peak's fixed-
    width summit window, not its raw chromStart/chromEnd. This matters for
    peaks with multiple summits (e.g. MACS2 --call-summits): sub-peaks
    sharing one broad raw span but with summits far enough apart that
    their own 300bp windows wouldn't overlap are allowed to seed separate
    unified peaks, instead of automatically collapsing to whichever
    summit has the strongest signal. It's also what actually guarantees
    the final unified set is non-overlapping -- two raw peaks that don't
    overlap each other can still have summits close enough that their
    fixed-width windows would, and window-based clustering catches that
    case where raw-span-based clustering would not.
    """
    peaks_df = compute_windows(peaks_df, window, genome_sizes)
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

        # Find every still-active peak whose fixed-width summit window
        # overlaps the seed's own summit window (this always includes the
        # seed itself, since its window is still in the tree).
        hits = list(tree.overlap(seed["window_start"], seed["window_end"]))
        if not hits:
            # Seed's own window was zero-length or otherwise missing from
            # the tree; fall back to treating it as a singleton cluster.
            overlapping_ids = [gid]
        else:
            overlapping_ids = [iv.data for iv in hits]

        start, end = int(seed["window_start"]), int(seed["window_end"])
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


def run_qc(peaks_df, mapping, qc_dir, min_pct_overlap_other_file, min_median_other_files):
    """Computes per-file QC metrics from an alignment pass, writes them to
    qc_dir, and returns (summary_df, distribution_df, auto_excluded_files)."""
    qc_input = peaks_df[["file_name", "global_id"]].copy()
    qc_input["unified_id"] = qc_input["global_id"].map(mapping)
    qc_input = compute_metrics(qc_input[["file_name", "unified_id"]])

    summary = build_summary(qc_input)
    distribution = build_distribution(qc_input)

    os.makedirs(qc_dir, exist_ok=True)
    summary_path = os.path.join(qc_dir, "file_quality_summary.tsv")
    dist_path = os.path.join(qc_dir, "file_overlap_distribution.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    distribution.to_csv(dist_path, sep="\t", index=False)

    auto_excluded = flag_low_quality_files(
        summary, min_pct_overlap_other_file, min_median_other_files
    )
    return summary, distribution, auto_excluded, summary_path, dist_path


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
    parser.add_argument(
        "--exclude", nargs="+", default=[],
        help=(
            "File basenames or glob patterns to drop from --input/"
            "--coord-input before pooling (e.g. flagged low-quality files "
            "from qc_peaks.py's flagged_files.txt). Matched against the "
            "basename of each resolved input file."
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
    parser.add_argument(
        "--min-pct-overlap-other-file", type=float, default=75.0,
        help=(
            "QC threshold: a file is auto-excluded from the final unified "
            "set if the percent of its peaks whose unified window is also "
            "occupied by a peak from another file falls below this value "
            "(default: 75)."
        ),
    )
    parser.add_argument(
        "--min-median-other-files", type=float, default=2.0,
        help=(
            "QC threshold: a file is auto-excluded if the median number of "
            "other files supporting its peaks falls below this value "
            "(default: 2)."
        ),
    )
    parser.add_argument(
        "--no-qc", action="store_true",
        help="Skip QC scoring and auto-exclusion entirely; align once on all resolved input files.",
    )
    parser.add_argument(
        "--qc-report-only", action="store_true",
        help="Compute and write QC metrics, but do not auto-exclude any files from the final unified set.",
    )
    args = parser.parse_args()

    def resolve(patterns):
        files = []
        for pattern in patterns:
            matches = sorted(glob.glob(pattern))
            files.extend(matches if matches else [pattern])
        return [f for f in files if os.path.isfile(f)]

    exclude_names = set()
    for pattern in args.exclude:
        exclude_names.add(os.path.basename(pattern))

    def apply_exclude(files):
        kept = []
        for f in files:
            base = os.path.basename(f)
            if any(fnmatch.fnmatch(base, pat) or base == pat for pat in exclude_names):
                print(f"Excluding {base} (matched --exclude)")
                continue
            kept.append(f)
        return kept

    input_files = resolve(args.input)
    if args.exclude:
        input_files = apply_exclude(input_files)
    if not input_files:
        sys.exit("No input files found for --input (after applying --exclude).")

    coord_files = resolve(args.coord_input) if args.coord_input else []
    if args.exclude and coord_files:
        coord_files = apply_exclude(coord_files)
    if args.coord_input and not coord_files:
        sys.exit("No input files found for --coord-input (after applying --exclude).")

    genome_sizes = GENOME_SIZES[args.genome]

    print(f"Loading {len(input_files)} scored peak file(s)...")
    if coord_files:
        print(f"Loading {len(coord_files)} coordinate-only peak file(s) (lower-priority tier)...")
    peaks_df = build_peak_pool(input_files, coord_files)
    print(f"Total peaks pooled: {len(peaks_df)}")

    unified_df, mapping = run_iterative_selection(peaks_df, args.window, genome_sizes)
    print(f"Unified peaks produced: {len(unified_df)}")

    if not args.no_qc:
        qc_dir = os.path.join(args.outdir, "qc")
        summary, distribution, auto_excluded, summary_path, dist_path = run_qc(
            peaks_df, mapping, qc_dir,
            args.min_pct_overlap_other_file, args.min_median_other_files,
        )
        print()
        print(summary.to_string(index=False))
        print(f"\nQC summary:      {summary_path}")
        print(f"QC distribution: {dist_path}")

        if auto_excluded:
            print(
                f"\n{len(auto_excluded)} file(s) failed QC thresholds "
                f"(pct_overlap_other_file < {args.min_pct_overlap_other_file} "
                f"or median_other_files < {args.min_median_other_files}):"
            )
            for f in auto_excluded:
                print(f"  {f}")
            excluded_path = os.path.join(qc_dir, "excluded_files.txt")
            with open(excluded_path, "w") as fh:
                fh.write("\n".join(auto_excluded) + "\n")
            print(f"Excluded file list: {excluded_path}")

            if args.qc_report_only:
                print("--qc-report-only set: keeping these files in the final unified set.")
            else:
                remaining_input = [
                    f for f in input_files if os.path.basename(f) not in auto_excluded
                ]
                remaining_coord = [
                    f for f in coord_files if os.path.basename(f) not in auto_excluded
                ]
                if not remaining_input:
                    sys.exit(
                        "All --input files failed the QC thresholds -- nothing left to "
                        "align. Loosen --min-pct-overlap-other-file/--min-median-other-files "
                        "or rerun with --qc-report-only."
                    )
                print(
                    f"\nRe-running alignment on the {len(remaining_input)} file(s) that "
                    f"passed QC ({len(remaining_coord)} coordinate-only)..."
                )
                peaks_df = build_peak_pool(remaining_input, remaining_coord)
                print(f"Total peaks pooled (post-QC): {len(peaks_df)}")
                unified_df, mapping = run_iterative_selection(peaks_df, args.window, genome_sizes)
                print(f"Unified peaks produced (post-QC): {len(unified_df)}")
        else:
            print("\nAll files passed QC thresholds -- nothing excluded.")

    unified_path, combined_path = write_outputs(peaks_df, unified_df, mapping, args.outdir)
    print(f"\nUnified peak set:      {unified_path}")
    print(f"Combined mapping:      {combined_path}")
    print(f"Per-file mapped beds:  {os.path.join(args.outdir, 'mapping')}/")


if __name__ == "__main__":
    main()
