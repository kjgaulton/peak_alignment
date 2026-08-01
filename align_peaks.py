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
any file that falls below any of three thresholds:

  - --min-pct-overlap-other-file (default 75): percent of a file's peaks
    whose unified window is also occupied by a peak from another file.
  - --min-median-other-files (default 2): median number of other files
    backing up a given peak from this file.
  - --min-max-other-files (default: half the number of files entering
    this QC pass): number of other files backing up this file's single
    BEST-supported peak. Catches a file with no peak reproduced across a
    meaningful fraction of the cohort, even if its other metrics pass.

Use --qc-report-only to compute and write these metrics without excluding
anything, or --no-qc to skip QC entirely.

Blacklist filtering
--------------------
If --blacklist-bed is given, any unified peak (from the final, post-QC
set) that overlaps a region in that BED file is dropped entirely, and any
original peak that had mapped to it is dropped from the mapping outputs
too. This does not bundle ENCODE blacklist region data itself -- supply
the BED file yourself, e.g. the official ENCODE blacklist v2 lists from
https://github.com/Boyle-Lab/Blacklist/tree/master/lists (see
blacklist.py for exact download commands).

Biosample annotation and filtering
------------------------------------
If --metadata-file is given (a delimited file mapping file accession ->
biosample_id -- see metadata.py for the exact matching rules), the final
unified_peaks.bed gets 'biosamples'/'n_biosamples' columns: the
deduplicated biosample IDs of every assay/file overlapping that peak, and
how many there are. Any unified peak backed by fewer than
--min-biosamples (default 2) distinct biosamples is then dropped
entirely -- by default this excludes peaks found in only one biosample.
This filtering only happens if --metadata-file is given.

Output
------
- <outdir>/unified_peaks.bed        Final non-overlapping 300bp consensus peaks
  (with seed/provenance columns, including overlapping_files and, if
  --metadata-file is given, biosamples/n_biosamples). Reflects the
  QC-filtered file set unless --no-qc or --qc-report-only is given, and
  excludes any blacklisted peaks if --blacklist-bed was given, and any
  peaks below --min-biosamples if --metadata-file was given.
- <outdir>/mapping/<file>.mapped.bed  Per input file, one row per original
  peak in its original order: a headerless BED4 of chrom, start, end,
  unified_id -- the coordinates of the unified peak that original peak was
  collapsed into (original peak's own coordinates/score/etc. are dropped).
- <outdir>/mapping/all_peaks_mapped.tsv  Same mapping for all input files
  combined, with a file_name column added.
- <outdir>/qc/file_quality_summary.tsv, <outdir>/qc/file_overlap_distribution.tsv,
  <outdir>/qc/excluded_files.txt  QC metrics from the initial (pre-exclusion)
  pass, and the list of files that were auto-excluded (unless --no-qc).
- <outdir>/blacklist_removed_peaks.bed  Unified peaks that were dropped for
  overlapping --blacklist-bed (only written if --blacklist-bed is given).

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
import traceback

import pandas as pd
from intervaltree import IntervalTree

from chrom_sizes import GENOME_SIZES
from qc_metrics import compute_metrics, build_summary, build_distribution, flag_low_quality_files
from blacklist import load_blacklist, filter_unified_peaks
from metadata import load_metadata

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
    df["file_stem"] = os.path.splitext(os.path.basename(path))[0]
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
    return df[COMMON_COLS + ["file_index", "file_name", "file_stem", "orig_line_index"]]


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
    df["file_stem"] = os.path.splitext(os.path.basename(path))[0]
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
    return df[COMMON_COLS + ["file_index", "file_name", "file_stem", "orig_line_index"]]


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
    chroms = sorted(peaks_df["chrom"].unique())
    print(f"Building interval trees for {len(chroms)} chromosome(s)...", flush=True)
    for chrom, sub in peaks_df.groupby("chrom"):
        tree = IntervalTree()
        for row in sub.itertuples():
            if row.window_end > row.window_start:
                tree[row.window_start:row.window_end] = row.global_id
        trees[chrom] = tree
        print(f"  {chrom}: {len(sub)} peak(s) indexed", flush=True)
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

    total = len(order)
    progress_step = max(total // 20, 100000)
    print(f"Running iterative seed selection over {total} peak(s)...", flush=True)

    for processed, gid in enumerate(order):
        if processed and processed % progress_step == 0:
            print(
                f"  ...processed {processed}/{total} candidates, "
                f"{unified_id} unified peak(s) so far",
                flush=True,
            )
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

        # Every distinct assay (file, extension stripped) with a peak
        # merged into this unified window, including the seed's own file.
        overlapping_files = sorted({
            by_id.loc[oid, "file_stem"] for oid in overlapping_ids
        })

        unified_rows.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "unified_id": uid,
            "seed_signalValue": seed["signalValue"],
            "seed_summit": int(seed["summit"]),
            "seed_file": seed["file_stem"],
            "n_peaks_merged": len(overlapping_ids),
            "overlapping_files": ",".join(overlapping_files),
            "n_assays": len(overlapping_files),
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


def annotate_biosamples(unified_df, metadata_map):
    """Adds 'biosamples' (deduplicated, comma-separated biosample_id for
    every file listed in a unified peak's overlapping_files) and
    'n_biosamples' (count of those) columns. Returns (annotated_df,
    sorted_list_of_unmatched_files).
    """
    unmatched = set()

    def lookup(files_str):
        biosamples = []
        for file_stem in files_str.split(","):
            if not file_stem:
                continue
            biosample = metadata_map.get(file_stem)
            if biosample is None:
                unmatched.add(file_stem)
                continue
            if biosample not in biosamples:
                biosamples.append(biosample)
        return ",".join(sorted(biosamples))

    unified_df = unified_df.copy()
    unified_df["biosamples"] = unified_df["overlapping_files"].map(lookup)
    unified_df["n_biosamples"] = unified_df["biosamples"].map(
        lambda s: len(s.split(",")) if s else 0
    )
    return unified_df, sorted(unmatched)


def filter_by_assay_count(unified_df, mapping, min_assays):
    """Drops any unified peak supported by fewer than min_assays distinct
    assays (files), based on the 'n_assays' column -- independent of
    biosample/metadata annotation, catching e.g. pseudo-replication where
    a peak only shows up in a single assay run more than once. Returns
    (kept_df, kept_mapping, removed_df), mirroring
    blacklist.filter_unified_peaks."""
    is_low = unified_df["n_assays"] < min_assays
    removed_df = unified_df.loc[is_low].reset_index(drop=True)
    kept_df = unified_df.loc[~is_low].reset_index(drop=True)

    removed_ids = set(removed_df["unified_id"])
    kept_mapping = {gid: uid for gid, uid in mapping.items() if uid not in removed_ids}
    return kept_df, kept_mapping, removed_df


def run_qc(
    peaks_df, mapping, qc_dir,
    min_pct_overlap_other_file, min_median_other_files, min_max_other_files,
):
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
        summary, min_pct_overlap_other_file, min_median_other_files, min_max_other_files
    )
    return summary, distribution, auto_excluded, summary_path, dist_path


def _explain_write_error(exc, outdir):
    """Called from an `except OSError` block around an output write. Gives
    a targeted message for the two most common real-world causes (NFS
    permission mismatch, disk full) instead of a bare traceback -- both of
    which have bitten this pipeline in practice on large remote runs."""
    if isinstance(exc, PermissionError):
        sys.exit(
            f"\nPermissionError writing output: {exc}\n\n"
            f"This is almost always a Docker + NFS root-squash mismatch, not a "
            f"bug: the container is writing as root, but '{outdir}' (or a "
            f"subdirectory already inside it) is on an NFS mount that maps "
            f"root to an unprivileged user with no write access.\n\n"
            f"Fix: make sure --outdir already exists and is owned by you, and "
            f"run the container as your own user/group, e.g.:\n"
            f"    mkdir -p {outdir}\n"
            f"    docker run --rm --user $(id -u):$(id -g) -v ... peak-alignment ...\n"
        )
    sys.exit(
        f"\nError writing output to {outdir}: {exc}\n\n"
        f"If this says 'No space left on device', the output disk/volume is "
        f"full -- check with `df -h {outdir}` and free space or point "
        f"--outdir somewhere with more room. Otherwise this may indicate an "
        f"unhealthy mounted volume/NFS share."
    )


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
    parser.add_argument(
        "--blacklist-bed", default=None,
        help=(
            "Optional BED file of blacklist regions (e.g. ENCODE blacklist "
            "v2 -- see blacklist.py for official download links). Any "
            "final unified peak overlapping a region in this file is "
            "dropped, along with any original peaks mapped to it."
        ),
    )
    parser.add_argument(
        "--metadata-file", default=None,
        help=(
            "Optional delimited file mapping file accession -> biosample_id "
            "(see metadata.py for the expected/matched column names). If "
            "given, adds 'biosamples'/'n_biosamples' columns to "
            "unified_peaks.bed: the deduplicated biosamples of every assay "
            "overlapping that peak, and how many there are."
        ),
    )
    parser.add_argument(
        "--min-assays", type=int, default=1,
        help=(
            "Drop any final unified peak supported by fewer than this many "
            "distinct assays/files (default: 1 -- no-op, since every "
            "unified peak has at least one supporting file; raise this to "
            "e.g. 2 to exclude peaks found in only a single assay). Based "
            "on the 'n_assays' count, independent of --metadata-file."
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
        "--min-max-other-files", type=float, default=None,
        help=(
            "QC threshold: a file is auto-excluded if even its single "
            "best-supported peak (max_other_files) is backed by fewer "
            "than this many other files -- i.e. it has no peak reproduced "
            "across a meaningful fraction of the cohort. Default: half "
            "the number of files entering this QC pass."
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

    print("=" * 70, flush=True)
    print("Starting align_peaks.py run", flush=True)
    print(f"  outdir: {args.outdir}", flush=True)
    print("=" * 70, flush=True)

    # Load and validate --blacklist-bed / --metadata-file upfront, before any
    # expensive computation, so a bad path or a malformed file fails in a
    # fraction of a second instead of after a full alignment (+QC) pass.
    blacklist_trees = None
    if args.blacklist_bed:
        if not os.path.isfile(args.blacklist_bed):
            sys.exit(
                f"--blacklist-bed file not found: {args.blacklist_bed}\n\n"
                f"If you're running in Docker, remember this path is looked up "
                f"INSIDE the container -- a host path only exists there if you "
                f"mounted it with -v. Either copy the blacklist file into an "
                f"already-mounted directory (e.g. alongside your peaks under "
                f"/data) or add another -v mount for wherever it lives, then "
                f"point --blacklist-bed at the container-side path."
            )
        try:
            blacklist_trees = load_blacklist(args.blacklist_bed)
        except (FileNotFoundError, ValueError, pd.errors.ParserError) as exc:
            sys.exit(f"Could not load --blacklist-bed {args.blacklist_bed}: {exc}")

    metadata_map = None
    if args.metadata_file:
        if not os.path.isfile(args.metadata_file):
            sys.exit(
                f"--metadata-file file not found: {args.metadata_file}\n\n"
                f"If you're running in Docker, remember this path is looked up "
                f"INSIDE the container -- a host path only exists there if you "
                f"mounted it with -v. Either copy the metadata file into an "
                f"already-mounted directory (e.g. alongside your peaks under "
                f"/data) or add another -v mount for wherever it lives, then "
                f"point --metadata-file at the container-side path."
            )
        try:
            metadata_map = load_metadata(args.metadata_file)
        except (FileNotFoundError, ValueError, pd.errors.ParserError) as exc:
            sys.exit(f"Could not load --metadata-file {args.metadata_file}: {exc}")

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

    print(f"Loading {len(input_files)} scored peak file(s)...", flush=True)
    if coord_files:
        print(f"Loading {len(coord_files)} coordinate-only peak file(s) (lower-priority tier)...", flush=True)
    peaks_df = build_peak_pool(input_files, coord_files)
    print(f"Total peaks pooled: {len(peaks_df)}", flush=True)

    unified_df, mapping = run_iterative_selection(peaks_df, args.window, genome_sizes)
    print(f"Unified peaks produced: {len(unified_df)}", flush=True)

    if not args.no_qc:
        n_files_in_pass = len(input_files) + len(coord_files)
        min_max_other_files = (
            args.min_max_other_files
            if args.min_max_other_files is not None
            else n_files_in_pass / 2
        )
        qc_dir = os.path.join(args.outdir, "qc")
        try:
            summary, distribution, auto_excluded, summary_path, dist_path = run_qc(
                peaks_df, mapping, qc_dir,
                args.min_pct_overlap_other_file, args.min_median_other_files,
                min_max_other_files,
            )
        except OSError as exc:
            _explain_write_error(exc, args.outdir)
        print()
        print(summary.to_string(index=False))
        print(f"\nQC summary:      {summary_path}")
        print(f"QC distribution: {dist_path}")

        if auto_excluded:
            print(
                f"\n{len(auto_excluded)} file(s) failed QC thresholds "
                f"(pct_overlap_other_file < {args.min_pct_overlap_other_file}, "
                f"median_other_files < {args.min_median_other_files}, or "
                f"max_other_files < {min_max_other_files:g} "
                f"[{'explicit' if args.min_max_other_files is not None else 'auto: half of ' + str(n_files_in_pass) + ' files'}]):"
            )
            for f in auto_excluded:
                print(f"  {f}")
            excluded_path = os.path.join(qc_dir, "excluded_files.txt")
            try:
                with open(excluded_path, "w") as fh:
                    fh.write("\n".join(auto_excluded) + "\n")
            except OSError as exc:
                _explain_write_error(exc, args.outdir)
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

    if blacklist_trees is not None:
        print(f"\nApplying blacklist filter: {args.blacklist_bed}")
        unified_df, mapping, removed_df = filter_unified_peaks(unified_df, mapping, blacklist_trees)
        peaks_df = peaks_df[peaks_df["global_id"].isin(mapping.keys())].reset_index(drop=True)
        print(f"Removed {len(removed_df)} unified peak(s) overlapping the blacklist.")
        if len(removed_df):
            try:
                os.makedirs(args.outdir, exist_ok=True)
                removed_path = os.path.join(args.outdir, "blacklist_removed_peaks.bed")
                removed_df.sort_values(["chrom", "start"]).to_csv(
                    removed_path, sep="\t", header=True, index=False
                )
            except OSError as exc:
                _explain_write_error(exc, args.outdir)
            print(f"Removed peaks written to: {removed_path}")
        print(f"Unified peaks remaining: {len(unified_df)}")
        if unified_df.empty:
            sys.exit(
                "Every unified peak was removed by --blacklist-bed -- nothing "
                "left to write. Check that the blacklist file and your peaks "
                "are on the same genome build, and that it isn't unexpectedly "
                "broad relative to your peak set."
            )

    if metadata_map is not None:
        unified_df, unmatched = annotate_biosamples(unified_df, metadata_map)
        if unmatched:
            preview = ", ".join(unmatched[:10]) + ("..." if len(unmatched) > 10 else "")
            print(
                f"\nWarning: {len(unmatched)} file(s) referenced in overlapping_files "
                f"had no match in --metadata-file (left out of biosamples): {preview}"
            )

    print(f"\nApplying assay-count filter: min {args.min_assays} assay(s)")
    unified_df, mapping, removed_assay_df = filter_by_assay_count(
        unified_df, mapping, args.min_assays
    )
    peaks_df = peaks_df[peaks_df["global_id"].isin(mapping.keys())].reset_index(drop=True)
    print(
        f"Removed {len(removed_assay_df)} unified peak(s) found in fewer than "
        f"{args.min_assays} assay(s)."
    )
    if len(removed_assay_df):
        try:
            os.makedirs(args.outdir, exist_ok=True)
            removed_assay_path = os.path.join(args.outdir, "low_assay_count_removed_peaks.bed")
            removed_assay_df.sort_values(["chrom", "start"]).to_csv(
                removed_assay_path, sep="\t", header=True, index=False
            )
        except OSError as exc:
            _explain_write_error(exc, args.outdir)
        print(f"Removed peaks written to: {removed_assay_path}")
    print(f"Unified peaks remaining: {len(unified_df)}")
    if unified_df.empty:
        sys.exit(
            "Every unified peak was removed by --min-assays -- nothing left "
            "to write. Lower --min-assays."
        )

    try:
        unified_path, combined_path = write_outputs(peaks_df, unified_df, mapping, args.outdir)
    except OSError as exc:
        _explain_write_error(exc, args.outdir)
    print(f"\nUnified peak set:      {unified_path}")
    print(f"Combined mapping:      {combined_path}")
    print(f"Per-file mapped beds:  {os.path.join(args.outdir, 'mapping')}/")
    print("=" * 70, flush=True)
    print("Run completed successfully.", flush=True)
    print("=" * 70, flush=True)


if __name__ == "__main__":
    try:
        main()
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception:
        print("\n" + "=" * 70, file=sys.stderr, flush=True)
        print("FATAL: align_peaks.py crashed with an unhandled exception.",
              file=sys.stderr, flush=True)
        print("=" * 70, file=sys.stderr, flush=True)
        traceback.print_exc()
        sys.stderr.flush()
        sys.exit(1)
