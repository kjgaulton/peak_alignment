# Peak Alignment Pipeline

Builds a single unified, non-overlapping peak set across multiple MACS2
`narrowPeak` files, maps every original peak back to the unified peak it
was collapsed into, and scores each input file's quality so low-quality
files can be automatically excluded from the final set.

## Algorithm

Peaks come in two priority tiers:

- Tier 0 ("scored") — narrowPeak entries, ranked by `signalValue`.
- Tier 1 ("coordinate-only") — plain BED entries with no score or summit
  (see `--coord-input` below), ranked by peak width, used only as a
  fallback when a real signal isn't available.

1. Pool every peak from every input file (both tiers) together.
2. Process peaks in priority order — every tier 0 peak (strongest signal
   first) before any tier 1 peak (widest first). While peaks remain:
   - Take the next unconsumed peak in that order as the seed.
   - Build a fixed-width window (default 300bp) centered on the seed's
     summit (narrowPeak: `chromStart` + summit offset, column 10;
     coordinate-only: interval midpoint), clipped to chromosome boundaries.
   - Remove the seed and every other pooled peak whose original interval
     overlaps it, regardless of tier. All of them are mapped to the new
     unified window.
3. Repeat until no peaks remain. The windows created in step 2 are the
   unified peak set, and by construction no two of them overlap.

Because tier 0 is always processed first, a coordinate-only peak that
overlaps any scored peak gets absorbed into that scored peak's window and
never gets a turn to seed its own — it only seeds a new unified window in
places no scored peak occupies.

## Quality control and auto-exclusion

After that first alignment pass, every original peak is mapped to a
unified window, which tells us which *other* input files also placed a
peak in that same window — a direct measure of reproducibility.
`align_peaks.py` scores every input file on that basis by default and
re-runs the alignment excluding any file that falls below either
threshold:

- **pct_overlap_other_file** (`--min-pct-overlap-other-file`, default
  **75**) — percent of a file's peaks whose unified window is also
  occupied by a peak from another file.
- **median_other_files** (`--min-median-other-files`, default **2**) —
  median number of other files backing up a given peak from this file.

A file is excluded if it fails *either* threshold. Use `--qc-report-only`
to compute and write the metrics without excluding anything, or `--no-qc`
to skip QC scoring entirely and just align once on all resolved input
files.

## Install

```bash
pip install -r requirements.txt
```

## Usage (local Python)

```bash
python3 align_peaks.py \
    --input peaks/*.narrowPeak \
    --outdir results \
    --window 300 \
    --genome hg38
```

With a lower-priority tier of coordinate-only BED files mixed in, and
custom QC thresholds:

```bash
python3 align_peaks.py \
    --input peaks/*.narrowPeak \
    --coord-input extra_peaks/*.bed \
    --outdir results \
    --window 300 \
    --genome hg38 \
    --min-pct-overlap-other-file 75 \
    --min-median-other-files 2
```

## Usage (Docker)

Build the image once:

```bash
docker build -t peak-alignment .
```

Run it, mounting your peaks directory to `/data` and an output directory to
`/results`. Create the output directory first and run as your own
user/group with `--user` — otherwise, if `/results` is on an NFS mount
(common on HPC/lab filesystems), root-squash will make writes from the
container's root user fail with a `PermissionError`:

```bash
mkdir -p /path/to/output

docker run --rm \
    --user $(id -u):$(id -g) \
    -v /path/to/peaks:/data \
    -v /path/to/output:/results \
    peak-alignment \
    --input "/data/*.narrowPeak" \
    --outdir /results \
    --window 300 \
    --genome hg38
```

Quote the `--input`/`--coord-input` globs so they're expanded inside the
container, not by your local shell.

## Arguments

- `--input`: one or more narrowPeak files or glob patterns (tier 0, scored).
- `--coord-input`: optional coordinate-only BED files or glob patterns
  (tier 1, lower priority — see Algorithm above).
- `--exclude`: file basenames or glob patterns to manually drop from
  `--input`/`--coord-input` before pooling (independent of QC
  auto-exclusion).
- `--min-pct-overlap-other-file`: QC threshold, default 75 (see above).
- `--min-median-other-files`: QC threshold, default 2 (see above).
- `--no-qc`: skip QC scoring and auto-exclusion entirely.
- `--qc-report-only`: compute/write QC metrics but don't exclude anything.
- `--outdir`: output directory (created if missing).
- `--window`: fixed window width in bp (default 300).
- `--genome`: `hg38`, `hg19`, or `mm10` — used only to clip windows at
  chromosome ends.

## Output

- `results/unified_peaks.bed` — the final consensus peak set, reflecting
  the QC-filtered file set (unless `--no-qc`/`--qc-report-only`). Columns:
  `chrom, start, end, unified_id, seed_tier, seed_signalValue, seed_width,
  seed_file, seed_peak_name, n_peaks_merged`. `seed_tier` is `0` if the
  window came from a scored peak, `1` if it came from a coordinate-only
  peak.
- `results/mapping/<original_file>.mapped.bed` — one per input file, one row
  per original peak in its original order. Headerless BED4: `chrom, start,
  end, unified_id` — the coordinates of the unified peak that original peak
  maps to (its own original coordinates/score/etc. are not included). Only
  covers files that passed QC (or all files if `--no-qc`/`--qc-report-only`).
- `results/mapping/all_peaks_mapped.tsv` — same mapping for all covered
  input files combined: `file_name, chrom, start, end, unified_id`.
- `results/qc/file_quality_summary.tsv` — one row per file **from the
  initial, pre-exclusion pass**: `file_name, n_peaks, n_unified_peaks,
  n_supported_by_other, pct_overlap_other_file, mean_other_files,
  median_other_files, max_other_files`. Sorted worst first. Not written if
  `--no-qc`.
- `results/qc/file_overlap_distribution.tsv` — long format: `file_name,
  n_other_files, n_peaks, pct_of_file_peaks`.
- `results/qc/excluded_files.txt` — one file_name per line, for files that
  failed a QC threshold (whether or not they were actually excluded from
  the final set — see `--qc-report-only`).

## Notes

- `--input` files must be standard 10-column narrowPeak (chrom, start, end,
  name, score, strand, signalValue, pValue, qValue, summit_offset). If a
  peak has no called summit (`summit_offset == -1`), the interval midpoint
  is used instead.
- `--coord-input` files are plain BED with just chrom, start, end (a 4th
  name column is used if present). They carry no score or summit, so they
  are ranked within their own tier by peak width and centered on the
  interval midpoint. They are always subordinate to `--input` peaks: a
  coordinate-only peak only seeds a unified window in a region no scored
  peak's original interval touches.
- Chromosome sizes for clipping live in `chrom_sizes.py`; add a build there
  if you need one beyond hg38/hg19/mm10.
- QC is a single pass: files are scored once against the initial full-data
  alignment, then removed all at once. It isn't recursive (removing a
  file doesn't trigger re-scoring and possibly excluding others).
- If all input files fail the QC thresholds, the script exits with an
  error rather than producing an empty unified set — loosen the
  thresholds or use `--qc-report-only`.
- The default thresholds assume a reasonably sized cohort. With very few
  input files, `median_other_files` can't exceed `n_files - 1`, so small
  file sets (e.g. 2-4 files, or files with several unique/singleton
  peaks) may fail the default `--min-median-other-files 2` even when
  everything looks fine — lower the threshold or use `--qc-report-only`
  in that case.

## Standalone QC re-analysis (qc_peaks.py)

`align_peaks.py` computes and applies QC automatically, so you normally
don't need to run this separately. `qc_peaks.py` is useful when you want
to recompute or re-threshold metrics from an *existing*
`all_peaks_mapped.tsv` without re-running the alignment (e.g. an older
result set, or experimenting with different thresholds):

```bash
python3 qc_peaks.py \
    --mapping-tsv results/mapping/all_peaks_mapped.tsv \
    --outdir results/qc_manual \
    --flag-below 50
```

This writes the same `file_quality_summary.tsv` / `file_overlap_distribution.tsv`,
plus (if `--flag-below` is given) a `flagged_files.txt` you can feed
straight into `align_peaks.py --exclude`.

The Docker image's entrypoint is `align_peaks.py`; to run `qc_peaks.py`
standalone, override the entrypoint:

```bash
docker run --rm \
    --user $(id -u):$(id -g) \
    --entrypoint python3 \
    -v /path/to/output:/results \
    peak-alignment \
    qc_peaks.py \
    --mapping-tsv /results/mapping/all_peaks_mapped.tsv \
    --outdir /results/qc_manual \
    --flag-below 50
```
