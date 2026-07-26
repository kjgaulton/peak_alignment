# Peak Alignment Pipeline

Builds a single unified, non-overlapping peak set across multiple MACS2
`narrowPeak` files, and maps every original peak back to the unified peak
it was collapsed into.

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

With a lower-priority tier of coordinate-only BED files mixed in:

```bash
python3 align_peaks.py \
    --input peaks/*.narrowPeak \
    --coord-input extra_peaks/*.bed \
    --outdir results \
    --window 300 \
    --genome hg38
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

The image's entrypoint is `align_peaks.py`; to run `qc_peaks.py` instead,
override the entrypoint (same `--user`/volume conventions apply):

```bash
docker run --rm \
    --user $(id -u):$(id -g) \
    --entrypoint python3 \
    -v /path/to/output:/results \
    peak-alignment \
    qc_peaks.py \
    --mapping-tsv /results/mapping/all_peaks_mapped.tsv \
    --outdir /results/qc \
    --flag-below 50
```

- `--input`: one or more narrowPeak files or glob patterns (tier 0, scored).
- `--coord-input`: optional coordinate-only BED files or glob patterns
  (tier 1, lower priority — see Algorithm above).
- `--exclude`: file basenames or glob patterns to drop from `--input`/
  `--coord-input` before pooling — see the QC workflow below.
- `--outdir`: output directory (created if missing).
- `--window`: fixed window width in bp (default 300).
- `--genome`: `hg38`, `hg19`, or `mm10` — used only to clip windows at
  chromosome ends.

## Output

- `results/unified_peaks.bed` — the final consensus peak set. Columns:
  `chrom, start, end, unified_id, seed_tier, seed_signalValue, seed_width,
  seed_file, seed_peak_name, n_peaks_merged`. `seed_tier` is `0` if the
  window came from a scored peak, `1` if it came from a coordinate-only
  peak.
- `results/mapping/<original_file>.mapped.bed` — one per input file, one row
  per original peak in its original order. Headerless BED4: `chrom, start,
  end, unified_id` — the coordinates of the unified peak that original peak
  maps to (its own original coordinates/score/etc. are not included).
- `results/mapping/all_peaks_mapped.tsv` — same mapping for all input files
  combined: `file_name, chrom, start, end, unified_id`.

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

## Quality metrics per input file (qc_peaks.py)

Once you've run `align_peaks.py`, `qc_peaks.py` scores each original input
file by how reproducible its peaks are, using the mapping it already
produced (no need to recompute overlaps):

```bash
python3 qc_peaks.py \
    --mapping-tsv results/mapping/all_peaks_mapped.tsv \
    --outdir results/qc \
    --flag-below 50
```

For every original peak, we already know which unified window it landed
in and, from that, which *other* input files also contributed a peak to
that same window. That gives two per-file metrics:

- **pct_overlap_other_file** — percent of a file's peaks whose unified
  window also contains a peak from at least one other file. A file
  scoring low here is contributing peaks nobody else supports, which is a
  reasonable proxy for a poor-quality replicate.
- **the distribution of "how many other files back up each peak"** — for
  a given file, how many of its peaks are supported by 0 other files, how
  many by 1, by 2, and so on (`file_overlap_distribution.tsv`).

Output:

- `results/qc/file_quality_summary.tsv` — one row per file: `file_name,
  n_peaks, n_unified_peaks, n_supported_by_other, pct_overlap_other_file,
  mean_other_files, median_other_files, max_other_files`. Sorted worst
  (lowest `pct_overlap_other_file`) first.
- `results/qc/file_overlap_distribution.tsv` — long format: `file_name,
  n_other_files, n_peaks, pct_of_file_peaks`.
- `results/qc/flagged_files.txt` — only written if `--flag-below` is
  given: one file_name per line, for files under that percentage.

### Regenerating the unified set without low-quality files

Review the summary table (and decide on a threshold — `--flag-below` is
just a suggestion, not an automatic cutoff), then rerun `align_peaks.py`
excluding the files you don't want:

```bash
python3 align_peaks.py \
    --input peaks/*.narrowPeak \
    --exclude $(cat results/qc/flagged_files.txt) \
    --outdir results_v2 \
    --window 300 --genome hg38
```

`--exclude` also accepts glob patterns and works the same way for
`--coord-input`.
