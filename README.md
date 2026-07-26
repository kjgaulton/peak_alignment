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

- `--input`: one or more narrowPeak files or glob patterns (tier 0, scored).
- `--coord-input`: optional coordinate-only BED files or glob patterns
  (tier 1, lower priority — see Algorithm above).
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
- `results/mapping/<original_file>.mapped.bed` — one per input file: every
  original peak plus the unified peak it was assigned to
  (`unified_id, unified_chrom, unified_start, unified_end`), along with its
  own `tier` and `width`.
- `results/mapping/all_peaks_mapped.tsv` — same mapping for all input files
  combined, with a `file_name` column.

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
