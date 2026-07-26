# Peak Alignment Pipeline

Builds a single unified, non-overlapping peak set across multiple MACS2
`narrowPeak` files, and maps every original peak back to the unified peak
it was collapsed into.

## Algorithm

1. Pool every peak from every input file together.
2. While peaks remain in the pool:
   - Take the peak with the strongest signal (`signalValue`, narrowPeak
     column 7).
   - Build a fixed-width window (default 300bp) centered on that peak's
     summit (`chromStart` + summit offset, narrowPeak column 10), clipped
     to chromosome boundaries.
   - Remove that peak and every other pooled peak whose original interval
     overlaps it. All of them are mapped to the new unified window.
3. Repeat until no peaks remain. The windows created in step 2 are the
   unified peak set, and by construction no two of them overlap.

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

## Usage (Docker)

Build the image once:

```bash
docker build -t peak-alignment .
```

Run it, mounting your peaks directory to `/data` and an output directory to
`/results`:

```bash
docker run --rm \
    -v /path/to/peaks:/data \
    -v /path/to/output:/results \
    peak-alignment \
    --input "/data/*.narrowPeak" \
    --outdir /results \
    --window 300 \
    --genome hg38
```

Quote the `--input` glob so it's expanded inside the container, not by your
local shell.

- `--input`: one or more narrowPeak files or glob patterns.
- `--outdir`: output directory (created if missing).
- `--window`: fixed window width in bp (default 300).
- `--genome`: `hg38`, `hg19`, or `mm10` — used only to clip windows at
  chromosome ends.

## Output

- `results/unified_peaks.bed` — the final consensus peak set. Columns:
  `chrom, start, end, unified_id, seed_signalValue, seed_file,
  seed_peak_name, n_peaks_merged`.
- `results/mapping/<original_file>.mapped.bed` — one per input file: every
  original peak plus the unified peak it was assigned to
  (`unified_id, unified_chrom, unified_start, unified_end`).
- `results/mapping/all_peaks_mapped.tsv` — same mapping for all input files
  combined, with a `file_name` column.

## Notes

- Input files must be standard 10-column narrowPeak (chrom, start, end,
  name, score, strand, signalValue, pValue, qValue, summit_offset). If a
  peak has no called summit (`summit_offset == -1`), the interval midpoint
  is used instead.
- Chromosome sizes for clipping live in `chrom_sizes.py`; add a build there
  if you need one beyond hg38/hg19/mm10.
