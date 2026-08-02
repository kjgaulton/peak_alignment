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

1. Pool every peak from every input file (both tiers) together. For each
   peak, precompute the fixed-width window (default 300bp) it would get if
   chosen as a seed: centered on its own summit (narrowPeak: `chromStart` +
   summit offset, column 10; coordinate-only: interval midpoint), clipped
   to chromosome boundaries.
2. Process peaks in priority order — every tier 0 peak (strongest signal
   first) before any tier 1 peak (widest first). While peaks remain:
   - Take the next unconsumed peak in that order as the seed; its
     precomputed window is the unified window.
   - Remove the seed and every other pooled peak whose own precomputed
     window overlaps the seed's window, regardless of tier. All of them
     are mapped to the new unified window.
3. Repeat until no peaks remain. The windows chosen in step 2 are the
   unified peak set, and by construction no two of them overlap.

Clustering is judged on each peak's own fixed-width summit window, not its
raw `chromStart`/`chromEnd` — that's what actually guarantees the final
set is non-overlapping (two raw peaks that don't touch can still have
summits close enough that their windows would overlap). It's also how
multi-summit peaks are handled: if a caller like MACS2 `--call-summits`
reports several summits sharing one broad raw span, each summit is just
another peak in the pool, and two of them only collapse into the same
unified peak if their own 300bp windows actually overlap. Summits far
enough apart within the same broad call can each seed their own separate
unified peak; summits close together still collapse to whichever has the
strongest signal, exactly like any other overlap.

Because tier 0 is always processed first, a coordinate-only peak whose
window overlaps a scored peak's window gets absorbed into that scored
peak's window and never gets a turn to seed its own — it only seeds a new
unified window in places no scored peak's window touches.

## Quality control and auto-exclusion

After that first alignment pass, every original peak is mapped to a
unified window, which tells us which *other* input files also placed a
peak in that same window — a direct measure of reproducibility.
`align_peaks.py` scores every input file on that basis by default and
re-runs the alignment excluding any file that falls below any of three
thresholds:

- **pct_overlap_other_file** (`--min-pct-overlap-other-file`, default
  **75**) — percent of a file's peaks whose unified window is also
  occupied by a peak from another file.
- **median_other_files** (`--min-median-other-files`, default **2**) —
  median number of other files backing up a given peak from this file.
- **max_other_files** (`--min-max-other-files`, default **half the
  number of files entering this QC pass**) — number of other files
  backing up this file's single *best-supported* peak. This catches a
  file with no peak reproduced across a meaningful fraction of the
  cohort even when it happens to pass the other two — e.g. a file whose
  peaks are all consistently backed by exactly 2 other files will pass
  `median_other_files` but, in a 6-file cohort, still gets excluded here
  since 2 is short of the default bar of 3.

A file is excluded if it fails *any* of the three thresholds. Use
`--qc-report-only` to compute and write the metrics without excluding
anything, or `--no-qc` to skip QC scoring entirely and just align once on
all resolved input files.

## Blacklist filtering

If `--blacklist-bed` is given, any final unified peak (post-QC) that
overlaps a region in that BED file is dropped entirely, along with any
original peak mapped to it. This does *not* bundle ENCODE blacklist region
data — supply the file yourself. Official ENCODE blacklist v2 lists
(hg38/hg19/mm10) are published by the Boyle lab:

```bash
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
```

```bash
python3 align_peaks.py \
    --input peaks/*.narrowPeak \
    --outdir results \
    --blacklist-bed hg38-blacklist.v2.bed
```

## Assay-count filtering

Every unified peak gets an `n_assays` column: the count of distinct
assays/files in `overlapping_files` supporting it. Any peak backed by
fewer than `--min-assays` (default **2**) distinct assays is dropped
from the final set and written to `low_assay_count_removed_peaks.bed`
— by default this excludes peaks found in only a single assay. This
filter doesn't depend on `--metadata-file` — it's based purely on the
`overlapping_files` count. Pass `--min-assays 1` to disable it and keep
everything (`n_assays` is always >= 1, so 1 is a no-op threshold).

```bash
python3 align_peaks.py \
    --input peaks/*.narrowPeak \
    --outdir results \
    --min-assays 2
```

## Biosample annotation

If `--metadata-file` is given, `unified_peaks.bed` gets `biosamples`/
`n_biosamples` columns: the deduplicated biosample IDs of every
assay/file overlapping that peak, and how many there are. The metadata
file is a delimited text file (tab or comma, auto-detected) with a header
row containing a file-identifier column and a biosample column — see
`metadata.py` for the exact column names matched (e.g.
`file_accession`/`accession`/`file_name` and `biosample_id`/`biosample`).
The file identifier is matched against peak file names with their
extension stripped, so it works whether or not your metadata file's
identifier column includes one.

This is annotation only — it doesn't filter anything on its own (unlike
`--min-assays` above). It's there so you can see, per peak, whether
multiple *assays* supporting it actually come from the same underlying
biosample (a sign of pseudo-replication rather than true reproducibility)
and filter on `n_biosamples` yourself downstream if you want to.

```bash
python3 align_peaks.py \
    --input peaks/*.narrowPeak \
    --outdir results \
    --metadata-file metadata.txt
```

Any file referenced in `overlapping_files` that isn't found in the
metadata is left out of `biosamples`/`n_biosamples` for that peak and
reported as a warning at the end of the run.

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
    --min-median-other-files 2 \
    --min-max-other-files 3
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
container, not by your local shell. If using `--blacklist-bed`, make sure
that file lives somewhere already mounted (e.g. alongside your peaks under
`/data`), since blacklist data isn't bundled in the image:

```bash
docker run --rm \
    --user $(id -u):$(id -g) \
    -v /path/to/peaks:/data \
    -v /path/to/output:/results \
    peak-alignment \
    --input "/data/*.narrowPeak" \
    --blacklist-bed /data/hg38-blacklist.v2.bed \
    --outdir /results \
    --window 300 \
    --genome hg38
```

## Arguments

- `--input`: one or more narrowPeak files or glob patterns (tier 0, scored).
- `--coord-input`: optional coordinate-only BED files or glob patterns
  (tier 1, lower priority — see Algorithm above).
- `--exclude`: file basenames or glob patterns to manually drop from
  `--input`/`--coord-input` before pooling (independent of QC
  auto-exclusion).
- `--min-pct-overlap-other-file`: QC threshold, default 75 (see above).
- `--min-median-other-files`: QC threshold, default 2 (see above).
- `--min-max-other-files`: QC threshold, default half the number of files
  entering the QC pass (see above).
- `--no-qc`: skip QC scoring and auto-exclusion entirely.
- `--qc-report-only`: compute/write QC metrics but don't exclude anything.
- `--blacklist-bed`: BED file of regions to drop final unified peaks
  overlapping (e.g. ENCODE blacklist) — see Blacklist filtering above.
- `--metadata-file`: file mapping file accession -> biosample_id, used to
  add `biosamples`/`n_biosamples` columns — see Biosample annotation
  above.
- `--min-assays`: drop unified peaks found in fewer than this many
  distinct assays/files, default 2 (excludes single-assay peaks). See
  Assay-count filtering above.
- `--outdir`: output directory (created if missing).
- `--window`: fixed window width in bp (default 300).
- `--genome`: `hg38`, `hg19`, or `mm10` — used only to clip windows at
  chromosome ends.

## Output

- `results/unified_peaks.bed` — the final consensus peak set, reflecting
  the QC-filtered file set (unless `--no-qc`/`--qc-report-only`). Columns:
  `chrom, start, end, unified_id, seed_signalValue, seed_summit, seed_file,
  n_peaks_merged, overlapping_files, n_assays[, biosamples, n_biosamples]`.
  - `seed_summit`: the genomic position of the seed peak's summit (what
    the unified window is centered on).
  - `seed_file`: the file that seeded this peak (its own summit/window
    became the unified window), with its extension stripped.
  - `overlapping_files`: comma-separated, extension-stripped names of
    every file with a peak that merged into this unified window
    (includes `seed_file`).
  - `n_assays`: count of `overlapping_files`; peaks below `--min-assays`
    have already been dropped by this point.
  - `biosamples`/`n_biosamples`: comma-separated, deduplicated biosample
    IDs for every file in `overlapping_files` (and their count), looked
    up via `--metadata-file`. Only present if `--metadata-file` was given.
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
- `results/blacklist_removed_peaks.bed` — unified peaks dropped for
  overlapping `--blacklist-bed` (only written if given, and if at least
  one peak was removed).
- `results/low_assay_count_removed_peaks.bed` — unified peaks dropped for
  having fewer than `--min-assays` distinct assays/files (only written if
  at least one peak was removed; not written at all with
  `--min-assays 1`, since nothing gets removed in that case).

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
  input files, `median_other_files`/`max_other_files` can't exceed
  `n_files - 1`, so small file sets (e.g. 2-4 files, or files with
  several unique/singleton peaks) may fail the default
  `--min-median-other-files 2` or the auto `--min-max-other-files`
  (n_files / 2) even when everything looks fine — lower the thresholds
  or use `--qc-report-only` in that case.
- On very large runs (tens of millions of pooled peaks), the underlying
  `intervaltree` package can, rarely, fail to locate an interval it just
  reported via an overlap query when removing it (an internal AVL
  rebalancing edge case, not specific to this pipeline). This used to
  surface as a crash (`KeyError`/`ValueError` from `tree.remove()`);
  it's now handled by discarding such intervals safely and filtering
  seed-selection results against already-consumed peaks, so it can no
  longer crash the run or double-count a peak into two unified windows.
  If it happens, you'll see a one-line `Note: N stale interval-tree
  hit(s) were encountered and safely ignored` in the log — informational
  only, does not affect correctness of the output.

## Running in the background / troubleshooting silent failures

For large cohorts (100+ files, tens of millions of pooled peaks), a run can
take long enough that you'll want to launch it in the background over SSH
and check on it later. A few things make that reliable instead of a run
that dies with no visible cause:

- **Redirect output to a log file and check the exit code, don't just
  background it.** `command &` alone loses stdout/stderr and the exit
  status once you disconnect. Use `nohup` (or `screen`/`tmux`) and always
  redirect:
  ```
  nohup python3 align_peaks.py --input ... --outdir ... > align.log 2>&1 &
  echo $! > align.pid
  ```
  When it's done, check `echo $?` (if run in `screen`/`tmux`) or
  `tail align.log` — the run now always ends with either
  `Run completed successfully.` or a `FATAL:` banner with a full Python
  traceback, so the log unambiguously tells you which happened. A log that
  just stops mid-line with neither banner means the *process* was killed
  from outside Python (see OOM note below), not that the script silently
  gave up.
- **Progress logging.** The run now prints per-chromosome interval-tree
  build progress and periodic `...processed N/total candidates` lines
  during seed selection, so a truncated log still shows roughly how far
  the run got before it stopped.
- **Unbuffered output.** `align_peaks.py`'s print statements flush
  immediately (the Docker image also sets `PYTHONUNBUFFERED=1`), so
  `tail -f align.log` reflects real progress rather than batches of
  buffered output appearing all at once.
- **Out-of-memory kills leave no traceback at all** — the OS SIGKILLs the
  process before Python gets a chance to run any exception handler, so the
  log will simply stop without a `FATAL:` banner. If you're running under
  Docker, `docker inspect <container> --format '{{.State.ExitCode}}'`
  showing `137` confirms an OOM/SIGKILL. Under SLURM, check `sacct -j
  <jobid> --format=MaxRSS,State,ExitCode` or the scheduler's own OOM log.
  The fix is more memory or fewer files per run (e.g. pre-filtering with
  `qc_peaks.py`), not a script change.
- **Disk-full and other I/O errors during output writing** are now caught
  explicitly (previously only permission errors were), so running out of
  space on `--outdir`'s volume produces a clear
  `No space left on device`-style message pointing at the offending path
  instead of a bare crash.
- If you do hit a `FATAL:` traceback and it's not one of the above, that's
  a real bug — the full traceback in the log has everything needed to
  diagnose it, so please share it rather than just the summary line.

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
    --min-pct-overlap-other-file 75 \
    --min-median-other-files 2 \
    --min-max-other-files 3
```

Uses the same three thresholds and defaults as `align_peaks.py`
(`--min-max-other-files` auto-computes as half the number of distinct
files in the mapping table if omitted). Writes the same
`file_quality_summary.tsv` / `file_overlap_distribution.tsv`, plus a
`flagged_files.txt` (if anything is flagged) you can feed straight into
`align_peaks.py --exclude`.

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
    --outdir /results/qc_manual
```

## Unified-peak-set QC plots (plot_qc.py)

`qc_peaks.py` scores individual input *files*; `plot_qc.py` instead looks
at the unified peak set as a whole, using the same mapping output:

```bash
python3 plot_qc.py \
    --mapping-tsv results/mapping/all_peaks_mapped.tsv \
    --outdir results/qc_plots
```

Produces two plots (plus the underlying TSVs):

- **`peak_support_distribution.png`** — a histogram of how many distinct
  input files support each unified peak (i.e. placed at least one
  original peak in that window). Shows whether the unified set skews
  toward singleton/low-reproducibility peaks or a reproducible core.
- **`peaks_retained_by_min_support.png`** — the number of unified peaks
  that would remain if you required at least K supporting files, for
  every K from 1 up to the total number of input files. A step plot that
  answers "how much of the peak set survives at each reproducibility
  bar," computed directly from the existing mapping output (no
  re-running the alignment on file subsets).

Requires `matplotlib` (in `requirements.txt`). Via Docker, override the
entrypoint the same way as `qc_peaks.py`:

```bash
docker run --rm \
    --user $(id -u):$(id -g) \
    --entrypoint python3 \
    -v /path/to/output:/results \
    peak-alignment \
    plot_qc.py \
    --mapping-tsv /results/mapping/all_peaks_mapped.tsv \
    --outdir /results/qc_plots
```
