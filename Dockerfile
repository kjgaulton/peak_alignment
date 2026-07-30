FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY align_peaks.py chrom_sizes.py qc_peaks.py qc_metrics.py blacklist.py plot_qc.py metadata.py ./

# /data is where input narrowPeak files should be mounted, /results is where
# output gets written. Both are just conventions -- override with -v. If
# using --blacklist-bed or --metadata-file, mount that file's directory too
# (e.g. under /data), since neither is bundled in the image -- see
# blacklist.py / metadata.py.
RUN mkdir -p /data /results

ENTRYPOINT ["python3", "align_peaks.py"]
CMD ["--help"]
