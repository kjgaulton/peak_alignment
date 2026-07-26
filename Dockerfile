FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY align_peaks.py chrom_sizes.py ./

# /data is where input narrowPeak files should be mounted, /results is where
# output gets written. Both are just conventions -- override with -v.
RUN mkdir -p /data /results

ENTRYPOINT ["python3", "align_peaks.py"]
CMD ["--help"]
