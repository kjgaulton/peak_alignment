"""
metadata.py

Loads a simple file-to-biosample mapping (e.g. a 4DN/ENCODE-style
metadata.txt) used to annotate the unified peak set with which biosamples
are represented at each peak.

Expected format: a delimited text file (tab or comma; auto-detected) with
a header row containing at least a file-identifier column and a
biosample column. Column names are matched case- and
whitespace/punctuation-insensitively against common variants:
  - file identifier: "file_accession", "file accession", "accession",
    "file", "file_name", "filename"
  - biosample: "biosample_id", "biosample id", "biosample",
    "biosample_term_name", "biosample_name"

The file identifier is matched against peak file names with their
extension stripped (e.g. "4DNFIX5FL62U.bed" -> "4DNFIX5FL62U"), matching
how file names are reported elsewhere in this pipeline's output. Both the
raw value and the extension-stripped value are indexed, so it works
whether or not the metadata file's identifier column includes an
extension.
"""
import csv
import os

import pandas as pd

FILE_COL_CANDIDATES = ["file_accession", "accession", "file", "file_name", "filename"]
BIOSAMPLE_COL_CANDIDATES = [
    "biosample_id", "biosample", "biosample_term_name", "biosample_name",
]


def _normalize(col):
    return col.strip().lower().replace(" ", "_").replace("-", "_")


def _find_column(columns, candidates):
    normalized = {_normalize(c): c for c in columns}
    for cand in candidates:
        if cand in normalized:
            return normalized[cand]
    return None


def load_metadata(path):
    """Returns a dict mapping file identifier (raw and extension-stripped)
    -> biosample_id."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Metadata file not found: {path}")

    with open(path, "r", newline="") as fh:
        sample = fh.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        sep = dialect.delimiter
    except csv.Error:
        sep = "\t"

    df = pd.read_csv(path, sep=sep, dtype=str)
    file_col = _find_column(df.columns, FILE_COL_CANDIDATES)
    biosample_col = _find_column(df.columns, BIOSAMPLE_COL_CANDIDATES)
    if file_col is None or biosample_col is None:
        raise ValueError(
            f"Could not find file-identifier/biosample columns in {path}. "
            f"Found columns: {list(df.columns)}. Expected something like "
            f"'file_accession' and 'biosample_id'."
        )

    mapping = {}
    for _, row in df.iterrows():
        file_id = str(row[file_col]).strip()
        biosample = str(row[biosample_col]).strip()
        if not file_id or file_id.lower() == "nan":
            continue
        mapping[file_id] = biosample
        mapping[os.path.splitext(file_id)[0]] = biosample
    return mapping
