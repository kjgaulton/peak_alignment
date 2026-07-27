"""
blacklist.py

Removes unified peaks that overlap ENCODE (or any other) blacklist
regions -- genomic regions known to produce artifactual high signal
across many unrelated experiments (repeats, centromeres/telomeres,
high-mappability artifacts, etc.) and that are standard practice to
exclude from downstream peak analysis.

This module doesn't bundle blacklist region data itself -- ENCODE
blacklists are maintained externally (see below) and change over time, so
the standard approach (matching e.g. nf-core/atacseq) is to point at a
BED file you supply via --blacklist-bed rather than have this pipeline
guess at a version. Official ENCODE blacklist v2 BED files (hg38, hg19,
mm10, etc.) are published by the Boyle lab:

    https://github.com/Boyle-Lab/Blacklist/tree/master/lists

e.g. for hg38:
    wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
    gunzip hg38-blacklist.v2.bed.gz
"""
import os

import pandas as pd
from intervaltree import IntervalTree


def load_blacklist(path):
    """Reads a blacklist BED file (chrom, start, end, ...extra columns
    ignored) into a dict of chrom -> IntervalTree of blacklisted
    intervals."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Blacklist BED file not found: {path}")
    df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2])
    df.columns = ["chrom", "start", "end"]

    trees = {}
    for chrom, sub in df.groupby("chrom"):
        tree = IntervalTree()
        for row in sub.itertuples():
            if row.end > row.start:
                tree[row.start:row.end] = True
        trees[chrom] = tree
    return trees


def filter_unified_peaks(unified_df, mapping, blacklist_trees):
    """Removes any unified peak that overlaps a blacklisted interval on
    its chromosome. Returns (kept_unified_df, kept_mapping, removed_df).

    kept_mapping drops any original-peak -> unified_id entries whose
    unified peak was removed, so downstream mapping outputs never point
    at a blacklisted/removed unified peak.
    """
    # Grouped, itertuples-based scan instead of DataFrame.apply(axis=1) --
    # apply is a per-row Python call and gets noticeably slow once the
    # unified set reaches the millions of rows seen on real datasets.
    blacklisted_flags = pd.Series(False, index=unified_df.index)
    for chrom, sub in unified_df.groupby("chrom"):
        tree = blacklist_trees.get(chrom)
        if tree is None:
            continue
        for row in sub.itertuples():
            if tree.overlap(row.start, row.end):
                blacklisted_flags.at[row.Index] = True

    removed_df = unified_df.loc[blacklisted_flags].reset_index(drop=True)
    kept_df = unified_df.loc[~blacklisted_flags].reset_index(drop=True)

    removed_ids = set(removed_df["unified_id"])
    kept_mapping = {
        gid: uid for gid, uid in mapping.items() if uid not in removed_ids
    }
    return kept_df, kept_mapping, removed_df
