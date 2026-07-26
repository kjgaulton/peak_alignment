"""
qc_metrics.py

Shared peak-file quality metrics, based on how peaks from a given input
file land in a unified peak set. Used both by align_peaks.py (to compute
metrics inline and auto-exclude low-quality files during a run) and by
qc_peaks.py (to recompute/report metrics standalone from an existing
mapping table).

All functions operate on a dataframe with (at least) one row per original
peak and columns "file_name" and "unified_id".
"""
import pandas as pd


def compute_metrics(df):
    """Adds n_other_files / supported_by_other columns.

    n_other_files for a given peak = (number of distinct files with at
    least one peak in that peak's unified window) - 1. This is constant
    across every row belonging to the same file at the same unified_id,
    even if that file contributed more than one of its own peaks there.
    """
    n_files_at_unified = (
        df.groupby("unified_id")["file_name"].nunique().rename("n_files_at_unified")
    )
    df = df.merge(n_files_at_unified, on="unified_id", how="left")
    df["n_other_files"] = df["n_files_at_unified"] - 1
    df["supported_by_other"] = df["n_other_files"] >= 1
    return df


def build_summary(df):
    """One row per file_name: n_peaks, n_unified_peaks,
    n_supported_by_other, pct_overlap_other_file, mean/median/max
    n_other_files. Sorted worst (lowest pct_overlap_other_file) first."""
    summary = df.groupby("file_name").agg(
        n_peaks=("unified_id", "size"),
        n_unified_peaks=("unified_id", "nunique"),
        n_supported_by_other=("supported_by_other", "sum"),
        mean_other_files=("n_other_files", "mean"),
        median_other_files=("n_other_files", "median"),
        max_other_files=("n_other_files", "max"),
    ).reset_index()
    summary["pct_overlap_other_file"] = (
        100 * summary["n_supported_by_other"] / summary["n_peaks"]
    )
    summary = summary[
        ["file_name", "n_peaks", "n_unified_peaks", "n_supported_by_other",
         "pct_overlap_other_file", "mean_other_files", "median_other_files",
         "max_other_files"]
    ]
    return summary.sort_values("pct_overlap_other_file").reset_index(drop=True)


def build_distribution(df):
    """Long format: file_name, n_other_files, n_peaks, pct_of_file_peaks."""
    dist = (
        df.groupby(["file_name", "n_other_files"])
        .size()
        .rename("n_peaks")
        .reset_index()
    )
    dist["pct_of_file_peaks"] = (
        100 * dist["n_peaks"] / dist.groupby("file_name")["n_peaks"].transform("sum")
    )
    return dist.sort_values(["file_name", "n_other_files"]).reset_index(drop=True)


def flag_low_quality_files(summary, min_pct_overlap_other_file, min_median_other_files):
    """Returns the list of file_names in `summary` that fail either
    threshold: pct_overlap_other_file below min_pct_overlap_other_file, or
    median_other_files below min_median_other_files."""
    fails = (
        (summary["pct_overlap_other_file"] < min_pct_overlap_other_file)
        | (summary["median_other_files"] < min_median_other_files)
    )
    return summary.loc[fails, "file_name"].tolist()
