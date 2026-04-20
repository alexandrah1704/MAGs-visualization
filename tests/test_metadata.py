from __future__ import annotations
import pandas as pd

# Bitte Import anpassen
from mags_visualization.heatmap import process_metadata_column


def test_process_metadata_column_boolean() -> None:
    meta = pd.Series(["yes", "no", "yes", None], index=["s1", "s2", "s3", "s4"])
    per_sample, colors, title = process_metadata_column(meta, "infection")

    assert title == "infection"
    assert "Unknown infection" in per_sample.values()
    assert "yes" in per_sample.values()
    assert "no" in per_sample.values()
    assert len(colors) >= 2


def test_process_metadata_column_numeric() -> None:
    meta = pd.Series([10, 12, 15, None], index=["s1", "s2", "s3", "s4"])
    per_sample, colors, title = process_metadata_column(meta, "temperature", meta_bin_width=5)

    assert "temperature" in title
    assert "Unknown temperature" in per_sample.values()
    assert len(colors) >= 1


def test_process_metadata_column_categorical() -> None:
    meta = pd.Series(["sunny", "rainy", "sunny"], index=["s1", "s2", "s3"])
    per_sample, colors, title = process_metadata_column(meta, "weather")

    assert title == "weather"
    assert set(per_sample.values()) == {"sunny", "rainy"}
    assert len(colors) == 2


def test_process_metadata_column_empty_strings_unknown() -> None:
    meta = pd.Series(["", "sunny", None], index=["s1", "s2", "s3"])
    per_sample, _, title = process_metadata_column(meta, "weather")

    assert title == "weather"
    assert per_sample["s1"] == "Unknown weather"
    assert per_sample["s3"] == "Unknown weather"