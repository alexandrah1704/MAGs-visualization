import re
import pandas as pd

GENOME_PATTERN = re.compile(r"([A-Za-z0-9_]+?)(?:\.[A-Za-z0-9]+)?(_bin_\d+)", re.I)

def normalize_genome_id(s: str) -> str:
    """
    Universal normalizing for MAG/Genome-IDs.
    Endformat: <id>_bin_<n>_fasta:
      - srr24759597_bin_73064_fasta
      - err9966616_bin_10_fasta
      - group_2_bin_200650_fasta
    """
    s = str(s).strip()

    s = re.sub(r".*[\\/]", "", s)

    s = re.sub(r"_Count$", "", s, flags=re.I)

    m = GENOME_PATTERN.search(s)
    if m:
        s = m.group(1) + m.group(2)
    else:
        s = s.replace(".", "_").replace("-", "_")
        s = re.sub(r"(?i)(\.fa(sta)?|\.fna|\.faa|\.gbff|\.gz)$", "", s)

    s = re.sub(r"_+", "_", s).strip("_")

    if not s.lower().endswith("_fasta"):
        s = s + "_fasta"

    return s.lower()


def normalize_genome_series(s: pd.Series) -> pd.Series:
    """Pandas Series normalizing"""
    return s.astype(str).map(normalize_genome_id)