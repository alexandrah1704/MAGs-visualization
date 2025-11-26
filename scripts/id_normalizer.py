import re
import pandas as pd

GENOME_PATTERN = re.compile(r"(SRR\d+_bin_\d+)", re.I)

def normalize_genome_id(s: str) -> str:
    """
    Universal normalizing for MAG/Genome-IDs.
    end-format: srrXXXXXXX_bin_YYYYY_fasta
    """
    s = str(s).strip()

    # remove path
    s = re.sub(r".*[\\/]", "", s)

    # BAKTA remove _count
    s = re.sub(r"_Count$", "", s, flags=re.I)

    # QUAST "SRR...fasta_SRR...fasta" -> SRR..._bin_...
    m = GENOME_PATTERN.search(s)
    if m:
        s = m.group(1)
    else:
        s = s.replace(".", "_").replace("-", "_")
        s = re.sub(r"(?i)(\.fa(sta)?|\.fna|\.faa|\.gbff|\.gz)$", "", s)
    
    s = re.sub(r"_+", "_", s).strip("_")

    if not s.lower().endswith("_fasta"):
        s = s + "_fasta"

    return s.lower()


def normalize_genome_series(s: pd.Series) -> pd.Series:
    """Panda series"""
    return s.astype(str).map(normalize_genome_id)