"""
Generate data/experimental/experimental.tsv from data/experimental/AP_RADD_merged.tsv.

Output columns:
    chr, start, end, RADD, AP, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, A, T, C, G

    AP          = mean(BW1, BW2, BW3)
    vipper      = mean vIP double-strand score across all valid trinucleotide positions
    7.x columns = number of trinucleotide positions in the region whose vIP
                  double-strand score is <= that threshold
    A/T/C/G     = nucleotide counts across the full region
"""
from viper import Viper
import polars as pl
from pyfaidx import Fasta
from pathlib import Path
from tqdm import tqdm

HOME_DIR = Path(__file__).resolve().parent.parent
INPUT_FILE  = HOME_DIR / "data" / "experimental" / "AP_RADD_merged.tsv"
OUTPUT_FILE = HOME_DIR / "data" / "experimental" / "experimental.tsv"
FASTA_DIR   = HOME_DIR / "hg38" / "chromosomes"

MOTIF_LEN  = 3
THRESHOLDS = ["7.2", "7.3", "7.4", "7.5", "7.6", "7.7", "7.8"]

# Precompute vIP double-strand scores for all 64 valid trinucleotides
_viper = Viper()
VIP_CACHE: dict[str, float] = {
    a + b + c: float(_viper.compute_double_strand_score(a + b + c))
    for a in "ACGT" for b in "ACGT" for c in "ACGT"
}


def process_region(seq: str) -> dict:
    """Return vIP well counts and nucleotide counts for a region sequence."""

    counts = {t: 0 for t in THRESHOLDS}
    for i in range(len(seq) - MOTIF_LEN + 1):
        score = VIP_CACHE.get(seq[i : i + MOTIF_LEN])
        if score is None:
            continue  # skip positions with N or other ambiguous bases
        for t in THRESHOLDS:
            if score <= float(t):
                counts[t] += 1
    vipper = None if "N" in seq else _viper.compute_double_strand_score(seq)
    return {"vipper": vipper, "counts": counts}


data = pl.read_csv(INPUT_FILE, separator="\t", schema_overrides={"chr": pl.Utf8})

all_results = []

for chrom in sorted(data["chr"].unique().to_list()):
    fasta = Fasta(FASTA_DIR / f"{chrom}.fa")
    chr_data = data.filter(pl.col("chr") == chrom)

    for row in tqdm(chr_data.iter_rows(named=True), total=len(chr_data), desc=chrom):
        start, end = row["start"], row["end"]
        seq = fasta[chrom][start:end].seq.upper()

        result = process_region(seq)
        ap = (row["BW1"] + row["BW2"] + row["BW3"]) / 3

        all_results.append({
            "chr":    chrom,
            "start":  start,
            "end":    end,
            "RADD":   row["RADD"],
            "AP":     ap,
            "vipper": result["vipper"],
            **result["counts"],
            "A": seq.count("A"),
            "T": seq.count("T"),
            "C": seq.count("C"),
            "G": seq.count("G"),
        })

if __name__ == "__main__":
    df = pl.DataFrame(all_results)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(OUTPUT_FILE, separator="\t")
    print(f"Written {len(df)} rows to {OUTPUT_FILE}")