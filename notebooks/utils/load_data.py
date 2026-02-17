import polars as pl
from pathlib import Path


def get_kmer_df(hg38_dir: Path | str, context_size):
    if isinstance(hg38_dir, str):
        hg38_dir = Path(hg38_dir)
    return pl.scan_csv(
            hg38_dir / f"kmer_counts_{context_size}N.csv",
            separator=",",
        ).with_columns(
            pl.exclude("context") / pl.exclude("context").sum(),
            pl.sum_horizontal(pl.exclude("context").alias("total")),
        ).with_columns(frequencies=pl.col("total") / pl.col("total").sum()).collect()



def load_context_count(context_size, dbs, base_dir):
    return (
        pl.concat(
            [
                pl.read_csv(
                    Path(base_dir)
                    / db_name
                    / f"context_{context_size}N"
                    / "countext_count.csv",
                    separator="\t",
                )
                .unique()
                .with_columns(pl.lit(db_key).alias("db"))
                .drop("cancer_type", strict=False)
                for db_key, db_name in dbs.items()
            ]
        )
        .filter(~pl.col("small_context").str.contains("N"))
        .drop_nulls()
    )
