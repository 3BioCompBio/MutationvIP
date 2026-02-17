"""
Usage: merge_big_wig.py <template.bedgraph[.gz]> <out.tsv> <rep1.bw> [rep2.bw ...]

Merges one or more BigWig files onto a bedgraph template grid.
Each BigWig value is computed as a weighted average over each template region.
Output columns: chr, start, end, RADD, BW1, BW2, ...
"""
import sys
import gzip
import polars as pl
import pyBigWig


def get_region_value(bw, chrom, start, end):
    try:
        intervals = bw.intervals(chrom, start, end)
    except RuntimeError:
        return 0.0
    if not intervals:
        return 0.0
    region_size = end - start
    region_value = 0.0
    for start_interval, end_interval, value in intervals:
        inner = min(end_interval, end) - max(start_interval, start)
        region_value += value * inner / (end_interval - start_interval)
    return region_value / region_size


def load_template(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        first = f.readline()
    has_header = not first.split("\t")[1].lstrip("-").isdigit()
    return pl.read_csv(
        path,
        separator="\t",
        has_header=has_header,
        new_columns=["chr", "start", "end", "RADD"],
        schema={"chr": pl.Utf8, "start": pl.Int64, "end": pl.Int64, "RADD": pl.Float64},
    )


def get_bw_regions(bw, template):
    values = []
    for chrom, start, end, _ in template.iter_rows():
        values.append(get_region_value(bw, chrom.strip('chr'), start, end))
    return values


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    template_path = sys.argv[1]
    output_path = sys.argv[2]
    bw_paths = sys.argv[3:]

    template = load_template(template_path)

    new_cols = []
    for i, bw_path in enumerate(bw_paths, start=1):
        print(f"Processing {bw_path} ...", flush=True)
        bw = pyBigWig.open(bw_path)
        new_cols.append(pl.Series(f"BW{i}", get_bw_regions(bw, template)))
        bw.close()

    result = template.with_columns(new_cols)
    result.write_csv(output_path, separator="\t")
    print(f"Written to {output_path}")