"""
Usage: liftover_bw.py <chain.gz> <input.bw> <output.bw>

Lifts over a BigWig file using a chain file (via pyliftover + pyBigWig).
Overlapping intervals produced by the liftover are resolved by taking the mean.
"""
import sys
from collections import defaultdict
import pyBigWig
from pyliftover import LiftOver


def liftover_bw(chain_path, input_path, output_path):
    lo = LiftOver(chain_path)
    bw_in = pyBigWig.open(input_path)

    # Liftover all intervals, collect by chrom
    lifted = defaultdict(list)  # chrom → [(start, end, value), ...]
    for chrom, size in bw_in.chroms().items():
        intervals = bw_in.intervals(chrom)
        if not intervals:
            continue
        for start, end, value in intervals:
            new_start = lo.convert_coordinate(chrom, start)
            new_end   = lo.convert_coordinate(chrom, end - 1)  # end is exclusive
            if not new_start or not new_end:
                continue
            new_chrom, ns, *_ = new_start[0]
            _, ne, *_          = new_end[0]
            ne += 1  # restore exclusive end
            if ns >= ne:
                continue
            lifted[new_chrom].append((ns, ne, value))

    bw_in.close()

    # Sort and merge overlapping intervals (mean of overlapping values)
    chrom_sizes = {}
    merged = {}
    for chrom, intervals in lifted.items():
        intervals.sort()
        stack = []
        for start, end, value in intervals:
            if stack and start < stack[-1][1]:  # overlap
                ps, pe, pv, pc = stack[-1]
                overlap_end = max(pe, end)
                # weighted mean by length
                total_len = pe - ps + end - start
                new_val = (pv * (pe - ps) + value * (end - start)) / total_len
                stack[-1] = (ps, overlap_end, new_val, pc + 1)
            else:
                stack.append((start, end, value, 1))
        merged[chrom] = [(s, e, v) for s, e, v, _ in stack]
        chrom_sizes[chrom] = max(e for _, e, _ in merged[chrom])

    # Write output BigWig
    bw_out = pyBigWig.open(output_path, "w")
    bw_out.addHeader(sorted(chrom_sizes.items()))
    for chrom in sorted(merged):
        for start, end, value in merged[chrom]:
            bw_out.addEntries([chrom], [start], ends=[end], values=[value])
    bw_out.close()

    print(f"Written to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    liftover_bw(sys.argv[1], sys.argv[2], sys.argv[3])
