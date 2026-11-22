"""Summarise the long-format CSV files of scripts/validate_behaviour.jl.

For every (exp, cond, snr, metric) prints the mean over field views x noise seeds and its
standard error over those units. Views whose name starts with "b." (super-Alfvenic runs)
are reported as a separate group.

    python3 scripts/summarize_behaviour.py validation_output/behaviour/*.csv [--metrics AM,dAM]
"""
import csv
import math
import sys
from collections import defaultdict


def main(argv):
    metrics = None
    files = []
    for a in argv:
        if a.startswith("--metrics="):
            metrics = set(a.split("=", 1)[1].split(","))
        else:
            files.append(a)
    groups = defaultdict(list)
    order = []
    for path in files:
        with open(path) as f:
            for row in csv.DictReader(f):
                if metrics and row["metric"] not in metrics:
                    continue
                grp = "superA" if row["view"].startswith("b.") else "subA"
                key = (row["exp"], grp, row["cond"], row["snr"], row["metric"])
                if key not in groups:
                    order.append(key)
                try:
                    v = float(row["value"])
                except ValueError:
                    continue
                if math.isfinite(v):
                    groups[key].append(v)
    print("exp\tgroup\tcond\tsnr\tmetric\tmean\tse\tn")
    for key in order:
        vals = groups[key]
        n = len(vals)
        if n == 0:
            continue
        m = sum(vals) / n
        se = math.sqrt(sum((x - m) ** 2 for x in vals) / (n - 1) / n) if n > 1 else float("nan")
        print("\t".join(key) + f"\t{m:.4f}\t{se:.4f}\t{n}")


if __name__ == "__main__":
    main(sys.argv[1:])
