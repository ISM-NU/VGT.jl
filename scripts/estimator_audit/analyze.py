"""Summarise the estimator-audit grid.

Rule fixed before the results were read: a paired difference between two
estimators counts as real only if the per-seed median errors differ in the same
direction in a two-sided sign test with p < 0.05 (9/10 or 10/10 seeds);
otherwise it is "within noise".

    python3 scripts/estimator_audit/analyze.py validation_output/estimator_audit
"""
import csv
import glob
import math
import statistics as st
import sys
from collections import defaultdict

d = sys.argv[1]
rows = []
for f in glob.glob(f"{d}/*_errors.csv"):
    rows += list(csv.DictReader(open(f)))
cal = []
for f in glob.glob(f"{d}/*_calib.csv"):
    cal += list(csv.DictReader(open(f)))

ESTS = ["vm", "gauss", "ms2", "ms4", "vmu", "st1", "st2", "st2ms2", "hg100", "hg100n", "hg36"]


def fl(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except ValueError:
        return None


def sign_p(k, n):
    if n == 0:
        return 1.0
    m = min(k, n - k)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(m + 1)) / 2 ** n)


g = defaultdict(dict)  # (scen, ref, region) -> est -> seed -> row
for r in rows:
    g[(r["scenario"], r["reference"], r["region"])].setdefault(r["estimator"], {})[int(r["seed"])] = r


def verdict(a, b):
    """Paired per-seed medians of a minus b; negative = a better."""
    seeds = sorted(set(a) & set(b))
    dl = [fl(a[s]["med"]) - fl(b[s]["med"]) for s in seeds if fl(a[s]["med"]) is not None and fl(b[s]["med"]) is not None]
    if not dl:
        return "n/a"
    k = sum(x < 0 for x in dl)
    nz = sum(x != 0 for x in dl)
    p = sign_p(k, nz)
    tag = ("better" if k > nz - k else "worse") if p < 0.05 else "~"
    return f"{st.median(dl):+.2f}({k}/{nz}){tag}"


regions = sys.argv[2].split(",") if len(sys.argv) > 2 else ["clean", "window_clean", "all_signal", "smoothed_all_signal", "all", "smoothed_all"]
ests = sys.argv[3].split(",") if len(sys.argv) > 3 else ["vm", "gauss", "ms4", "vmu", "ms2", "st2"]
print("cells: median over seeds of per-seed median error [deg] / median p90; deltas: median paired delta (seeds lower/nonzero) verdict")
for key in sorted(g):
    scen, ref, region = key
    if region not in regions:
        continue
    e = g[key]
    nseeds = len(e.get("vm", {}))
    line = f"{scen:20s} {ref:10s} {region:20s} n={nseeds:2d} "
    for est in ests:
        if est not in e:
            continue
        meds = [fl(r["med"]) for r in e[est].values() if fl(r["med"]) is not None]
        p90 = [fl(r["p90"]) for r in e[est].values() if fl(r["p90"]) is not None]
        sup = [fl(r["support"]) for r in e[est].values() if fl(r["support"]) is not None]
        if meds:
            line += f"| {est} {st.median(meds):.2f}/{st.median(p90):.1f} s{st.median(sup):.2f} "
    print(line)
    for est in ests:
        if est in ("vm",) or est not in e:
            continue
        print(f"{'':54s} {est:6s} vs vm {verdict(e[est], e['vm']):24s}" +
              (f" vs gauss {verdict(e[est], e['gauss'])}" if "gauss" in e and est != "gauss" else ""))

print()
print("calibration: mean Spearman(error, confidence) over seeds (negative = higher confidence, lower error)")
c = defaultdict(list)
dec = defaultdict(list)
for r in cal:
    v = fl(r["spearman"])
    if v is not None:
        c[(r["suite"], r["map"], r["estimator"])].append(v)
    dec[(r["suite"], r["map"], r["estimator"])].append([fl(r[f"d{k}"]) for k in range(1, 11)])
for key in sorted(c):
    if key[2] not in ests:
        continue
    ds = dec[key]
    lo = st.median([x[0] for x in ds if x[0] is not None])
    hi = st.median([x[9] for x in ds if x[9] is not None])
    print(f"{key[0]:4s} {key[1]:30s} {key[2]:7s} rho {st.mean(c[key]):+.2f} (n={len(c[key])})  median err lowest-conf decile {lo:6.2f}  highest {hi:5.2f}")
