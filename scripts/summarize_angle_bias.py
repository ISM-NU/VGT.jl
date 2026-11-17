"""Summary of scripts/validate_angle_bias.jl: harmonic fits of the angle bias per stencil, and paired differences
between stencils on the same cubes (realisation noise cancels).

Model per curve: b(t) = c0 + A4 sin 4t + B4 cos 4t + A8 sin 8t, t = true field angle (0, 5, ..., 85 deg).
GRF: fitted per seed, mean and standard error over seeds. CATS (one snapshot per run): ordinary least squares over
the 18 angles; the error is from the fit residuals (sigma_res^2 (X'X)^-1), which is honest when the residual is
noise and conservative when it is unmodelled shape. Paired rows fit the difference of two stencils' curves.

    python3 scripts/summarize_angle_bias.py <grf out>/grf_runs.csv [<grf out 2>/grf_runs.csv ...] <cats out>/cats_runs.csv
"""
import sys
from collections import defaultdict

import numpy as np

MODES = ("central_checked", "sobel_checked", "isotropic_checked")
PAIRS = (("central_checked", "isotropic_checked"), ("sobel_checked", "isotropic_checked"))


GRF_FIELDS = ("tag", "beta", "zeta", "r", "dn", "mode", "seed", "angle_deg", "bias_deg", "am", "rbar2", "npix")
CATS_FIELDS = ("run", "snr", "mode", "angle_deg", "bias_deg", "bias_err", "am", "nblocks")


def read_rows(path, fields):
    """Rows as dicts; NUL bytes (a write that failed on a full disk leaves a hole) and malformed or header lines are
    skipped, so a damaged file loses rows instead of breaking the summary (incomplete seeds are then dropped)."""
    out = []
    for line in open(path, errors="replace").read().replace("\x00", "").splitlines():
        parts = line.strip().split(",")
        if len(parts) != len(fields) or parts[0] in ("tag", "run"):
            continue
        try:
            float(parts[-1]); float(parts[-3])
        except ValueError:
            continue
        out.append(dict(zip(fields, parts)))
    return out


def design(t_deg):
    t = np.radians(np.asarray(t_deg, float))
    return np.stack([np.ones_like(t), np.sin(4 * t), np.cos(4 * t), np.sin(8 * t)], 1)


def ols(t, b):
    X = design(t)
    coef, *_ = np.linalg.lstsq(X, b, rcond=None)
    res = b - X @ coef
    dof = max(len(b) - X.shape[1], 1)
    cov = (res @ res / dof) * np.linalg.inv(X.T @ X)
    return coef, np.sqrt(np.diag(cov))


def fmt(c, e):
    return "  ".join(f"{n} {c[i]:+.3f}±{e[i]:.3f}" for i, n in enumerate(("c0", "A4", "B4", "A8")))


def grf(paths):
    """Several CSVs may be given (seeds split over jobs); a seed's rows are keyed by seed number, so a seed present in
    two files counts once, and a seed with an incomplete angle scan (a job stopped mid-seed) is dropped."""
    d = defaultdict(dict)   # (tag, r, dn) -> mode -> seed -> {angle: bias}
    am = defaultdict(list)
    rb = defaultdict(list)
    rows = [row for path in paths for row in read_rows(path, GRF_FIELDS)]
    for row in rows:
        key = (row["tag"], row["r"], row["dn"])
        d[key].setdefault(row["mode"], defaultdict(dict))[int(row["seed"])][int(row["angle_deg"])] = float(row["bias_deg"])
        am[(key, row["mode"])].append(float(row["am"]))
        rb[(key, row["mode"])].append(float(row["rbar2"]))
    print("GRF (deg; mean ± s.e. over seeds)")
    for key in d:
        print(f"  {key[0]} r={key[1]} dn={key[2]}")
        curves = {}
        for mode in MODES:
            seeds = {s: a for s, a in d[key][mode].items() if len(a) == 18}
            H = []
            for s, ang in seeds.items():
                t = sorted(ang)
                H.append(ols(t, np.array([ang[a] for a in t]))[0])
                curves.setdefault(mode, {})[s] = (t, np.array([ang[a] for a in t]))
            H = np.array(H)
            print(f"    {mode:18s} {fmt(H.mean(0), H.std(0, ddof=1) / np.sqrt(len(H)))}  AM {np.mean(am[(key, mode)]):.3f}"
                  f"  rbar2 {np.mean(rb[(key, mode)]):.3f}  seeds {len(H)}")
        for a, b in PAIRS:
            H = np.array([ols(curves[a][s][0], curves[a][s][1] - curves[b][s][1])[0] for s in curves[a] if s in curves[b]])
            if not len(H):
                continue
            print(f"    {a[:7]} - {b[:9]:9s}     {fmt(H.mean(0), H.std(0, ddof=1) / np.sqrt(len(H)))}")


def cats(path):
    d = defaultdict(dict)
    am = defaultdict(list)
    for row in read_rows(path, CATS_FIELDS):
        key = (row["run"], row["snr"])
        d[key].setdefault(row["mode"], {})[int(row["angle_deg"])] = float(row["bias_deg"])
        am[(key, row["mode"])].append(float(row["am"]))
    print("CATS (deg; OLS over angles, residual-based errors; c0 is the snapshot's own offset, not grid bias)")
    for key in d:
        print(f"  {key[0]} S/N {key[1]}")
        for mode in MODES:
            t = sorted(d[key][mode])
            c, e = ols(t, np.array([d[key][mode][a] for a in t]))
            print(f"    {mode:18s} {fmt(c, e)}  AM {np.mean(am[(key, mode)]):.3f}")
        for a, b in PAIRS:
            t = sorted(set(d[key][a]) & set(d[key][b]))
            c, e = ols(t, np.array([d[key][a][x] - d[key][b][x] for x in t]))
            print(f"    {a[:7]} - {b[:9]:9s}     {fmt(c, e)}")


if __name__ == "__main__":
    g = [p for p in sys.argv[1:] if "grf" in p]
    g and grf(g)
    for p in sys.argv[1:]:
        if "grf" not in p:
            cats(p)
