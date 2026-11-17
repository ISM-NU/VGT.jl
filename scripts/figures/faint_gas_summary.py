"""Summary of scripts/validate_faint_gas.jl: paired AM differences of the faint-gas VGT variants against plain VGT on
CATS, the per-window AM against the gradient S/N rho_g, and an example map. Applies the pass rule declared in the
Julia script's header.

Errors: hierarchical paired bootstrap, cubes (run x line of sight) resampled, then 32 px blocks within each drawn cube;
seeds are averaged per block first (they share the MHD field). At 120' the blocks of one cube are strongly correlated,
so there the error is set by the cube resampling.

    python3 scripts/figures/faint_gas_summary.py <out dir of validate_faint_gas.jl> LO,HI
LO,HI: the block-median rho_g range of the target data's faint layer (computed as in validate_faint_gas.jl).
Writes summary.txt and faint_gas.png into that folder.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

VARS = ["amp2", "pre4", "amp2pre4"]


def boot(diff, nboot=2000, seed=5):
    """diff: DataFrame with cube, block, d, w. Returns mean and bootstrap std (cubes, then blocks within cubes)."""
    rng = np.random.default_rng(seed)
    diff = diff[np.isfinite(diff["d"]) & (diff["w"] > 0)]
    if diff.empty:
        return np.nan, np.nan
    groups = [(g["d"].to_numpy(), g["w"].to_numpy()) for _, g in diff.groupby("cube")]
    tot = sum((d * w).sum() for d, w in groups) / sum(w.sum() for _, w in groups)
    out = np.empty(nboot)
    for b in range(nboot):
        num = den = 0.0
        for i in rng.integers(0, len(groups), len(groups)):
            d, w = groups[i]
            j = rng.integers(0, d.size, d.size)
            num += (d[j] * w[j]).sum()
            den += w[j].sum()
        out[b] = num / den
    return tot, out.std()


def main(out):
    b = pd.read_csv(os.path.join(out, "blocks.csv"))
    b["cube"] = b["run"] + "_l" + b["los"].astype(str)
    b["scale"] = b["scale"].map({15: 10, 180: 120})   # FWHM in px -> the equivalent scale in arcmin on a 40" grid
    b = b[b["npix"] > 0]
    agg = b.groupby(["cube", "snr", "layer", "scale", "block", "variant"]).agg(am=("am", "mean"), w=("npix", "sum")).reset_index()
    piv = agg.pivot_table(index=["cube", "snr", "layer", "scale", "block"], columns="variant", values=["am", "w"]).reset_index()
    lines = []
    p = lines.append
    w = pd.read_csv(os.path.join(out, "windows.csv"))
    w["cube"] = w["run"] + "_l" + w["los"].astype(str)
    wp = w[w["variant"] == "plain"]
    p("rho_g (block median per window; median over windows and cubes, seed 1) by peak S/N:")
    matched = []
    for snr, g in wp[wp["seed"] == 1].groupby("snr"):
        rw = g[g["vcen_std"] < -1]["rho_g"].median()
        ra = g["rho_g"].median()
        p(f"  S/N {snr:5.1f}: wing windows {rw:.2f}, all windows {ra:.2f}")
        if TARGET[0] <= rw <= TARGET[1]:
            matched.append(snr)
    p(f"S/N levels whose wing rho_g lies in the target range {TARGET}: {matched}")
    p("")
    p("AM (mean over cubes, seeds, blocks) and paired dAM = variant - plain +- bootstrap std; cubes with dAM > 0:")
    res = {}
    for layer in ("wing", "all"):
        for scale in (10, 120):
            for snr in sorted(piv["snr"].unique()):
                g = piv[(piv["layer"] == layer) & (piv["scale"] == scale) & (piv["snr"] == snr)]
                if g.empty:
                    continue
                base = (g[("am", "plain")] * g[("w", "plain")]).sum() / g[("w", "plain")].sum()
                nb = int((g[("w", "plain")] > 0).sum())
                row = f"  {layer:4s} {scale:3d}' S/N {snr:5.1f}  plain AM {base:.3f} ({nb} blocks)"
                for v in VARS:
                    d = pd.DataFrame({"cube": g["cube"], "d": g[("am", v)] - g[("am", "plain")], "w": g[("w", "plain")]})
                    m, e = boot(d)
                    d = d[np.isfinite(d["d"]) & (d["w"] > 0)]
                    pos = sum((dd["d"] * dd["w"]).sum() > 0 for _, dd in d.groupby("cube"))
                    res[(layer, scale, snr, v)] = (m, e, pos, d["cube"].nunique())
                    row += f" | {v} {m:+.4f} +- {e:.4f} ({pos}/{d['cube'].nunique()})"
                p(row)
    p("")
    p("Pooled over S/N 8, 12, 16 (cubes x S/N as units):")
    for layer in ("wing", "all"):
        for scale in (10, 120):
            g = piv[(piv["layer"] == layer) & (piv["scale"] == scale) & piv["snr"].isin([8.0, 12.0, 16.0])]
            row = f"  {layer:4s} {scale:3d}'"
            for v in VARS:
                d = pd.DataFrame({"cube": g["cube"] + "_" + g["snr"].astype(str), "d": g[("am", v)] - g[("am", "plain")],
                                  "w": g[("w", "plain")]})
                m, e = boot(d)
                res[(layer, scale, "pool", v)] = (m, e)
                row += f" | {v} {m:+.4f} +- {e:.4f}"
            p(row)
    p("")
    p("Pass rule (declared in validate_faint_gas.jl):")
    for v in ("amp2", "pre4"):
        c1 = []
        nan4 = (np.nan, np.nan, 0, 0)
        for key in [("pool",)] + [(s,) for s in matched]:
            m10, e10 = res.get(("wing", 10, key[0], v), nan4)[:2]
            m120, e120 = res.get(("wing", 120, key[0], v), nan4)[:2]
            c1.append(m10 > 2 * e10 and (not np.isfinite(m120) or m120 > -e120))
        c2 = all(not np.isfinite(res.get((layer, sc, hi, v), nan4)[0]) or res[(layer, sc, hi, v)][0] >= -0.005
                 for layer in ("wing", "all") for sc in (10, 120) for hi in (23.0, 46.0))
        c3 = all(res.get(("wing", 10, s, v), nan4)[2] >= 6 for s in sorted(set([8.0, 12.0, 16.0] + matched)))
        p(f"  (missing cells, e.g. wing at 120' where a window has no coverage, count as untestable)")
        p(f"  {v}: (1) faint gain {all(c1)}  (2) no harm at S/N 23, 46 {c2}  (3) sign in >= 6/8 cubes {c3}  -> "
          f"{'PASS' if all(c1) and c2 and c3 else 'fail'}")
    p("")
    p("Per-window AM at 10' (window truth) by rho_g bin, mean over windows/cubes/seeds:")
    bins = [0, 0.5, 0.8, 1.2, 1.7, 2.5, 4, 100]
    w["rbin"] = pd.cut(w["rho_g"], bins)
    t = w.groupby(["rbin", "variant"], observed=True)["am10_window"].mean().unstack()
    n = w[w["variant"] == "plain"].groupby("rbin", observed=True).size()
    p(pd.concat([t, n.rename("n")], axis=1).to_string(float_format=lambda x: f"{x:.3f}"))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(out, "summary.txt"), "w").write(txt + "\n")
    figure(out, w, bins)


def figure(out, w, bins):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.io import fits
    maps = sorted(glob.glob(os.path.join(out, "maps", "*snr12_wing.fits")))
    fig = plt.figure(figsize=(14, 7.5))
    ax0 = fig.add_subplot(2, 3, 1)
    cen = 0.5 * (np.array(bins[:-1]) + np.array(bins[1:]))
    cen[-1] = 5
    for v, col in zip(["plain"] + VARS, ["0.2", "#1f77b4", "#d62728", "#9467bd"]):
        g = w[w["variant"] == v].groupby("rbin", observed=True)["am10_window"].mean()
        ax0.plot([cen[i] for i, iv in enumerate(pd.IntervalIndex.from_breaks(bins)) if iv in g.index], g.values,
                 "o-", color=col, label=v)
    ax0.axvspan(*TARGET, color="0.9", zorder=0, label="target faint layer")
    ax0.set_xlabel("gradient S/N rho_g (block median)")
    ax0.set_ylabel("window AM at 10' vs window truth")
    ax0.legend(fontsize=8)
    if maps:
        with fits.open(maps[0]) as f:
            col = f["COLUMN"].data
            tr = f["TRUE10"].data
            for i, name in enumerate(["PLAIN10", "AMP210", "PRE410"]):
                ax = fig.add_subplot(2, 3, 2 + i if i < 2 else 4)
                d = np.degrees(np.angle(np.exp(2j * (f[name].data - tr))) / 2)
                im = ax.imshow(d.T, origin="lower", cmap="twilight_shifted", vmin=-90, vmax=90)
                am = np.nanmean(np.cos(np.radians(2 * d)))
                ax.set_title(f"{name[:-2].lower()}: psi - true at 10', AM {am:.3f}", fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(im, ax=ax, shrink=0.8, label="deg")
            ax = fig.add_subplot(2, 3, 5)
            ax.imshow(col.T, origin="lower", cmap="gray")
            ax.set_title("wing-layer column (noisy)", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(os.path.basename(maps[0]).replace(".fits", "") + " (S/N 12, seed 1)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "faint_gas.png"), dpi=90)


if __name__ == "__main__":
    TARGET = tuple(float(x) for x in sys.argv[2].split(","))
    main(sys.argv[1])
