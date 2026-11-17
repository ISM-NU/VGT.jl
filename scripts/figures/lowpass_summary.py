"""Summary of the lp variant of scripts/validate_faint_gas.jl (beam-scaled low-pass before the gradients): paired AM
difference lp - plain on CATS with the 4 sims as units, and the rule declared in the Julia script's header.

Errors: paired bootstrap over sims (lines of sight 2 and 3 pooled within a sim), then 32 px blocks within each drawn
sim; seeds are averaged per block first. Reuses boot() of faint_gas_summary.py with the sim as the unit.

    python3 scripts/figures/lowpass_summary.py <out dir of validate_faint_gas.jl --variants plain,lp>
Writes lowpass_summary.txt into that folder.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from faint_gas_summary import boot  # noqa: E402


def main(out):
    b = pd.read_csv(os.path.join(out, "blocks.csv"))
    b = b[b["npix"] > 0]
    b["unit"] = b["run"]
    b["blk"] = "l" + b["los"].astype(str) + "_" + b["block"].astype(str)
    b["scale"] = b["scale"].map({15: 10, 180: 120})
    agg = b.groupby(["unit", "snr", "layer", "scale", "blk", "variant"]).agg(am=("am", "mean"), w=("npix", "sum")).reset_index()
    piv = agg.pivot_table(index=["unit", "snr", "layer", "scale", "blk"], columns="variant", values=["am", "w"]).reset_index()
    lines = []
    p = lines.append
    p("lp - plain, paired dAM +- bootstrap std (sims, then blocks); sims with dAM > 0")
    res = {}
    for layer in ("wing", "all"):
        for sc in (10, 120):
            for snr in sorted(piv["snr"].unique()):
                g = piv[(piv["layer"] == layer) & (piv["scale"] == sc) & (piv["snr"] == snr)]
                g = g[np.isfinite(g[("am", "lp")]) & np.isfinite(g[("am", "plain")])]
                if g.empty:
                    continue
                d = pd.DataFrame({"cube": g["unit"], "d": g[("am", "lp")] - g[("am", "plain")], "w": g[("w", "plain")]})
                m, e = boot(d)
                pos = sum((dd["d"] * dd["w"]).sum() > 0 for _, dd in d.groupby("cube"))
                res[(layer, sc, snr)] = (m, e)
                plain = np.average(g[("am", "plain")], weights=g[("w", "plain")])
                p(f"  {layer:4s} {sc:3d}' S/N {snr:5.1f}: plain AM {plain:.3f}, dAM {m:+.4f} +- {e:.4f} "
                  f"({pos}/{d['cube'].nunique()}, {len(g)} blocks)")
            g = piv[(piv["layer"] == layer) & (piv["scale"] == sc)]
            g = g[np.isfinite(g[("am", "lp")]) & np.isfinite(g[("am", "plain")])]
            if not g.empty:
                d = pd.DataFrame({"cube": g["unit"] + "_" + g["snr"].astype(str),
                                  "d": g[("am", "lp")] - g[("am", "plain")], "w": g[("w", "plain")]})
                m, e = boot(d)
                res[(layer, sc, "pool")] = (m, e)
                p(f"  {layer:4s} {sc:3d}' pooled over S/N (sims x S/N as units): dAM {m:+.4f} +- {e:.4f}")
    ok = True
    for layer in ("wing", "all"):
        for snr in (23.0, 46.0):
            m, e = res.get((layer, 10, snr), (np.nan, np.nan))
            ok &= np.isfinite(m) and (m >= -0.005 or abs(m) <= 2 * e)
        m, e = res.get((layer, 10, "pool"), (np.nan, np.nan))
        ok &= np.isfinite(m) and m >= -2 * e
    p(f"Rule (declared in validate_faint_gas.jl): lp is {'FREE' if ok else 'COSTLY'} on clean data")
    w = pd.read_csv(os.path.join(out, "windows.csv"))
    bins = [0, 0.5, 0.8, 1.2, 1.7, 2.5, 4, 100]
    w["rbin"] = pd.cut(w["rho_g"], bins)
    t = w.groupby(["rbin", "variant"], observed=True)["am10_window"].mean().unstack()
    p("Per-window AM at 10' by rho_g bin (descriptive):")
    p(t.to_string(float_format=lambda x: f"{x:.3f}"))
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(out, "lowpass_summary.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
