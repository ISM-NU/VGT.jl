"""Tables, the declared selection rules and figures from scripts/validate_resolution_uncertainty.jl outputs.

Merges every <root>/*/{tf_runs,tf_fit,sweep_runs,mc_runs}.csv (one directory per part and S/N), applies rules R1,
R2 and the adoption rule of docs/validation-resolution-uncertainty.md, and writes <root>/summary.txt,
<root>/fig_transfer.png, <root>/fig_accuracy.png, <root>/fig_calibration.png.

    python3 scripts/summarize_resolution_uncertainty.py --root validation_output/resolution
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PROD = (20, 4, 0)


def load(root, name):
    fs = sorted(glob.glob(os.path.join(root, "*", name)))
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def tie_tol(se):
    return max(0.005, 2 * se)


def rules(sw, out):
    fcols = [c for c in sw.columns if c.startswith("am_F")]
    g = sw.groupby(["field", "snr", "dn", "ker", "rb"])
    agg = g.agg(fwhm=("fwhm_px", "first"), n=("seed", "count"),
                **{c: (c, "mean") for c in ["am_native"] + fcols},
                **{c + "_se": (c, lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0)
                   for c in ["am_native"] + fcols}).reset_index()
    lines = []
    for (field, snr), d in agg.groupby(["field", "snr"]):
        d = d.copy()
        best = d.loc[d.am_native.idxmax()]
        tol = tie_tol(best.am_native_se)
        tied = d[d.am_native >= best.am_native - tol]
        pick = tied.loc[tied.fwhm.idxmin()]
        prod = d[(d.dn == PROD[0]) & (d.ker == PROD[1]) & (d.rb == PROD[2])]
        pa = prod.am_native.iloc[0] if len(prod) else np.nan
        lines.append(f"R1 {field} r={snr}: best dn{int(best.dn)} ker{int(best.ker)} rb{int(best.rb)} "
                     f"AM {best.am_native:.4f} (FWHM {best.fwhm:.1f}); chosen dn{int(pick.dn)} ker{int(pick.ker)} "
                     f"rb{int(pick.rb)} AM {pick.am_native:.4f} FWHM {pick.fwhm:.1f}; production AM {pa:.4f}")
        for c in fcols:
            F = float(c[4:])
            e = d[d.fwhm <= F + 1e-9]
            if e.empty:
                lines.append(f"R2 {field} r={snr} F={F}: no eligible config")
                continue
            b = e.loc[e[c].idxmax()]
            tied = e[e[c] >= b[c] - tie_tol(b[c + "_se"])]
            if ((tied.dn == PROD[0]) & (tied.ker == PROD[1]) & (tied.rb == PROD[2])).any():
                p = tied[(tied.dn == PROD[0]) & (tied.ker == PROD[1]) & (tied.rb == PROD[2])].iloc[0]
            else:
                p = tied.sort_values(["dn", "fwhm"]).iloc[0]
            pv = prod[c].iloc[0] if len(prod) and prod.fwhm.iloc[0] <= F + 1e-9 else np.nan
            lines.append(f"R2 {field} r={snr} F={F}: best dn{int(b.dn)} ker{int(b.ker)} rb{int(b.rb)} {b[c]:.4f}; "
                         f"chosen dn{int(p.dn)} ker{int(p.ker)} rb{int(p.rb)} {p[c]:.4f}; production {pv:.4f}; "
                         f"chosen - production {p[c] - pv:+.4f}")
    agg.to_csv(os.path.join(out, "sweep_summary.csv"), index=False)
    return lines, agg


def fig_transfer(tr, fit, out):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    d = tr[(tr.ker == 4) & (tr.rb == 0) & (tr.snr == np.inf)]
    for dn, g in d.groupby("dn"):
        m = g.groupby("k").T.mean()
        ax[0].plot(m.index, m.values, "o-", ms=3, label=f"dn {dn}")
    ax[0].axhline(0.5, color="0.6", lw=0.8)
    ax[0].set(xscale="log", xlabel="k of the orientation modulation (rad/px)", ylabel="recovered / true amplitude T",
              title="Transfer function, ker 4, no noise")
    ax[0].legend(fontsize=7)
    for ker, mk in ((1, "o"), (2, "s"), (4, "^")):
        f = fit[(fit.ker == ker) & (fit.rb == 0) & (fit.snr == np.inf)].sort_values("dn")
        ax[1].plot(f.dn, f.fwhm_fit_px, mk + "-", label=f"measured ker {ker}")
        ax[1].plot(f.dn, f.fwhm_formula_px, mk + ":", color=ax[1].lines[-1].get_color(), label=f"formula ker {ker}")
    ax[1].set(xlabel="dn (px)", ylabel="FWHM (px)", title="Effective FWHM, no noise")
    ax[1].legend(fontsize=7)
    for snr, g in fit[(fit.ker == 4) & (fit.rb == 0)].groupby("snr"):
        g = g.sort_values("dn")
        ax[2].plot(g.dn, g.fwhm_fit_px, "o-", label=f"r = {snr}")
    ax[2].set(xlabel="dn (px)", ylabel="FWHM (px)", title="Effective FWHM against S/N, ker 4")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_transfer.png"), dpi=130)


def fig_accuracy(agg, out):
    fields = sorted(agg.field.unique())
    snrs = sorted(agg.snr.unique(), reverse=True)
    fig, ax = plt.subplots(len(fields), len(snrs), figsize=(3.6 * len(snrs), 3.2 * len(fields)), squeeze=False)
    for i, f in enumerate(fields):
        for j, s in enumerate(snrs):
            d = agg[(agg.field == f) & (agg.snr == s)]
            a = ax[i, j]
            for ker, mk in ((1, "o"), (2, "s"), (4, "^")):
                e = d[(d.ker == ker) & (d.rb == 0)].sort_values("dn")
                a.errorbar(e.fwhm, e.am_native, yerr=e.am_native_se, fmt=mk + "-", ms=4, label=f"ker {ker}")
                for _, r in e.iterrows():
                    a.annotate(str(int(r.dn)), (r.fwhm, r.am_native), fontsize=6)
            e = d[d.rb == 1]
            a.plot(e.fwhm, e.am_native, "kx", label="rbar2 weight")
            p = d[(d.dn == 20) & (d.ker == 4) & (d.rb == 0)]
            a.plot(p.fwhm, p.am_native, "r*", ms=10, label="production")
            a.set(title=f"{f}, r = {s}", xlabel="effective FWHM (px)", ylabel="AM vs pixel truth")
            if i == 0 and j == 0:
                a.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_accuracy.png"), dpi=130)


def fig_calibration(mc, out):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for (field, snr, dn, scale), g in mc.groupby(["field", "snr", "dn", "scale"]):
        if dn != 20:
            continue
        bs = np.mean([[float(x) for x in s.split(";")] for s in g.bin_sigma_deg], axis=0)
        be = np.mean([[float(x) for x in s.split(";")] for s in g.bin_rms_err_deg], axis=0)
        a = ax[0] if scale == 1.0 else ax[1]
        a.plot(bs, be, "o-", label=f"{field} r={snr}")
    for a, s in zip(ax, (1.0, 0.5)):
        lim = a.get_xlim()[1]
        a.plot([0, lim], [0, lim], "k:", lw=0.8)
        a.set(xlabel="predicted sigma_psi (deg, bin median)", ylabel="rms actual error (deg)",
              title=f"Calibration, dn 20, scale {s}")
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_calibration.png"), dpi=130)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    lines = []
    tr, fit = load(a.root, "tf_runs.csv"), load(a.root, "tf_fit.csv")
    if len(fit):
        lines.append("Effective FWHM (px), rb 0:")
        lines.append(fit[fit.rb == 0].sort_values(["snr", "ker", "dn"]).to_string(index=False))
        fig_transfer(tr, fit, a.root)
    sw = load(a.root, "sweep_runs.csv")
    if len(sw):
        rl, agg = rules(sw, a.root)
        lines += rl
        fig_accuracy(agg, a.root)
    mc = load(a.root, "mc_runs.csv")
    if len(mc):
        g = mc.groupby(["field", "snr", "dn", "scale"]).agg(
            cov1=("cov1", "mean"), cov1_lo=("cov1_lo", "min"), cov1_hi=("cov1_hi", "max"), cov2=("cov2", "mean"),
            rms_z=("rms_z", "mean"), sigma=("median_sigma_deg", "mean"), rms_err=("rms_err_deg", "mean"),
            mean_err=("mean_err_deg", "mean"), offset=("median_offset_deg", "mean"), n=("seed", "count"))
        lines.append("Uncertainty coverage (target 0.68 / 0.95, rms z 1):")
        lines.append(g.reset_index().to_string(index=False))
        fig_calibration(mc, a.root)
    txt = "\n".join(lines)
    open(os.path.join(a.root, "summary.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
