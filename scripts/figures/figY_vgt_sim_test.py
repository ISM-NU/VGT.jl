"""figY: VGT on CATS MHD simulations against the true plane-of-sky field.

    python figY_vgt_sim_test.py <vda_sim dir> <out.png> [--degraded snr23]

Rows: strong field (b1p.032, M_A 0.7, M_s 4) and weak field (b.1p.032, M_A 7). Columns: integrated
intensity with the true field as segments | VGT field | angle error. Segments are drawn every
STEP px, each the axial mean over the STEP x STEP block (true field: from psi_true; VGT: from the
VGT map), so both are sampled identically. Error statistics are over all pixels, as in
scripts/score_vgt_vs_true.jl (AM = <cos 2 dpsi>, median |dpsi|).
Default input: the clean cubes (vda_sim/<run>/ppv.fits, vgt_none). --degraded snr23 uses
vda_sim/degraded/<run>/snr23 (3 px beam, white noise at peak S/N 23, dn 20, VGT noise threshold
= the noise), the closest to the LMC data.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

COBALT, AMBER, INK, INK2 = "#1f5fae", "#e3a21a", "#14202e", "#4b5867"
RUNS = [("b1p.032", "Strong field (Alfvén Mach 0.7, sonic Mach 4)"),
        ("b.1p.032", "Weak field (Alfvén Mach 7, sonic Mach 4)")]
STEP, LEN = 16, 13


def rd(p):
    return np.asarray(fits.getdata(p), dtype=float)


def headless(d):
    return (d + np.pi / 2) % np.pi - np.pi / 2


def block_axial(psi, n):
    """Axial mean angle of every n x n block (circular mean of 2 psi); NaN if no finite pixel."""
    ny, nx = (psi.shape[0] // n) * n, (psi.shape[1] // n) * n
    a = psi[:ny, :nx].reshape(ny // n, n, nx // n, n)
    c = np.nansum(np.cos(2 * a), axis=(1, 3))
    s = np.nansum(np.sin(2 * a), axis=(1, 3))
    ok = np.isfinite(a).any(axis=(1, 3))
    return np.where(ok, 0.5 * np.arctan2(s, c), np.nan)


def segments(ax, ang, color, lw=1.5):
    ny, nx = ang.shape
    yy, xx = np.mgrid[:ny, :nx]
    x0, y0 = (xx + 0.5) * STEP, (yy + 0.5) * STEP
    dx, dy = 0.5 * LEN * np.cos(ang), 0.5 * LEN * np.sin(ang)
    ok = np.isfinite(ang)
    lines = np.stack([np.stack([x0 - dx, y0 - dy], -1), np.stack([x0 + dx, y0 + dy], -1)], 2)[ok]
    from matplotlib.collections import LineCollection
    ax.add_collection(LineCollection(lines, colors=color, linewidths=lw, capstyle="round", zorder=3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root"); ap.add_argument("out"); ap.add_argument("--degraded", default="")
    a = ap.parse_args()
    fig, axs = plt.subplots(2, 3, figsize=(10.8, 7.0), dpi=200, constrained_layout=True)
    for r, (run, label) in enumerate(RUNS):
        if a.degraded:
            d = os.path.join(a.root, "degraded", run, a.degraded)
            psi_v = rd(os.path.join(d, "vgt_none_dn20", "psi_cube_none.fits"))
        else:
            d = os.path.join(a.root, run)
            psi_v = rd(os.path.join(d, "vgt_none", "psi_cube_none.fits"))
        psi_t = rd(os.path.join(d, "psi_true.fits"))
        im = rd(os.path.join(d, "ppv.fits")).sum(0)          # integrated intensity (code units)
        err = np.degrees(np.abs(headless(psi_v - psi_t)))
        ok = np.isfinite(err)
        am = float(np.mean(np.cos(2 * np.radians(err[ok]))))   # cos is even, so |d| is enough
        med = float(np.median(err[ok]))
        print(f"{run} {'degraded ' + a.degraded if a.degraded else 'clean'}: AM {am:.3f}  median {med:.1f} deg  pixels {ok.sum()}", flush=True)
        n = im.shape[0]
        lo, hi = np.percentile(im, [1, 99.5])
        bt, bv = block_axial(psi_t, STEP), block_axial(psi_v, STEP)
        for c in range(3):
            ax = axs[r, c]
            ax.set_xlim(0, n); ax.set_ylim(0, n); ax.set_aspect("equal")
            if c < 2:
                m = ax.imshow(im, origin="lower", cmap="Greys", vmin=lo, vmax=hi, extent=(0, n, 0, n), zorder=0)
            if c == 0:
                segments(ax, bt, "#c0392b")
                ax.set_title("true field on integrated intensity", fontsize=11, color=INK)
            elif c == 1:
                segments(ax, bv, COBALT)
                ax.set_title("VGT field (same sampling)", fontsize=11, color=INK)
            else:
                e = ax.imshow(err, origin="lower", cmap="magma_r", vmin=0, vmax=45, extent=(0, n, 0, n))
                ax.set_title(f"angle error: median {med:.0f}°, mean cos 2Δ = {am:.2f}", fontsize=11, color=INK)
            ax.set_xticks([0, 128, 256]); ax.set_yticks([0, 128, 256])
            if c == 0:
                ax.set_ylabel(f"{label}\ny (pixels)", fontsize=10, color=INK)
            else:
                ax.set_yticklabels([])
            if r == 1:
                ax.set_xlabel("x (pixels)", fontsize=10)
            else:
                ax.set_xticklabels([])
            ax.tick_params(labelsize=9)
        cb = fig.colorbar(m, ax=axs[r, :2], fraction=0.025, pad=0.01, location="right")
        cb.set_label("integrated intensity (code units)", fontsize=9); cb.ax.tick_params(labelsize=8)
        cb2 = fig.colorbar(e, ax=axs[r, 2], fraction=0.05, pad=0.02)
        cb2.set_label("|VGT − true| (degrees)", fontsize=9); cb2.ax.tick_params(labelsize=8)
    what = ("degraded to the LMC data: 3 px beam, noise at peak S/N 23, VGT window 20 px" if a.degraded
            else "noise-free cube, VGT window 20 px")
    fig.suptitle(f"VGT against the true field in MHD simulations ({what})", fontsize=11.5, color=INK2)
    fig.savefig(a.out, facecolor="white")


if __name__ == "__main__":
    main()
