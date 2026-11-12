"""figAA: J2 maps of the LMC cube (64/32/16 px blocks) next to two CATS simulations.

    python figAA_lmc_j2.py <magnetization dir> <figures dir> <out.png>

LMC: lmc_data_b{64,32,16}.fits (HDU J2; recipe of validation_output/magnetization/run_lmc_j2.sh:
smooth 1.3 px, 3 sigma threshold, 20 px = 20 arcsec per pixel). Isotropic-null levels are the
medians printed in lmc_j2.out. Simulations (figures dir, from validation_output/figures/job1.sh):
transonic CATS runs b1p.32 (strong) and b.1p.32 (weak) at the LMC beam (3 px) and peak S/N 23,
same recipe, 32 px blocks, null from 20 surrogates (sim_j2_*.out in fig_prep.out).
Each group shares one colour scale (LMC panels: J2 0-0.4; simulations: 0-0.8).
"""
import re
import sys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

INK, INK2 = "#14202e", "#4b5867"
PX_DEG = 20 / 3600.0


def nulls(path):
    """{(source substring, block): null median} from the gradient_magnetization.jl log."""
    out, cur = {}, None
    for line in open(path):
        m = re.search(r"(\S+) centroid smooth [\d.]+ (?:highpass [\d.]+ )?block (\d+):.*median J2 ([\d.]+)", line)
        if m:
            cur = (m.group(1), int(m.group(2)), float(m.group(3)))
        m = re.search(r"isotropic null.*median J2 ([\d.]+)", line)
        if m and cur:
            out[(cur[0], cur[1])] = (cur[2], float(m.group(1)))
    return out


def j2map(path):
    with fits.open(path) as h:
        return np.asarray(h["J2"].data, dtype=float)


def main():
    mag, figs, outp = sys.argv[1:4]
    nl = nulls(os.path.join(mag, "lmc_j2.out"))
    ns = nulls(os.path.join(figs, "fig_prep.out"))
    fig = plt.figure(figsize=(12.5, 6.2), dpi=200)
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.5, wspace=0.3, left=0.06, right=0.93, top=0.86, bottom=0.07)
    # LMC row
    key = lambda d, b: [k for k in d if k[1] == b and ("crop_y1_x113.fits" in k[0] and "layers" not in k[0] or "sim" in k[0] or "calib3" in k[0])]
    top, bot = [], []
    for c, b in enumerate((64, 32, 16)):
        ax = fig.add_subplot(gs[0, c]); top.append(ax)
        j2 = j2map(os.path.join(mag, f"lmc_data_b{b}.fits"))
        k = [k for k in nl if k[1] == b and "rohsa/lmc2048_n16_smooth10" in k[0]][0]
        med, nul = nl[k]
        ext = (0, j2.shape[1] * b * PX_DEG, 0, j2.shape[0] * b * PX_DEG)
        im = ax.imshow(j2, origin="lower", cmap="viridis", vmin=0, vmax=0.4, extent=ext)
        ax.set_title(f"LMC, {b} px blocks ({b * 20 / 60:.0f}′)\nmedian J2 {med:.3f}, null {nul:.3f}  (ratio {med / nul:.1f})", fontsize=10.5, color=INK)
        ax.set_xlabel("x (degrees)", fontsize=9); ax.set_ylabel("y (degrees)", fontsize=9); ax.tick_params(labelsize=8)
    cb = fig.colorbar(im, ax=top, fraction=0.02, pad=0.01); cb.set_label("J2", fontsize=9); cb.ax.tick_params(labelsize=8)
    # simulations row
    sims = [("b1p.32", "strong field (transonic)"), ("b.1p.32", "weak field (transonic)")]
    for c, (run, lab) in enumerate(sims):
        ax = fig.add_subplot(gs[1, c]); bot.append(ax)
        j2 = j2map(os.path.join(figs, f"sim_j2_{run}_b32.fits"))
        k = [k for k in ns if k[1] == 32 and f"/{run}/" in k[0]][0]
        med, nul = ns[k]
        im2 = ax.imshow(j2, origin="lower", cmap="viridis", vmin=0, vmax=0.8, extent=(0, 256, 0, 256))
        ax.set_title(f"CATS {lab}, 32 px blocks\nmedian J2 {med:.2f}, null {nul:.2f}  (ratio {med / nul:.1f})", fontsize=10.5, color=INK)
        ax.set_xlabel("x (pixels)", fontsize=9); ax.set_ylabel("y (pixels)", fontsize=9); ax.tick_params(labelsize=8)
    cb2 = fig.colorbar(im2, ax=bot, fraction=0.03, pad=0.02); cb2.set_label("J2", fontsize=9); cb2.ax.tick_params(labelsize=8)
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off"); ax.set_position([0.70, 0.07, 0.25, 0.35])
    ax.text(0, 0.95, "J2 = 1: gradients all along one direction\nJ2 = 0: no preferred direction\n\n"
            "Simulations: LMC beam (3 px),\nS/N 23, same recipe, own colour\nscale, 8 x 8 blocks of 32 px.\n\nThe LMC sits only a few times above its\nisotropic null; the strong-field run is\n"
            "about four times above its own.", va="top", fontsize=10, color=INK2, transform=ax.transAxes)
    fig.suptitle("Anisotropy of velocity-centroid gradients (J2) per block: LMC H I against simulations", fontsize=12, color=INK)
    fig.savefig(outp, facecolor="white")
    for b in (64, 32, 16):
        k = [k for k in nl if k[1] == b and "rohsa/lmc2048_n16_smooth10" in k[0]][0]
        print("LMC", b, nl[k], "ratio", round(nl[k][0] / nl[k][1], 1))
    for k, v in ns.items():
        print("sim", k, v, "ratio", round(v[0] / v[1], 1))


main()
