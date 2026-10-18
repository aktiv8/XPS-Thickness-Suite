"""
Regenerate the plots used in the documentation, straight from the code (so they cannot disagree with the app).

    python docs/make_figures.py

Writes docs/images/topofactor-curves.png. Needs matplotlib only (no display).
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import topofactor as TF  # noqa: E402

IMAGES = os.path.join(HERE, "images")
COLOURS = {"cylinder": "#1f77b4", "sphere": "#d62728", "sphere_shadowed": "#2ca02c"}
NAMES = {"cylinder": "Cylinder (fibre)", "sphere": "Sphere, X-rays reach all of it",
         "sphere_shadowed": "Sphere, X-ray shadowed"}


def topofactor_curves(path):
    t = np.linspace(0.005, 3.2, 700)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4), dpi=130)

    # (a) the three shapes at B = 1, exactly as the app applies them (thin limit below 0.1, Eq. 6 above)
    for key in ("cylinder", "sphere", "sphere_shadowed"):
        y = [TF.topofactor(key, v, 1.0, "eq6").tf if v <= 3.0 else np.nan for v in t]
        ax1.plot(t, y, color=COLOURS[key], lw=2, label=NAMES[key])
        ax1.axhline(TF.SHAPES[key].constant, color=COLOURS[key], lw=1, ls=":")
    ax1.axvspan(0, TF.T_THIN, color="0.9", zorder=0)
    ax1.text(TF.T_THIN / 2, 0.97, "thin-film\nlimit", ha="center", va="top", fontsize=8, color="0.35")
    ax1.set_xlim(0, 3.2)
    ax1.set_ylim(0.4, 1.0)
    ax1.set_xlabel("planar (Thickogram) thickness t / L")
    ax1.set_ylabel("Topofactor  T / t")
    ax1.set_title("Eq. (6) at B = 1  (dotted: the paper's constant averages)", fontsize=10)
    ax1.legend(loc="lower right", fontsize=8.5, frameon=False)
    ax1.grid(alpha=0.3)

    # (b) sensitivity of the cylinder factor to B
    for b, ls in ((0.7, "--"), (1.0, "-"), (1.4, "-.")):
        y = [TF.topofactor("cylinder", v, b, "eq6").tf if v <= 3.0 else np.nan for v in t]
        ax2.plot(t, y, color=COLOURS["cylinder"], lw=2, ls=ls, label=f"B = {b:g}")
    ax2.axvspan(0, TF.T_THIN, color="0.9", zorder=0)
    ax2.set_xlim(0, 3.2)
    ax2.set_ylim(0.4, 1.0)
    ax2.set_xlabel("planar (Thickogram) thickness t / L")
    ax2.set_title("Cylinder: effect of B = (Eo/Es)^0.75", fontsize=10)
    ax2.legend(loc="lower right", fontsize=8.5, frameon=False)
    ax2.grid(alpha=0.3)

    fig.suptitle("Topofactors for fibres and particles (after Shard, Wang and Spencer, 2009)", fontsize=11.5,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(IMAGES, exist_ok=True)
    out = os.path.join(IMAGES, "topofactor-curves.png")
    topofactor_curves(out)
    print("wrote", out)
