#!/usr/bin/env python3
"""
figure1_composite.py
Composite Figure 1:
  (a) Harvard Forest 6-panel flux profiles (top, full width)
  (b) Yale Myers Forest single-panel profile (bottom-left, half width)
  (c) IMG_5926_edited.jpg – DinoLift photo (bottom-right top quarter)
  (d) IMG_6437.jpg – arborist climbing photo (bottom-right bottom quarter)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
from statsmodels.nonparametric.smoothers_lowess import lowess
import re, os, warnings
warnings.filterwarnings("ignore")

base_dir = "/sessions/tender-zen-fermi/mnt/whole_tree_flux"

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.linewidth": 0.4,
    "xtick.major.width": 0.4,
    "ytick.major.width": 0.4,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
})

component_colors = {"stem": "#8B4513", "branch": "#4682B4", "leaf": "#2E8B57"}
component_labels = {"stem": "Stem", "branch": "Branch", "leaf": "Leaf"}

# ── Data: Harvard Forest ─────────────────────────────────────────────────────
dat = pd.read_csv(os.path.join(base_dir, "data processing", "goFlux_reprocessing",
                                "results", "canopy_flux_goFlux_compiled.csv"))

# Remove Nyssa sylvatica tag 3
dat = dat[~((dat["Species"] == "bg") & (dat["Tree_Tag"] == 3))].copy()
dat["Component"] = dat["Type"].replace({"leaf (shaded)": "leaf"})
dat = dat.dropna(subset=["CH4_best.flux", "Height_m"])

species_lookup = {"bg": "Nyssa sylvatica", "rm": "Acer rubrum",
                  "ro": "Quercus rubra", "hem": "Tsuga canadensis"}
dat["Species_full"] = dat["Species"].map(species_lookup)
dat["Tree_label"] = dat["Species_full"] + " (" + dat["Site"] + ")"

# Order: group by site then alphabetically
order = (dat[["Tree_label", "Site", "Species_full"]]
         .drop_duplicates()
         .sort_values(["Site", "Species_full"])["Tree_label"].tolist())

# ── Data: Yale Myers Forest ──────────────────────────────────────────────────
ymf = pd.read_csv(os.path.join(base_dir, "data processing", "goFlux_reprocessing",
                                "ymf_black_oak", "results",
                                "ymf_black_oak_flux_compiled.csv"))
ymf["Height_num"] = ymf["Height_m"].apply(
    lambda x: float(re.match(r"^[\d.]+", str(x)).group()) if re.match(r"^[\d.]+", str(x)) else np.nan
)
ymf = ymf.dropna(subset=["Height_num", "CH4_best.flux"])

# ── Photos ───────────────────────────────────────────────────────────────────
img_c = Image.open(os.path.join(base_dir, "IMG_5926_edited.jpg"))
img_d = Image.open(os.path.join(base_dir, "IMG_6437.jpg"))

# ── Figure layout ────────────────────────────────────────────────────────────
# Using nested GridSpec:
#   Top 2/3: 2x3 grid for Harvard Forest panels (a)
#   Bottom 1/3: left half = Yale (b), right half split into top (c) + bottom (d)

fig = plt.figure(figsize=(7.5, 10))

# Outer grid: top (Harvard) and bottom (Yale + photos)
outer = gridspec.GridSpec(2, 1, height_ratios=[1.8, 1.2], hspace=0.32, figure=fig)

# Top: Harvard Forest 2x3
gs_hf = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=outer[0],
                                          hspace=0.45, wspace=0.45)

# Bottom: left = Yale, right = photos stacked
gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1],
                                           wspace=0.12)
# Bottom-right: two photos stacked
gs_photos = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_bot[0, 1],
                                              hspace=0.08)

# ── Helper: LOESS + confidence band ─────────────────────────────────────────
def plot_loess(ax, x, y, color, frac=0.6):
    """Plot LOESS smooth with shaded CI band."""
    if len(x) < 4:
        return
    sorted_idx = np.argsort(x)
    xs, ys = x[sorted_idx], y[sorted_idx]
    result = lowess(ys, xs, frac=frac, return_sorted=True)
    xl, yl = result[:, 0], result[:, 1]
    # Bootstrap-like residual band
    residuals = np.interp(xs, xl, yl) - ys
    se = np.std(residuals)
    ax.fill_betweenx(xl, yl - 1.5 * se, yl + 1.5 * se,
                      color=color, alpha=0.12, linewidth=0)
    ax.plot(yl, xl, color=color, linewidth=0.8)

# ── Panel (a): Harvard Forest subplots ───────────────────────────────────────
hf_axes = []
for i, label in enumerate(order):
    row, col = divmod(i, 3)
    ax = fig.add_subplot(gs_hf[row, col])
    hf_axes.append(ax)

    sub = dat[dat["Tree_label"] == label]

    # Zero line
    ax.axvline(0, color="black", linestyle="--", linewidth=0.4)

    # Plot each component
    for comp in ["stem", "branch", "leaf"]:
        csub = sub[sub["Component"] == comp]
        if len(csub) == 0:
            continue
        ax.scatter(csub["CH4_best.flux"].values, csub["Height_m"].values,
                   c=component_colors[comp], s=18, alpha=0.8,
                   label=component_labels[comp], zorder=3, edgecolors="none")

    # LOESS on stem only
    stem = sub[sub["Component"] == "stem"]
    if len(stem) >= 4:
        plot_loess(ax, stem["Height_m"].values, stem["CH4_best.flux"].values,
                   component_colors["stem"])

    # Title: italic species
    ax.set_title(label, fontsize=8, fontstyle="italic")
    ax.set_ylabel("Height (m)" if col == 0 else "", fontsize=8)
    ax.set_xlabel(r"CH$_4$ flux (nmol m$^{-2}$ s$^{-1}$)" if row == 1 else "",
                  fontsize=8)
    ax.tick_params(labelsize=7)

# Panel tag "a" on first subplot
hf_axes[0].text(-0.25, 1.08, "a", transform=hf_axes[0].transAxes,
                fontsize=14, fontweight="bold", va="top")

# Legend on bottom of Harvard block
handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=component_colors[c],
                       markersize=6, label=component_labels[c])
           for c in ["stem", "branch", "leaf"]]
# Place legend in bottom-right of last Harvard panel
hf_axes[-1].legend(handles=handles, loc="lower right", frameon=False, fontsize=8,
                    borderpad=0.3, handletextpad=0.3)

# ── Panel (b): Yale Myers Forest ─────────────────────────────────────────────
ax_ymf = fig.add_subplot(gs_bot[0, 0])
ax_ymf.axvline(0, color="black", linestyle="--", linewidth=0.4)
ax_ymf.scatter(ymf["CH4_best.flux"].values, ymf["Height_num"].values,
               c=component_colors["stem"], s=25, alpha=0.8, zorder=3,
               edgecolors="none")
plot_loess(ax_ymf, ymf["Height_num"].values, ymf["CH4_best.flux"].values,
           component_colors["stem"], frac=0.7)
ax_ymf.set_title("Quercus velutina\n(Yale Myers Forest)",
                  fontsize=9, fontstyle="italic")
ax_ymf.set_ylabel("Height (m)", fontsize=9)
ax_ymf.set_xlabel(r"CH$_4$ flux (nmol m$^{-2}$ s$^{-1}$)", fontsize=9)
ax_ymf.tick_params(labelsize=8)
ax_ymf.text(-0.18, 1.05, "b", transform=ax_ymf.transAxes,
            fontsize=14, fontweight="bold", va="top")

# ── Panel (c): DinoLift photo ────────────────────────────────────────────────
ax_c = fig.add_subplot(gs_photos[0, 0])
ax_c.imshow(img_c, aspect="auto")
ax_c.set_axis_off()
ax_c.text(0.02, 0.98, "c", transform=ax_c.transAxes,
          fontsize=14, fontweight="bold", va="top", color="white",
          bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.5, lw=0))

# ── Panel (d): Arborist photo ────────────────────────────────────────────────
ax_d = fig.add_subplot(gs_photos[1, 0])
ax_d.imshow(img_d, aspect="auto")
ax_d.set_axis_off()
ax_d.text(0.02, 0.98, "d", transform=ax_d.transAxes,
          fontsize=14, fontweight="bold", va="top", color="white",
          bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.5, lw=0))

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir = os.path.join(base_dir, "figures")
os.makedirs(out_dir, exist_ok=True)

fig.savefig(os.path.join(out_dir, "figure1_composite.pdf"),
            bbox_inches="tight", dpi=300)
fig.savefig(os.path.join(out_dir, "figure1_composite.png"),
            bbox_inches="tight", dpi=300)

print("Saved figure1_composite.pdf/.png")
