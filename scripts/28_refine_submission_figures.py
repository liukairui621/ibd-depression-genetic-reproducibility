#!/usr/bin/env python3
"""Rebuild the submission figures from frozen derived tables.

This script changes layout only. It does not recompute statistics or alter the
underlying estimates. Legends use dedicated exterior regions so they cannot
cover plotted data, and panel labels use a common coordinate rule.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from PIL import Image


REPO = Path(__file__).resolve().parents[1]
REPO_FIGURES = REPO / "figures"
MAIN_OUT = REPO_FIGURES
SUPP_OUT = REPO_FIGURES

FONT = "Arial"
TEXT = "#202124"
GREY = "#8A9299"
TEAL = "#007C83"
CORAL = "#C84B31"
BLUE = "#3F77A8"
GOLD = "#E8B94B"

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.labelweight": "bold",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "text.color": TEXT,
    "axes.labelcolor": TEXT,
    "axes.edgecolor": TEXT,
    "xtick.color": TEXT,
    "ytick.color": TEXT,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.8, length=3)
    ax.grid(False)


def panel_label(ax: plt.Axes, label: str, x: float = -0.20) -> None:
    ax.text(
        x, 1.08, label, transform=ax.transAxes, ha="left", va="bottom",
        fontsize=13, fontweight="bold", color=TEXT, clip_on=False,
    )


def copy_to_repository_figures(source: Path) -> None:
    destination = REPO_FIGURES / source.name
    if source.resolve() != destination.resolve():
        shutil.copy2(source, destination)


def save_main_figure(fig: plt.Figure, stem: str) -> None:
    MAIN_OUT.mkdir(parents=True, exist_ok=True)
    REPO_FIGURES.mkdir(parents=True, exist_ok=True)
    pdf = MAIN_OUT / f"Figure_{stem}.pdf"
    png = MAIN_OUT / f"Figure_{stem}.png"
    tif = MAIN_OUT / f"Fig{stem}.tif"
    fig.savefig(pdf, bbox_inches=None)
    fig.savefig(png, dpi=300, bbox_inches=None, facecolor="white")
    fig.savefig(
        tif, dpi=300, bbox_inches=None, facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    # Matplotlib may retain an opaque alpha channel in TIFF output. Flatten it
    # explicitly to standard RGB for journal image preflight.
    with Image.open(tif) as image:
        rgb = Image.new("RGB", image.size, "white")
        if image.mode == "RGBA":
            rgb.paste(image, mask=image.getchannel("A"))
        else:
            rgb.paste(image.convert("RGB"))
        rgb.save(tif, compression="tiff_lzw", dpi=(300, 300))
    for source in (pdf, png, tif):
        copy_to_repository_figures(source)


def build_figure_1() -> None:
    derived = REPO / "results" / "derived"
    rg = read_tsv(derived / "global_core" / "rg_summary.tsv")
    primary = rg[
        ((rg.trait1 == "MDD_Howard2019") & (rg.trait2 == "IBD_FinnGen_R12"))
        | ((rg.trait1 == "DEP_FinnGen_R12") & (rg.trait2 == "IBD_deLange2017"))
    ].copy()
    primary["lower"] = primary.rg - 1.96 * primary.se
    primary["upper"] = primary.rg + 1.96 * primary.se
    primary["label"] = np.where(
        primary.trait1.eq("MDD_Howard2019"),
        "Howard DEP ×\nFinnGen IBD",
        "FinnGen DEP ×\nde Lange IBD",
    )
    primary = primary.set_index("trait1").loc[
        ["MDD_Howard2019", "DEP_FinnGen_R12"]
    ].reset_index()

    local = read_tsv(derived / "local" / "bivar_q05_primary.tsv")
    same = read_tsv(derived / "local" / "same_block_eligible_both.tsv")

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 7.0))
    ax_a, ax_b, ax_c, ax_d = axes.flat
    fig.subplots_adjust(left=0.20, right=0.97, top=0.93, bottom=0.14, wspace=0.78, hspace=0.82)

    # A: reciprocal genome-wide correlations.
    y = np.arange(len(primary))[::-1]
    ax_a.axvline(0, color=GREY, linewidth=0.8)
    ax_a.errorbar(
        primary.rg, y,
        xerr=[primary.rg - primary.lower, primary.upper - primary.rg],
        fmt="o", markersize=6, markerfacecolor="white", markeredgecolor=TEXT,
        markeredgewidth=1.2, ecolor="#5F6368", elinewidth=1.1, capsize=3,
    )
    ax_a.set_yticks(y, primary.label)
    ax_a.set_xlim(-0.01, 0.21)
    ax_a.set_ylim(-0.6, 1.6)
    ax_a.set_xlabel("Genetic correlation (r$_g$)")
    ax_a.set_title("Reciprocal genome-wide correlation", loc="center", pad=22)
    ax_a.text(0.5, 1.02, "Estimates with 95% confidence intervals", transform=ax_a.transAxes,
              ha="center", va="bottom", fontsize=8.5, color="#555B61")
    panel_label(ax_a, "A", x=-0.24)
    clean_axis(ax_a)

    # B: pair-specific local correlations.
    local = local.copy()
    local["pair"] = np.where(local.pair_id.eq("HowardMDD__FinnGenIBD"), "Howard", "FinnGen")
    local["label"] = local.apply(
        lambda r: f"Locus {int(r.locus)}\nchr{int(r.chr)}:{r.start/1e6:.1f}–{r.stop/1e6:.1f} Mb", axis=1
    )
    local = local.sort_values(["pair", "locus"], ascending=[True, True]).reset_index(drop=True)
    y = np.arange(len(local))[::-1]
    ax_b.axvline(0, color=GREY, linewidth=0.8)
    for i, row in local.iterrows():
        color = TEAL if row["pair"] == "Howard" else CORAL
        ax_b.errorbar(
            row.rho, y[i], xerr=[[row.rho - row["rho.lower"]], [row["rho.upper"] - row.rho]],
            fmt="o", markersize=6, color=color, ecolor=color, elinewidth=1.1, capsize=3,
        )
    ax_b.set_yticks(y, local.label)
    ax_b.set_xlim(-1.1, 1.1)
    ax_b.set_xlabel("Local correlation (ρ)")
    ax_b.set_title("Pair-specific signals", loc="center", pad=22)
    ax_b.text(0.5, 1.02, f"{len(local)} q<0.05; no same-block replication", transform=ax_b.transAxes,
              ha="center", va="bottom", fontsize=8.5, color="#555B61")
    panel_label(ax_b, "B")
    clean_axis(ax_b)

    # C: blocks estimable in both reciprocal pairs.
    y = np.arange(len(same))[::-1]
    labels = [
        f"Locus {int(r.locus)}\nchr{int(r.chr)}:{r.start/1e6:.1f}–{r.stop/1e6:.1f} Mb"
        for _, r in same.iterrows()
    ]
    ax_c.axvline(0, color=GREY, linewidth=0.8)
    for i, row in same.iterrows():
        ax_c.errorbar(
            row.pair1_rho, y[i] + 0.08,
            xerr=[[row.pair1_rho - row.pair1_rho_lower], [row.pair1_rho_upper - row.pair1_rho]],
            fmt="o", markersize=5.5, color=TEAL, ecolor=TEAL, elinewidth=1.05, capsize=3,
        )
        ax_c.errorbar(
            row.pair2_rho, y[i] - 0.08,
            xerr=[[row.pair2_rho - row.pair2_rho_lower], [row.pair2_rho_upper - row.pair2_rho]],
            fmt="o", markersize=5.5, color=CORAL, ecolor=CORAL, elinewidth=1.05, capsize=3,
        )
    ax_c.set_yticks(y, labels)
    ax_c.set_xlim(-1.0, 0.85)
    ax_c.set_xlabel("Local genetic correlation (ρ; 95% CI)")
    ax_c.set_title("Blocks estimable in both pairs", loc="center", pad=22)
    ax_c.text(0.5, 1.02, f"{len(same)} blocks; no same-direction q<0.05", transform=ax_c.transAxes,
              ha="center", va="bottom", fontsize=8.5, color="#555B61")
    panel_label(ax_c, "C")
    clean_axis(ax_c)

    # D: evidence ladder.
    layers = ["Genome-wide LDSC", "Local LAVA", "Replicated PLACO", "Fine-mapping / QTL"]
    details = [
        "2/2 reciprocal pairs\npositive (P<0.025)",
        "0 same-block\ncovariance components",
        "3 non-MHC\npleiotropic blocks",
        "None established\nacross four GWAS",
    ]
    colors = ["#2E7D32", "#B3261E", "#2E7D32", "#6A5ACD"]
    y = np.arange(len(layers))[::-1]
    ax_d.barh(y, [0.26] * 4, left=0.02, height=0.55, color=colors)
    for yi, detail in zip(y, details):
        ax_d.text(0.34, yi, detail, ha="left", va="center", fontsize=8.7)
    ax_d.set_yticks(y, layers)
    ax_d.set_xlim(0, 1.12)
    ax_d.set_xticks([])
    ax_d.set_title("Evidence by scale", loc="center", pad=22)
    ax_d.text(0.5, 1.02, "Prespecified decisions at each analytical scale", transform=ax_d.transAxes,
              ha="center", va="bottom", fontsize=8.5, color="#555B61")
    panel_label(ax_d, "D")
    for spine in ax_d.spines.values():
        spine.set_visible(False)
    ax_d.tick_params(axis="y", length=0)

    legend = [
        Line2D([0], [0], color=TEAL, marker="o", linewidth=1.2, label="Howard DEP × FinnGen IBD"),
        Line2D([0], [0], color=CORAL, marker="o", linewidth=1.2, label="FinnGen DEP × de Lange IBD"),
    ]
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=2,
               frameon=False, columnspacing=1.8, handlelength=2.0)
    save_main_figure(fig, "1")
    plt.close(fig)


def build_s1() -> None:
    ext = REPO / "extensions" / "shared_locus_finemap_qtl" / "results"
    locus = read_tsv(ext / "locus_finemap_summary.tsv")
    evidence = read_tsv(ext / "candidate_evidence_tiers.tsv")
    pqtl = read_tsv(ext / "pqtl_signal_family_summary.tsv")

    traits = ["Howard_DEP", "FinnGen_DEP", "deLange_IBD", "FinnGen_IBD"]
    trait_labels = ["Howard DEP", "FinnGen DEP", "de Lange IBD", "FinnGen IBD"]
    colors = [BLUE, "#72B7B2", "#E45756", GOLD]

    fig = plt.figure(figsize=(10.0, 7.3))
    outer = fig.add_gridspec(
        2, 1, left=0.10, right=0.94, top=0.92, bottom=0.09,
        height_ratios=[1.0, 1.08], hspace=0.62,
    )
    top = outer[0].subgridspec(
        1, 4, width_ratios=[2.0, 1.35, 2.55, 0.20], wspace=0.48,
    )
    bottom = outer[1].subgridspec(
        1, 2, width_ratios=[5.3, 1.55], wspace=0.22,
    )
    ax_a = fig.add_subplot(top[0, 0])
    ax_a_leg = fig.add_subplot(top[0, 1])
    ax_b = fig.add_subplot(top[0, 2])
    ax_b_cbar = fig.add_subplot(top[0, 3])
    ax_c = fig.add_subplot(bottom[0, 0])
    ax_c_leg = fig.add_subplot(bottom[0, 1])

    # A: number of credible sets.
    x = np.arange(len(locus))
    width = 0.18
    for i, (trait, label, color) in enumerate(zip(traits, trait_labels, colors)):
        ax_a.bar(x + (i - 1.5) * width, locus[f"{trait}_n_cs"], width, color=color, label=label)
    ax_a.set_xticks(x, ["Block 154", "Block 464", "Block 1671"], rotation=18, ha="right")
    ax_a.set_ylabel("Number of 95% credible sets")
    ax_a.set_ylim(0, 2.2)
    ax_a.set_title("Regional fine-mapping", loc="center", pad=14)
    panel_label(ax_a, "A", x=-0.28)
    clean_axis(ax_a)
    ax_a_leg.axis("off")
    ax_a_leg.legend(
        handles=[Patch(facecolor=c, label=l) for c, l in zip(colors, trait_labels)],
        loc="center left", frameon=False, borderaxespad=0, handlelength=1.2,
    )

    # B: maximum SuSiE colocalisation H4.
    genes = [g for g in ["GPR25", "MST1", "FADS1", "TMEM258"] if g in set(evidence.gene_symbol)]
    matrix = np.full((len(genes), len(traits)), np.nan)
    for i, gene in enumerate(genes):
        selected = evidence[evidence.gene_symbol.eq(gene)]
        for j, trait in enumerate(traits):
            if len(selected) and pd.notna(selected[trait].iloc[0]):
                matrix[i, j] = selected[trait].iloc[0]
    image = ax_b.imshow(matrix, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    ax_b.set_xticks(range(len(traits)), ["Howard\nDEP", "FinnGen\nDEP", "de Lange\nIBD", "FinnGen\nIBD"])
    ax_b.set_yticks(range(len(genes)), genes)
    ax_b.set_title("Maximum SuSiE coloc PP.H4", loc="center", pad=14)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            label = "NE" if np.isnan(value) else f"{value:.2f}"
            ax_b.text(j, i, label, ha="center", va="center", fontsize=8.5,
                      color="white" if np.isfinite(value) and value < 0.45 else TEXT)
    fig.colorbar(image, cax=ax_b_cbar)
    panel_label(ax_b, "B", x=-0.18)

    # C: distinct MST1 pQTL signals.
    plot_p = pqtl[(pqtl.max_ibd_h4.ge(0.8)) | (pqtl.max_depression_h4.ge(0.8))].copy()
    plot_p = plot_p.sort_values(["max_ibd_h4", "max_depression_h4"], ascending=False).head(6)
    plot_p = plot_p.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_p))
    ax_c.barh(y + 0.18, plot_p.max_ibd_h4.fillna(0), 0.34, color="#E45756")
    ax_c.barh(y - 0.18, plot_p.max_depression_h4.fillna(0), 0.34, color=BLUE)
    ax_c.set_yticks(y, [x if isinstance(x, str) and x else "unlabelled" for x in plot_p.pqtl_variant_rsids])
    ax_c.set_xlim(0, 1.02)
    ax_c.axvline(0.8, color="#555555", linestyle="--", linewidth=0.9)
    ax_c.set_xlabel("Maximum Open Targets PP.H4")
    ax_c.set_title("MST1 plasma-pQTL signals", loc="center", pad=14)
    panel_label(ax_c, "C", x=-0.10)
    clean_axis(ax_c)
    ax_c_leg.axis("off")
    ax_c_leg.legend(
        handles=[Patch(facecolor="#E45756", label="IBD GWAS"), Patch(facecolor=BLUE, label="Depression GWAS")],
        loc="center left", bbox_to_anchor=(0.0, 0.5), frameon=False,
    )

    SUPP_OUT.mkdir(parents=True, exist_ok=True)
    pdf = SUPP_OUT / "S1_Fig.pdf"
    png = SUPP_OUT / "S1_Fig.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300, facecolor="white")
    copy_to_repository_figures(pdf)
    copy_to_repository_figures(png)
    plt.close(fig)


def build_s2() -> None:
    path = REPO / "extensions" / "annotation_stratified_gcov" / "results" / "PRIMARY_ANNOTATION_DECISION.tsv"
    data = read_tsv(path).sort_values("order").reset_index(drop=True)
    labels = {
        "Coding_UCSC": "Coding (UCSC)",
        "Conserved_LindbladToh": "Conserved (Lindblad-Toh)",
        "DHS_Trynka": "DNase hypersensitive sites",
        "Enhancer_Andersson": "Enhancer (Andersson)",
        "H3K27ac_Hnisz": "H3K27ac (Hnisz)",
        "H3K4me1_Trynka": "H3K4me1 (Trynka)",
        "Promoter_UCSC": "Promoter (UCSC)",
        "SuperEnhancer_Hnisz": "Super-enhancer (Hnisz)",
        "TFBS_ENCODE": "TFBS (ENCODE)",
        "Repressed_Hoffman": "Repressed (Hoffman)",
    }
    y = np.arange(len(data))[::-1]
    disc_lo = data.covariance_discovery - 1.96 * data.covariance_se_discovery
    disc_hi = data.covariance_discovery + 1.96 * data.covariance_se_discovery
    repl_lo = data.covariance_replication - 1.96 * data.covariance_se_replication
    repl_hi = data.covariance_replication + 1.96 * data.covariance_se_replication
    low = min(disc_lo.min(), repl_lo.min())
    high = max(disc_hi.max(), repl_hi.max())
    pad = (high - low) * 0.08

    fig, ax = plt.subplots(figsize=(7.5, 5.9))
    fig.subplots_adjust(left=0.35, right=0.78, top=0.86, bottom=0.14)
    ax.axvline(0, color=GREY, linewidth=0.8, linestyle="--")
    ax.errorbar(
        data.covariance_discovery, y + 0.12,
        xerr=[data.covariance_discovery - disc_lo, disc_hi - data.covariance_discovery],
        fmt="o", markersize=5, color="#0072B2", ecolor="#0072B2", elinewidth=1.0, capsize=2.5,
    )
    ax.errorbar(
        data.covariance_replication, y - 0.12,
        xerr=[data.covariance_replication - repl_lo, repl_hi - data.covariance_replication],
        fmt="^", markersize=5.5, color="#D55E00", ecolor="#D55E00", elinewidth=1.0, capsize=2.5,
    )
    ax.set_yticks(y, [labels.get(x, x.replace("_", " ")) for x in data.annotation])
    ax.set_xlim(low - pad, high + pad)
    ax.set_xlabel("")
    ax.set_title("Annotation-specific genetic covariance (95% CI)", loc="center", pad=18)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="#0072B2", linestyle="none", label="Discovery"),
            Line2D([0], [0], marker="^", color="#D55E00", linestyle="none", label="Reciprocal replication"),
        ],
        loc="upper left", bbox_to_anchor=(1.03, 1.0), frameon=False, borderaxespad=0,
    )

    pdf = SUPP_OUT / "S2_Fig.pdf"
    png = SUPP_OUT / "S2_Fig.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300, facecolor="white")
    copy_to_repository_figures(pdf)
    copy_to_repository_figures(png)
    plt.close(fig)


def main() -> None:
    global MAIN_OUT, SUPP_OUT
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--package-root", type=Path,
        help="Optional PLOS package root; figures are always copied to repository/figures.",
    )
    args = parser.parse_args()
    if args.package_root:
        package = args.package_root.resolve()
        MAIN_OUT = package / "02_Figures"
        SUPP_OUT = package / "03_Supporting_Information"
    build_figure_1()
    build_s1()
    build_s2()
    print("Rebuilt Fig 1, S1 Fig, and S2 Fig from frozen derived tables.")


if __name__ == "__main__":
    main()
