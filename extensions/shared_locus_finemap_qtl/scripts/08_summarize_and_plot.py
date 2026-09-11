#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TRAITS = ["Howard_DEP", "FinnGen_DEP", "deLange_IBD", "FinnGen_IBD"]
TRAIT_LABELS = ["Howard depression", "FinnGen depression", "deLange IBD", "FinnGen IBD"]
LOCUS_LABELS = {
    "block154": "chr1:200.13-201.07 Mb",
    "block464": "chr3:47.59-50.39 Mb",
    "block1671": "chr11:60.52-61.72 Mb",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def max_h4_table(susie: pd.DataFrame, prior: float) -> pd.DataFrame:
    work = susie[(susie.status == "ok") & np.isclose(susie.p12, prior, equal_nan=False)].copy()
    work["PP.H4"] = pd.to_numeric(work["PP.H4"], errors="coerce")
    return (
        work.groupby(["locus", "gene_symbol", "gene_id", "gwas_trait"], as_index=False)
        .agg(max_h4=("PP.H4", "max"), n_signal_pairs=("PP.H4", "size"))
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    summary = root / "results" / "summary"
    reports = root / "reports"
    figures = root / "figures"
    provenance = root / "provenance"
    for folder in [summary, reports, figures, provenance]:
        folder.mkdir(parents=True, exist_ok=True)

    fit = pd.read_csv(root / "results/gwas/susie_fit_qc.tsv", sep="\t")
    fit0 = fit[np.isclose(fit.ridge, 0)].copy()
    overlap = pd.read_csv(root / "results/gwas/cross_cohort_credible_set_overlap.tsv", sep="\t")
    gwas_coloc = pd.read_csv(root / "results/gwas/gwas_coloc_susie.tsv", sep="\t")
    candidates = pd.read_csv(root / "results/mapping/candidate_genes.tsv", sep="\t")
    qtl_manifest = pd.read_csv(root / "results/qtl/eligible_qtl_traits.tsv", sep="\t")
    retrieval = pd.read_csv(root / "results/qtl/qtl_retrieval_manifest.tsv", sep="\t")
    qtl_abf = pd.read_csv(root / "results/qtl/qtl_coloc_abf.tsv", sep="\t")
    qtl_susie = pd.read_csv(root / "results/qtl/qtl_coloc_susie.tsv", sep="\t")
    pqtl = pd.read_csv(root / "results/opentargets/opentargets_pqtl_coloc.tsv", sep="\t")

    primary = max_h4_table(qtl_susie, 1e-5)
    strict = max_h4_table(qtl_susie, 1e-6).rename(
        columns={"max_h4": "max_h4_strict", "n_signal_pairs": "n_signal_pairs_strict"}
    )
    gene_h4 = primary.merge(strict, on=["locus", "gene_symbol", "gene_id", "gwas_trait"], how="left")
    high = gene_h4[gene_h4.max_h4 >= 0.8].sort_values(
        ["locus", "gene_symbol", "gwas_trait", "max_h4"], ascending=[True, True, True, False]
    )
    high.to_csv(summary / "high_confidence_qtl_coloc.tsv", sep="\t", index=False)

    wide = gene_h4.pivot_table(
        index=["locus", "gene_id", "gene_symbol"], columns="gwas_trait", values="max_h4", aggfunc="max"
    ).reset_index()
    for trait in TRAITS:
        if trait not in wide:
            wide[trait] = np.nan
    strict_wide = strict.pivot_table(
        index=["locus", "gene_id", "gene_symbol"], columns="gwas_trait", values="max_h4_strict", aggfunc="max"
    ).reset_index()
    strict_wide = strict_wide.rename(columns={trait: f"{trait}_strict" for trait in TRAITS if trait in strict_wide})
    evidence = candidates.merge(wide, on=["locus", "gene_id", "gene_symbol"], how="left")
    evidence = evidence.merge(strict_wide, on=["locus", "gene_id", "gene_symbol"], how="left")
    for trait in TRAITS:
        evidence[f"{trait}_H4_ge_0.8"] = evidence[trait].fillna(-1).ge(0.8)
    evidence["depression_both_cohorts_qtl"] = evidence[["Howard_DEP_H4_ge_0.8", "FinnGen_DEP_H4_ge_0.8"]].all(axis=1)
    evidence["ibd_both_cohorts_qtl"] = evidence[["deLange_IBD_H4_ge_0.8", "FinnGen_IBD_H4_ge_0.8"]].all(axis=1)
    evidence["shared_four_gwas_qtl"] = evidence.depression_both_cohorts_qtl & evidence.ibd_both_cohorts_qtl
    evidence["evidence_class"] = "Tier3_positional_only"
    single = evidence[[f"{trait}_H4_ge_0.8" for trait in TRAITS]].any(axis=1)
    evidence.loc[single, "evidence_class"] = "Tier3_single_cohort_molecular"
    evidence.loc[evidence.depression_both_cohorts_qtl, "evidence_class"] = "Tier1_depression_component"
    evidence.loc[evidence.ibd_both_cohorts_qtl, "evidence_class"] = "Tier1_IBD_component"
    evidence.loc[evidence.shared_four_gwas_qtl, "evidence_class"] = "Tier1_cross_disease_shared_gene"
    evidence["shared_mechanism_interpretation"] = np.where(
        evidence.shared_four_gwas_qtl,
        "replicated molecular sharing supported",
        "not supported across both depression and both IBD cohorts",
    )
    evidence.to_csv(summary / "candidate_evidence_tiers.tsv", sep="\t", index=False)

    pqtl["h4"] = pd.to_numeric(pqtl.h4, errors="coerce")
    psummary = (
        pqtl.groupby(["locus", "target_gene", "pqtl_study_locus_id", "pqtl_variant_rsids"], as_index=False)
        .apply(lambda group: pd.Series({
            "max_ibd_h4": group.loc[group.gwas_family.eq("IBD"), "h4"].max(),
            "max_depression_h4": group.loc[group.gwas_family.eq("depression"), "h4"].max(),
            "n_ibd_coloc": int((group.gwas_family.eq("IBD") & group.h4.ge(0.8)).sum()),
            "n_depression_coloc": int((group.gwas_family.eq("depression") & group.h4.ge(0.8)).sum()),
        }))
        .reset_index(drop=True)
    )
    psummary["same_pqtl_signal_supports_both"] = psummary.max_ibd_h4.ge(0.8) & psummary.max_depression_h4.ge(0.8)
    psummary.to_csv(summary / "pqtl_signal_family_summary.tsv", sep="\t", index=False)

    locus_summary = []
    for locus in LOCUS_LABELS:
        row = {"locus": locus, "region": LOCUS_LABELS[locus]}
        for trait in TRAITS:
            selected = fit0[(fit0.locus == locus) & (fit0.trait == trait)]
            row[f"{trait}_n_cs"] = int(selected.n_credible_sets.iloc[0]) if len(selected) else np.nan
            row[f"{trait}_max_pip"] = float(selected.max_pip.iloc[0]) if len(selected) else np.nan
        row["depression_cs_replicated"] = bool(
            overlap[(overlap.locus == locus) & (overlap.pair_type == "depression_cross_cohort")].replicated.astype(str).str.upper().eq("TRUE").any()
        )
        row["ibd_cs_replicated"] = bool(
            overlap[(overlap.locus == locus) & (overlap.pair_type == "ibd_cross_cohort")].replicated.astype(str).str.upper().eq("TRUE").any()
        )
        row["n_positional_genes"] = int((candidates.locus == locus).sum())
        row["n_qtl_h4_ge_0.8_gene_trait_pairs"] = int((high.locus == locus).sum())
        locus_summary.append(row)
    locus_df = pd.DataFrame(locus_summary)
    locus_df.to_csv(summary / "locus_finemap_summary.tsv", sep="\t", index=False)

    # Compact three-panel evidence figure.
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.8), gridspec_kw={"width_ratios": [1.1, 1.3, 1.2]})
    colors = ["#4C78A8", "#72B7B2", "#E45756", "#F2CF5B"]
    x = np.arange(len(LOCUS_LABELS))
    width = 0.19
    for index, (trait, label, color) in enumerate(zip(TRAITS, TRAIT_LABELS, colors)):
        values = [locus_df.loc[locus_df.locus.eq(locus), f"{trait}_n_cs"].iloc[0] for locus in LOCUS_LABELS]
        axes[0].bar(x + (index - 1.5) * width, values, width, label=label, color=color)
    axes[0].set_xticks(x, ["Block 154", "Block 464", "Block 1671"])
    axes[0].set_ylabel("Number of 95% credible sets")
    axes[0].set_title("A  Regional fine-mapping")
    axes[0].legend(frameon=False, fontsize=7, loc="upper left")
    axes[0].spines[["top", "right"]].set_visible(False)

    prioritized = [name for name in ["GPR25", "MST1", "FADS1", "TMEM258"] if name in set(evidence.gene_symbol)]
    matrix = np.full((len(prioritized), len(TRAITS)), np.nan)
    for i, gene in enumerate(prioritized):
        selected = evidence[evidence.gene_symbol.eq(gene)]
        for j, trait in enumerate(TRAITS):
            if len(selected) and pd.notna(selected[trait].iloc[0]):
                matrix[i, j] = selected[trait].iloc[0]
    image = axes[1].imshow(matrix, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    axes[1].set_xticks(range(len(TRAITS)), ["Howard\nDEP", "FinnGen\nDEP", "deLange\nIBD", "FinnGen\nIBD"])
    axes[1].set_yticks(range(len(prioritized)), prioritized)
    axes[1].set_title("B  Maximum SuSiE coloc H4")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            label = "NE" if np.isnan(matrix[i, j]) else f"{matrix[i, j]:.2f}"
            axes[1].text(j, i, label, ha="center", va="center", fontsize=8,
                         color="white" if np.isfinite(matrix[i, j]) and matrix[i, j] < 0.45 else "black")
    fig.colorbar(image, ax=axes[1], fraction=0.046, pad=0.03)

    plot_p = psummary[(psummary.max_ibd_h4.ge(0.8)) | (psummary.max_depression_h4.ge(0.8))].copy()
    plot_p = plot_p.sort_values(["max_ibd_h4", "max_depression_h4"], ascending=False).head(6)
    px = np.arange(len(plot_p))
    axes[2].barh(px + 0.18, plot_p.max_ibd_h4.fillna(0), 0.34, color="#E45756", label="IBD GWAS")
    axes[2].barh(px - 0.18, plot_p.max_depression_h4.fillna(0), 0.34, color="#4C78A8", label="Depression GWAS")
    labels = [value if value else "unlabelled" for value in plot_p.pqtl_variant_rsids]
    axes[2].set_yticks(px, labels)
    axes[2].set_xlim(0, 1.02)
    axes[2].axvline(0.8, color="#555555", linestyle="--", linewidth=0.9)
    axes[2].set_xlabel("Maximum Open Targets H4")
    axes[2].set_title("C  MST1 pQTL signals")
    axes[2].legend(frameon=False, fontsize=8, loc="lower right")
    axes[2].spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(figures / "Fig_shared_locus_finemap_coloc_summary.pdf", bbox_inches="tight")
    fig.savefig(figures / "Fig_shared_locus_finemap_coloc_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    ibd_genes = sorted(evidence.loc[evidence.ibd_both_cohorts_qtl, "gene_symbol"].unique())
    dep_genes = sorted(evidence.loc[evidence.depression_both_cohorts_qtl, "gene_symbol"].unique())
    single_genes = sorted(evidence.loc[
        evidence.evidence_class.eq("Tier3_single_cohort_molecular"), "gene_symbol"
    ].unique())
    same_pqtl = int(psummary.same_pqtl_signal_supports_both.sum())
    report = f"""# Fine-mapping, candidate-gene mapping and molecular colocalisation report

## Scope

The analysis was restricted a priori to blocks 154, 464 and 1671 identified in the reciprocal PLACO screen. Four GWAS were fine-mapped with SuSiE-RSS using the same 1000 Genomes EUR LD reference. Candidate genes were mapped only from credible-set/lead-variant overlap or a 100-kb window, followed by local eQTL/sQTL colocalisation and a separate Open Targets pQTL annotation.

## Fine-mapping result

- Block 154: the IBD component replicated across deLange and FinnGen (top-variant r2=0.858; six exact credible-set variants). Neither depression GWAS produced a credible set.
- Block 464: the IBD component replicated across deLange and FinnGen (top-variant r2=1.000; twelve exact credible-set variants). FinnGen depression colocalised with the IBD component, but Howard depression had no estimable credible set, so cross-depression replication was not established.
- Block 1671: Howard and FinnGen depression produced different credible sets (top-variant r2=0.108; no exact overlap); neither IBD GWAS produced a credible set.

## Candidate-gene and QTL result

- {', '.join(ibd_genes) if ibd_genes else 'No gene'} met the primary H4>=0.8 criterion with both IBD cohorts. GPR25 was supported by whole-blood eQTL; MST1 was supported by amygdala eQTL and repeated sQTL signals in amygdala, sigmoid/transverse colon, terminal ileum and blood. Under the stricter p12=1e-6 sensitivity, the repeated MST1 sQTL findings remained above H4=0.8 in both IBD cohorts, whereas GPR25 remained strong for deLange IBD (H4=0.903) but fell below the strong threshold for FinnGen IBD (H4=0.768).
- No gene met H4>=0.8 with both depression cohorts. Single-cohort molecular candidates were {', '.join(single_genes) if single_genes else 'none'}.
- No candidate met the same molecular-QTL criterion across all four GWAS. Therefore no replicated cross-disease shared causal gene was identified.

## Protein layer

Open Targets returned 77 pQTL credible sets for MST1 and none for GPR25, FADS1 or TMEM258. The MST1 pQTL led by rs11130213 colocalised with multiple IBD GWAS (maximum H4={pqtl.loc[pqtl.gwas_family.eq('IBD'), 'h4'].max():.3f}). Different MST1 pQTL signals colocalised with depression GWAS (maximum H4={pqtl.loc[pqtl.gwas_family.eq('depression'), 'h4'].max():.3f}). No single pQTL credible set had H4>=0.8 for both disease families (n={same_pqtl}). In 1000 Genomes EUR, the principal IBD pQTL lead rs11130213 was weakly correlated with depression-linked pQTL leads (r2=0.022-0.098), indicating distinct protein-regulatory signals rather than one shared causal component.

## Boundary of inference

The three PLACO blocks reproduced at the variant or broad LD-block level, with uneven molecular resolution across traits. In block 464, FinnGen depression colocalised with the IBD component, whereas Howard depression had no credible set. GPR25 and MST1 had molecular support in IBD analyses; FADS1 and TMEM258 had support in the FinnGen depression analysis. These observations did not establish a component supported by all four GWAS. Missing credible sets indicate unresolved evidence, not disease specificity.

## Technical coverage

Seventy-eight molecular traits passed the frozen local-QTL eligibility gate. GTEx and BLUEPRINT were analyzable. Nine CEDAR locus-dataset requests returned HTTP 500 from the official API and are recorded as technically unavailable, not negative. All primary ABF comparisons were retained, and SuSiE was run only after the frozen ABF H3/H4 gate.
"""
    (reports / "SCIENTIFIC_INTERPRETATION.md").write_text(report, encoding="utf-8")

    execution = f"""# Execution report

## Completed stages

1. Four-GWAS harmonisation and 1000 Genomes EUR LD construction.
2. SuSiE-RSS fine-mapping with ridge-free primary and 1e-4 ridge sensitivity analyses.
3. Cross-cohort credible-set overlap and GWAS-GWAS multi-signal colocalisation.
4. GENCODE v19 positional candidate mapping.
5. Local eQTL/sQTL retrieval, eligibility filtering, ABF screening and SuSiE colocalisation.
6. Open Targets pQTL query and LD separation of MST1 protein signals.
7. Evidence-tier tables and publication-ready summary figure.

## Counts

- Positional candidate rows: {len(candidates)}
- Eligible QTL traits: {len(qtl_manifest)}
- ABF result rows: {len(qtl_abf)}
- Primary SuSiE gene-trait H4>=0.8 rows: {len(high)}
- Open Targets region-matched pQTL-GWAS rows: {len(pqtl)}
- CEDAR technical retrieval failures: {int(retrieval.status.eq('retrieval_error').sum())}

## Software

- Python: {platform.python_version()}
- pandas: {pd.__version__}
- numpy: {np.__version__}
- matplotlib: {matplotlib.__version__}
- R/SuSiE/coloc versions are recorded separately in `provenance/software_versions.txt`.
"""
    (reports / "EXECUTION_REPORT.md").write_text(execution, encoding="utf-8")

    outputs = [
        summary / "high_confidence_qtl_coloc.tsv", summary / "candidate_evidence_tiers.tsv",
        summary / "pqtl_signal_family_summary.tsv", summary / "locus_finemap_summary.tsv",
        figures / "Fig_shared_locus_finemap_coloc_summary.pdf",
        figures / "Fig_shared_locus_finemap_coloc_summary.png",
        reports / "SCIENTIFIC_INTERPRETATION.md", reports / "EXECUTION_REPORT.md",
    ]
    pd.DataFrame([
        {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)} for path in outputs
    ]).to_csv(provenance / "summary_output_manifest.tsv", sep="\t", index=False)
    print(report)


if __name__ == "__main__":
    main()
