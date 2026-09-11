#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    checks = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    gfit = pd.read_csv(root / "results/gwas/susie_fit_qc.tsv", sep="\t")
    check("all_gwas_susie_converged", bool(gfit.converged.all()), f"{int(gfit.converged.sum())}/{len(gfit)}")

    overlap = pd.read_csv(root / "results/gwas/cross_cohort_credible_set_overlap.tsv", sep="\t")
    rep = overlap.replicated.astype(str).str.upper().eq("TRUE")
    expected_ibd = set(overlap.loc[rep & overlap.pair_type.eq("ibd_cross_cohort"), "locus"])
    expected_dep = set(overlap.loc[rep & overlap.pair_type.eq("depression_cross_cohort"), "locus"])
    check("replicated_ibd_components", expected_ibd == {"block154", "block464"}, str(sorted(expected_ibd)))
    check("no_replicated_depression_component", not expected_dep, str(sorted(expected_dep)))

    retrieval = pd.read_csv(root / "results/qtl/qtl_retrieval_manifest.tsv", sep="\t")
    cedar_errors = retrieval.dataset_label.str.startswith("CEDAR") & retrieval.status.eq("retrieval_error")
    check("cedar_failures_classified_technical", int(cedar_errors.sum()) == 9, f"{int(cedar_errors.sum())}/9")
    check("non_cedar_retrieval_complete", not retrieval.loc[~retrieval.dataset_label.str.startswith("CEDAR"), "status"].eq("retrieval_error").any(), "GTEx and BLUEPRINT have no retrieval errors")

    eligible = pd.read_csv(root / "results/qtl/eligible_qtl_traits.tsv", sep="\t")
    check("qtl_gate_count", len(eligible) == 78, str(len(eligible)))
    check("qtl_gate_minimum_snps", bool(eligible.n_aligned.ge(200).all()), f"minimum={int(eligible.n_aligned.min())}")
    check("qtl_gate_cis_signal", bool(eligible.min_qtl_p.lt(1e-5).all()), f"maximum min-P={eligible.min_qtl_p.max():.3g}")

    abf = pd.read_csv(root / "results/qtl/qtl_coloc_abf.tsv", sep="\t")
    check("complete_abf_grid", len(abf) == len(eligible) * 4 * 2, f"{len(abf)} vs {len(eligible) * 8}")
    check("all_abf_runs_ok", bool(abf.status.eq("ok").all()), abf.status.value_counts().to_dict().__str__())

    qfit = pd.read_csv(root / "results/qtl/qtl_susie_fit_qc.tsv", sep="\t")
    check("all_screened_qtl_susie_converged", bool(qfit.qtl_converged.all()), f"{int(qfit.qtl_converged.sum())}/{len(qfit)}")
    check("all_screened_gwas_susie_converged", bool(qfit.gwas_converged.all()), f"{int(qfit.gwas_converged.sum())}/{len(qfit)}")

    evidence = pd.read_csv(root / "results/summary/candidate_evidence_tiers.tsv", sep="\t")
    ibd = set(evidence.loc[evidence.ibd_both_cohorts_qtl, "gene_symbol"])
    dep = set(evidence.loc[evidence.depression_both_cohorts_qtl, "gene_symbol"])
    shared = set(evidence.loc[evidence.shared_four_gwas_qtl, "gene_symbol"])
    single = set(evidence.loc[evidence.evidence_class.eq("Tier3_single_cohort_molecular"), "gene_symbol"])
    check("ibd_component_genes", ibd == {"GPR25", "MST1"}, str(sorted(ibd)))
    check("no_cross_cohort_depression_qtl_gene", not dep, str(sorted(dep)))
    check("no_four_gwas_shared_gene", not shared, str(sorted(shared)))
    check("single_cohort_depression_candidates", single == {"FADS1", "TMEM258"}, str(sorted(single)))

    psummary = pd.read_csv(root / "results/summary/pqtl_signal_family_summary.tsv", sep="\t")
    check("no_single_pqtl_signal_supports_both_families", not psummary.same_pqtl_signal_supports_both.any(), f"n={int(psummary.same_pqtl_signal_supports_both.sum())}")
    pqtl = pd.read_csv(root / "results/opentargets/opentargets_pqtl_coloc.tsv", sep="\t")
    check("mst1_pqtl_ibd_support", pqtl.loc[pqtl.gwas_family.eq("IBD"), "h4"].max() >= 0.8, f"max H4={pqtl.loc[pqtl.gwas_family.eq('IBD'), 'h4'].max():.6f}")
    check("mst1_pqtl_depression_support_distinct", pqtl.loc[pqtl.gwas_family.eq("depression"), "h4"].max() >= 0.8, f"max H4={pqtl.loc[pqtl.gwas_family.eq('depression'), 'h4'].max():.6f}")

    ld_path = root / "results/opentargets/ld_checks/MST1_selected_pairwise_ld.tsv"
    ld = pd.read_csv(ld_path, sep=r"\s+", header=None, names=["chr_a", "bp_a", "snp_a", "chr_b", "bp_b", "snp_b", "r2"])
    depression_leads = {"rs6774202", "rs111812549", "rs113465792"}
    dep_ld = ld[ld.snp_b.isin(depression_leads)]
    check("mst1_ibd_vs_depression_pqtl_signals_low_ld", len(dep_ld) == 3 and dep_ld.r2.max() < 0.1, f"r2 range={dep_ld.r2.min():.6f}-{dep_ld.r2.max():.6f}")

    frame = pd.DataFrame(checks)
    out = root / "provenance/FINAL_VALIDATION.tsv"
    frame.to_csv(out, sep="\t", index=False)
    passed = frame.status.eq("PASS").all()
    report = "# Final validation\n\n" + "\n".join(
        f"- {row.status}: {row.check} ({row.detail})" for row in frame.itertuples(index=False)
    ) + f"\n\nOverall: {'PASS' if passed else 'FAIL'}\n"
    report_path = root / "reports/FINAL_VALIDATION.md"
    report_path.write_text(report, encoding="utf-8")
    (root / "provenance/FINAL_VALIDATION.sha256").write_text(
        f"{sha256(out)}  {out}\n{sha256(report_path)}  {report_path}\n", encoding="utf-8"
    )
    print(report)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
