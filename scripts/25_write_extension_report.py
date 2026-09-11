#!/usr/bin/env python3
import csv
import math
from pathlib import Path


ROOT = Path("/root/IBD/20_Reproducibility_Ladder")
GLOBAL = ROOT / "results" / "global_extension"
ATTR = GLOBAL / "attribution"
REPORT = ROOT / "reports" / "EXTENSION_UKB_DEPRESSION_CD_UC_REPORT.md"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def f(value, digits=3):
    number = float(value)
    if number == 0:
        return "0"
    if abs(number) < 0.001:
        return f"{number:.2e}"
    return f"{number:.{digits}f}"


h2 = read_tsv(GLOBAL / "h2_summary.tsv")
rg = read_tsv(GLOBAL / "rg_ukb_definition_matrix.tsv")
ratio = read_tsv(ATTR / "dispersion_ratio.tsv")[0]
decision = read_tsv(ATTR / "instability_attribution_decision.tsv")[0]
omnibus = read_tsv(ATTR / "omnibus_tests.tsv")
interactions = read_tsv(ATTR / "interaction_tests.tsv")
contrasts = read_tsv(ATTR / "sampling_covariance_aware_contrasts.tsv")
leave_one = read_tsv(ATTR / "leave_one_definition.tsv")
leave_one_subtype = read_tsv(ATTR / "leave_one_subtype.tsv")
cd_uc_direction = read_tsv(ATTR / "cd_uc_direction_by_cohort.tsv")

definition_order = [
    "UKB_LifetimeMDD",
    "UKB_MDDRecur",
    "UKB_GPpsy",
    "UKB_ICD10Dep",
]
ibd_order = [
    "IBD_deLange2017",
    "CD_deLange2017",
    "UC_deLange2017",
    "IBD_FinnGen_R12",
    "CD_FinnGen_R12",
    "UC_FinnGen_R12",
]
lookup = {(row["trait1"], row["trait2"]): row for row in rg}
significant_primary = [
    row for row in rg if float(row["fdr_primary_extension"]) < 0.05
]
minimum_rg = min(rg, key=lambda row: float(row["rg"]))
maximum_rg = max(rg, key=lambda row: float(row["rg"]))
direction_lookup = {row["ibd_cohort"]: row for row in cd_uc_direction}

lines = [
    "# UKB depression-definition and CD/UC extension",
    "",
    "## Question and design",
    "",
    (
        "This extension tests whether variability in IBD-depression genetic "
        "correlation is more closely associated with depression phenotype "
        "definition or IBD subtype. Four UK Biobank depression definitions "
        "were crossed with IBD, CD, and UC from de Lange 2017 and FinnGen R12. "
        "The UKB definitions share participants and are phenotype-sensitivity "
        "analyses, not independent replications."
    ),
    "",
    "## Layer-0 QC",
    "",
    "| Trait | h2_obs | SE | z | Intercept | Mean chi2 | h2 z>=4 |",
    "|---|---:|---:|---:|---:|---:|:---:|",
]
for row in h2:
    lines.append(
        "| {trait} | {h2} | {se} | {z} | {intercept} | {mean} | {gate} |".format(
            trait=row["trait"],
            h2=f(row["h2_obs"]),
            se=f(row["h2_se"]),
            z=f(row["h2_z"], 2),
            intercept=f(row["intercept"], 4),
            mean=f(row["mean_chisq"], 3),
            gate="yes" if row["passes_h2_z4"] == "1" else "no",
        )
    )

lines.extend(
    [
        "",
        "## Crossed UKB definition by IBD phenotype matrix",
        "",
        "Cells are rg (SE); a dagger marks |genetic-covariance intercept z| >= 2.",
        "",
        "| Depression definition | "
        + " | ".join(ibd_order)
        + " |",
        "|---|" + "|".join(["---:"] * len(ibd_order)) + "|",
    ]
)
for definition in definition_order:
    cells = []
    for ibd_trait in ibd_order:
        row = lookup[(definition, ibd_trait)]
        flag = "†" if row["gcov_intercept_flag"] == "1" else ""
        cells.append(f"{f(row['rg'])} ({f(row['se'])}){flag}")
    lines.append(f"| {definition} | " + " | ".join(cells) + " |")

lines.extend(
    [
        "",
        (
            f"All 24 estimates were positive, ranging from "
            f"{f(minimum_rg['rg'])} "
            f"({minimum_rg['trait1']} x {minimum_rg['trait2']}) to "
            f"{f(maximum_rg['rg'])} "
            f"({maximum_rg['trait1']} x {maximum_rg['trait2']}); "
            f"{len(significant_primary)}/24 passed BH-FDR < 0.05. "
            "These counts describe signal detection and are not used to "
            "attribute instability because precision differs across cells."
        ),
    ]
)

lines.extend(
    [
        "",
        "## Attribution of instability",
        "",
        f"Classification: **{decision['classification']}**.",
        "",
        decision["reason"],
        "",
        (
            "The adjusted between-definition SD was "
            f"{f(ratio['depression_definition_sd'])}; the adjusted "
            f"between-subtype SD was {f(ratio['ibd_subtype_sd'])}. "
            f"Their ratio was {f(ratio['ratio'])} "
            f"(simulation 95% CI {f(ratio['ratio_ci_low'])}-"
            f"{f(ratio['ratio_ci_high'])}; "
            f"Pr[ratio>1]={f(ratio['probability_ratio_gt_1'])})."
        ),
        "",
        (
            "The point estimate therefore leaned toward greater variation "
            "across depression definitions, but the interval was wide and "
            "crossed 1. Excluding the broad GP-consultation definition moved "
            f"the ratio to {f(decision['leave_gp_ratio'])}, so this tendency "
            "was not leave-one-definition robust."
        ),
        "",
        (
            "IBD subtype showed an omnibus difference "
            f"(P={f(next(row['p'] for row in omnibus if row['test'] == 'ibd_subtype'))}), "
            "whereas depression definition did not "
            f"(P={f(next(row['p'] for row in omnibus if row['test'] == 'depression_definition'))}). "
            "However, CD exceeded UC for "
            f"{direction_lookup['deLange2017']['n_cd_gt_uc']}/4 UKB "
            "definitions in de Lange, while CD was lower than UC for "
            f"{direction_lookup['FinnGen_R12']['n_cd_lt_uc']}/4 in FinnGen. "
            "This opposite cross-cohort direction and the absence of "
            "FDR-significant CD-versus-UC contrasts prevent a subtype-dominant "
            "interpretation."
        ),
        "",
        (
            "Taken together, the extension supports a consistently positive "
            "IBD-depression genetic-correlation direction, but it does not "
            "support assigning the observed magnitude instability primarily "
            "to either depression definition or IBD subtype."
        ),
        "",
        "### Omnibus and interaction tests",
        "",
        "| Test | Statistic | df | P |",
        "|---|---:|---:|---:|",
    ]
)
for row in omnibus + interactions:
    lines.append(
        f"| {row['test']} | {f(row['statistic'])} | "
        f"{row['df']} | {f(row['p'])} |"
    )

lines.extend(
    [
        "",
        "### Leave-one-definition sensitivity",
        "",
        "| Excluded UKB definition | Definition SD | Subtype SD | Ratio |",
        "|---|---:|---:|---:|",
    ]
)
for row in leave_one:
    lines.append(
        f"| {row['excluded_definition']} | "
        f"{f(row['depression_definition_sd'])} | "
        f"{f(row['ibd_subtype_sd'])} | {f(row['ratio'])} |"
    )

lines.extend(
    [
        "",
        "### Leave-one-subtype sensitivity",
        "",
        "| Excluded subtype | Definition SD | Retained-subtype SD | Ratio |",
        "|---|---:|---:|---:|",
    ]
)
for row in leave_one_subtype:
    lines.append(
        f"| {row['excluded_subtype']} | "
        f"{f(row['depression_definition_sd'])} | "
        f"{f(row['ibd_subtype_sd'])} | {f(row['ratio'])} |"
    )

lines.extend(
    [
        "",
        "### CD-versus-UC direction across UKB definitions",
        "",
        "| IBD cohort | CD > UC | CD < UC | Prespecified direction rule |",
        "|---|---:|---:|:---:|",
    ]
)
for row in cd_uc_direction:
    lines.append(
        f"| {row['ibd_cohort']} | {row['n_cd_gt_uc']} | "
        f"{row['n_cd_lt_uc']} | "
        f"{'pass' if row['passes_three_of_four'] == 'TRUE' else 'fail'} |"
    )

significant_definition = [
    row
    for row in contrasts
    if row["axis"] == "depression_definition"
    and float(row["q_within_axis"]) < 0.05
]
significant_subtype = [
    row
    for row in contrasts
    if row["axis"] == "ibd_subtype_CD_vs_UC"
    and float(row["q_within_axis"]) < 0.05
]
gcov_flags = [row for row in rg if row["gcov_intercept_flag"] == "1"]

lines.extend(
    [
        "",
        "## Contrast summary",
        "",
        (
            f"Sampling-covariance-aware contrasts identified "
            f"{len(significant_definition)} FDR-significant pairwise "
            f"differences among UKB depression definitions and "
            f"{len(significant_subtype)} FDR-significant CD-versus-UC "
            f"differences. {len(gcov_flags)} of 24 primary extension cells "
            "had an absolute genetic-covariance intercept z-score of at least "
            "2 and are flagged rather than removed."
        ),
        "",
        "## Interpretation boundaries",
        "",
        (
            "This analysis attributes variability in observed LDSC rg "
            "estimates. It does not establish which disease causes the other, "
            "and it does not partition biological mechanisms. The UKB "
            "definitions are correlated because they use overlapping "
            "participants. The multivariable LDSC sampling covariance matrix "
            "was therefore used for inferential contrasts. Significance "
            "counts alone were not used to classify stability."
        ),
        "",
        "## Reproducibility files",
        "",
        "- Prespecified protocol: "
        "`provenance/extension_prespecified_protocol.md`",
        "- Pairwise LDSC output: `results/global_extension/`",
        "- GenomicSEM S and V matrices: "
        "`results/global_extension/genomicsem/`",
        "- Attribution models and contrasts: "
        "`results/global_extension/attribution/`",
        "- Figure panels: `figures/extension/`",
        "- Final checksums: `provenance/extension_final_sha256.tsv`",
    ]
)

REPORT.parent.mkdir(parents=True, exist_ok=True)
REPORT.write_text("\n".join(lines) + "\n")
print(REPORT)
