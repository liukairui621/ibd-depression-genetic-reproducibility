#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260823_annotation_stratified_gcov"
config <- read.delim(file.path(root, "config", "primary_annotations.tsv"), check.names = FALSE)
discovery <- read.delim(file.path(root, "results", "discovery", "annotation_covariance_all.tsv"), check.names = FALSE)
replication <- read.delim(file.path(root, "results", "replication", "annotation_covariance_all.tsv"), check.names = FALSE)

pick <- function(x, suffix) {
  x <- merge(config, x, by = "annotation", all.x = TRUE, sort = FALSE)
  x <- x[order(x$order), ]
  keep <- c(
    "annotation", "annotation_prop", "covariance", "covariance_se",
    "covariance_z", "covariance_p_two_sided", "covariance_enrichment",
    "h2_trait_1", "h2_trait_2", "correlation_descriptive",
    "estimable_positive_h2"
  )
  x <- x[, keep]
  names(x)[-1] <- paste0(names(x)[-1], "_", suffix)
  x
}

summary <- merge(
  config,
  merge(pick(discovery, "discovery"), pick(replication, "replication"), by = "annotation"),
  by = "annotation",
  sort = FALSE
)
summary <- summary[order(summary$order), ]
summary$discovery_fdr_bh <- p.adjust(summary$covariance_p_two_sided_discovery, method = "BH")
summary$same_direction <- sign(summary$covariance_discovery) == sign(summary$covariance_replication)
summary$replication_p_directional <- ifelse(
  summary$same_direction,
  summary$covariance_p_two_sided_replication / 2,
  1 - summary$covariance_p_two_sided_replication / 2
)
summary$positive_h2_both_pairs <- (
  summary$estimable_positive_h2_discovery &
    summary$estimable_positive_h2_replication
)
summary$replicated_annotation_covariance <- (
  summary$discovery_fdr_bh < 0.05 &
    summary$same_direction &
    summary$replication_p_directional < 0.05 &
    summary$positive_h2_both_pairs
)
summary$replicated_concentration <- (
  summary$replicated_annotation_covariance &
    summary$covariance_enrichment_discovery > 1 &
    summary$covariance_enrichment_replication > 1
)
summary$decision <- ifelse(
  summary$replicated_concentration,
  "replicated_concentration",
  ifelse(
    summary$replicated_annotation_covariance,
    "replicated_covariance_not_enriched_in_both",
    ifelse(
      !summary$positive_h2_both_pairs,
      "not_estimable_positive_h2",
      ifelse(
        summary$discovery_fdr_bh < 0.05 & !summary$same_direction,
        "discovery_fdr_opposite_replication_direction",
        ifelse(
          summary$discovery_fdr_bh < 0.05,
          "discovery_only",
          "not_discovery_fdr"
        )
      )
    )
  )
)

write.table(
  summary,
  file.path(root, "results", "PRIMARY_ANNOTATION_DECISION.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

base_discovery <- discovery[discovery$annotation == "base", ]
base_replication <- replication[replication$annotation == "base", ]
n_rep <- sum(summary$replicated_annotation_covariance, na.rm = TRUE)
n_conc <- sum(summary$replicated_concentration, na.rm = TRUE)

report <- c(
  "# Annotation-stratified genetic covariance report",
  "",
  "## Frozen decision",
  "",
  sprintf("- Replicated annotation-level covariance: %d of 10 annotations.", n_rep),
  sprintf("- Replicated concentration above SNP proportion: %d of 10 annotations.", n_conc),
  "- Primary regression used HapMap3 weights excluding the extended MHC.",
  "",
  "## Genome-wide base annotation",
  "",
  sprintf(
    "- Discovery (deLange IBD x FinnGen depression): covariance %.6g (SE %.6g, P %.4g).",
    base_discovery$covariance, base_discovery$covariance_se,
    base_discovery$covariance_p_two_sided
  ),
  sprintf(
    "- Replication (FinnGen IBD x Howard public no-23andMe depression): covariance %.6g (SE %.6g, P %.4g).",
    base_replication$covariance, base_replication$covariance_se,
    base_replication$covariance_p_two_sided
  ),
  "",
  "## Annotation decisions",
  ""
)
for (i in seq_len(nrow(summary))) {
  report <- c(
    report,
    sprintf(
      "- %s: discovery covariance %.6g (P %.4g, BH-FDR %.4g); replication covariance %.6g (directional P %.4g); %s.",
      summary$annotation[[i]],
      summary$covariance_discovery[[i]],
      summary$covariance_p_two_sided_discovery[[i]],
      summary$discovery_fdr_bh[[i]],
      summary$covariance_replication[[i]],
      summary$replication_p_directional[[i]],
      summary$decision[[i]]
    )
  )
}
report <- c(
  report,
  "",
  "## Interpretation boundary",
  "",
  if (n_rep == 0L) {
    "No frozen annotation carried reproducible covariance across the reciprocal cohort pairs. Under the prespecified stopping rule, this extension does not support a stable functional concentration and does not trigger further shared-gene or shared-pathway searches."
  } else {
    "At least one frozen annotation carried reproducible covariance across reciprocal cohort pairs. This localizes covariance to a broad genomic context but does not identify a shared causal pathway, cell type, gene or treatment target."
  },
  "",
  "The prior three-locus fine-mapping result remains the local adjudication: no cross-disease shared causal gene was identified."
)
writeLines(report, file.path(root, "reports", "ANNOTATION_STRATIFIED_GCOV_REPORT.md"))

pdf(file.path(root, "figures", "annotation_covariance_forest.pdf"), width = 9, height = 7)
op <- par(mar = c(5, 12, 2, 1))
y <- rev(seq_len(nrow(summary)))
lim <- range(c(
  summary$covariance_discovery + c(-1, 1) * 1.96 * summary$covariance_se_discovery,
  summary$covariance_replication + c(-1, 1) * 1.96 * summary$covariance_se_replication
), finite = TRUE)
plot(
  NA,
  xlim = lim,
  ylim = c(0.5, nrow(summary) + 0.5),
  yaxt = "n",
  ylab = "",
  xlab = "Annotation-specific genetic covariance (95% CI)",
  bty = "n"
)
axis(2, at = y, labels = summary$annotation, las = 1, tick = FALSE)
abline(v = 0, lty = 2, col = "grey60")
arrows(
  summary$covariance_discovery - 1.96 * summary$covariance_se_discovery,
  y + 0.12,
  summary$covariance_discovery + 1.96 * summary$covariance_se_discovery,
  y + 0.12,
  angle = 90,
  code = 3,
  length = 0.03,
  col = "#0072B2"
)
points(summary$covariance_discovery, y + 0.12, pch = 16, col = "#0072B2")
arrows(
  summary$covariance_replication - 1.96 * summary$covariance_se_replication,
  y - 0.12,
  summary$covariance_replication + 1.96 * summary$covariance_se_replication,
  y - 0.12,
  angle = 90,
  code = 3,
  length = 0.03,
  col = "#D55E00"
)
points(summary$covariance_replication, y - 0.12, pch = 17, col = "#D55E00")
legend(
  "topright",
  legend = c("Discovery", "Reciprocal replication"),
  pch = c(16, 17),
  col = c("#0072B2", "#D55E00"),
  bty = "n"
)
par(op)
dev.off()

png(file.path(root, "figures", "annotation_covariance_forest.png"), width = 1800, height = 1400, res = 200)
op <- par(mar = c(5, 12, 2, 1))
plot(
  NA,
  xlim = lim,
  ylim = c(0.5, nrow(summary) + 0.5),
  yaxt = "n",
  ylab = "",
  xlab = "Annotation-specific genetic covariance (95% CI)",
  bty = "n"
)
axis(2, at = y, labels = summary$annotation, las = 1, tick = FALSE)
abline(v = 0, lty = 2, col = "grey60")
arrows(summary$covariance_discovery - 1.96 * summary$covariance_se_discovery, y + 0.12, summary$covariance_discovery + 1.96 * summary$covariance_se_discovery, y + 0.12, angle = 90, code = 3, length = 0.03, col = "#0072B2")
points(summary$covariance_discovery, y + 0.12, pch = 16, col = "#0072B2")
arrows(summary$covariance_replication - 1.96 * summary$covariance_se_replication, y - 0.12, summary$covariance_replication + 1.96 * summary$covariance_se_replication, y - 0.12, angle = 90, code = 3, length = 0.03, col = "#D55E00")
points(summary$covariance_replication, y - 0.12, pch = 17, col = "#D55E00")
legend("topright", legend = c("Discovery", "Reciprocal replication"), pch = c(16, 17), col = c("#0072B2", "#D55E00"), bty = "n")
par(op)
dev.off()

