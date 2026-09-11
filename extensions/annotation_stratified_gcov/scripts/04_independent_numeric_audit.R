#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260823_annotation_stratified_gcov"
x <- read.delim(file.path(root, "results", "PRIMARY_ANNOTATION_DECISION.tsv"), check.names = FALSE)

stopifnot(nrow(x) == 10L)
stopifnot(identical(x$order, seq_len(10L)))
stopifnot(all(is.finite(x$covariance_discovery)))
stopifnot(all(is.finite(x$covariance_replication)))

p_discovery_recalc <- 2 * pnorm(
  abs(x$covariance_discovery / x$covariance_se_discovery),
  lower.tail = FALSE
)
p_replication_recalc <- 2 * pnorm(
  abs(x$covariance_replication / x$covariance_se_replication),
  lower.tail = FALSE
)
q_recalc <- p.adjust(p_discovery_recalc, method = "BH")
same_direction_recalc <- sign(x$covariance_discovery) == sign(x$covariance_replication)
p_directional_recalc <- ifelse(
  same_direction_recalc,
  p_replication_recalc / 2,
  1 - p_replication_recalc / 2
)
replicated_recalc <- (
  q_recalc < 0.05 &
    same_direction_recalc &
    p_directional_recalc < 0.05 &
    x$positive_h2_both_pairs
)

tol <- 1e-12
stopifnot(max(abs(p_discovery_recalc - x$covariance_p_two_sided_discovery)) < tol)
stopifnot(max(abs(p_replication_recalc - x$covariance_p_two_sided_replication)) < tol)
stopifnot(max(abs(q_recalc - x$discovery_fdr_bh)) < tol)
stopifnot(identical(same_direction_recalc, x$same_direction))
stopifnot(max(abs(p_directional_recalc - x$replication_p_directional)) < tol)
stopifnot(identical(replicated_recalc, x$replicated_annotation_covariance))

lines <- c(
  "INDEPENDENT NUMERIC AUDIT: PASS",
  paste0("n_frozen_annotations=", nrow(x)),
  paste0("n_discovery_fdr_lt_0.05=", sum(q_recalc < 0.05)),
  paste0("n_replicated_annotation_covariance=", sum(replicated_recalc)),
  paste0(
    "replicated_annotations=",
    paste(x$annotation[replicated_recalc], collapse = ",")
  ),
  paste0("tolerance=", tol),
  "checks=p_values,BH_FDR,direction,directional_replication_p,decision_boolean"
)
writeLines(lines, file.path(root, "reports", "INDEPENDENT_NUMERIC_AUDIT.txt"))

