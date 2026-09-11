#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(
  "/root/IBD/20_Reproducibility_Ladder/software/R_lib",
  .libPaths()
))
suppressPackageStartupMessages(library(GenomicSEM))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: 01_run_sldsc_pair.R PAIR_ID TRAIT1 TRAIT2 NAME1 NAME2")
}

pair_id <- args[[1]]
traits <- args[2:3]
trait_names <- args[4:5]
root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260823_annotation_stratified_gcov"
out_dir <- file.path(root, "results", pair_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.exists(traits)))

rds_file <- file.path(out_dir, "sldsc_full_output.rds")
if (file.exists(rds_file)) {
  message("Loading completed s_ldsc RDS: ", rds_file)
  result <- readRDS(rds_file)
} else {
  result <- s_ldsc(
    traits = traits,
    sample.prev = c(NA_real_, NA_real_),
    population.prev = c(NA_real_, NA_real_),
    ld = "/root/IBD/00_RawData/Reference/SLDSC/baselineLD.",
    wld = paste0(
      "/root/IBD/00_RawData/Reference/LD_EUR/",
      "1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
    ),
    frq = paste0(
      "/root/IBD/00_RawData/Reference/LD_EUR/",
      "1000G_Phase3_frq/1000G.EUR.QC."
    ),
    trait.names = trait_names,
    n.blocks = 200,
    ldsc.log = file.path(out_dir, pair_id),
    exclude_cont = TRUE
  )
  saveRDS(result, rds_file)
}

annotation_names_raw <- names(result$S)
if (is.null(annotation_names_raw) || length(annotation_names_raw) == 0L) {
  stop("s_ldsc returned unnamed annotations")
}
annotation_names <- sub("L2$", "", annotation_names_raw)

prop <- result$Prop$Prop
if (length(prop) != length(annotation_names)) {
  stop("Annotation proportion vector does not match S output")
}

extract_row <- function(i) {
  s <- result$S[[i]]
  v <- result$V[[i]]
  st <- result$S_Tau[[i]]
  vt <- result$V_Tau[[i]]
  cov_est <- s[2, 1]
  cov_se <- sqrt(v[2, 2])
  tau_est <- st[2, 1]
  tau_se <- sqrt(vt[2, 2])
  h2_1 <- s[1, 1]
  h2_2 <- s[2, 2]
  corr <- if (is.finite(h2_1) && is.finite(h2_2) && h2_1 > 0 && h2_2 > 0) {
    cov_est / sqrt(h2_1 * h2_2)
  } else {
    NA_real_
  }
  data.frame(
    pair_id = pair_id,
    trait_1 = trait_names[[1]],
    trait_2 = trait_names[[2]],
    annotation = annotation_names[[i]],
    annotation_prop = prop[[i]],
    covariance = cov_est,
    covariance_se = cov_se,
    covariance_z = cov_est / cov_se,
    covariance_p_two_sided = 2 * pnorm(abs(cov_est / cov_se), lower.tail = FALSE),
    tau = tau_est,
    tau_se = tau_se,
    h2_trait_1 = h2_1,
    h2_trait_2 = h2_2,
    correlation_descriptive = corr,
    estimable_positive_h2 = is.finite(corr)
  )
}

all_annotations <- do.call(rbind, lapply(seq_along(annotation_names), extract_row))
base_cov <- all_annotations$covariance[all_annotations$annotation == "base"]
if (length(base_cov) != 1L) {
  stop("Expected exactly one base annotation")
}
all_annotations$covariance_enrichment <- (
  all_annotations$covariance / base_cov
) / all_annotations$annotation_prop

write.table(
  all_annotations,
  file.path(out_dir, "annotation_covariance_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  result$I,
  file.path(out_dir, "cross_trait_intercept.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
