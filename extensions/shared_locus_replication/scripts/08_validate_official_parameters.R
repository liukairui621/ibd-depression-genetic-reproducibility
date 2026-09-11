#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication"
source(file.path(root, "software/PLACO/PLACO_v0.2.0.R"))
streamed <- fread(file.path(root, "results/raw/pair_harmonization_stats.tsv"))
pairs <- c("External_Howard_deLange", "FinnGen_DEP_IBD")
out <- vector("list", length(pairs))

for (i in seq_along(pairs)) {
  target_pair <- pairs[i]
  message("Validating official estimators for ", target_pair)
  path <- file.path(root, "results/raw", paste0(target_pair, ".harmonized.tsv.gz"))
  dat <- fread(cmd = paste("zcat", shQuote(path)), select = c("Z_DEP", "Z_IBD", "P_DEP", "P_IBD"), showProgress = FALSE)
  z <- as.matrix(dat[, .(Z_DEP, Z_IBD)])
  p <- as.matrix(dat[, .(P_DEP, P_IBD)])
  official_var <- var.placo(z, p, p.threshold = 1e-4)
  official_cor <- cor.pearson(z, p, p.threshold = 1e-4, returnMatrix = FALSE)
  ref <- streamed[streamed[["pair_id"]] == target_pair]
  if (nrow(ref) != 1) stop("Expected one streamed parameter row for ", target_pair)
  diffs <- c(
    abs(official_var[1] - ref$var_z_dep),
    abs(official_var[2] - ref$var_z_ibd),
    abs(official_cor - ref$cor_z)
  )
  out[[i]] <- data.frame(
    pair_id = target_pair,
    streamed_var_z_dep = ref$var_z_dep,
    official_var_z_dep = official_var[1],
    abs_diff_var_z_dep = diffs[1],
    streamed_var_z_ibd = ref$var_z_ibd,
    official_var_z_ibd = official_var[2],
    abs_diff_var_z_ibd = diffs[2],
    streamed_cor_z = ref$cor_z,
    official_cor_z = official_cor,
    abs_diff_cor_z = diffs[3],
    tolerance = 1e-10,
    status = if (all(diffs < 1e-10)) "PASS" else "FAIL"
  )
  rm(dat, z, p)
  gc()
}

out <- rbindlist(out)
fwrite(out, file.path(root, "reports/OFFICIAL_PARAMETER_VALIDATION.tsv"), sep = "\t")
if (any(out$status != "PASS")) stop("Official parameter validation failed")
capture.output(sessionInfo(), file = file.path(root, "provenance/R_session_parameter_validation.txt"))
