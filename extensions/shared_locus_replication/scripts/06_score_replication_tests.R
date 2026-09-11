#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication"
source(file.path(root, "software/PLACO/PLACO_v0.2.0.R"))

score <- function(z_dep, z_ibd, method, var_dep, var_ibd, corz) {
  if (!is.finite(z_dep) || !is.finite(z_ibd)) return(NA_real_)
  varz <- c(var_dep, var_ibd)
  ans <- if (method == "placo_plus") {
    placo.plus(c(z_dep, z_ibd), VarZ = varz, CorZ = corz)$p.placo.plus
  } else {
    placo(c(z_dep, z_ibd), VarZ = varz)$p.placo
  }
  pmin(1, pmax(0, ans))
}

rep_path <- file.path(root, "results/loci/replication_tests_input.tsv")
rep <- fread(rep_path)
if (nrow(rep)) {
  rep[, VALIDATION_P_PLACO := vapply(seq_len(.N), function(i) {
    if (VALIDATION_FOUND[i] != "TRUE") return(NA_real_)
    score(VALIDATION_Z_DEP[i], VALIDATION_Z_IBD[i], VALIDATION_METHOD[i],
          VALIDATION_VAR_Z_DEP[i], VALIDATION_VAR_Z_IBD[i], VALIDATION_COR_Z[i])
  }, numeric(1))]
} else {
  rep[, VALIDATION_P_PLACO := numeric()]
}
fwrite(rep, file.path(root, "results/loci/replication_tests_scored.tsv"), sep = "\t", na = "")

support_path <- file.path(root, "results/loci/four_pair_support_input.tsv")
support <- fread(support_path)
if (nrow(support)) {
  support[, P_PLACO := vapply(seq_len(.N), function(i) {
    if (FOUND[i] != "TRUE") return(NA_real_)
    score(Z_DEP[i], Z_IBD[i], METHOD[i], VAR_Z_DEP[i], VAR_Z_IBD[i], COR_Z[i])
  }, numeric(1))]
} else {
  support[, P_PLACO := numeric()]
}
fwrite(support, file.path(root, "results/loci/four_pair_lead_snp_support.tsv"), sep = "\t", na = "")
capture.output(sessionInfo(), file = file.path(root, "provenance/R_session_replication_scoring.txt"))
