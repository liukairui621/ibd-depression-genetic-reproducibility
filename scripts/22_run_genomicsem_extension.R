#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(
  "/root/IBD/20_Reproducibility_Ladder/software/R_lib",
  .libPaths()
))
suppressPackageStartupMessages(library(GenomicSEM))

root <- "/root/IBD/20_Reproducibility_Ladder"
munged <- file.path(root, "data", "munged")
out <- file.path(root, "results", "global_extension", "genomicsem")
ld_dir <- file.path(root, "data", "genomicsem_reference", "ld")
wld_dir <- file.path(root, "data", "genomicsem_reference", "wld")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wld_dir, recursive = TRUE, showWarnings = FALSE)

ld_source <- "/root/IBD/00_RawData/Reference/LD_EUR/LDscore"
wld_source <- paste0(
  "/root/IBD/00_RawData/Reference/LD_EUR/",
  "1000G_Phase3_weights_hm3_no_MHC"
)
for (chromosome in 1:22) {
  file.symlink(
    file.path(ld_source, paste0("LDscore.", chromosome, ".l2.ldscore.gz")),
    file.path(ld_dir, paste0(chromosome, ".l2.ldscore.gz"))
  )
  file.symlink(
    file.path(ld_source, paste0("LDscore.", chromosome, ".l2.M_5_50")),
    file.path(ld_dir, paste0(chromosome, ".l2.M_5_50"))
  )
  file.symlink(
    file.path(
      wld_source,
      paste0("weights.hm3_noMHC.", chromosome, ".l2.ldscore.gz")
    ),
    file.path(wld_dir, paste0(chromosome, ".l2.ldscore.gz"))
  )
}

trait_names <- c(
  "MDD_Howard2019",
  "DEP_FinnGen_R12",
  "UKB_LifetimeMDD",
  "UKB_MDDRecur",
  "UKB_GPpsy",
  "UKB_ICD10Dep",
  "IBD_deLange2017",
  "CD_deLange2017",
  "UC_deLange2017",
  "IBD_FinnGen_R12",
  "CD_FinnGen_R12",
  "UC_FinnGen_R12"
)
traits <- file.path(munged, paste0(trait_names, ".sumstats.gz"))
stopifnot(all(file.exists(traits)))

result <- ldsc(
  traits = traits,
  sample.prev = rep(NA_real_, length(traits)),
  population.prev = rep(NA_real_, length(traits)),
  ld = ld_dir,
  wld = wld_dir,
  trait.names = trait_names,
  sep_weights = TRUE,
  chr = 22,
  n.blocks = 200,
  ldsc.log = file.path(out, "multivariable"),
  stand = TRUE,
  select = FALSE,
  chisq.max = NA
)
dimnames(result$S) <- list(trait_names, trait_names)
dimnames(result$S_Stand) <- list(trait_names, trait_names)
dimnames(result$I) <- list(trait_names, trait_names)

saveRDS(result, file.path(out, "genomicsem_ldsc.rds"))

write.table(
  result$S,
  file.path(out, "genetic_covariance_S.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  result$S_Stand,
  file.path(out, "genetic_correlation_S_Stand.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  result$I,
  file.path(out, "intercept_I.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  result$V,
  file.path(out, "sampling_covariance_V.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
write.table(
  result$V_Stand,
  file.path(out, "sampling_covariance_V_Stand.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

parameter_index <- data.frame()
counter <- 0L
for (column in seq_along(trait_names)) {
  for (row in column:length(trait_names)) {
    counter <- counter + 1L
    parameter_index <- rbind(
      parameter_index,
      data.frame(
        parameter_index = counter,
        trait1 = trait_names[column],
        trait2 = trait_names[row],
        is_variance = as.integer(column == row),
        estimate_covariance = result$S[row, column],
        estimate_standardized = result$S_Stand[row, column],
        se_covariance = sqrt(result$V[counter, counter]),
        se_standardized = sqrt(result$V_Stand[counter, counter])
      )
    )
  }
}
stopifnot(nrow(parameter_index) == nrow(result$V))
write.table(
  parameter_index,
  file.path(out, "parameter_index.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out, "sessionInfo.txt")
)
description <- packageDescription("GenomicSEM")
writeLines(
  c(
    paste0("GenomicSEM_version=", description$Version),
    paste0("GenomicSEM_RemoteSha=", description$RemoteSha),
    paste0("n_blocks=200"),
    paste0("stand=TRUE"),
    paste0("sep_weights=TRUE")
  ),
  file.path(out, "method_versions.txt")
)
