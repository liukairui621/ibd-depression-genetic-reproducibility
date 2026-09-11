#!/usr/bin/env Rscript

base <- "/root/IBD/20_Reproducibility_Ladder"
genomicsem <- file.path(base, "results", "global_extension", "genomicsem")
out <- file.path(base, "results", "global_core")

result <- readRDS(file.path(genomicsem, "genomicsem_ldsc.rds"))
index <- read.delim(file.path(genomicsem, "parameter_index.tsv"))
rg <- read.delim(file.path(out, "rg_summary.tsv"))

find_parameter <- function(trait_a, trait_b) {
  hit <- which(
    (index$trait1 == trait_a & index$trait2 == trait_b) |
      (index$trait1 == trait_b & index$trait2 == trait_a)
  )
  stopifnot(length(hit) == 1L)
  index$parameter_index[hit]
}

i1 <- find_parameter("MDD_Howard2019", "IBD_FinnGen_R12")
i2 <- find_parameter("DEP_FinnGen_R12", "IBD_deLange2017")
covariance <- result$V_Stand[i1, i2]
sampling_correlation <- covariance / sqrt(
  result$V_Stand[i1, i1] * result$V_Stand[i2, i2]
)

primary <- rg[
  (rg$trait1 == "MDD_Howard2019" & rg$trait2 == "IBD_FinnGen_R12") |
    (rg$trait1 == "DEP_FinnGen_R12" & rg$trait2 == "IBD_deLange2017"),
]
stopifnot(nrow(primary) == 2L)

diagnostic <- data.frame(
  pair1_parameter_index = i1,
  pair2_parameter_index = i2,
  standardized_sampling_covariance = covariance,
  standardized_sampling_correlation = sampling_correlation,
  independence_assumption_supported = FALSE,
  synthesis = "report_reciprocal_estimates_separately_no_fixed_effect_pooling"
)
write.table(
  diagnostic,
  file.path(out, "primary_pair_dependence_diagnostic.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

primary$synthesis <- "reciprocal_pair_reported_separately"
write.table(
  primary,
  file.path(out, "primary_pair_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(sprintf("sampling covariance = %.12g\n", covariance))
cat(sprintf("sampling correlation = %.12g\n", sampling_correlation))
