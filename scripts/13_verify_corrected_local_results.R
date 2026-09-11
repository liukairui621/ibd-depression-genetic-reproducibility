#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
repo <- normalizePath(file.path(dirname(script_arg), ".."), mustWork = TRUE)
local_dir <- file.path(repo, "results", "derived", "local")
infile <- file.path(local_dir, "bivar_all_primary.tsv")
samefile <- file.path(local_dir, "same_block_eligible_both.tsv")
outfile <- file.path(local_dir, "local_independent_verification.tsv")

x <- read.delim(infile, check.names = FALSE)
stopifnot(nrow(x) == 110L)
x$q_recomputed_R <- p.adjust(x$p, method = "BH")
max_delta <- max(abs(x$q_recomputed_R - x$q_global_primary), na.rm = TRUE)
stopifnot(max_delta < 1e-12)

same <- read.delim(samefile, check.names = FALSE)
stopifnot(nrow(same) == 3L)

out <- data.frame(
  check = c(
    "eligible_bivariate_tests",
    "maximum_absolute_BH_q_difference",
    "q_less_than_0.05_count",
    "blocks_estimable_in_both_pairs",
    "replicated_same_direction_q_less_than_0.05"
  ),
  observed = c(
    nrow(x),
    format(max_delta, scientific = TRUE, digits = 6),
    sum(x$q_recomputed_R < 0.05, na.rm = TRUE),
    nrow(same),
    sum(same$replicated_same_direction_q05 %in% TRUE, na.rm = TRUE)
  ),
  expected = c(110, "<1e-12", 4, 3, 0),
  pass = c(
    nrow(x) == 110L,
    max_delta < 1e-12,
    sum(x$q_recomputed_R < 0.05, na.rm = TRUE) == 4L,
    nrow(same) == 3L,
    sum(same$replicated_same_direction_q05 %in% TRUE, na.rm = TRUE) == 0L
  )
)

write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(dirname(outfile), "local_independent_verification_sessionInfo.txt"))
cat("Independent local verification passed; output:", outfile, "\n")
