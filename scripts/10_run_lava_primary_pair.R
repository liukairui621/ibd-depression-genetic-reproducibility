#!/usr/bin/env Rscript

set.seed(42)
suppressPackageStartupMessages({
  library(LAVA)
  library(data.table)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 10_run_lava_primary_pair.R <pair_id>")
}

pair_id <- args[[1]]
root <- "/root/IBD/20_Reproducibility_Ladder"
sumstats_dir <- file.path(root, "data", "munged_core")
out_dir <- file.path(root, "results", "local", pair_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pair_specs <- list(
  HowardMDD__FinnGenIBD = data.frame(
    phenotype = c("MDD_Howard2019", "IBD_FinnGen_R12"),
    cases = c(170756, 10960),
    controls = c(329443, 489388),
    filename = file.path(
      sumstats_dir,
      c("MDD_Howard2019.sumstats.gz", "IBD_FinnGen_R12.sumstats.gz")
    )
  ),
  FinnGenDEP__deLangeIBD = data.frame(
    phenotype = c("DEP_FinnGen_R12", "IBD_deLange2017"),
    cases = c(59333, 25042),
    controls = c(434831, 34915),
    filename = file.path(
      sumstats_dir,
      c("DEP_FinnGen_R12.sumstats.gz", "IBD_deLange2017.sumstats.gz")
    )
  )
)

if (!pair_id %in% names(pair_specs)) {
  stop(sprintf("Unknown pair_id: %s", pair_id))
}

info <- pair_specs[[pair_id]]
phenos <- info$phenotype
info_path <- file.path(out_dir, "lava_input.info")
write.table(info, info_path, sep = "\t", row.names = FALSE, quote = FALSE)

overlap <- diag(2)
rownames(overlap) <- phenos
colnames(overlap) <- phenos
overlap_path <- file.path(out_dir, "sample_overlap.txt")
write.table(overlap, overlap_path, sep = "\t", quote = FALSE)

plink_ref <- paste0(
  "/root/IBD/00_RawData/Reference/LD_EUR/",
  "1000G_EUR_Phase3_plink/1000G.EUR.QC.ALL"
)
blocks_path <- "/root/IBD/00_RawData/Reference/LAVA/LAVA_blocks_fixed.txt"
blocks <- read.loci(blocks_path)
n_primary_traits <- 4L
univ_threshold <- 0.05 / (nrow(blocks) * n_primary_traits)
workers <- as.integer(Sys.getenv("LAVA_WORKERS", "6"))

cat(sprintf("Pair: %s\n", pair_id))
cat(sprintf("Workers: %d\n", workers))
cat(sprintf("Loci: %d\n", nrow(blocks)))
cat(sprintf("Univariate threshold: %.12g\n", univ_threshold))

inp <- process.input(
  input.info.file = info_path,
  sample.overlap.file = overlap_path,
  ref.prefix = plink_ref,
  phenos = phenos
)

run_univ_locus <- function(i) {
  loc <- tryCatch(
    process.locus(locus = blocks[i, ], input = inp, phenos = phenos),
    error = function(e) NULL
  )
  if (is.null(loc)) return(NULL)

  result <- tryCatch(run.univ(loc), error = function(e) NULL)
  if (is.null(result)) return(NULL)

  result <- as.data.frame(result)
  result$locus <- blocks$LOC[i]
  result$chr <- blocks$CHR[i]
  result$start <- blocks$START[i]
  result$stop <- blocks$STOP[i]
  result
}

univ_list <- mclapply(
  seq_len(nrow(blocks)),
  run_univ_locus,
  mc.cores = workers,
  mc.preschedule = TRUE,
  mc.set.seed = TRUE
)
univ <- rbindlist(univ_list, fill = TRUE)
univ$pair_id <- pair_id
univ$univ_threshold <- univ_threshold
fwrite(univ, file.path(out_dir, "univ_all.tsv"), sep = "\t", na = "NA")

qualified <- univ[!is.na(p) & p < univ_threshold]
loci_1 <- unique(qualified[phen == phenos[[1]], locus])
loci_2 <- unique(qualified[phen == phenos[[2]], locus])
eligible_loci <- intersect(loci_1, loci_2)

cat(sprintf(
  "Qualified loci: %s=%d; %s=%d; eligible=%d\n",
  phenos[[1]], length(loci_1), phenos[[2]], length(loci_2),
  length(eligible_loci)
))

run_bivar_locus <- function(locus_id) {
  block <- blocks[blocks$LOC == locus_id, ]
  loc <- tryCatch(
    process.locus(locus = block, input = inp, phenos = phenos),
    error = function(e) NULL
  )
  if (is.null(loc)) return(NULL)

  result <- tryCatch(
    run.bivar(loc, phenos = phenos),
    error = function(e) NULL
  )
  if (is.null(result)) return(NULL)

  result <- as.data.frame(result)
  result$locus <- locus_id
  result$chr <- block$CHR
  result$start <- block$START
  result$stop <- block$STOP
  result$pair_id <- pair_id
  result
}

if (length(eligible_loci) > 0L) {
  bivar_list <- mclapply(
    eligible_loci,
    run_bivar_locus,
    mc.cores = workers,
    mc.preschedule = TRUE,
    mc.set.seed = TRUE
  )
  bivar <- rbindlist(bivar_list, fill = TRUE)
} else {
  bivar <- data.table()
}
fwrite(bivar, file.path(out_dir, "bivar_eligible.tsv"), sep = "\t", na = "NA")

writeLines(
  c(
    sprintf("pair_id=%s", pair_id),
    sprintf("n_blocks=%d", nrow(blocks)),
    sprintf("n_primary_traits=%d", n_primary_traits),
    sprintf("univ_threshold=%.12g", univ_threshold),
    sprintf("n_univ_rows=%d", nrow(univ)),
    sprintf("n_eligible_bivar=%d", nrow(bivar)),
    sprintf("workers=%d", workers),
    sprintf("seed=%d", 42)
  ),
  file.path(out_dir, "run_summary.txt")
)

sink(file.path(out_dir, "sessionInfo.txt"))
sessionInfo()
sink()
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(out_dir, "DONE"))
