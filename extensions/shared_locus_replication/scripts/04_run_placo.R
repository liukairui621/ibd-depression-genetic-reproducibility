#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication"
source(file.path(root, "software/PLACO/PLACO_v0.2.0.R"))
stats <- fread(file.path(root, "results/raw/pair_harmonization_stats.tsv"))
workers <- min(8L, parallel::detectCores(logical = FALSE))

score_variant <- function(z_dep, z_ibd, method, varz, corz) {
  tryCatch({
    if (method == "placo_plus") {
      ans <- placo.plus(c(z_dep, z_ibd), VarZ = varz, CorZ = corz)
      c(stat = ans$T.placo.plus, p = ans$p.placo.plus, error = 0)
    } else {
      ans <- placo(c(z_dep, z_ibd), VarZ = varz)
      c(stat = ans$T.placo, p = ans$p.placo, error = 0)
    }
  }, error = function(e) {
    c(stat = z_dep * z_ibd, p = NA_real_, error = 1)
  })
}

summaries <- vector("list", nrow(stats))
for (i in seq_len(nrow(stats))) {
  pair_id <- stats$pair_id[i]
  method <- stats$method[i]
  varz <- c(stats$var_z_dep[i], stats$var_z_ibd[i])
  corz <- stats$cor_z[i]
  infile <- file.path(root, "results/raw", paste0(pair_id, ".integration_candidates.tsv.gz"))
  message("Scoring ", pair_id, " with ", method)
  dat <- fread(cmd = paste("zcat", shQuote(infile)), showProgress = FALSE)

  scored <- mclapply(seq_len(nrow(dat)), function(j) {
    score_variant(dat$Z_DEP[j], dat$Z_IBD[j], method, varz, corz)
  }, mc.cores = workers, mc.preschedule = TRUE)
  scored <- do.call(rbind, scored)
  dat[, `:=`(
    T_PLACO = scored[, "stat"],
    P_PLACO = pmin(1, pmax(0, scored[, "p"])),
    INTEGRATION_ERROR = as.integer(scored[, "error"])
  )]
  n_error <- sum(dat$INTEGRATION_ERROR != 0 | !is.finite(dat$P_PLACO))
  if (n_error > 0) warning(pair_id, " had ", n_error, " integration errors")
  dat <- dat[INTEGRATION_ERROR == 0 & is.finite(P_PLACO) & P_PLACO <= 1.000001e-4]
  setorder(dat, P_PLACO)
  dat[, Q_BH_GENOMEWIDE := p.adjust(P_PLACO, method = "BH", n = stats$retained[i])]
  dat[, `:=`(
    GWS = P_PLACO < 5e-8,
    EXPLORATORY_P1E6 = P_PLACO < 1e-6,
    METHOD = method,
    VAR_Z_DEP = varz[1],
    VAR_Z_IBD = varz[2],
    COR_Z = corz
  )]
  outfile <- file.path(root, "results/raw", paste0(pair_id, ".PLACO.tsv.gz"))
  fwrite(dat, outfile, sep = "\t", compress = "gzip")
  summaries[[i]] <- data.frame(
    pair_id = pair_id, method = method, workers = workers,
    exact_p_le_1e4 = nrow(dat), integration_errors = n_error,
    p_lt_1e6 = sum(dat$P_PLACO < 1e-6),
    p_lt_5e8 = sum(dat$P_PLACO < 5e-8),
    minimum_p = if (nrow(dat)) min(dat$P_PLACO) else NA_real_
  )
}

fwrite(rbindlist(summaries), file.path(root, "results/raw/placo_run_summary.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(root, "provenance/R_session_placo.txt"))
