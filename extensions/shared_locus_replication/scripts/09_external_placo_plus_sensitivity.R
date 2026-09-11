#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

root <- "/root/IBD/20_Reproducibility_Ladder/extensions/20260817_shared_locus_replication"
source(file.path(root, "software/PLACO/PLACO_v0.2.0.R"))
stats <- fread(file.path(root, "results/raw/pair_harmonization_stats.tsv"))
ref <- stats[pair_id == "External_Howard_deLange"]
varz <- c(ref$var_z_dep, ref$var_z_ibd)
corz <- ref$cor_z
workers <- min(8L, parallel::detectCores(logical = FALSE))

input <- file.path(root, "results/raw/External_Howard_deLange.integration_candidates.tsv.gz")
dat <- fread(cmd = paste("zcat", shQuote(input)), showProgress = FALSE)
scores <- mclapply(seq_len(nrow(dat)), function(i) {
  placo.plus(c(dat$Z_DEP[i], dat$Z_IBD[i]), VarZ = varz, CorZ = corz)$p.placo.plus
}, mc.cores = workers, mc.preschedule = TRUE)
dat[, P_PLACO_PLUS := pmin(1, pmax(0, unlist(scores)))]
dat[, Q_BH_GENOMEWIDE := p.adjust(P_PLACO_PLUS, method = "BH", n = ref$retained)]
setorder(dat, P_PLACO_PLUS)
fwrite(dat, file.path(root, "results/raw/External_Howard_deLange.PLACOplus_sensitivity.tsv.gz"), sep = "\t", compress = "gzip")

annotation <- fread(file.path(root, "results/loci/External_Howard_deLange.variants_p1e6.tsv"),
                    select = c("SNP", "LOC", "MHC"))
gws <- dat[P_PLACO_PLUS < 5e-8]
gws <- annotation[gws, on = "SNP"]
gws_loci <- unique(gws[!is.na(LOC), LOC])

support <- fread(file.path(root, "results/loci/four_pair_lead_snp_support.tsv"))
external_support <- support[PAIR_ID == "External_Howard_deLange" & FOUND == TRUE]
external_support <- dat[external_support, on = c(SNP = "LEAD_SNP")]
external_support[, P_PLACO_PLUS_TARGETED := P_PLACO_PLUS]
fwrite(external_support, file.path(root, "reports/PRIMARY_LEADS_EXTERNAL_PLACO_PLUS.tsv"), sep = "\t", na = "")

rep <- fread(file.path(root, "results/loci/cross_cohort_replication.tsv"))
lookup <- setNames(external_support$P_PLACO_PLUS_TARGETED, external_support$SNP)
rep[, SENS_DISCOVERY_P := ifelse(discovery_pair == "External_Howard_deLange", lookup[LEAD_SNP], as.numeric(DISCOVERY_P_PLACO))]
rep[, SENS_VALIDATION_P := ifelse(validation_pair == "External_Howard_deLange", lookup[LEAD_SNP], as.numeric(VALIDATION_P_PLACO))]
rep[, SENS_REPLICATED := MHC != TRUE &
      is.finite(SENS_DISCOVERY_P) & SENS_DISCOVERY_P < 5e-8 &
      is.finite(SENS_VALIDATION_P) & SENS_VALIDATION_P < as.numeric(VALIDATION_BONF_THRESHOLD) &
      BOTH_VALIDATION_MARGINAL_P_LT_0.05 == TRUE & DIRECTION_CONCORDANT == TRUE]
fwrite(rep, file.path(root, "reports/REPLICATION_WITH_EXTERNAL_PLACO_PLUS.tsv"), sep = "\t", na = "")

summary <- data.table(
  metric = c("external_original_gws_variants", "external_placo_plus_gws_variants",
             "external_placo_plus_gws_loci", "frozen_replicated_rows",
             "replicated_rows_with_external_placo_plus", "replicated_blocks_with_external_placo_plus"),
  value = c(sum(fread(file.path(root, "results/raw/External_Howard_deLange.PLACO.tsv.gz"))$P_PLACO < 5e-8),
            nrow(gws), length(gws_loci), sum(rep$PRIMARY_REPLICATED == TRUE),
            sum(rep$SENS_REPLICATED == TRUE), uniqueN(rep[SENS_REPLICATED == TRUE, LOC]))
)
fwrite(summary, file.path(root, "reports/PLACO_PLUS_SENSITIVITY_SUMMARY.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(root, "provenance/R_session_external_placo_plus.txt"))
