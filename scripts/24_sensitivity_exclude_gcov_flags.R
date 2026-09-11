#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(42)
suppressPackageStartupMessages(library(Matrix))

base <- "/root/IBD/20_Reproducibility_Ladder"
global <- file.path(base, "results", "global_extension")
genomicsem <- file.path(global, "genomicsem")
out <- file.path(global, "sensitivity_exclude_gcov_flags")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

result <- readRDS(file.path(genomicsem, "genomicsem_ldsc.rds"))
index <- read.delim(file.path(genomicsem, "parameter_index.tsv"))
pairwise <- read.delim(file.path(global, "rg_cross_disease.tsv"))

ukb_definitions <- c(
  "UKB_LifetimeMDD", "UKB_MDDRecur", "UKB_GPpsy", "UKB_ICD10Dep"
)
ibd_traits <- c(
  "IBD_deLange2017", "CD_deLange2017", "UC_deLange2017",
  "IBD_FinnGen_R12", "CD_FinnGen_R12", "UC_FinnGen_R12"
)

find_parameter <- function(trait_a, trait_b) {
  hit <- which(
    (index$trait1 == trait_a & index$trait2 == trait_b) |
      (index$trait1 == trait_b & index$trait2 == trait_a)
  )
  stopifnot(length(hit) == 1L)
  index$parameter_index[hit]
}

cells <- expand.grid(
  depression_definition = ukb_definitions,
  ibd_trait = ibd_traits,
  KEEP.OUT.ATTRS = FALSE
)
cells$ibd_cohort <- ifelse(
  grepl("deLange", cells$ibd_trait), "deLange2017", "FinnGen_R12"
)
cells$ibd_subtype <- sub("_.*", "", cells$ibd_trait)
cells$parameter_index <- mapply(
  find_parameter, cells$depression_definition, cells$ibd_trait
)
cells$rg_genomicsem <- mapply(
  function(a, b) result$S_Stand[b, a],
  match(cells$depression_definition, colnames(result$S_Stand)),
  match(cells$ibd_trait, rownames(result$S_Stand))
)

pair_key <- paste(pairwise$trait1, pairwise$trait2, sep = "__")
cell_key <- paste(cells$depression_definition, cells$ibd_trait, sep = "__")
match_index <- match(cell_key, pair_key)
stopifnot(!anyNA(match_index))
cells$rg_pairwise <- pairwise$rg[match_index]
cells$se_pairwise <- pairwise$se[match_index]
cells$p_pairwise <- pairwise$p[match_index]
cells$gcov_intercept <- pairwise$gcov_intercept[match_index]
cells$gcov_intercept_se <- pairwise$gcov_intercept_se[match_index]
cells$gcov_intercept_z <- pairwise$gcov_intercept_z[match_index]
cells$gcov_intercept_flag <- pairwise$gcov_intercept_flag[match_index]

flagged <- cells[cells$gcov_intercept_flag == 1, , drop = FALSE]
retained <- cells[cells$gcov_intercept_flag == 0, , drop = FALSE]
stopifnot(nrow(flagged) == 6L, nrow(retained) == 18L)
retained$fdr_sensitivity_18 <- p.adjust(retained$p_pairwise, method = "BH")

retained$depression_definition <- factor(
  retained$depression_definition, levels = ukb_definitions
)
retained$ibd_subtype <- factor(
  retained$ibd_subtype, levels = c("IBD", "CD", "UC")
)
retained$ibd_cohort <- factor(
  retained$ibd_cohort, levels = c("deLange2017", "FinnGen_R12")
)

V <- result$V_Stand[
  retained$parameter_index, retained$parameter_index, drop = FALSE
]
V <- (V + t(V)) / 2
eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
used_near_pd <- FALSE
if (min(eigenvalues) <= 1e-10 || rcond(V) < 1e-10) {
  V <- as.matrix(nearPD(V, corr = FALSE, keepDiag = TRUE)$mat)
  used_near_pd <- TRUE
}
W <- solve(V)

fit_gls <- function(formula, data, weights) {
  X <- model.matrix(formula, data)
  if (qr(X)$rank != ncol(X)) stop("Rank-deficient sensitivity model")
  information <- t(X) %*% weights %*% X
  covariance <- solve(information)
  beta <- as.vector(covariance %*% t(X) %*% weights %*% data$rg_genomicsem)
  residual <- data$rg_genomicsem - as.vector(X %*% beta)
  list(
    formula = formula, X = X, beta = beta, covariance = covariance,
    q = as.numeric(t(residual) %*% weights %*% residual)
  )
}

additive <- fit_gls(
  rg_genomicsem ~ depression_definition + ibd_subtype + ibd_cohort,
  retained, W
)
dep_interaction <- fit_gls(
  rg_genomicsem ~ depression_definition * ibd_subtype + ibd_cohort,
  retained, W
)
cohort_interaction <- fit_gls(
  rg_genomicsem ~ depression_definition + ibd_subtype * ibd_cohort,
  retained, W
)

nested_test <- function(reduced, full, label) {
  statistic <- max(0, reduced$q - full$q)
  df <- length(full$beta) - length(reduced$beta)
  data.frame(
    test = label, statistic = statistic, df = df,
    p = pchisq(statistic, df = df, lower.tail = FALSE)
  )
}
interaction_tests <- rbind(
  nested_test(additive, dep_interaction, "depression_definition_x_ibd_subtype"),
  nested_test(additive, cohort_interaction, "ibd_subtype_x_ibd_cohort")
)

wald_block <- function(fit, pattern, label) {
  positions <- grep(pattern, colnames(fit$X))
  beta <- fit$beta[positions]
  covariance <- fit$covariance[positions, positions, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% solve(covariance) %*% beta)
  data.frame(
    test = label, statistic = statistic, df = length(positions),
    p = pchisq(statistic, df = length(positions), lower.tail = FALSE)
  )
}
omnibus_tests <- rbind(
  wald_block(additive, "^depression_definition", "depression_definition"),
  wald_block(additive, "^ibd_subtype", "ibd_subtype"),
  wald_block(additive, "^ibd_cohort", "ibd_cohort")
)

reference_grid <- expand.grid(
  depression_definition = factor(
    ukb_definitions, levels = levels(retained$depression_definition)
  ),
  ibd_subtype = factor(c("IBD", "CD", "UC"), levels = levels(retained$ibd_subtype)),
  ibd_cohort = factor(
    c("deLange2017", "FinnGen_R12"), levels = levels(retained$ibd_cohort)
  ),
  KEEP.OUT.ATTRS = FALSE
)
X_grid <- model.matrix(delete.response(terms(additive$formula)), reference_grid)
marginal_rows <- list()
for (definition in ukb_definitions) {
  contrast <- colMeans(
    X_grid[reference_grid$depression_definition == definition, , drop = FALSE]
  )
  estimate <- sum(contrast * additive$beta)
  se <- sqrt(as.numeric(t(contrast) %*% additive$covariance %*% contrast))
  marginal_rows[[length(marginal_rows) + 1L]] <- data.frame(
    axis = "depression_definition", level = definition,
    estimate = estimate, se = se,
    ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se
  )
}
for (subtype in c("IBD", "CD", "UC")) {
  contrast <- colMeans(
    X_grid[reference_grid$ibd_subtype == subtype, , drop = FALSE]
  )
  estimate <- sum(contrast * additive$beta)
  se <- sqrt(as.numeric(t(contrast) %*% additive$covariance %*% contrast))
  marginal_rows[[length(marginal_rows) + 1L]] <- data.frame(
    axis = "ibd_subtype", level = subtype,
    estimate = estimate, se = se,
    ci_low = estimate - 1.96 * se, ci_high = estimate + 1.96 * se
  )
}
marginals <- do.call(rbind, marginal_rows)

summary_row <- data.frame(
  n_total = nrow(cells),
  n_flagged_removed = nrow(flagged),
  n_retained = nrow(retained),
  n_positive = sum(retained$rg_pairwise > 0),
  n_nominal_p_lt_0_05 = sum(retained$p_pairwise < 0.05),
  n_bh_fdr_lt_0_05 = sum(retained$fdr_sensitivity_18 < 0.05),
  rg_min = min(retained$rg_pairwise),
  rg_max = max(retained$rg_pairwise),
  near_pd_used = used_near_pd,
  min_original_v_eigenvalue = min(eigenvalues)
)

write.table(flagged, file.path(out, "excluded_flagged_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(retained, file.path(out, "retained_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_row, file.path(out, "sensitivity_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(marginals, file.path(out, "adjusted_marginal_rg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(omnibus_tests, file.path(out, "omnibus_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(interaction_tests, file.path(out, "interaction_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(additive = additive, dep_interaction = dep_interaction, cohort_interaction = cohort_interaction), file.path(out, "models.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

print(summary_row)
print(omnibus_tests)
print(interaction_tests)
