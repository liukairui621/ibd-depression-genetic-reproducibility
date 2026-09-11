#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(42)

root <- "/root/IBD/20_Reproducibility_Ladder"
global <- file.path(root, "results", "global_extension")
genomicsem <- file.path(global, "genomicsem")
out <- file.path(global, "attribution")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

result <- readRDS(file.path(genomicsem, "genomicsem_ldsc.rds"))
index <- read.delim(file.path(genomicsem, "parameter_index.tsv"))
pairwise <- read.delim(file.path(global, "rg_cross_disease.tsv"))

ukb_definitions <- c(
  "UKB_LifetimeMDD",
  "UKB_MDDRecur",
  "UKB_GPpsy",
  "UKB_ICD10Dep"
)
ibd_traits <- c(
  "IBD_deLange2017",
  "CD_deLange2017",
  "UC_deLange2017",
  "IBD_FinnGen_R12",
  "CD_FinnGen_R12",
  "UC_FinnGen_R12"
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
  grepl("deLange", cells$ibd_trait),
  "deLange2017",
  "FinnGen_R12"
)
cells$ibd_subtype <- sub("_.*", "", cells$ibd_trait)
cells$parameter_index <- mapply(
  find_parameter,
  cells$depression_definition,
  cells$ibd_trait
)
cells$rg <- mapply(
  function(a, b) result$S_Stand[b, a],
  match(cells$depression_definition, colnames(result$S_Stand)),
  match(cells$ibd_trait, rownames(result$S_Stand))
)
cells$depression_definition <- factor(
  cells$depression_definition,
  levels = ukb_definitions
)
cells$ibd_subtype <- factor(
  cells$ibd_subtype,
  levels = c("IBD", "CD", "UC")
)
cells$ibd_cohort <- factor(
  cells$ibd_cohort,
  levels = c("deLange2017", "FinnGen_R12")
)

V <- result$V_Stand[cells$parameter_index, cells$parameter_index, drop = FALSE]
V <- (V + t(V)) / 2
eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
used_near_pd <- FALSE
if (min(eigenvalues) <= 1e-10 || rcond(V) < 1e-10) {
  V <- as.matrix(Matrix::nearPD(V, corr = FALSE, keepDiag = TRUE)$mat)
  used_near_pd <- TRUE
}
W <- solve(V)

fit_gls <- function(formula, data, weights) {
  X <- model.matrix(formula, data)
  information <- t(X) %*% weights %*% X
  covariance <- solve(information)
  beta <- covariance %*% t(X) %*% weights %*% data$rg
  fitted <- as.vector(X %*% beta)
  residual <- data$rg - fitted
  q <- as.numeric(t(residual) %*% weights %*% residual)
  list(
    formula = formula,
    X = X,
    beta = as.vector(beta),
    covariance = covariance,
    fitted = fitted,
    residual = residual,
    q = q,
    df = nrow(X) - ncol(X)
  )
}

additive <- fit_gls(
  rg ~ depression_definition + ibd_subtype + ibd_cohort,
  cells,
  W
)
dep_interaction <- fit_gls(
  rg ~ depression_definition * ibd_subtype + ibd_cohort,
  cells,
  W
)
cohort_interaction <- fit_gls(
  rg ~ depression_definition + ibd_subtype * ibd_cohort,
  cells,
  W
)

nested_test <- function(reduced, full, label) {
  statistic <- max(0, reduced$q - full$q)
  df <- length(full$beta) - length(reduced$beta)
  data.frame(
    test = label,
    statistic = statistic,
    df = df,
    p = pchisq(statistic, df = df, lower.tail = FALSE)
  )
}

interaction_tests <- rbind(
  nested_test(
    additive,
    dep_interaction,
    "depression_definition_x_ibd_subtype"
  ),
  nested_test(
    additive,
    cohort_interaction,
    "ibd_subtype_x_ibd_cohort"
  )
)

wald_block <- function(fit, columns, label) {
  positions <- which(colnames(fit$X) %in% columns)
  beta <- fit$beta[positions]
  covariance <- fit$covariance[positions, positions, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% solve(covariance) %*% beta)
  data.frame(
    test = label,
    statistic = statistic,
    df = length(positions),
    p = pchisq(statistic, df = length(positions), lower.tail = FALSE)
  )
}

coefficient_names <- colnames(additive$X)
omnibus_tests <- rbind(
  wald_block(
    additive,
    grep("^depression_definition", coefficient_names, value = TRUE),
    "depression_definition"
  ),
  wald_block(
    additive,
    grep("^ibd_subtype", coefficient_names, value = TRUE),
    "ibd_subtype"
  ),
  wald_block(
    additive,
    grep("^ibd_cohort", coefficient_names, value = TRUE),
    "ibd_cohort"
  )
)

reference_grid <- expand.grid(
  depression_definition = factor(
    ukb_definitions,
    levels = levels(cells$depression_definition)
  ),
  ibd_subtype = factor(
    c("IBD", "CD", "UC"),
    levels = levels(cells$ibd_subtype)
  ),
  ibd_cohort = factor(
    c("deLange2017", "FinnGen_R12"),
    levels = levels(cells$ibd_cohort)
  ),
  KEEP.OUT.ATTRS = FALSE
)
X_grid <- model.matrix(
  delete.response(terms(additive$formula)),
  reference_grid
)
marginal_rows <- list()
for (definition in ukb_definitions) {
  selector <- reference_grid$depression_definition == definition
  contrast <- colMeans(X_grid[selector, , drop = FALSE])
  estimate <- sum(contrast * additive$beta)
  se <- sqrt(as.numeric(t(contrast) %*% additive$covariance %*% contrast))
  marginal_rows[[length(marginal_rows) + 1L]] <- data.frame(
    axis = "depression_definition",
    level = definition,
    estimate = estimate,
    se = se,
    ci_low = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se
  )
}
for (subtype in c("IBD", "CD", "UC")) {
  selector <- reference_grid$ibd_subtype == subtype
  contrast <- colMeans(X_grid[selector, , drop = FALSE])
  estimate <- sum(contrast * additive$beta)
  se <- sqrt(as.numeric(t(contrast) %*% additive$covariance %*% contrast))
  marginal_rows[[length(marginal_rows) + 1L]] <- data.frame(
    axis = "ibd_subtype",
    level = subtype,
    estimate = estimate,
    se = se,
    ci_low = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se
  )
}
marginal_means <- do.call(rbind, marginal_rows)

draws <- MASS::mvrnorm(
  n = 50000,
  mu = additive$beta,
  Sigma = additive$covariance
)
depression_contrasts <- lapply(
  ukb_definitions,
  function(definition) {
    colMeans(
      X_grid[
        reference_grid$depression_definition == definition,
        ,
        drop = FALSE
      ]
    )
  }
)
subtype_contrasts <- lapply(
  c("IBD", "CD", "UC"),
  function(subtype) {
    colMeans(
      X_grid[reference_grid$ibd_subtype == subtype, , drop = FALSE]
    )
  }
)
depression_draws <- sapply(
  depression_contrasts,
  function(contrast) draws %*% contrast
)
subtype_draws <- sapply(
  subtype_contrasts,
  function(contrast) draws %*% contrast
)
depression_sd <- apply(depression_draws, 1, sd)
subtype_sd <- apply(subtype_draws, 1, sd)
dispersion_ratio <- depression_sd / subtype_sd

observed_dep_sd <- sd(
  marginal_means$estimate[
    marginal_means$axis == "depression_definition"
  ]
)
observed_subtype_sd <- sd(
  marginal_means$estimate[marginal_means$axis == "ibd_subtype"]
)
ratio_summary <- data.frame(
  depression_definition_sd = observed_dep_sd,
  ibd_subtype_sd = observed_subtype_sd,
  ratio = observed_dep_sd / observed_subtype_sd,
  ratio_ci_low = unname(quantile(dispersion_ratio, 0.025)),
  ratio_ci_high = unname(quantile(dispersion_ratio, 0.975)),
  probability_ratio_gt_1 = mean(dispersion_ratio > 1),
  simulation_draws = 50000,
  seed = 42
)

empirical_rows <- list()
for (ibd_trait in ibd_traits) {
  values <- cells$rg[cells$ibd_trait == ibd_trait]
  empirical_rows[[length(empirical_rows) + 1L]] <- data.frame(
    axis = "depression_definition",
    stratum = ibd_trait,
    n_levels = length(values),
    sd = sd(values),
    mad = mad(values, constant = 1),
    range = diff(range(values))
  )
}
for (definition in ukb_definitions) {
  for (cohort in c("deLange2017", "FinnGen_R12")) {
    values <- cells$rg[
      cells$depression_definition == definition &
        cells$ibd_cohort == cohort
    ]
    empirical_rows[[length(empirical_rows) + 1L]] <- data.frame(
      axis = "ibd_subtype",
      stratum = paste(definition, cohort, sep = "__"),
      n_levels = length(values),
      sd = sd(values),
      mad = mad(values, constant = 1),
      range = diff(range(values))
    )
  }
}
empirical_dispersion <- do.call(rbind, empirical_rows)

leave_one_definition <- lapply(
  ukb_definitions,
  function(excluded) {
    retained <- setdiff(ukb_definitions, excluded)
    dep_values <- marginal_means$estimate[
      marginal_means$axis == "depression_definition" &
        marginal_means$level %in% retained
    ]
    data.frame(
      excluded_definition = excluded,
      depression_definition_sd = sd(dep_values),
      ibd_subtype_sd = observed_subtype_sd,
      ratio = sd(dep_values) / observed_subtype_sd
    )
  }
)
leave_one_definition <- do.call(rbind, leave_one_definition)

leave_one_subtype <- lapply(
  c("IBD", "CD", "UC"),
  function(excluded) {
    retained <- setdiff(c("IBD", "CD", "UC"), excluded)
    subtype_values <- marginal_means$estimate[
      marginal_means$axis == "ibd_subtype" &
        marginal_means$level %in% retained
    ]
    data.frame(
      excluded_subtype = excluded,
      depression_definition_sd = observed_dep_sd,
      ibd_subtype_sd = sd(subtype_values),
      ratio = observed_dep_sd / sd(subtype_values)
    )
  }
)
leave_one_subtype <- do.call(rbind, leave_one_subtype)

contrast_rows <- list()
append_contrast <- function(axis, stratum, level1, level2, parameter1, parameter2) {
  difference <- index$estimate_standardized[
    index$parameter_index == parameter1
  ] - index$estimate_standardized[index$parameter_index == parameter2]
  variance <- result$V_Stand[parameter1, parameter1] +
    result$V_Stand[parameter2, parameter2] -
    2 * result$V_Stand[parameter1, parameter2]
  se <- sqrt(max(variance, 0))
  z <- difference / se
  contrast_rows[[length(contrast_rows) + 1L]] <<- data.frame(
    axis = axis,
    stratum = stratum,
    level1 = level1,
    level2 = level2,
    difference = difference,
    se = se,
    z = z,
    p = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}

for (ibd_trait in ibd_traits) {
  pairs <- combn(ukb_definitions, 2, simplify = FALSE)
  for (pair in pairs) {
    append_contrast(
      "depression_definition",
      ibd_trait,
      pair[1],
      pair[2],
      find_parameter(pair[1], ibd_trait),
      find_parameter(pair[2], ibd_trait)
    )
  }
}
for (definition in c("MDD_Howard2019", "DEP_FinnGen_R12", ukb_definitions)) {
  for (cohort in c("deLange2017", "FinnGen_R12")) {
    suffix <- if (cohort == "deLange2017") "deLange2017" else "FinnGen_R12"
    cd_trait <- paste0("CD_", suffix)
    uc_trait <- paste0("UC_", suffix)
    append_contrast(
      "ibd_subtype_CD_vs_UC",
      paste(definition, cohort, sep = "__"),
      "CD",
      "UC",
      find_parameter(definition, cd_trait),
      find_parameter(definition, uc_trait)
    )
  }
}
contrasts <- do.call(rbind, contrast_rows)
contrasts$q_within_axis <- ave(
  contrasts$p,
  contrasts$axis,
  FUN = function(values) p.adjust(values, method = "BH")
)

ukb_cd_uc <- contrasts[
  contrasts$axis == "ibd_subtype_CD_vs_UC" &
    grepl("^UKB_", contrasts$stratum),
]
cd_uc_direction_by_cohort <- do.call(
  rbind,
  lapply(
    c("deLange2017", "FinnGen_R12"),
    function(cohort) {
      values <- ukb_cd_uc$difference[
        grepl(paste0("__", cohort, "$"), ukb_cd_uc$stratum)
      ]
      n_positive <- sum(values > 0)
      n_negative <- sum(values < 0)
      modal_direction <- if (
        n_positive >= 3
      ) {
        "CD_gt_UC"
      } else if (
        n_negative >= 3
      ) {
        "CD_lt_UC"
      } else {
        "inconsistent"
      }
      data.frame(
        ibd_cohort = cohort,
        n_definitions = length(values),
        n_cd_gt_uc = n_positive,
        n_cd_lt_uc = n_negative,
        modal_direction = modal_direction,
        passes_three_of_four = modal_direction != "inconsistent"
      )
    }
  )
)
cd_uc_direction_consistent <- (
  all(cd_uc_direction_by_cohort$passes_three_of_four) &&
    length(unique(cd_uc_direction_by_cohort$modal_direction)) == 1
)

interaction_dep_p <- interaction_tests$p[
  interaction_tests$test == "depression_definition_x_ibd_subtype"
]
interaction_cohort_p <- interaction_tests$p[
  interaction_tests$test == "ibd_subtype_x_ibd_cohort"
]
leave_gp_ratio <- leave_one_definition$ratio[
  leave_one_definition$excluded_definition == "UKB_GPpsy"
]
leave_ibd_ratio <- leave_one_subtype$ratio[
  leave_one_subtype$excluded_subtype == "IBD"
]

classification <- "mixed_or_indeterminate"
reason <- paste(
  "The 95% interval for the adjusted dispersion ratio includes 1,",
  "or interaction/leave-one-out sensitivity prevents a stable attribution."
)
if (
  ratio_summary$ratio_ci_low > 1 &&
    leave_gp_ratio > 1 &&
    interaction_dep_p >= 0.05
) {
  classification <- "depression_definition_dominant"
  reason <- paste(
    "Adjusted depression-definition dispersion exceeded subtype dispersion,",
    "the 95% interval excluded 1, and the result persisted after excluding",
    "the broad GP consultation phenotype."
  )
} else if (
  ratio_summary$ratio_ci_high < 1 &&
    interaction_cohort_p >= 0.05 &&
    cd_uc_direction_consistent &&
    leave_ibd_ratio < 1
) {
  classification <- "ibd_subtype_dominant"
  reason <- paste(
    "Adjusted IBD-subtype dispersion exceeded depression-definition",
    "dispersion, the 95% interval excluded 1, and subtype effects did not",
    "materially differ by IBD cohort. CD-versus-UC direction was consistent",
    "across at least three of four UKB definitions in both IBD cohorts, and",
    "the result persisted when the combined IBD phenotype was excluded."
  )
}

decision <- data.frame(
  classification = classification,
  reason = reason,
  near_pd_used = used_near_pd,
  min_original_v_eigenvalue = min(eigenvalues),
  depression_x_subtype_p = interaction_dep_p,
  subtype_x_cohort_p = interaction_cohort_p,
  leave_gp_ratio = leave_gp_ratio,
  leave_ibd_ratio = leave_ibd_ratio,
  cd_uc_direction_consistent = cd_uc_direction_consistent,
  cd_uc_deLange_direction = cd_uc_direction_by_cohort$modal_direction[
    cd_uc_direction_by_cohort$ibd_cohort == "deLange2017"
  ],
  cd_uc_FinnGen_direction = cd_uc_direction_by_cohort$modal_direction[
    cd_uc_direction_by_cohort$ibd_cohort == "FinnGen_R12"
  ]
)

write.table(
  cells,
  file.path(out, "ukb_crossed_cells.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  marginal_means,
  file.path(out, "adjusted_marginal_rg.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  ratio_summary,
  file.path(out, "dispersion_ratio.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  empirical_dispersion,
  file.path(out, "empirical_dispersion.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  leave_one_definition,
  file.path(out, "leave_one_definition.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  leave_one_subtype,
  file.path(out, "leave_one_subtype.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  contrasts,
  file.path(out, "sampling_covariance_aware_contrasts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  cd_uc_direction_by_cohort,
  file.path(out, "cd_uc_direction_by_cohort.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  omnibus_tests,
  file.path(out, "omnibus_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  interaction_tests,
  file.path(out, "interaction_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  decision,
  file.path(out, "instability_attribution_decision.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(
  list(
    additive = additive,
    dep_interaction = dep_interaction,
    cohort_interaction = cohort_interaction,
    dispersion_ratio_draws = dispersion_ratio,
    decision = decision
  ),
  file.path(out, "attribution_models.rds")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out, "sessionInfo.txt")
)
