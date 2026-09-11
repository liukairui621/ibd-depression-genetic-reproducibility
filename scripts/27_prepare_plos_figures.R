#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("scripts/27_prepare_plos_figures.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
corrected_repo <- file.path(repo_root, "results", "derived")
out_dir <- Sys.getenv("PLOS_FIGURE_OUT", file.path(repo_root, "figures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

global_core <- file.path(corrected_repo, "global_core")
local_dir <- file.path(corrected_repo, "local")
global_ext <- file.path(corrected_repo, "global_extension")
attrib_dir <- file.path(global_ext, "attribution")

theme_set(
  theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      text = element_text(colour = "#202124"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, colour = "#4D5156"),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(7, 9, 7, 7)
    )
)

pair_colours <- c(
  "Howard public depression x FinnGen IBD" = "#007C83",
  "FinnGen depression x de Lange IBD" = "#C84B31"
)

# Figure 1: prespecified reproducibility ladder.
rg <- fread(file.path(global_core, "rg_summary.tsv"))
primary <- rg[
  (trait1 == "MDD_Howard2019" & trait2 == "IBD_FinnGen_R12") |
    (trait1 == "DEP_FinnGen_R12" & trait2 == "IBD_deLange2017")
]
primary[, label := fifelse(
  trait1 == "MDD_Howard2019",
  "Howard public depression x FinnGen IBD",
  "FinnGen depression x de Lange IBD"
)]
primary[, `:=`(lower = rg - 1.96 * se, upper = rg + 1.96 * se)]

forest_data <- primary[, .(
  label, estimate = rg, lower, upper, p, qc = gcov_intercept_z
)]
forest_data[, label := factor(
  label,
  levels = rev(c(
    "Howard public depression x FinnGen IBD",
    "FinnGen depression x de Lange IBD"
  ))
)]

panel_1a <- ggplot(forest_data, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, colour = "#9AA0A6", linewidth = 0.4) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.16, colour = "#5F6368", linewidth = 0.55
  ) +
  geom_point(
    aes(shape = is.na(qc)), size = 2.5, fill = "#FFFFFF",
    colour = "#202124", stroke = 0.8
  ) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 18), guide = "none") +
  labs(
    title = "Reciprocal genome-wide correlation",
    subtitle = "Points show rg with 95% confidence intervals",
    x = "Genetic correlation (rg)", y = NULL
  ) +
  coord_cartesian(xlim = c(0, 0.20))

local_sig <- fread(file.path(local_dir, "bivar_q05_primary.tsv"))
local_sig[, pair_label := fifelse(
  pair_id == "HowardMDD__FinnGenIBD",
  "Howard public depression x FinnGen IBD",
  "FinnGen depression x de Lange IBD"
)]
local_sig[, locus_label := sprintf(
  "Locus %s, chr%s:%0.1f-%0.1f Mb",
  locus, chr, start / 1e6, stop / 1e6
)]
local_sig[, locus_label := factor(locus_label, levels = rev(locus_label))]

panel_1b <- ggplot(local_sig, aes(x = rho, y = locus_label, colour = pair_label)) +
  geom_vline(xintercept = 0, colour = "#9AA0A6", linewidth = 0.4) +
  geom_errorbarh(
    aes(xmin = rho.lower, xmax = rho.upper),
    height = 0.16, linewidth = 0.65
  ) +
  geom_point(size = 2.7) +
  scale_colour_manual(values = pair_colours) +
  labs(
    title = "Pair-specific signals",
    subtitle = sprintf("%d q<0.05; no replication", nrow(local_sig)),
    x = "Local correlation (rho)", y = NULL
  ) +
  theme(legend.position = "none")

same <- fread(file.path(local_dir, "same_block_eligible_both.tsv"))
same_long <- rbind(
  same[, .(
    locus, chr, start, stop,
    pair_label = "Howard public depression x FinnGen IBD",
    rho = pair1_rho, lower = pair1_rho_lower,
    upper = pair1_rho_upper, q = pair1_q
  )],
  same[, .(
    locus, chr, start, stop,
    pair_label = "FinnGen depression x de Lange IBD",
    rho = pair2_rho, lower = pair2_rho_lower,
    upper = pair2_rho_upper, q = pair2_q
  )]
)
same_long[, locus_label := sprintf(
  "Locus %s, chr%s:%0.1f-%0.1f Mb",
  locus, chr, start / 1e6, stop / 1e6
)]
same_long[, locus_label := factor(locus_label, levels = rev(unique(locus_label)))]

if (nrow(same_long) > 0L) {
  panel_1c <- ggplot(
    same_long, aes(x = rho, y = locus_label, colour = pair_label)
  ) +
    geom_vline(xintercept = 0, colour = "#9AA0A6", linewidth = 0.4) +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper), height = 0.16, linewidth = 0.65
    ) +
    geom_point(size = 2.7) +
    scale_colour_manual(values = pair_colours) +
    labs(
      title = "Blocks estimable in both pairs",
      subtitle = sprintf("%d blocks; no same-direction q<0.05", nrow(same)),
      x = "Local genetic correlation (rho; 95% CI)", y = NULL
    )
} else {
  panel_1c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No block was eligible in both pairs") +
    labs(title = "Blocks eligible in both reciprocal pairs", x = NULL, y = NULL) +
    theme_void()
}

ladder <- data.table(
  layer = factor(
    c("Genome-wide LDSC", "Local LAVA", "Replicated PLACO", "Fine-mapping / QTL"),
    levels = rev(c("Genome-wide LDSC", "Local LAVA", "Replicated PLACO", "Fine-mapping / QTL"))
  ),
  status = c("Reproduced", "Not reproduced", "Reproduced", "Adjudicated"),
  detail = c(
    "2/2 reciprocal pairs\npositive (P<0.025)",
    "0 same-block\ncovariance components",
    "3 non-MHC\npleiotropic blocks",
    "None established\nacross four GWAS"
  )
)
status_colours <- c(
  "Reproduced" = "#2E7D32",
  "Not reproduced" = "#B3261E",
  "Adjudicated" = "#6A5ACD"
)
panel_1d <- ggplot(ladder, aes(x = 1, y = layer, fill = status)) +
  geom_tile(width = 0.24, height = 0.65, colour = "#FFFFFF", linewidth = 0.5) +
  geom_text(aes(x = 1.18, label = detail), hjust = 0, size = 3.0, colour = "#202124") +
  scale_fill_manual(values = status_colours) +
  scale_x_continuous(limits = c(0.86, 2.25), expand = c(0, 0)) +
  labs(
    title = "Evidence by scale",
    subtitle = "Separate frozen decision rules",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line = element_blank(), legend.position = "none"
  )

figure_1 <- (panel_1a | panel_1b) / (panel_1c | panel_1d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

fig1_width <- 7.5
fig1_height <- fig1_width * 8.5 / 12
ggsave(file.path(out_dir, "Figure_1.pdf"), figure_1, width = fig1_width, height = fig1_height, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_1.png"), figure_1, width = fig1_width, height = fig1_height, dpi = 300)
ggsave(
  file.path(out_dir, "Fig1.tif"), figure_1, width = fig1_width, height = fig1_height,
  dpi = 300, device = "tiff", compression = "lzw"
)

# Figure 2: global phenotype-definition extension.
ext_rg <- fread(file.path(global_ext, "rg_ukb_definition_matrix.tsv"))
marginal <- fread(file.path(attrib_dir, "adjusted_marginal_rg.tsv"))

definition_labels <- c(
  UKB_LifetimeMDD = "Lifetime MDD",
  UKB_MDDRecur = "Recurrent MDD",
  UKB_GPpsy = "GP consultation",
  UKB_ICD10Dep = "ICD-10 depression"
)
ibd_labels <- c(
  IBD_deLange2017 = "IBD\nde Lange",
  CD_deLange2017 = "CD\nde Lange",
  UC_deLange2017 = "UC\nde Lange",
  IBD_FinnGen_R12 = "IBD\nFinnGen",
  CD_FinnGen_R12 = "CD\nFinnGen",
  UC_FinnGen_R12 = "UC\nFinnGen"
)
axis_colours <- c(depression_definition = "#0072B2", ibd_subtype = "#D55E00")

ext_rg[, depression_label := factor(
  definition_labels[trait1], levels = rev(unname(definition_labels))
)]
ext_rg[, ibd_label := factor(ibd_labels[trait2], levels = unname(ibd_labels))]
ext_rg[, label := sprintf(
  "%.2f%s%s", rg,
  ifelse(fdr_primary_extension < 0.05, "#", ""),
  ifelse(gcov_intercept_flag == 1, "*", "")
)]

panel_2a <- ggplot(ext_rg, aes(x = ibd_label, y = depression_label, fill = rg)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
    limits = c(-max(abs(ext_rg$rg)), max(abs(ext_rg$rg))), oob = squish
  ) +
  labs(
    title = "Genetic correlation across phenotype definitions",
    subtitle = "# BH-FDR < 0.05; * |genetic-covariance intercept z| >= 2",
    x = NULL, y = NULL, fill = "rg"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")

marginal[, display_level := fifelse(
  axis == "depression_definition", definition_labels[level], level
)]
marginal[, display_level := factor(
  display_level,
  levels = rev(c(unname(definition_labels), "IBD", "CD", "UC"))
)]
panel_2b <- ggplot(
  marginal,
  aes(x = estimate, y = display_level, color = axis, xmin = ci_low, xmax = ci_high)
) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.65) +
  geom_point(size = 2.2) +
  scale_color_manual(
    values = axis_colours,
    labels = c(
      depression_definition = "Depression definition",
      ibd_subtype = "IBD phenotype"
    )
  ) +
  labs(
    title = "Adjusted marginal correlations",
    x = "Adjusted genetic correlation", y = NULL,
    colour = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.margin = margin(0, 0, 0, 0)
  )

figure_2 <- panel_2a / panel_2b +
  plot_layout(heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

fig2_width <- 7.5
fig2_height <- fig2_width * 8.0 / 11.5
ggsave(file.path(out_dir, "Figure_2.pdf"), figure_2, width = fig2_width, height = fig2_height, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure_2.png"), figure_2, width = fig2_width, height = fig2_height, dpi = 300)
ggsave(
  file.path(out_dir, "Fig2.tif"), figure_2, width = fig2_width, height = fig2_height,
  dpi = 300, device = "tiff", compression = "lzw"
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "figure_sessionInfo.txt"))

cat("Created PLOS ONE submission figures in", out_dir, "\n")
