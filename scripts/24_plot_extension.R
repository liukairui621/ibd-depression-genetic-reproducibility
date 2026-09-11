#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

root <- "/root/IBD/20_Reproducibility_Ladder"
global <- file.path(root, "results", "global_extension")
attribution <- file.path(global, "attribution")
figures <- file.path(root, "figures", "extension")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

rg <- read.delim(file.path(global, "rg_ukb_definition_matrix.tsv"))
marginal <- read.delim(file.path(attribution, "adjusted_marginal_rg.tsv"))
dispersion <- read.delim(file.path(attribution, "empirical_dispersion.tsv"))
ratio <- read.delim(file.path(attribution, "dispersion_ratio.tsv"))

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
axis_colors <- c(
  depression_definition = "#0072B2",
  ibd_subtype = "#D55E00"
)

theme_manuscript <- function() {
  theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold"),
      plot.margin = margin(7, 10, 7, 7)
    )
}

rg$depression_label <- factor(
  definition_labels[rg$trait1],
  levels = rev(unname(definition_labels))
)
rg$ibd_label <- factor(
  ibd_labels[rg$trait2],
  levels = unname(ibd_labels)
)
rg$label <- sprintf("%.2f", rg$rg)

heatmap <- ggplot(rg, aes(x = ibd_label, y = depression_label, fill = rg)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-max(abs(rg$rg)), max(abs(rg$rg))),
    oob = squish
  ) +
  labs(x = NULL, y = NULL, fill = expression(r[g])) +
  theme_manuscript() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )
ggsave(
  file.path(figures, "FigExt1A_UKB_definition_by_IBD_subtype_rg_heatmap.pdf"),
  heatmap,
  width = 7.2,
  height = 3.5,
  device = cairo_pdf
)

marginal$display_level <- ifelse(
  marginal$axis == "depression_definition",
  definition_labels[marginal$level],
  marginal$level
)
marginal$axis_label <- ifelse(
  marginal$axis == "depression_definition",
  "UKB depression definition",
  "IBD subtype"
)
marginal$display_level <- factor(
  marginal$display_level,
  levels = rev(c(
    unname(definition_labels),
    "IBD",
    "CD",
    "UC"
  ))
)

forest <- ggplot(
  marginal,
  aes(
    x = estimate,
    y = display_level,
    color = axis,
    xmin = ci_low,
    xmax = ci_high
  )
) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.65) +
  geom_point(size = 2.2) +
  facet_grid(axis_label ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = axis_colors) +
  labs(x = "Adjusted genetic correlation", y = NULL) +
  theme_manuscript() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0)
  )
ggsave(
  file.path(figures, "FigExt1B_adjusted_marginal_rg_forest.pdf"),
  forest,
  width = 5.4,
  height = 4.7,
  device = cairo_pdf
)

dispersion$axis_label <- factor(
  dispersion$axis,
  levels = c("depression_definition", "ibd_subtype"),
  labels = c("Depression definition", "IBD subtype")
)

dispersion_plot <- ggplot(
  dispersion,
  aes(x = axis_label, y = sd, color = axis)
) +
  geom_boxplot(
    aes(fill = axis),
    width = 0.45,
    alpha = 0.15,
    outlier.shape = NA,
    color = "grey40"
  ) +
  geom_jitter(width = 0.08, height = 0, size = 2.1) +
  scale_color_manual(values = axis_colors) +
  scale_fill_manual(values = axis_colors) +
  annotate(
    "label",
    x = 1.5,
    y = max(dispersion$sd) * 1.08,
    label = sprintf(
      "Adjusted SD ratio %.2f (95%% CI %.2f-%.2f)",
      ratio$ratio,
      ratio$ratio_ci_low,
      ratio$ratio_ci_high
    ),
    size = 3,
    linewidth = 0.2
  ) +
  expand_limits(y = max(dispersion$sd) * 1.22) +
  labs(x = NULL, y = expression("Within-stratum SD of " * r[g])) +
  theme_manuscript() +
  theme(legend.position = "none")
ggsave(
  file.path(figures, "FigExt1C_instability_dispersion_comparison.pdf"),
  dispersion_plot,
  width = 4.6,
  height = 4.1,
  device = cairo_pdf
)

writeLines(
  capture.output(sessionInfo()),
  file.path(figures, "sessionInfo.txt")
)
