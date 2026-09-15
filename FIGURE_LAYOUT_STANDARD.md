# Figure layout standard

Apply these checks to every multi-panel figure before a submission package is frozen.

1. Use one font family and one size hierarchy across all panels: 9 pt ticks, 10 pt axis labels, 11 pt panel titles, and 13 pt panel letters at final figure size.
2. Center every panel title over its plotting area. Place panel letters from a fixed axes-relative anchor so letters align across rows and columns.
3. Reserve explicit layout columns or margins for legends and color bars. Legends must never be placed on top of plotted data.
4. Wrap or shorten long category labels before reducing font size. Confirm that the leftmost and rightmost glyphs remain inside the exported canvas.
5. Set figure dimensions before plotting and reserve panel spacing explicitly. Do not rely on `tight_layout()` to solve crowded multi-panel layouts.
6. Export a 300-dpi PNG preview and inspect it at 100% zoom. Check panel titles, labels, legends, plot boundaries, and white-space balance.
7. Render the submission PDF/TIFF after the PNG passes. Confirm that vector and raster exports show the same geometry and labels.
8. Keep the source table, plotting script, session information, final PDF/TIFF, preview PNG, and visual-QA record together.

The submission layouts for Fig 1, S1 Fig, and S2 Fig are implemented in
`scripts/28_refine_submission_figures.py`.
