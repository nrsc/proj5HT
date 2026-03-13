#!/usr/bin/env Rscript

# 5HT1 vs 5HT2 postsynaptic signaling summary table
# - Outputs a data.frame
# - Optionally renders a publication-ready table with gt

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

make_5ht_signaling_table <- function() {
  tibble::tibble(
    Feature = c(
      "G-protein",
      "cAMP",
      "PIP2",
      "Ca2+",
      "GIRK",
      "M-current",
      "TRPC",
      "Membrane potential",
      "Input resistance",
      "Spike output"
    ),
    `5-HT1 (Gi/o)` = c(
      "Gi/o",
      "↓",
      "—",
      "↓ VGCC",
      "ON",
      "—",
      "—",
      "Hyperpolarize",
      "↓",
      "↓"
    ),
    `5-HT2 (Gq/11)` = c(
      "Gq/11",
      "indirect / mixed",
      "↓",
      "↑ via IP3",
      "—",
      "OFF",
      "ON",
      "Depolarize",
      "↑",
      "↑ / gain modulation"
    )
  )
}

# ---- main ----
df <- make_5ht_signaling_table()

# Print as plain text table for terminals / logs
print(df)

# Optional: render to HTML (or PNG via webshot) if gt is installed
if (requireNamespace("gt", quietly = TRUE)) {
  gt_tbl <- df |>
    gt::gt() |>
    gt::tab_header(
      title = "Postsynaptic signaling: 5-HT1 vs 5-HT2 receptor families"
    ) |>
    gt::cols_label(
      Feature = "Feature",
      `5-HT1 (Gi/o)` = "5-HT₁ (Gi/o)",
      `5-HT2 (Gq/11)` = "5-HT₂ (Gq/11)"
    ) |>
    gt::tab_options(
      table.font.size = gt::px(14),
      heading.title.font.size = gt::px(16)
    )

  # Save HTML next to where you run the script
  out_html <- "5HT1_vs_5HT2_postsyn_signaling_table.html"
  gt::gtsave(gt_tbl, out_html)
  message("Saved gt table to: ", out_html)
} else {
  message("Package 'gt' not installed; skipping formatted table export.")
  message("Install with: install.packages('gt')")
}
