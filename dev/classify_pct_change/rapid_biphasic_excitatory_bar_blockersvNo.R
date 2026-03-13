suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

## --- 1) Ensure blockers column exists in per_cell2 (join from df if needed) ---
if (!"blockers" %in% names(per_cell2)) {
  blockers_map <- df %>% dplyr::distinct(cell_name, blockers)
  per_cell2 <- per_cell2 %>% dplyr::left_join(blockers_map, by = "cell_name")
}

## --- 2) Make 2-level blockers bin (no_blockers vs has_blockers) ---
per_cell2 <- per_cell2 %>%
  dplyr::mutate(
    blockers_bin = dplyr::if_else(is.na(blockers) | blockers == "", "no_blockers", "has_blockers"),
    blockers_bin = factor(blockers_bin, levels = c("no_blockers", "has_blockers"))
  )

## --- 3) Filter to rapid biphasic with excitation then inhibition ---
## This supports either:
##  - response_type == "biphasic" AND biphasic_order == "exc_then_inh"
##  - (optional) rapid_pattern strings if you used those earlier
per_cell_exc_then_inh <- per_cell2 %>%
  dplyr::filter(response_type == "biphasic") %>%
  dplyr::filter(biphasic_order == "exc_then_inh") %>%
  # choose what magnitude you want to compare:
  # A) exc_mag = peak excitation above baseline (%)
  dplyr::filter(is.finite(exc_mag), exc_mag > 0)

## quick sanity check
print(table(per_cell_exc_then_inh$blockers_bin, useNA = "ifany"))

## --- 4) Bar graph: excitation magnitude by blockers bin (mean ± 95% CI) ---
exc_summ <- per_cell_exc_then_inh %>%
  dplyr::group_by(blockers_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_exc = mean(exc_mag, na.rm = TRUE),
    sd_exc   = sd(exc_mag, na.rm = TRUE),
    se_exc   = sd_exc / sqrt(n),
    ci95     = 1.96 * se_exc,
    .groups = "drop"
  )

p_bar <- ggplot(exc_summ, aes(x = blockers_bin, y = mean_exc)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_exc - ci95, ymax = mean_exc + ci95),
                width = 0.18, linewidth = 0.7) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.6, size = 3.8) +
  labs(
    x = NULL,
    y = "Excitation magnitude (% above baseline)",
    title = "Rapid biphasic (exc→inh): excitation magnitude vs blockers",
    subtitle = "Mean ± 95% CI; cells restricted to L23_IT/L5_ET/L5_IT upstream"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_bar)

## --- 5) Optional QC: show the distribution too (often more informative than bars) ---
p_dist <- ggplot(per_cell_exc_then_inh, aes(x = blockers_bin, y = exc_mag)) +
  geom_violin(trim = TRUE, alpha = 0.55) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1) +
  labs(
    x = NULL,
    y = "Excitation magnitude (% above baseline)",
    title = "QC: excitation magnitude distribution (rapid biphasic exc→inh)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_dist)
