library(dplyr)
library(ggplot2)

prep_for_split_violins <- function(per_cell, baseline = 100, gap = 5) {

  stopifnot(all(c("assigned_subclass", "response_call", "response_peak") %in% names(per_cell)))

  per_cell %>%
    dplyr::filter(
      !is.na(assigned_subclass),
      !is.na(response_call),
      !is.na(response_peak)
    ) %>%
    dplyr::mutate(
      # Infer direction from response_call (much safer than trusting a dropped column)
      direction = dplyr::case_when(
        grepl("inhibition", response_call) ~ "inhibition",
        grepl("excitation", response_call) ~ "excitation",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(direction)) %>%
    dplyr::mutate(
      # True delta from baseline (e.g. 90 -> -10 ; 120 -> +20)
      delta = response_peak - baseline,
      mag   = abs(delta),

      # Build split geometry with a gap around 0 so violins don't touch
      y_plot = dplyr::case_when(
        direction == "excitation" ~  (mag + gap),
        direction == "inhibition" ~ -(mag + gap),
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::filter(!is.na(y_plot))
}

per_cell_split <- prep_for_split_violins(per_cell_plot, baseline = 100, gap = 5)


plot_split_violin_baseline_clean <- function(per_cell_split,
                                             gap = 5,
                                             title = "Per-cell response magnitude (Δ from baseline, split)") {

  cols <- c(inhibition = "#B2182B", excitation = "#2166AC")

  ggplot(per_cell_split,
         aes(x = assigned_subclass, y = y_plot, fill = direction)) +

    # Baseline reference (0 = no change)
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, color = "black") +

    # Split violins
    geom_violin(
      color = "black",
      linewidth = 0.35,
      trim = TRUE
    ) +

    # Light QC points (optional – comment out if you want ultra-clean slides)
    geom_point(
      aes(color = direction),
      position = position_jitter(width = 0.10, height = 0),
      size = 1.2,
      alpha = 0.45,
      show.legend = FALSE
    ) +

    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +

    # Axis shows TRUE delta from baseline (gap removed from labels)
    scale_y_continuous(
      breaks = function(lims) {
        b <- pretty(lims, n = 7)
        # avoid ticks inside the artificial gap
        b[abs(b) >= gap]
      },
      labels = function(v) {
        signv <- sign(v)
        true_delta <- pmax(abs(v) - gap, 0) * signv
        sprintf("%+g", round(true_delta, 0))
      }
    ) +

    coord_flip() +

    labs(
      x = NULL,
      y = "Δ from baseline (%)   (0 = no change; −10 = 10% inhibition; +10 = 10% excitation)",
      fill = "Direction",
      title = title
    ) +

    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}


p_split <- plot_split_violin_baseline_clean(per_cell_split, gap = 5,
                                            title = "Per-cell puff responses (baseline-centered, split)")
print(p_split)

