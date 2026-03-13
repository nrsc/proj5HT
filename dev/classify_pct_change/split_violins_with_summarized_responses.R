library(dplyr)
library(ggplot2)
library(scales)

summarise_puff_responses <- function(df,
                                     baseline_window = c(-5, 0),
                                     response_window = c(0, Inf),
                                     baseline_min_hz = 0.5,
                                     keep_subclasses = c("L23_IT","L5_ET","L5_IT")) {

  # Fix subclass naming (merge L3c → L23_IT)
  df2 <- df %>%
    dplyr::mutate(
      assigned_subclass = dplyr::if_else(assigned_subclass == "L3c", "L23_IT", assigned_subclass)
    ) %>%
    dplyr::filter(assigned_subclass %in% keep_subclasses)

  # ---- Baseline: mean rate in last 5 seconds before puff ----
  base_tbl <- df2 %>%
    dplyr::filter(time >= baseline_window[1], time < baseline_window[2]) %>%
    dplyr::group_by(cell_name, assigned_subclass) %>%
    dplyr::summarise(
      baseline_rate = mean(instRate, na.rm = TRUE),
      baseline_n    = dplyr::n(),
      .groups = "drop"
    )

  # ---- Response window: after puff onset ----
  resp_tbl <- df2 %>%
    dplyr::filter(time >= response_window[1], time <= response_window[2]) %>%
    dplyr::group_by(cell_name, assigned_subclass) %>%
    dplyr::summarise(
      response_min_hz = suppressWarnings(min(instRate, na.rm = TRUE)),
      response_max_hz = suppressWarnings(max(instRate, na.rm = TRUE)),
      response_n      = dplyr::n(),
      .groups = "drop"
    )

  per_cell <- base_tbl %>%
    dplyr::inner_join(resp_tbl, by = c("cell_name","assigned_subclass"))

  # ---- QC + percent-of-baseline metrics ----
  per_cell <- per_cell %>%
    dplyr::mutate(
      valid = baseline_rate >= baseline_min_hz & baseline_n >= 10 & response_n >= 10,

      # Percent of baseline (100 = no change, 90 = 10% inhibition, 0 = full silence)
      pc_min = 100 * response_min_hz / baseline_rate,
      pc_max = 100 * response_max_hz / baseline_rate,

      # Direction based on dominant effect
      direction = dplyr::case_when(
        pc_min < 100 & pc_max <= 100 ~ "inhibition",
        pc_max > 100 & pc_min >= 100 ~ "excitation",
        (100 - pc_min) > (pc_max - 100) ~ "inhibition",
        (pc_max - 100) > (100 - pc_min) ~ "excitation",
        TRUE ~ "no_change"
      ),

      # Peak response (percent-of-baseline scale, centered at 100)
      response_peak = dplyr::if_else(
        direction == "inhibition", pc_min, pc_max
      ),

      magnitude = abs(response_peak - 100)
    )

  per_cell
}

per_cell_plot <- summarise_puff_responses(df)

summary(per_cell_plot$response_peak)
table(per_cell_plot$direction)

prep_for_split_violins <- function(per_cell, baseline = 100, gap = 5) {

  stopifnot("direction" %in% names(per_cell))

  per_cell %>%
    dplyr::filter(valid) %>%
    dplyr::filter(direction %in% c("excitation","inhibition")) %>%
    dplyr::mutate(
      side = direction,

      # shift slightly apart so they are two separate violins
      subclass_side = paste(assigned_subclass, side, sep = "_"),

      # signed value around baseline
      value = dplyr::if_else(
        direction == "excitation",
        response_peak - baseline + gap,
        response_peak - baseline - gap
      )
    )
}

per_cell_split <- prep_for_split_violins(per_cell_plot)

plot_split_violins <- function(per_cell_split,
                               title = "Per-cell puff responses (baseline = 100%)") {

  ggplot(per_cell_split,
         aes(x = assigned_subclass,
             y = value,
             fill = side)) +

    geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +

    # zero line = baseline
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +

    scale_fill_manual(
      values = c(
        excitation = "#2166AC",   # blue
        inhibition = "#B2182B"    # red
      ),
      name = "Direction"
    ) +

    scale_y_continuous(
      labels = function(x) paste0(round(x + 100), "%")
    ) +

    labs(
      x = "Subclass",
      y = "Response (% of baseline)",
      title = title
    ) +

    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

p_split <- plot_split_violins(per_cell_split)
print(p_split)


per_cell_plot <- per_cell_plot %>%
  dplyr::mutate(
    hits_zero = response_min_hz <= 0.05
  )

table(per_cell_plot$hits_zero, per_cell_plot$direction)



