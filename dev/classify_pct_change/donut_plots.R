# res$per_cell must have: cell_name, response_call
# If you also have delta_pct or similar in per_cell, great; not required for this plot.

make_response_donut <- function(per_cell,
                                label_min_frac = 0.04,  # only label wedges >= 4%
                                title = "Distribution of Puff Responses") {

  stopifnot("response_call" %in% names(per_cell))

  levels_order <- c(
    "strong_inhibition",
    "moderate_inhibition",
    "mild_inhibition",
    "no_change",
    "mild_excitation",
    "moderate_excitation",
    "strong_excitation"
  )

  response_cols <- c(
    strong_inhibition   = "#08306B",
    moderate_inhibition = "#2171B5",
    mild_inhibition     = "#6BAED6",
    no_change           = "grey80",
    mild_excitation     = "#FC9272",
    moderate_excitation = "#FB6A4A",
    strong_excitation   = "#CB181D"
  )

  pie_df <- per_cell %>%
    dplyr::filter(!is.na(response_call)) %>%
    dplyr::mutate(response_call = factor(response_call, levels = levels_order)) %>%
    dplyr::count(response_call, name = "n") %>%
    dplyr::mutate(
      frac = n / sum(n),
      pct  = scales::percent(frac, accuracy = 0.1),
      # show wedge label only if slice is "big enough"
      wedge_label = ifelse(frac >= label_min_frac, pct, ""),
      # legend label always has both n and percent
      legend_label = paste0(as.character(response_call), " — ", n, " (", pct, ")")
    ) %>%
    dplyr::arrange(response_call) %>%
    dplyr::mutate(ypos = cumsum(frac) - 0.5 * frac)

  ggplot(pie_df, aes(x = 2, y = frac, fill = response_call)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +                           # donut hole
    #geom_text(aes(y = ypos), size = 3) +
    scale_fill_manual(
      values = response_cols,
      drop = FALSE,
      labels = pie_df$legend_label
    ) +
    labs(title = title, fill = "Response") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# Use it:
p_donut <- make_response_donut(res$per_cell, label_min_frac = 0.04)
print(p_donut)

make_response_donut_by_subclass <- function(per_cell,
                                            label_min_frac = 0.06,   # only label wedges >= 6%
                                            title = "Distribution of Puff Responses by Subclass") {

  stopifnot("response_call" %in% names(per_cell))
  stopifnot("assigned_subclass" %in% names(per_cell))

  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(scales)
  })

  # Order of response classes
  levels_order <- c(
    "strong_inhibition",
    "moderate_inhibition",
    "mild_inhibition",
    "no_change",
    "mild_excitation",
    "moderate_excitation",
    "strong_excitation"
  )

  # Color palette
  response_cols <- c(
    strong_inhibition   = "#08306B",
    moderate_inhibition = "#2171B5",
    mild_inhibition     = "#6BAED6",
    no_change           = "grey80",
    mild_excitation     = "#FC9272",
    moderate_excitation = "#FB6A4A",
    strong_excitation   = "#CB181D"
  )

  # -----------------------------
  # Build per-subclass pie data
  # -----------------------------
  pie_df <- per_cell %>%
    dplyr::filter(!is.na(response_call), !is.na(assigned_subclass)) %>%
    dplyr::mutate(response_call = factor(response_call, levels = levels_order)) %>%
    dplyr::count(assigned_subclass, response_call, name = "n") %>%
    dplyr::group_by(assigned_subclass) %>%
    dplyr::mutate(
      frac = n / sum(n),
      pct  = scales::percent(frac, accuracy = 0.1),
      wedge_label = ifelse(frac >= label_min_frac, pct, ""),
      ypos = cumsum(frac) - 0.5 * frac,
      n_total = sum(n)
    ) %>%
    dplyr::ungroup()

  # Make a nicer facet label with N per subclass
  pie_df <- pie_df %>%
    dplyr::mutate(
      subclass_label = paste0(assigned_subclass, "  (N = ", n_total, ")")
    )

  # -----------------------------
  # Plot
  # -----------------------------
  ggplot(pie_df, aes(x = 2, y = frac, fill = response_call)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +   # donut hole

    # Only label large enough wedges
    geom_text(
      aes(y = ypos, label = wedge_label),
      size = 2.8,
      color = "black"
    ) +

    scale_fill_manual(
      values = response_cols,
      drop = FALSE
    ) +

    facet_wrap(~ subclass_label) +

    labs(
      title = title,
      fill = "Response class"
    ) +

    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "right"
    )
}


p_subclass <- make_response_donut_by_subclass(res$per_cell, label_min_frac = 0.06)
print(p_subclass)

p_subclass <- make_response_donut_by_subclass(
  res$per_cell_sub,
  label_min_frac = 0.06,
  title = "Puff responses (L23_IT, L5_ET, L5_IT only)"
)
print(p_subclass)

