# ============================================================
# Puff response plotting toolkit (talk slides + QC)
# - Blue = excitation, Red = inhibition
# - Handles subclass fixups (L3c -> L23_IT)
# - Filters to L23_IT / L5_ET / L5_IT
# - Multiple visualization options
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

# ----------------------------
# 0) Helpers: subclass handling
# ----------------------------

ensure_subclass <- function(per_cell, raw_df = NULL, cell_col = "cell_name") {
  stopifnot(is.data.frame(per_cell))
  stopifnot(cell_col %in% names(per_cell))

  # If already present, return
  if ("assigned_subclass" %in% names(per_cell)) return(per_cell)

  if (is.null(raw_df)) {
    stop("per_cell has no 'assigned_subclass'. Provide raw_df (your long df) to join it.")
  }
  stopifnot(is.data.frame(raw_df))
  stopifnot(cell_col %in% names(raw_df))
  stopifnot("assigned_subclass" %in% names(raw_df))

  map_sub <- raw_df %>%
    dplyr::select(dplyr::all_of(c(cell_col, "assigned_subclass"))) %>%
    dplyr::distinct()

  per_cell %>%
    dplyr::left_join(map_sub, by = setNames(cell_col, cell_col))
}

prep_subclasses <- function(per_cell,
                            keep = c("L23_IT","L5_ET","L5_IT"),
                            recode_map = c("L3c" = "L23_IT")) {
  stopifnot(is.data.frame(per_cell))
  stopifnot("assigned_subclass" %in% names(per_cell))

  per_cell %>%
    dplyr::mutate(
      assigned_subclass = dplyr::if_else(
        assigned_subclass %in% names(recode_map),
        unname(recode_map[assigned_subclass]),
        assigned_subclass
      ),
      assigned_subclass = as.character(assigned_subclass)
    ) %>%
    dplyr::filter(assigned_subclass %in% keep) %>%
    dplyr::mutate(
      assigned_subclass = factor(assigned_subclass, levels = keep)
    )
}

# ----------------------------
# 1) Shared settings: response levels + colors
# ----------------------------

response_levels <- c(
  "strong_inhibition",
  "moderate_inhibition",
  "mild_inhibition",
  "no_change",
  "mild_excitation",
  "moderate_excitation",
  "strong_excitation"
)

# Conventional (requested): BLUE = excitation, RED = inhibition
# (Grey for no change)
response_cols <- c(
  strong_inhibition   = "#B2182B",
  moderate_inhibition = "#D6604D",
  mild_inhibition     = "#F4A582",
  no_change           = "grey80",
  mild_excitation     = "#92C5DE",
  moderate_excitation = "#4393C3",
  strong_excitation   = "#2166AC"
)

pretty_labels <- c(
  strong_inhibition   = "Strong inhibition",
  moderate_inhibition = "Moderate inhibition",
  mild_inhibition     = "Mild inhibition",
  no_change           = "No change",
  mild_excitation     = "Mild excitation",
  moderate_excitation = "Moderate excitation",
  strong_excitation   = "Strong excitation"
)

# ----------------------------
# 2) Prep plotting data
# ----------------------------

prep_plot_df <- function(per_cell) {
  stopifnot("response_call" %in% names(per_cell))
  stopifnot("assigned_subclass" %in% names(per_cell))

  per_cell %>%
    dplyr::filter(!is.na(response_call)) %>%
    dplyr::mutate(
      response_call = factor(response_call, levels = response_levels)
    )
}

# ----------------------------
# 3) Plot A: 100% stacked bars (best “at-a-glance” slide figure)
# ----------------------------

plot_stacked_100 <- function(per_cell,
                             title = "Puff response composition by subclass",
                             show_counts_in_strip = TRUE) {
  dat <- prep_plot_df(per_cell)

  # N per subclass (for facet titles / captions)
  n_sub <- dat %>% count(assigned_subclass, name = "N")

  p <- dat %>%
    count(assigned_subclass, response_call, name = "n") %>%
    group_by(assigned_subclass) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = assigned_subclass, y = frac, fill = response_call)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.3) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = response_cols, drop = FALSE, labels = pretty_labels) +
    labs(x = NULL, y = "Percent of cells", fill = "Response class", title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  if (show_counts_in_strip) {
    # Add N above bars
    p <- p +
      geom_text(
        data = n_sub,
        aes(x = assigned_subclass, y = 1.02, label = paste0("N=", N)),
        inherit.aes = FALSE,
        size = 3.6,
        vjust = 0
      ) +
      coord_cartesian(clip = "off")
  }

  p
}

# ----------------------------
# 4) Plot B: Heatmap of proportions (great for QC comparisons)
# ----------------------------

plot_heatmap_prop <- function(per_cell,
                              title = "Response proportions (heatmap)") {
  dat <- prep_plot_df(per_cell)

  hm <- dat %>%
    count(assigned_subclass, response_call, name = "n") %>%
    group_by(assigned_subclass) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup()

  ggplot(hm, aes(x = response_call, y = assigned_subclass, fill = frac)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = percent(frac, accuracy = 1)), size = 3.4) +
    scale_x_discrete(labels = pretty_labels) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = NULL, fill = "Percent", title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# ----------------------------
# 5) Plot C: Dot + CI (binomial) for each response class by subclass
#    More “stat-y”, good for internal QC
# ----------------------------

plot_prop_points <- function(per_cell,
                             title = "Response proportions with 95% CI") {
  dat <- prep_plot_df(per_cell)

  dfp <- dat %>%
    count(assigned_subclass, response_call, name = "n") %>%
    group_by(assigned_subclass) %>%
    mutate(
      N = sum(n),
      frac = n / N
    ) %>%
    ungroup()

  # Approximate binomial SE CI (Wilson would be nicer, but this is simple & stable)
  dfp <- dfp %>%
    mutate(
      se = sqrt(pmax(frac * (1 - frac) / pmax(N, 1), 0)),
      lo = pmax(frac - 1.96 * se, 0),
      hi = pmin(frac + 1.96 * se, 1)
    )

  ggplot(dfp, aes(x = response_call, y = frac, color = response_call)) +
    geom_point(size = 2.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.5) +
    facet_wrap(~assigned_subclass, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_x_discrete(labels = pretty_labels) +
    scale_color_manual(values = response_cols, drop = FALSE, guide = "none") +
    labs(x = NULL, y = "Percent of cells", title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# ----------------------------
# 6) Plot D: “Small multiples” donuts (if you still want donuts)
#    Fixes label placement by labeling OUTSIDE with leader lines.
# ----------------------------

plot_donuts_by_subclass <- function(per_cell,
                                    label_min_frac = 0.06,
                                    title = "Puff responses (donut, by subclass)") {
  dat <- prep_plot_df(per_cell)

  pie_df <- dat %>%
    count(assigned_subclass, response_call, name = "n") %>%
    group_by(assigned_subclass) %>%
    mutate(
      frac = n / sum(n),
      pct  = percent(frac, accuracy = 0.1),
      show_label = frac >= label_min_frac
    ) %>%
    arrange(assigned_subclass, response_call) %>%
    group_by(assigned_subclass) %>%
    mutate(
      ymax = cumsum(frac),
      ymin = dplyr::lag(ymax, default = 0),
      ymid = (ymax + ymin) / 2,
      # for outside labels:
      label = ifelse(show_label, pct, "")
    ) %>%
    ungroup()

  # Put labels just outside donut: compute x position and alignment
  pie_df <- pie_df %>%
    mutate(
      x = 2,
      x_label = 2.55,
      hjust = ifelse(ymid > 0.5, 1, 0) # crude but works ok for polar coords
    )

  ggplot(pie_df, aes(x = x, y = frac, fill = response_call)) +
    geom_col(color = "white", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.7, 2.9) +  # donut hole + space for outside labels
    facet_wrap(~assigned_subclass, nrow = 1) +
    scale_fill_manual(values = response_cols, drop = FALSE, labels = pretty_labels) +
    # outside labels
    # geom_text(
    #   data = subset(pie_df, show_label),
    #   aes(x = x_label, y = ymid, label = label),
    #   inherit.aes = FALSE,
    #   size = 3.2
    # ) +
    labs(title = title, fill = "Response class") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

# ----------------------------
# 7) QC plot E: per-cell magnitudes (directional)
#    Uses your per_cell columns: magnitude/inh_mag/exc_mag
# ----------------------------

plot_magnitude_violin <- function(per_cell,
                                  title = "Per-cell response magnitude by subclass",
                                  y_col_preference = c("magnitude", "magnitude_signed", "delta_pc", "exc_mag")) {

  stopifnot(is.data.frame(per_cell))
  stopifnot("assigned_subclass" %in% names(per_cell))

  # Choose a magnitude-like column that exists
  ycol <- y_col_preference[y_col_preference %in% names(per_cell)][1]
  if (is.na(ycol) || length(ycol) == 0) {
    stop("No magnitude column found. Tried: ", paste(y_col_preference, collapse = ", "))
  }

  dfm <- per_cell %>%
    dplyr::mutate(
      mag_plot = .data[[ycol]],
      # direction: use existing if present, else infer from sign if mag_plot is signed
      direction_plot = if ("direction" %in% names(per_cell)) as.character(direction) else ifelse(mag_plot >= 0, "excitation", "inhibition"),
      direction_plot = dplyr::case_when(
        direction_plot %in% c("excitation","Excitation","exc") ~ "excitation",
        direction_plot %in% c("inhibition","Inhibition","inh") ~ "inhibition",
        TRUE ~ direction_plot
      ),
      direction_plot = factor(direction_plot, levels = c("inhibition","excitation"))
    ) %>%
    dplyr::filter(!is.na(mag_plot))

  ggplot(dfm, aes(x = assigned_subclass, y = mag_plot, fill = direction_plot)) +
    geom_violin(trim = TRUE, alpha = 0.85) +
    geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.75) +
    scale_fill_manual(values = c(inhibition = "#B2182B", excitation = "#2166AC")) +
    labs(x = NULL, y = ycol, fill = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
}

# Use it:
if (any(c("magnitude","magnitude_signed","delta_pc","exc_mag") %in% names(per_cell_plot))) {
  p_mag <- plot_magnitude_violin(per_cell_plot, title = "Per-cell response magnitude (QC)")
  print(p_mag)
}


# ----------------------------
# 8) RUN SECTION: create per_cell_plot and draw plots
# ----------------------------

# Starting objects assumed:
# - res$per_cell exists
# - df exists ONLY if res$per_cell lacks assigned_subclass
per_cell_plot <- ensure_subclass(df_prepped$per_cell, raw_df = df) %>%
  prep_subclasses(keep = c("L23_IT","L5_ET","L5_IT"),
                  recode_map = c("L3c" = "L23_IT"))

# A) Recommended slide figure:
p_stack <- plot_stacked_100(per_cell2,
                            title = "Puff responses by subclass (L23_IT, L5_ET, L5_IT)")
print(p_stack)

# B) QC-friendly heatmap:
p_heat <- plot_heatmap_prop(per_cell_plot,
                            title = "Response proportions (heatmap)")
print(p_heat)

# C) Points + CI:
p_pts <- plot_prop_points(per_cell_plot,
                          title = "Response proportions with 95% CI")
print(p_pts)

# D) Donuts (optional):
p_donut <- plot_donuts_by_subclass(per_cell_plot,
                                   label_min_frac = 0.06,
                                   title = "Puff responses")
print(p_donut)

# E) Magnitude distribution (optional internal QC):
# (Only works if you have per_cell$magnitude or similar)
if ("magnitude" %in% names(per_cell_plot) || "exc_mag" %in% names(per_cell_plot)) {
  p_mag <- plot_magnitude_violin(per_cell_plot,
                                 title = "Per-cell response magnitude (QC)")
  print(p_mag)
}

# ----------------------------
# 9) QC loop: save per-cell trace figures for spot-checking (optional)
#    Requires your long df with time, percent_change, protocol, cell_name.
#    Also assumes your percent_change is already baseline-relative (100=baseline).
# ----------------------------

# qc_save_traces <- function(raw_df,
#                            per_cell,
#                            out_dir = "qc_traces",
#                            n_cells_per_subclass = 8,
#                            seed = 1) {
#   stopifnot(is.data.frame(raw_df), is.data.frame(per_cell))
#   req <- c("cell_name","time","percent_change","protocol")
#   if (!all(req %in% names(raw_df))) {
#     stop("raw_df must have: ", paste(req, collapse = ", "))
#   }
#   if (!all(c("cell_name","assigned_subclass","response_call") %in% names(per_cell))) {
#     stop("per_cell must have: cell_name, assigned_subclass, response_call")
#   }
#
#   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#   set.seed(seed)
#
#   # sample cells per subclass for QC
#   pick <- per_cell %>%
#     filter(!is.na(response_call)) %>%
#     group_by(assigned_subclass) %>%
#     slice_sample(n = min(n_cells_per_subclass, n()), replace = FALSE) %>%
#     ungroup()
#
#   for (i in seq_len(nrow(pick))) {
#     cn <- pick$cell_name[i]
#     sc <- as.character(pick$assigned_subclass[i])
#     rc <- as.character(pick$response_call[i])
#
#     d0 <- raw_df %>%
#       filter(cell_name == cn) %>%
#       filter(!tolower(protocol) %in% c("baseline"))  # remove baseline rows if present
#
#     if (nrow(d0) == 0) next
#
#     p <- ggplot(d0, aes(x = time, y = percent_change)) +
#       geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4) +
#       geom_line(linewidth = 0.6) +_
#

per_cell_plot <- per_cell_plot %>%
  dplyr::mutate(
    response_peak = dplyr::case_when(
      direction == "inhibition" ~ pc_min,
      direction == "excitation" ~ pc_max,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::filter(!is.na(response_peak))


plot_centered_violin <- function(per_cell,
                                 title = "Per-cell response magnitude (baseline-centered QC)") {

  stopifnot("assigned_subclass" %in% names(per_cell))
  stopifnot("response_peak" %in% names(per_cell))

  ggplot(per_cell, aes(x = assigned_subclass, y = response_peak)) +

    # Baseline reference
    geom_hline(yintercept = 100, linetype = 2, linewidth = 0.6, color = "black") +

    # Main violin (single distribution per subclass)
    geom_violin(
      fill = "grey85",
      color = "black",
      linewidth = 0.4,
      trim = TRUE
    ) +

    # Optional: light jittered points for QC (can comment out if too busy)
    geom_point(
      aes(color = direction),
      position = position_jitter(width = 0.12),
      size = 1.4,
      alpha = 0.6
    ) +

    # Blue = excitation, Red = inhibition (as requested)
    scale_color_manual(
      values = c(inhibition = "#B2182B", excitation = "#2166AC"),
      guide = "none"
    ) +

    coord_flip() +

    labs(
      x = NULL,
      y = "Percent of baseline (100 = no change)",
      title = title
    ) +

    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

p_center <- plot_centered_violin(per_cell_plot,
                                 title = "Per-cell puff responses (baseline-centered)")
print(p_center)

