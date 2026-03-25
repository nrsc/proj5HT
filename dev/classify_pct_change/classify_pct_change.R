# ============================================================
# Classify puff responses from percent_change traces
# - Filters out Baseline using protocol / Protocol
# - Interprets percent_change relative to 100:
#       90  -> 10% inhibition
#       110 -> 10% excitation
# - Produces per-cell classification into 7 buckets
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------------------------------------
# Main classification function
# ------------------------------------------------------------
classify_percent_change <- function(
    df,
    protocol_values_to_drop = c("baseline"),
    time_window = NULL,          # e.g. c(0, 30) to restrict to puff window
    min_points = 20,
    mild_thresh = 10,            # 10% change = mild
    moderate_thresh = 25,        # 25% change = moderate
    strong_thresh = 60           # 60% change = strong
) {

  stopifnot(is.data.frame(df))
  stopifnot("cell_name" %in% names(df))
  stopifnot("percent_change" %in% names(df))

  # ----------------------------
  # Identify protocol column
  # ----------------------------
  protocol_col <- NULL
  if ("Protocol" %in% names(df)) protocol_col <- "Protocol"
  if ("protocol" %in% names(df)) protocol_col <- "protocol"

  if (is.null(protocol_col)) {
    message("WARNING: No Protocol/protocol column found — baseline filtering skipped.")
  }

  # ----------------------------
  # Filter out baseline rows
  # ----------------------------
  df2 <- df

  if (!is.null(protocol_col)) {
    drop_set <- tolower(protocol_values_to_drop)

    df2 <- df2 %>%
      dplyr::filter(
        !is.na(.data[[protocol_col]]),
        !(tolower(.data[[protocol_col]]) %in% drop_set)
      )
  }

  # ----------------------------
  # Optional time window (post-0 etc.)
  # ----------------------------
  if (!is.null(time_window)) {
    stopifnot("time" %in% names(df2))
    stopifnot(is.numeric(time_window), length(time_window) == 2)

    tmin <- min(time_window)
    tmax <- max(time_window)

    df2 <- df2 %>%
      dplyr::filter(.data$time >= tmin, .data$time <= tmax)
  }

  # ----------------------------
  # Per-cell classification
  # ----------------------------
  per_cell <- df2 %>%
    group_by(cell_name) %>%
    summarise(
      n_points = sum(!is.na(percent_change)),

      pc_min = min(percent_change, na.rm = TRUE),
      pc_max = max(percent_change, na.rm = TRUE),

      # deviations relative to baseline (100)
      inh_mag = 100 - pc_min,   # e.g. 90  -> 10% inhibition
      exc_mag = pc_max - 100,   # e.g. 110 -> 10% excitation

      .groups = "drop"
    ) %>%
    mutate(
      valid = n_points >= min_points,

      # Decide dominant direction
      direction = case_when(
        !valid ~ "insufficient-data",
        is.infinite(pc_min) | is.infinite(pc_max) ~ "insufficient-data",
        inh_mag > exc_mag ~ "inhibition",
        exc_mag > inh_mag ~ "excitation",
        TRUE ~ "mixed/flat"
      ),

      magnitude = pmax(inh_mag, exc_mag, na.rm = TRUE),

      # Strength bins
      strength = case_when(
        direction == "insufficient-data" ~ "insufficient-data",
        magnitude >= strong_thresh   ~ "strong",
        magnitude >= moderate_thresh ~ "moderate",
        magnitude >= mild_thresh     ~ "mild",
        TRUE                          ~ "no_change"
      ),

      # Final 7-level call
      response_call = case_when(
        strength == "no_change" ~ "no_change",
        direction == "inhibition" ~ paste0(strength, "_inhibition"),
        direction == "excitation" ~ paste0(strength, "_excitation"),
        direction == "mixed/flat" ~ "mixed/flat",
        TRUE ~ direction
      ),

      # Numeric rank (useful later)
      response_rank = case_when(
        response_call == "strong_excitation"   ~  3L,
        response_call == "moderate_excitation" ~  2L,
        response_call == "mild_excitation"     ~  1L,
        response_call == "no_change"            ~  0L,
        response_call == "mild_inhibition"     ~ -1L,
        response_call == "moderate_inhibition" ~ -2L,
        response_call == "strong_inhibition"   ~ -3L,
        TRUE ~ NA_integer_
      )
    )

  list(
    df_filtered = df2,
    per_cell = per_cell
  )
}

# ============================================================
# Example usage on your data
# ============================================================

# Classify all non-baseline data (no time restriction yet)
res <- classify_percent_change(
  df,
  protocol_values_to_drop = c("baseline"),
  time_window = c(0, 100),   # only analyze 0–30 s after puff
  mild_thresh = 20,
  moderate_thresh = 35,
  strong_thresh = 70
)

file.edit("dev/build_response_classifier/puff_helpers.R")
file.edit("dev/build_response_classifier/puff_processing.R")
file.edit("dev/build_response_classifier/puff_plotting.R")
file.edit("dev/build_response_classifier/run_puff_classifier.R")

# ------------------------------------------------------------
# Check the inhibitory cell you gave
# ------------------------------------------------------------
res$per_cell %>%
  dplyr::filter(cell_name == "QF25.26.024.19.07.04")

# res$per_cell %>%
#   dplyr::filter(cell_name == "Q21.26.021.1A.02.01")

# You SHOULD see something like:
#   pc_min ~ 13–40
#   inh_mag huge (60–90)
#   response_call = "strong_inhibition"

# ------------------------------------------------------------
# Global sanity checks
# ------------------------------------------------------------

# Distribution of calls
res$per_cell %>%
  count(response_call, sort = TRUE)

# Most inhibitory cells
res$per_cell %>%
  arrange(pc_min) %>%
  select(cell_name, pc_min, pc_max, inh_mag, exc_mag, response_call) %>%
  head(20)

# Most excitatory cells
res$per_cell %>%
  arrange(desc(pc_max)) %>%
  select(cell_name, pc_min, pc_max, inh_mag, exc_mag, response_call) %>%
  head(20)

res$per_cell <- res$per_cell %>%
  left_join(
    df %>%
      distinct(cell_name, assigned_subclass),
    by = "cell_name"
  )

res$per_cell <- res$per_cell %>%
  dplyr::mutate(
    assigned_subclass = dplyr::recode(
      assigned_subclass,
      "L3c" = "L23_IT"
    )
  )

res$per_cell_sub <- res$per_cell %>%
  dplyr::filter(assigned_subclass %in% c("L23_IT", "L5_ET", "L5_IT"))




