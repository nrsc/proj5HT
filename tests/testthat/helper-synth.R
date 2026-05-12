# Synthetic trace builder used by the classifier tests.
#
# Each cell is a long-format dataframe with `cell_name`, `time`,
# `percent_change`, and a baseline-locked sampling rate. `kind` controls
# the shape of the response:
#   "flat"            — no change
#   "excitation"      — sustained rise above 100 in (0, 30)
#   "inhibition"      — sustained dip below 100 in (0, 30)
#   "biphasic_ei"     — rise first then dip
#   "biphasic_ie"     — dip first then rise
make_synth_cell <- function(cell_name,
                            kind = c("flat", "excitation", "inhibition",
                                     "biphasic_ei", "biphasic_ie"),
                            assigned_subclass = "L23_IT",
                            Species = "Macaque",
                            t_step = 0.5,
                            tmax = 50) {
  kind <- match.arg(kind)
  t <- seq(-10, tmax, by = t_step)
  pct <- rep(100, length(t))

  switch(
    kind,
    flat        = NULL,
    excitation  = { pct[t >= 0 & t <= 30] <- 180 },
    inhibition  = { pct[t >= 0 & t <= 30] <- 30 },
    biphasic_ei = {
      pct[t >= 0  & t <= 10] <- 180
      pct[t >  10 & t <= 30] <- 40
    },
    biphasic_ie = {
      pct[t >= 0  & t <= 10] <- 40
      pct[t >  10 & t <= 30] <- 180
    }
  )

  data.frame(
    cell_name         = cell_name,
    time              = t,
    percent_change    = pct,
    assigned_subclass = assigned_subclass,
    Species           = Species,
    x_bin             = floor(t / 2) * 2,
    stringsAsFactors  = FALSE
  )
}

make_synth_df <- function() {
  do.call(rbind, list(
    make_synth_cell("c_flat", "flat"),
    make_synth_cell("c_exc",  "excitation"),
    make_synth_cell("c_inh",  "inhibition", assigned_subclass = "L5_IT"),
    make_synth_cell("c_ei",   "biphasic_ei"),
    make_synth_cell("c_ie",   "biphasic_ie", assigned_subclass = "L5_ET")
  ))
}

# Cell that looks like a strong excitation but has a brief silent period at
# puff onset (spiking shuts off for ~1.5 s, then rebounds high). Used to
# verify the raw zero-spiking detector flips this from "excitation" to a
# biphasic call.
make_synth_silent_then_exc <- function(cell_name = "c_silent_exc",
                                       t_step = 0.5) {
  t   <- seq(-10, 50, by = t_step)
  pct <- rep(100, length(t))
  pct[t >= 0  & t <  1.5] <- 0      # brief silence at puff
  pct[t >= 1.5 & t <= 30] <- 200    # large excitation
  data.frame(
    cell_name         = cell_name,
    time              = t,
    percent_change    = pct,
    assigned_subclass = "L23_IT",
    Species           = "Macaque",
    x_bin             = floor(t / 2) * 2,
    stringsAsFactors  = FALSE
  )
}

# Pathological pattern: a single hard-zero sample at puff onset bracketed
# by spiking and immediate high-rate rebound. With 2 s bin grids this is
# easy for run-length detectors to miss (a one-sample run has duration 0).
make_synth_single_zero_rebound <- function(cell_name = "c_one_zero",
                                           bin_step  = 2) {
  t <- seq(-10, 50, by = bin_step)
  pct <- rep(110, length(t))           # baseline-ish spiking
  pct[t == 0] <- 0                     # one hard zero at puff
  pct[t > 0 & t <= 20] <- 280          # high rebound
  pct[t > 20] <- 150
  data.frame(
    cell_name         = cell_name,
    time              = t,
    percent_change    = pct,
    assigned_subclass = "L23_IT",
    Species           = "Macaque",
    x_bin             = floor(t / 2) * 2,
    stringsAsFactors  = FALSE
  )
}
