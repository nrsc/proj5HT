import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO

def _get_sweep_table_df(nwbfile):
    if not hasattr(nwbfile, "sweep_table") or nwbfile.sweep_table is None:
        return pd.DataFrame()

    st = nwbfile.sweep_table.to_dataframe().copy()
    if "sweep_number" not in st.columns:
        st = st.reset_index().rename(columns={"index": "sweep_number"})

    def _series_name(obj):
        try:
            return obj.name
        except Exception:
            return None

    # Common columns vary by producer; probe
    acq_candidates  = ["series", "acquisition", "acquisition_series", "response"]
    stim_candidates = ["stimulus", "stimulus_series"]

    acq_name = None
    for c in acq_candidates:
        if c in st.columns:
            acq_name = st[c].apply(_series_name)
            break
    if acq_name is None:
        acq_name = pd.Series([None] * len(st))
    st["acquisition_name"] = acq_name

    stim_name = None
    for c in stim_candidates:
        if c in st.columns:
            stim_name = st[c].apply(_series_name)
            break
    if stim_name is None:
        stim_name = pd.Series([None] * len(st))
    st["stimulus_name"] = stim_name

    return st


def _find_pulse_intervals(stim, sr,
                         amp_thresh=None,
                         min_dur_ms=50,
                         max_dur_ms=200,
                         merge_gap_ms=5):
    """
    Detect contiguous intervals in the stimulus vector where |stim| exceeds threshold.
    Returns list of (start_idx, stop_idx) inclusive, plus start/stop times in seconds.
    """

    stim = np.asarray(stim, dtype=np.float32)
    n = stim.size
    if n == 0:
        return []

    # choose threshold if not given (robust): 10x MAD above baseline noise
    if amp_thresh is None:
        med = np.median(stim)
        mad = np.median(np.abs(stim - med)) + 1e-12
        amp_thresh = float(10.0 * mad)

        # if stim is nearly all zeros, MAD ~0; make a small fallback
        if amp_thresh < 1e-9:
            amp_thresh = float(np.max(np.abs(stim)) * 0.2)

    above = np.abs(stim) >= amp_thresh
    if not np.any(above):
        return []

    # rising/falling edges
    idx = np.flatnonzero(above)
    # group contiguous indices
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.insert(idx[breaks + 1], 0, idx[0])
    stops  = np.append(idx[breaks], idx[-1])

    # merge nearby groups (small gaps)
    merge_gap = int(round((merge_gap_ms / 1000.0) * sr))
    merged = []
    for s, e in zip(starts, stops):
        if not merged:
            merged.append([int(s), int(e)])
        else:
            if s - merged[-1][1] <= merge_gap:
                merged[-1][1] = int(e)
            else:
                merged.append([int(s), int(e)])

    # duration filter
    min_dur = (min_dur_ms / 1000.0) * sr
    max_dur = (max_dur_ms / 1000.0) * sr

    out = []
    for s, e in merged:
        dur = (e - s + 1)
        if dur >= min_dur and dur <= max_dur:
            out.append((s, e))

    return out


def get_tp_testpulse_epochs(nwb_path,
                            sweep_numbers,
                            amp_thresh=None,
                            min_dur_ms=70,
                            max_dur_ms=140,
                            merge_gap_ms=5,
                            pad_ms=5):
    """
    For each sweep_number:
      - resolve stimulus + acquisition series via sweep_table
      - detect 100ms test pulse intervals from stimulus series
      - return an epoch-like table: start/stop times (s) relative to sweep

    Returns:
      dict:
        'pulses' : pd.DataFrame with one row per pulse interval
        'sweeps' : pd.DataFrame sweep metadata (names, rates, n_samples)
    """

    sweep_numbers = [int(x) for x in sweep_numbers]
    with NWBHDF5IO(nwb_path, "r") as io:
        nwbfile = io.read()
        st = _get_sweep_table_df(nwbfile)
        if st.empty:
            raise ValueError("NWB missing sweep_table; cannot map sweep_number to series.")

        st2 = st[st["sweep_number"].astype(int).isin(set(sweep_numbers))].copy()

        rows = []
        sweep_rows = []

        for _, sw in st2.iterrows():
            sn = int(sw["sweep_number"])
            acq_name = sw.get("acquisition_name", None)
            stim_name = sw.get("stimulus_name", None)

            if not acq_name or acq_name not in nwbfile.acquisition:
                continue
            if not stim_name or stim_name not in nwbfile.stimulus:
                # some files store stimulus elsewhere; if so, you can extend this lookup
                continue

            acq = nwbfile.acquisition[acq_name]
            stim = nwbfile.stimulus[stim_name]

            sr = float(stim.rate) if hasattr(stim, "rate") else float(acq.rate)
            stim_vec = np.asarray(stim.data[:], dtype=np.float32)

            intervals = _find_pulse_intervals(
                stim_vec, sr,
                amp_thresh=amp_thresh,
                min_dur_ms=min_dur_ms,
                max_dur_ms=max_dur_ms,
                merge_gap_ms=merge_gap_ms
            )

            pad = float(pad_ms) / 1000.0

            # sweep metadata
            sweep_rows.append({
                "sweep_number": sn,
                "acquisition_name": acq_name,
                "stimulus_name": stim_name,
                "sampling_rate": float(acq.rate),
                "stim_rate": sr,
                "starting_time": float(acq.starting_time),
                "n_samples": int(len(acq.data)),
                "unit": getattr(acq, "unit", None)
            })

            # pulse rows
            for (s, e) in intervals:
                t0 = s / sr
                t1 = (e + 1) / sr
                rows.append({
                    "sweep_number": sn,
                    "acquisition_name": acq_name,
                    "stimulus_name": stim_name,
                    "pulse_start_s": float(t0),
                    "pulse_stop_s": float(t1),
                    "pulse_start_s_padded": float(max(0.0, t0 - pad)),
                    "pulse_stop_s_padded": float(t1 + pad),
                    "pulse_dur_ms": float((t1 - t0) * 1000.0)
                })

        pulses_df = pd.DataFrame(rows)
        sweeps_df = pd.DataFrame(sweep_rows)

    return {"pulses": pulses_df, "sweeps": sweeps_df}


def make_tp_mask_for_sweep(v_len, sr, pulse_intervals_s, invert=False):
    """
    Build boolean mask for samples.
    pulse_intervals_s: array-like of (start_s, stop_s) relative to sweep.
    invert=False -> mask is True for "keep" (excluding pulse windows)
    invert=True  -> mask is True for "pulse windows"
    """
    mask = np.ones(v_len, dtype=bool)
    for (a, b) in pulse_intervals_s:
        i0 = int(np.floor(a * sr))
        i1 = int(np.ceil(b * sr))
        i0 = max(0, i0)
        i1 = min(v_len, i1)
        mask[i0:i1] = False
    return ~mask if invert else mask
