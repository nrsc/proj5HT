import numpy as np
from scipy.optimize import minimize

# --- Helper functions ---
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def predict_features(features_list, w_vec, w0, tau_rec):
    """Predict sodium availability before each spike"""
    a_pre = 1.0
    a_trace = []
    for k, feat in enumerate(features_list):
        a_trace.append(a_pre)
        F = np.array([
            feat['V_peak'],
            feat['width'],
            feat['dVdt_max'],
            0 if np.isnan(feat['ISI']) else feat['ISI']
        ])
        d_k = sigmoid(w0 + np.dot(w_vec, F))
        a_post = a_pre * (1 - d_k)
        dt = F[3] if F[3] > 0 else 0.01
        a_pre = 1 - (1 - a_post) * np.exp(-dt / tau_rec)
    return np.array(a_trace)

def loss_fn(params, FI_data):
    w0 = params[0]
    w_vec = params[1:5]
    tau_rec = params[5]

    total_loss = 0
    for step in FI_data:
        features_list = step['spike_features']
        if len(features_list) == 0:
            continue
        observed_dvdt = np.array([f['dVdt_max'] for f in features_list])
        predicted_a = predict_features(features_list, w_vec, w0, tau_rec)
        # simple linear mapping
        alpha, beta = 1.0, 0.0
        predicted_dvdt = alpha * predicted_a + beta
        total_loss += np.sum((predicted_dvdt - observed_dvdt) ** 2)
        # optional FI curve contribution
        observed_rate = len(features_list) / step['duration']
        predicted_rate = observed_rate * np.mean(predicted_a)
        total_loss += 0.1 * (predicted_rate - observed_rate) ** 2
    return total_loss

# --- FI data preparation ---
def prepare_FI_data(voltage_traces, stim_values, fs, spike_thresh=-20):
    """
    Convert column-major voltage traces to FI_data format.
    Robust to single-column matrices, empty traces, or no spikes.
    """
    FI_data = []

    # Convert 2D array to list of columns
    if isinstance(voltage_traces, np.ndarray) and voltage_traces.ndim == 2:
        traces_list = [voltage_traces[:, i] for i in range(voltage_traces.shape[1])]
    else:
        traces_list = [np.asarray(tr) for tr in voltage_traces]

    for i, voltage in enumerate(traces_list):
        current = stim_values[i]
        voltage = np.ravel(voltage)  # ensure 1D
        duration = len(voltage) / fs

        if len(voltage) < 2:
            FI_data.append({'current': current,
                            'duration': duration,
                            'spike_times': np.array([]),
                            'spike_features': []})
            continue

        # --- Spike detection ---
        dv = np.diff(voltage, prepend=voltage[0])
        spike_idx = np.where((voltage[1:] > spike_thresh) & (dv[1:] > 0))[0]
        spike_idx = spike_idx.astype(int)

        if len(spike_idx) > 0:
            # Remove duplicates within the same spike
            spike_idx = spike_idx[np.insert(np.diff(spike_idx) > 5, 0, True)]
        else:
            spike_idx = np.array([], dtype=int)

        # --- Extract features ---
        features_list = []
        dt = 1 / fs
        for k, idx in enumerate(spike_idx):
            win = slice(max(idx - 5, 0), min(idx + 50, len(voltage)))
            v_spike = voltage[win]

            V_peak = np.max(v_spike)
            dVdt = np.diff(v_spike) / dt
            dVdt_max = np.max(dVdt) if len(dVdt) > 0 else np.nan

            half_max = (V_peak + v_spike[0]) / 2
            above_half = np.where(v_spike >= half_max)[0]
            width = (above_half[-1] - above_half[0]) * dt if len(above_half) > 1 else 0

            ahp_window = slice(idx, min(idx + int(0.02 * fs), len(voltage)))
            AHP = np.min(voltage[ahp_window]) - v_spike[0] if len(voltage[ahp_window]) > 0 else np.nan

            ISI = (spike_idx[k + 1] - idx) / fs if k + 1 < len(spike_idx) else np.nan

            features_list.append({'V_peak': V_peak,
                                  'dVdt_max': dVdt_max,
                                  'width': width,
                                  'AHP': AHP,
                                  'ISI': ISI})

        FI_data.append({'current': current,
                        'duration': duration,
                        'spike_times': spike_idx / fs,
                        'spike_features': features_list})

    return FI_data

# --- Full fitting pipeline ---
def fit_sodium_from_traces(voltage_traces, stim_values, fs, spike_thresh=-20, initial_guess=None):
    """
    Full pipeline: generate FI_data from column-major traces,
    then fit w_vec, w0, tau_rec.
    """
    FI_data = prepare_FI_data(voltage_traces, stim_values, fs, spike_thresh)

    if initial_guess is None:
        initial_guess = np.array([-2, 0.01, -0.5, 0.001, -0.1, 0.05])

    bounds = [(-10, 10),
              (-5, 5), (-5, 5), (-5, 5), (-5, 5),
              (0.001, 1)]

    res = minimize(loss_fn, initial_guess, args=(FI_data,), bounds=bounds, method='L-BFGS-B')

    fitted_params = {
        'w0': res.x[0],
        'w_vec': res.x[1:5],
        'tau_rec': res.x[5],
        'success': res.success,
        'loss': res.fun
    }
    return fitted_params
