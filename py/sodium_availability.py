import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

def detect_spikes(voltage, fs, thresh=-20):
    """
    Detect spike times in a voltage trace.
    voltage: array of voltage in mV
    fs: sampling frequency in Hz
    thresh: threshold for spike detection (mV)
    
    Returns: spike indices (samples)
    """
    dv = np.diff(voltage, prepend=voltage[0])
    spike_idx = np.where((voltage[1:] > thresh) & (dv[1:] > 0))[0]
    # Remove duplicates (multiple points per spike)
    spike_idx = spike_idx[np.insert(np.diff(spike_idx) > 5, 0, True)]
    return spike_idx

def extract_spike_features(voltage, spike_idx, fs):
    """
    Extract spike features:
    - V_peak
    - dV/dt_max
    - width at half-max
    - AHP
    - ISI
    
    Returns: list of dicts
    """
    voltage = np.asarray(voltage)
    spike_idx = np.asarray(spike_idx, dtype=int) 
    
    features_list = []
    dt = 1/fs
    for i, idx in enumerate(spike_idx):
        # Define window around spike
        win = slice(max(idx-5,0), min(idx+50, len(voltage)))
        v_spike = voltage[win]
        t_spike = np.arange(len(v_spike))*dt

        V_peak = np.max(v_spike)
        dVdt = np.diff(v_spike)/dt
        dVdt_max = np.max(dVdt)
        
        # Spike width at half-height
        half_max = (V_peak + v_spike[0])/2
        above_half = np.where(v_spike >= half_max)[0]
        width = (above_half[-1]-above_half[0])*dt if len(above_half)>1 else 0
        
        # Afterhyperpolarization (min after spike within 20 ms)
        ahp_window = slice(idx, min(idx + int(0.02*fs), len(voltage)))
        AHP = np.min(voltage[ahp_window]) - v_spike[0]
        
        # ISI
        ISI = (spike_idx[i+1]-idx)/fs if i+1<len(spike_idx) else np.nan
        
        features_list.append({'V_peak': V_peak,
                              'dVdt_max': dVdt_max,
                              'width': width,
                              'AHP': AHP,
                              'ISI': ISI})
    return features_list



def sigmoid(x):
    return 1/(1+np.exp(-x))

def infer_sodium_availability(features_list, w_vec, w0=0, tau_rec=0.05):
    """
    Infer sodium availability from spike features
    features_list: list of dicts from extract_spike_features
    w_vec: weights for [V_peak, width, dVdt_max, ISI]
    w0: bias
    tau_rec: recovery time constant in seconds
    
    Returns: list of a_pre (availability before each spike)
    """
    a_pre = 1.0
    a_trace = []
    Vref = 0  # reference for normalization
    
    for k, feat in enumerate(features_list):
        a_trace.append(a_pre)
        F = np.array([feat['V_peak']-Vref, feat['width'], feat['dVdt_max'], 0 if np.isnan(feat['ISI']) else feat['ISI']])
        d_k = sigmoid(w0 + np.dot(w_vec, F))
        a_post = a_pre * (1 - d_k)  # depletion
        # recovery until next spike
        dt = F[3] if F[3]>0 else 0.01  # small dt if last spike
        a_pre = 1 - (1 - a_post) * np.exp(-dt / tau_rec)
    return np.array(a_trace)

def infer_sodium_long_trace(voltage_long, fs, fitted_params, spike_thresh=-20):
    """
    Infer sodium availability in a long trace using parameters from FI traces.
    
    Parameters
    ----------
    voltage_long : 1D array
        Voltage trace of the long recording.
    fs : float
        Sampling frequency in Hz.
    fitted_params : dict
        Output from fit_sodium_from_traces:
        {'w0': ..., 'w_vec': ..., 'tau_rec': ...}
    spike_thresh : float
        Threshold for spike detection in mV.
        
    Returns
    -------
    spike_times : array
        Times of spikes in seconds.
    a_spike : array
        Sodium availability before each spike.
    a_continuous : array
        Continuous sodium availability over entire trace (length = voltage_long).
    """
    
    # Extract fitted parameters
    w0 = fitted_params['w0']
    w_vec = fitted_params['w_vec']
    tau_rec = fitted_params['tau_rec']
    
    voltage = np.ravel(voltage_long)  # ensure 1D
    if len(voltage) < 2:
        return np.array([]), np.array([]), np.ones_like(voltage)
    
    # --- Spike detection ---
    dv = np.diff(voltage, prepend=voltage[0])
    spike_idx = np.where((voltage[1:] > spike_thresh) & (dv[1:] > 0))[0]
    spike_idx = spike_idx.astype(int)
    
    if len(spike_idx) > 0:
        spike_idx = spike_idx[np.insert(np.diff(spike_idx) > 5, 0, True)]
    else:
        spike_idx = np.array([], dtype=int)
    
    spike_times = spike_idx / fs
    
    # --- Extract spike features ---
    dt = 1/fs
    features_list = []
    for k, idx in enumerate(spike_idx):
        win = slice(max(idx-5,0), min(idx+50, len(voltage)))
        v_spike = voltage[win]
        V_peak = np.max(v_spike)
        dVdt = np.diff(v_spike)/dt
        dVdt_max = np.max(dVdt) if len(dVdt) > 0 else np.nan
        half_max = (V_peak + v_spike[0])/2
        above_half = np.where(v_spike >= half_max)[0]
        width = (above_half[-1]-above_half[0])*dt if len(above_half) > 1 else 0
        ahp_window = slice(idx, min(idx + int(0.02*fs), len(voltage)))
        AHP = np.min(voltage[ahp_window]) - v_spike[0] if len(voltage[ahp_window]) > 0 else np.nan
        ISI = (spike_idx[k+1]-idx)/fs if k+1 < len(spike_idx) else np.nan
        features_list.append({'V_peak': V_peak,
                              'dVdt_max': dVdt_max,
                              'width': width,
                              'AHP': AHP,
                              'ISI': ISI})
    
    # --- Predict sodium availability before each spike ---
    def predict_features(features_list, w_vec, w0, tau_rec):
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
            d_k = 1/(1+np.exp(-(w0 + np.dot(w_vec,F))))
            a_post = a_pre*(1-d_k)
            dt = F[3] if F[3] > 0 else 0.01
            a_pre = 1 - (1 - a_post)*np.exp(-dt/tau_rec)
        return np.array(a_trace)
    
    a_spike = predict_features(features_list, w_vec, w0, tau_rec)
    
    # --- Continuous interpolation ---
    a_continuous = np.ones_like(voltage)
    if len(spike_times) > 0:
        a_continuous = np.interp(np.arange(len(voltage))/fs, spike_times, a_spike, left=1.0, right=a_spike[-1])
    
    return spike_times, a_spike, a_continuous
