"""
sodium_tool.py
Estimate sodium channel availability dynamics (a(t)) from voltage traces.
Supports:
  • Fitting from FI + VI sweeps (fit_params_from_FI_VI)
  • Fitting from one long baseline sweep (fit_params_from_baseline)
  • Continuous baseline→drug inference (infer_availability_continuous_pair)
"""

import numpy as np
from scipy.optimize import minimize

# ============================================================
# Utilities
# ============================================================

def _to_1d(x):
    return np.ravel(np.asarray(x))

def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))

def _safe_spike_detect(v, fs, thresh=-20, refractory_samples=5):
    v = _to_1d(v)
    if v.size < 2:
        return np.array([], int)
    dv = np.diff(v, prepend=v[0])
    idx = np.where((v[1:] > thresh) & (dv[1:] > 0))[0].astype(int)
    if idx.size:
        keep = np.insert(np.diff(idx) > refractory_samples, 0, True)
        idx = idx[keep]
    return idx

def _extract_spike_features(v, fs, sidx):
    v = _to_1d(v); dt = 1.0 / fs
    feats = []
    for k, i in enumerate(sidx):
        win = slice(max(i-5,0), min(i+50, v.size))
        seg = v[win]
        Vp = float(np.max(seg))
        dVdt = np.diff(seg)/dt
        dVdt_max = float(np.max(dVdt)) if dVdt.size else np.nan
        half = 0.5*(Vp + seg[0])
        above = np.where(seg >= half)[0]
        width = float((above[-1]-above[0])*dt) if above.size>1 else 0.0
        ahp_win = slice(i, min(i+int(0.02*fs), v.size))
        AHP = float(np.min(v[ahp_win]) - seg[0]) if ahp_win.stop>ahp_win.start else np.nan
        ISI = float((sidx[k+1]-i)/fs) if k+1<len(sidx) else np.nan
        feats.append(dict(V_peak=Vp, dVdt_max=dVdt_max, width=width, AHP=AHP, ISI=ISI))
    return feats

def _stack_features_list(feats):
    arr=[]
    for f in feats:
        ISI = 0.0 if (np.isnan(f['ISI']) or f['ISI'] is None) else f['ISI']
        arr.append([f['V_peak'], f['width'], f['dVdt_max'], ISI])
    return np.asarray(arr,float) if arr else np.zeros((0,4),float)

def _compute_scaler(FI):
    X=[]
    for step in FI:
        F=_stack_features_list(step['spike_features'])
        if F.size: X.append(F)
    X=np.vstack(X) if X else np.zeros((0,4))
    med=np.nanmedian(X,axis=0) if X.size else np.zeros(4)
    mad=np.nanmedian(np.abs(X-med),axis=0) if X.size else np.ones(4)
    mad[mad==0]=1.0
    return med,mad

def _scale_feats(F,med,mad): return (F-med)/mad

# ============================================================
# Prepare FI and VI data
# ============================================================

def prepare_FI(voltage_traces, stim_values, fs, spike_thresh=-20):
    if isinstance(voltage_traces,np.ndarray) and voltage_traces.ndim==2:
        traces=[voltage_traces[:,i] for i in range(voltage_traces.shape[1])]
    else: traces=[np.asarray(tr) for tr in voltage_traces]
    FI=[]
    for i,v in enumerate(traces):
        v=_to_1d(v)
        sidx=_safe_spike_detect(v,fs,spike_thresh)
        feats=_extract_spike_features(v,fs,sidx)
        FI.append(dict(current=float(stim_values[i]),
                       duration=v.size/fs,
                       spike_times=sidx/fs,
                       spike_features=feats))
    return FI

def prepare_VI(voltage_traces, stim_values, fs, spike_thresh=-20):
    if isinstance(voltage_traces,np.ndarray) and voltage_traces.ndim==2:
        traces=[voltage_traces[:,i] for i in range(voltage_traces.shape[1])]
    else: traces=[np.asarray(tr) for tr in voltage_traces]
    rows=[]
    for i,v in enumerate(traces):
        sidx=_safe_spike_detect(v,fs,spike_thresh)
        feats=_extract_spike_features(v,fs,sidx)
        if feats:
            f0=feats[0]
            rows.append(dict(current=float(stim_values[i]),
                             V_peak=f0['V_peak'],
                             dVdt_max=f0['dVdt_max'],
                             width=f0['width']))
    return rows

# ============================================================
# Availability and fitting (FI + VI)
# ============================================================

def _availability_from_features(feats,w_vec,w0,tau,med,mad):
    a=1.0; seq=[]
    for f in feats:
        seq.append(a)
        F=_stack_features_list([f])[0]; Fz=_scale_feats(F,med,mad)
        d=_sigmoid(w0+np.dot(w_vec,Fz))
        a_post=a*(1-d)
        dt=0.01 if (np.isnan(F[3]) or F[3]<=0) else F[3]
        a=1-(1-a_post)*np.exp(-dt/tau)
    return np.asarray(seq)

def _loss(p,FI,VI,med,mad,lam_var=0.01,lam_vi=0.1):
    w0=p[0]; wv=p[1:5]; tau=p[5]
    loss=0.0
    for step in FI:
        f=step['spike_features']
        if not f: continue
        y=np.array([s['dVdt_max'] for s in f])
        a=_availability_from_features(f,wv,w0,tau,med,mad)
        A=np.vstack([a,np.ones_like(a)]).T
        alpha,beta=np.linalg.lstsq(A,y,rcond=None)[0]
        yhat=alpha*a+beta
        loss+=np.nansum((yhat-y)**2)+lam_var*np.var(a)
    if VI and lam_vi>0:
        rows=[r for r in VI if not np.isnan(r['dVdt_max'])]
        if len(rows)>2:
            cur=np.array([r['current'] for r in rows])
            dv=np.array([r['dVdt_max'] for r in rows])
            wid=np.array([r['width'] for r in rows])
            def _c(x,y):
                return 1-np.corrcoef(x,y)[0,1] if np.std(x)>0 and np.std(y)>0 else 1
            loss+=lam_vi*(_c(cur,dv)+_c(cur,-wid))
    return float(loss)

def fit_params_from_FI_VI(FI,VI=None,initial=None):
    med,mad=_compute_scaler(FI)
    if initial is None:
        initial=np.array([-1.5,0.5,-0.5,0.5,-0.5,0.05])
    bnds=[(-5,5),(-3,3),(-3,3),(-3,3),(-3,3),(0.003,0.5)]
    res=minimize(lambda p:_loss(p,FI,VI,med,mad),initial,method="L-BFGS-B",bounds=bnds)
    return dict(success=res.success,loss=res.fun,
                w0=res.x[0],w_vec=res.x[1:5],tau_rec=res.x[5],
                scaler_med=med,scaler_mad=mad)

# ============================================================
# Fit from one long baseline trace
# ============================================================

def _availability_from_features_single(feats,wv,w0,tau,med,mad):
    a=1.0; seq=[]
    for f in feats:
        seq.append(a)
        F=_stack_features_list([f])[0]; Fz=_scale_feats(F,med,mad)
        d=_sigmoid(w0+np.dot(wv,Fz))
        a_post=a*(1-d)
        dt=0.01 if (np.isnan(F[3]) or F[3]<=0) else F[3]
        a=1-(1-a_post)*np.exp(-dt/tau)
    return np.asarray(seq)

def fit_params_from_baseline(voltage_baseline,fs,spike_thresh=-20,
                             initial=None,lambda_var=0.02,lambda_l2=0.001):
    v=_to_1d(voltage_baseline)
    sidx=_safe_spike_detect(v,fs,spike_thresh)
    feats=_extract_spike_features(v,fs,sidx)
    if len(feats)<5:
        return dict(success=False,msg="Too few spikes",loss=np.nan)
    X=_stack_features_list(feats)
    med=np.nanmedian(X,0); mad=np.nanmedian(np.abs(X-med),0); mad[mad==0]=1
    y=np.array([f['dVdt_max'] for f in feats])
    if initial is None:
        initial=np.array([-1.5,0.6,-0.6,0.6,-0.6,0.05])
    bnds=[(-5,5),(-3,3),(-3,3),(-3,3),(-3,3),(0.003,0.5)]

    def loss(p):
        w0=p[0]; wv=p[1:5]; tau=p[5]
        a=_availability_from_features_single(feats,wv,w0,tau,med,mad)
        A=np.vstack([a,np.ones_like(a)]).T
        alpha,beta=np.linalg.lstsq(A,y,rcond=None)[0]
        yhat=alpha*a+beta
        data=np.sum((yhat-y)**2)
        var_pen=-lambda_var*np.var(a)
        l2=lambda_l2*(w0**2+np.sum(wv**2))
        return float(data + l2 - var_pen)

    res=minimize(loss,initial,method="L-BFGS-B",bounds=bnds)
    return dict(success=res.success,loss=res.fun,
                w0=res.x[0],w_vec=res.x[1:5],tau_rec=res.x[5],
                scaler_med=med,scaler_mad=mad)

# ============================================================
# Continuous baseline→drug inference
# ============================================================

def infer_availability_continuous_pair(v_base, v_drug, fitted,
                                       fs_base, fs_drug, spike_thresh=-20):
    """
    Continuous Na availability across Baseline -> Drug with *exact timing*.
    Uses per-segment sampling rates; no artificial gap; dt between spikes is
    taken from absolute spike times (seconds), not from segment-local ISI.
    """

    # --- Unpack params and scalers
    w0  = float(fitted['w0'])
    wv  = np.asarray(fitted['w_vec'], float)
    tau = float(fitted['tau_rec'])
    med = np.asarray(fitted['scaler_med'], float)
    mad = np.asarray(fitted['scaler_mad'], float)

    # --- Vectors and time axes
    v_base = np.ravel(v_base);  v_drug = np.ravel(v_drug)
    t_base = np.arange(v_base.size) / float(fs_base)
    t0_drug = t_base[-1] + (1.0 / float(fs_base)) if v_base.size > 0 else 0.0
    t_drug = t0_drug + (np.arange(v_drug.size) / float(fs_drug))

    # --- Spike detection per segment (own fs)
    sidx_b = _safe_spike_detect(v_base, fs_base, spike_thresh)
    sidx_d = _safe_spike_detect(v_drug, fs_drug, spike_thresh)

    # Absolute spike times (seconds)
    tspk_b = sidx_b / float(fs_base)
    tspk_d = t0_drug + (sidx_d / float(fs_drug))

    # --- Features per segment (use local fs for shape features)
    feats_b = _extract_spike_features(v_base, fs_base, sidx_b)
    feats_d = _extract_spike_features(v_drug, fs_drug, sidx_d)

    # --- Concatenate spikes by true time
    tspk_all = np.concatenate([tspk_b, tspk_d])
    src_flags = np.concatenate([np.zeros_like(tspk_b, dtype=int),  # 0=base
                                np.ones_like(tspk_d, dtype=int)])  # 1=drug
    # indices to sort by time
    order = np.argsort(tspk_all)
    tspk_all = tspk_all[order]
    src_flags = src_flags[order]

    # helper to fetch the matching feature dict in time order
    def get_feat(k_sorted):
        if src_flags[k_sorted] == 0:
            # baseline
            # where am I among baseline spikes?
            i = np.sum((order[:k_sorted] < len(tspk_b)) & (src_flags[order[:k_sorted]] == 0))
            return feats_b[i]
        else:
            # drug
            i = np.sum(src_flags[order[:k_sorted]] == 1)
            return feats_d[i]

    # --- Build full time vector for continuous rendering
    t_all = np.concatenate([t_base, t_drug])
    a_cont = np.ones_like(t_all, dtype=float)

    # --- Propagate availability using true dt between spikes
    a_prev = 1.0
    a_post = 1.0
    a_spike_all = []
    last_t = t_all[0] if t_all.size else 0.0

    for k in range(tspk_all.size):
        t_k = tspk_all[k]
        a_start = a_post if k > 0 else a_prev

        # fill recovery from last_t .. t_k
        seg = (t_all >= last_t) & (t_all < t_k)
        if np.any(seg):
            a_cont[seg] = 1.0 - (1.0 - a_start) * np.exp(-(t_all[seg] - last_t) / tau)

        # depletion at spike k using FEATS from correct segment
        f = get_feat(k)
        F = _stack_features_list([f])[0].copy()

        # overwrite ISI with the *true* dt (seconds) since previous spike
        if k == 0:
            dt_true = max(0.01, (tspk_all[0] - t_all[0]))  # small positive
        else:
            dt_true = max(0.001, tspk_all[k] - tspk_all[k-1])
        F[3] = dt_true

        Fz = _scale_feats(F, med, mad)
        d_k = _sigmoid(w0 + np.dot(wv, Fz))
        a_post = a_start * (1.0 - d_k)

        a_spike_all.append(a_start)
        last_t = t_k

    # tail segment to the end
    if t_all.size:
        seg = (t_all >= last_t)
        a_cont[seg] = 1.0 - (1.0 - a_post) * np.exp(-(t_all[seg] - last_t) / tau)

    # --- Split back to baseline / drug outputs with *correct* times
    mask_base = t_all <= (t_base[-1] if t_base.size else -np.inf)
    a_cont_base = a_cont[mask_base]
    a_cont_drug = a_cont[~mask_base]

    tspk_base = tspk_all[tspk_all <= (t_base[-1] if t_base.size else -np.inf)]
    tspk_drug = tspk_all[tspk_all >  (t_base[-1] if t_base.size else -np.inf)] - t0_drug
    a_spk_base = np.asarray(a_spike_all[:tspk_base.size])
    a_spk_drug = np.asarray(a_spike_all[tspk_base.size:])

    return dict(
        time_base=t_base,
        time_drug=t_drug - t0_drug,       # start drug at 0 s for plotting
        a_cont_base=a_cont_base,
        a_cont_drug=a_cont_drug,
        a_spike_base=a_spk_base,
        a_spike_drug=a_spk_drug,
        spike_times_base=tspk_base,
        spike_times_drug=tspk_drug
    )


# ============================================================
# Fit parameters automatically from first 5 s of drug trace
# ============================================================

def fit_params_from_drug_segment(v_drug, fs_drug, spike_thresh=-20,
                                 window_s=5.0,
                                 initial=None,
                                 lambda_var=0.02, lambda_l2=0.001):
    """
    Estimate (w0, w_vec, tau_rec) automatically from the first `window_s`
    seconds of the drug trace.

    Parameters
    ----------
    v_drug : 1D array
        Drug‐condition voltage trace.
    fs_drug : float
        Sampling rate of the drug trace (Hz).
    spike_thresh : float
        Spike detection threshold in mV.
    window_s : float
        Time window in seconds to analyze (default 5 s).
    initial, lambda_var, lambda_l2 :
        Same meaning as in fit_params_from_baseline.

    Returns
    -------
    dict with fitted parameters and scaling:
        {'w0', 'w_vec', 'tau_rec', 'scaler_med', 'scaler_mad', 'loss', 'success'}
    """

    v = _to_1d(v_drug)
    n_samp = int(window_s * fs_drug)
    if n_samp < 10:
        raise ValueError("Drug trace too short for 5-second window.")
    v_seg = v[:n_samp]

    sidx = _safe_spike_detect(v_seg, fs_drug, spike_thresh)
    feats = _extract_spike_features(v_seg, fs_drug, sidx)

    if len(feats) < 5:
        return dict(success=False, msg="Too few spikes in first 5 s", loss=np.nan)

    X = _stack_features_list(feats)
    med = np.nanmedian(X, axis=0)
    mad = np.nanmedian(np.abs(X - med), axis=0)
    mad[mad == 0] = 1.0

    y = np.array([f['dVdt_max'] for f in feats], float)

    if initial is None:
        initial = np.array([-1.5, 0.6, -0.6, 0.6, -0.6, 0.05], float)

    bounds = [(-5, 5), (-3, 3), (-3, 3), (-3, 3), (-3, 3), (0.003, 0.5)]

    def loss_fn(p):
        w0 = p[0]; wv = p[1:5]; tau = p[5]
        a = _availability_from_features_single(feats, wv, w0, tau, med, mad)
        A = np.vstack([a, np.ones_like(a)]).T
        alpha, beta = np.linalg.lstsq(A, y, rcond=None)[0]
        y_hat = alpha * a + beta
        data = np.sum((y_hat - y) ** 2)
        var_pen = -lambda_var * np.var(a)
        l2_pen = lambda_l2 * (w0 ** 2 + np.sum(wv ** 2))
        return float(data + l2_pen - var_pen)

    res = minimize(loss_fn, initial, method="L-BFGS-B", bounds=bounds)

    return dict(
        success=bool(res.success),
        loss=float(res.fun),
        w0=float(res.x[0]),
        w_vec=res.x[1:5].astype(float),
        tau_rec=float(res.x[5]),
        scaler_med=med.astype(float),
        scaler_mad=mad.astype(float),
    )

# ============================================================
# Inference on a single drug trace only
# ============================================================

def infer_availability_drug(v_drug, fitted, fs_drug, spike_thresh=-20):
    """
    Infer sodium availability on a single drug trace only (no baseline).

    Parameters
    ----------
    v_drug : 1D array
        Voltage trace during the drug condition.
    fitted : dict
        Fitted parameters (e.g., from fit_params_from_baseline or
        fit_params_from_drug_segment).
    fs_drug : float
        Sampling rate of the drug trace in Hz.
    spike_thresh : float
        Threshold for spike detection in mV.

    Returns
    -------
    dict with:
      time_drug : array of time (s)
      a_cont_drug : continuous sodium availability a(t)
      a_spike_drug : sodium availability before each spike
      spike_times_drug : spike times (s)
    """

    # Extract parameters
    w0 = fitted['w0']
    wv = np.asarray(fitted['w_vec'])
    tau = fitted['tau_rec']
    med = np.asarray(fitted['scaler_med'])
    mad = np.asarray(fitted['scaler_mad'])

    v = np.ravel(v_drug)
    fs = float(fs_drug)
    t = np.arange(len(v)) / fs

    # Spike detection and feature extraction
    sidx = _safe_spike_detect(v, fs, spike_thresh)
    feats = _extract_spike_features(v, fs, sidx)

    # Initialize availability state
    a_pre = 1.0
    a_post = 1.0
    a_spk = []
    a_cont = np.ones_like(t)
    last_t = 0.0

    # Iterate through spikes and propagate depletion/recovery
    for k, i in enumerate(sidx):
        t_k = i / fs
        a_start = a_post if k > 0 else a_pre
        mask = (t >= last_t) & (t < t_k)
        a_cont[mask] = 1 - (1 - a_start) * np.exp(-(t[mask] - last_t) / tau)

        F = _stack_features_list([feats[k]])[0]
        Fz = _scale_feats(F, med, mad)
        d = _sigmoid(w0 + np.dot(wv, Fz))
        a_post = a_start * (1 - d)
        a_spk.append(a_start)
        last_t = t_k

    # Tail recovery after the last spike
    mask = (t >= last_t)
    a_cont[mask] = 1 - (1 - a_post) * np.exp(-(t[mask] - last_t) / tau)

    return dict(
        time_drug=t,
        a_cont_drug=a_cont,
        a_spike_drug=np.array(a_spk),
        spike_times_drug=np.array(sidx) / fs,
    )

# sodium_tool.py

import numpy as np
from scipy.optimize import minimize
from scipy.signal import find_peaks


def detect_spikes(voltage, thresh=-20):
    """Detect spike indices from a voltage trace."""
    voltage = np.asarray(voltage)
    peaks, _ = find_peaks(voltage, height=thresh)
    return peaks


def availability_model(isi, w_vec, w0, tau_rec):
    """Compute sodium availability a(t) for a given ISI sequence."""
    a = np.zeros(len(isi) + 1)
    a[0] = 1.0
    for i in range(1, len(a)):
        dt = isi[i - 1]
        a[i] = a[i - 1] + (1 - a[i - 1]) * (1 - np.exp(-dt / tau_rec)) - w0 * w_vec * a[i - 1]
        a[i] = np.clip(a[i], 0, 1)
    return a


def loss_fn(params, isi, amplitudes):
    """Loss = RMSE between predicted availability and normalized amplitude proxy."""
    w_vec, w0, tau_rec = params
    if tau_rec <= 0 or w_vec < 0 or w0 < 0:
        return np.inf

    a = availability_model(isi, w_vec, w0, tau_rec)
    pred = a[:-1]
    n = min(len(pred), len(amplitudes))
    return np.sqrt(np.mean((pred[:n] - amplitudes[:n]) ** 2))


def fit_params_base_drug(v_base, v_drug, fs_base, fs_drug, spike_thresh=-20):
    """
    Fit sodium-availability parameters using spikes from baseline and drug traces only.
    Returns dictionary with w_vec, w0, tau_rec, scaler_med, scaler_std.
    """

    # --- Combine traces ---
    v_base = np.asarray(v_base).ravel()
    v_drug = np.asarray(v_drug).ravel()
    v_combined = np.concatenate([v_base, v_drug])
    fs = (fs_base + fs_drug) / 2.0

    # --- Spike detection ---
    spikes = detect_spikes(v_combined, thresh=spike_thresh)
    if len(spikes) < 3:
        raise ValueError("Too few spikes detected for fitting.")

    isi = np.diff(spikes) / fs

    # --- Simple amplitude proxy from ISIs (shorter ISI → less availability) ---
    amplitudes = np.exp(-isi / np.mean(isi))
    amplitudes /= np.max(amplitudes)

    # --- Fit parameters ---
    x0 = [0.1, 0.05, 100.0]  # [w_vec, w0, tau_rec]
    bounds = [(1e-5, 5.0), (1e-5, 1.0), (1.0, 2000.0)]

    result = minimize(loss_fn, x0, args=(isi, amplitudes),
                      bounds=bounds, method="L-BFGS-B")

    if not result.success:
        print("Warning: fitting did not converge:", result.message)

    w_vec, w0, tau_rec = result.x

    # --- Compute scale factors to satisfy infer_availability_continuous_pair() ---
    scaler_med = np.median(amplitudes)
    scaler_std = np.std(amplitudes)
    scaler_mad = np.median(np.abs(amplitudes - scaler_med))
    
    return {
    "w_vec": float(np.ravel(w_vec)[0]),
    "w0": float(np.ravel(w0)[0]),
    "tau_rec": float(np.ravel(tau_rec)[0]),
    "scaler_med": float(np.ravel(scaler_med)[0]),
    "scaler_std": float(np.ravel(scaler_std)[0]),
    "scaler_mad": float(np.ravel(scaler_mad)[0])
}
