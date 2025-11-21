# file: filter_and_plot.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import firwin, filtfilt

def filter_and_plot(signal, fs, cutoff=500.0, filter_type="lowpass"):
    """
    Zero-phase FIR filtering of a short voltage trace (e.g. 200 ms).
    Inputs:
        signal : list or numpy array (numeric vector)
        fs     : sampling rate (Hz)
    Returns:
        filtered : numpy array (filtered signal)
    """
    sig = np.array(signal, dtype=float)
    duration = len(sig) / fs
    t = np.arange(0, duration, 1/fs)

    # ---- filter design ----
    cutoff = 500.0        # Hz (adjustable)
    trans_width = 200.0   # Hz (wider => shorter filter)
    nyq = fs / 2
    N = int(np.ceil(4.0 * fs / trans_width))
    if N % 2 == 0:
        N += 1
    coef = firwin(N, cutoff/nyq, window='hamming')

    # ---- reflection padding ----
    padlen = max(3*len(coef), 100)
    left_pad = sig[1:padlen+1][::-1]
    right_pad = sig[-padlen-1:-1][::-1]
    sig_padded = np.concatenate([left_pad, sig, right_pad])

    # ---- zero-phase filtering ----
    filtered_padded = filtfilt(coef, [1.0], sig_padded, padlen=padlen)
    filtered = filtered_padded[padlen:-padlen]

    # ---- plotting ----
    plt.figure(figsize=(10, 6))

    # Full trace
    #plt.subplot(2, 1, 1)
    plt.plot(t*1000, sig, label='Raw', alpha=0.6)
    plt.plot(t*1000, filtered, label='Filtered', linewidth=2)
    plt.title('Voltage trace: Raw vs Filtered')
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    return filtered
