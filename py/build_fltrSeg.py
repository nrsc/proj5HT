import numpy as np
from scipy.signal import firwin, filtfilt
import matplotlib.pyplot as plt

# Example parameters (adjust to your data)
fs = 50000.0          # sampling rate in Hz (change to your fs)
cutoff = 500.0        # lowpass cutoff in Hz
trans_width = 200.0   # transition width in Hz (wider => shorter filter)
nyq = 0.5 * fs

# Design a linear-phase FIR (Hamming-windowed)
# Approximate length: N ~= 4 / (trans_width / fs)  (rule-of-thumb)
N = int(np.ceil(4.0 * fs / trans_width))
if N % 2 == 0:
    N += 1           # make odd for symmetric linear phase
coef = firwin(N, cutoff/nyq, window='hamming')


# Example signal (replace with your 200 ms voltage array)
sig = r.isi_seg
t = np.arange(0, len(sig), 1/fs)   # 200 ms
#sig = np.sin(2*np.pi*50*t) + 0.2*np.random.randn(len(t))  # demo signal




# Padding by reflection
padlen = max(3*len(coef), 100)   # pad len in samples (safety)
left_pad = sig[1:padlen+1][::-1]   # reflect left
right_pad = sig[-padlen-1:-1][::-1]# reflect right
sig_padded = np.concatenate([left_pad, sig, right_pad])

# Zero-phase filtering
filtered_padded = filtfilt(coef, [1.0], sig_padded, padlen=padlen)

# Remove padding
filtered = filtered_padded[padlen:-padlen]

# -----------------------------
# 5. Plot results
# -----------------------------
plt.figure(figsize=(5, 6))

# (a) Full signal comparison
#plt.subplot(2, 1, 1)
plt.plot(sig, label='Raw signal', alpha=0.6)
plt.plot(filtered, label='Filtered (zero-phase FIR)', linewidth=2)
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (a.u.)')
plt.legend()
plt.grid(True)

# (b) Zoom near end to check for edge artifacts
plt.subplot(2, 1, 2)
zoom_start = int(0.8 * len(t))
plt.plot(t[zoom_start:]*1000, sig[zoom_start:], label='Raw', alpha=0.6)
plt.plot(t[zoom_start:]*1000, filtered[zoom_start:], label='Filtered', linewidth=2)
plt.title('Zoom near end (check for artifacts)')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (a.u.)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()


