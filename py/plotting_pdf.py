import os
os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')

import numpy as np
import matplotlib
matplotlib.use('Agg', force=True)

from pynwb import NWBHDF5IO

def _get_sweep_from_sweep_table(nwbfile, sweep_number):
    st = getattr(nwbfile, 'sweep_table', None)
    if st is None:
        return None
    candidates = [sweep_number, sweep_number-1, sweep_number+1]
    for s in candidates:
        try:
            sw = st.get_sweep(int(s))
            return int(s), sw
        except Exception:
            pass
    return None

def _get_series_by_acquisition_key(nwbfile, key):
    acq = getattr(nwbfile, 'acquisition', None)
    if acq is None:
        return None
    if key in acq:
        return key, acq[key]
    return None

def _get_series_from_acquisition_patterns(nwbfile, sweep_number):
    acq = getattr(nwbfile, 'acquisition', None)
    if acq is None:
        return None
    keys = list(acq.keys())
    patterns = [
        f'data_{int(sweep_number):05d}_',   # matches your data_00090_AD0 style
        f'sweep_{sweep_number}',
        f'Sweep_{sweep_number}',
        f'sweep{sweep_number}',
        f'Sweep{sweep_number}',
        f'{sweep_number}'
    ]
    for k in keys:
        lk = k.lower()
        for p in patterns:
            if p.lower() in lk:
                return k, acq[k]
    return None

def _timeseries_to_xy(ts):
    data = np.asarray(ts.data[:], dtype=float)

    if ts.timestamps is not None:
        t = np.asarray(ts.timestamps[:], dtype=float)
    else:
        st = float(getattr(ts, 'starting_time', 0.0))
        rate = float(getattr(ts, 'rate', np.nan))
        if not np.isfinite(rate) or rate <= 0:
            raise ValueError('TimeSeries has no timestamps and invalid rate.')
        t = st + np.arange(data.shape[0], dtype=float) / rate

    units = getattr(ts, 'unit', '')
    return t, data, units

def _maybe_convert_units(y, units, y_unit_target):
    # Minimal helper: convert Volts->mV if requested
    if y_unit_target is None:
        return y, units
    u = (units or '').strip().lower()
    tgt = y_unit_target.strip().lower()
    if u in ['v', 'volt', 'volts'] and tgt in ['mv', 'millivolt', 'millivolts']:
        return y * 1e3, 'mV'
    return y, units

def _apply_tlim(t, y, tlim):
    if tlim is None:
        return t, y
    t0, t1 = float(tlim[0]), float(tlim[1])
    m = (t >= t0) & (t <= t1)
    return t[m], y[m]

def _void_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')

def _add_scalebar(ax,
                  dx=None, dy=None,
                  labelx=None, labely=None,
                  loc='lower right',
                  inset=0.08,
                  thickness=1.2,
                  text_offset=0.015):

    if dx is None and dy is None:
        return

    # Auto-label if labels not provided
    if labelx is None and dx is not None:
        labelx = f"{float(dx):g} s"
    if labely is None and dy is not None:
        # use axis label if present (e.g. 'mV'), otherwise omit units
        ylab = (ax.get_ylabel() or "").strip()
        unit = ylab if ylab else ""
        labely = f"{float(dy):g} {unit}".strip()

    # Axes limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xr = xlim[1] - xlim[0]
    yr = ylim[1] - ylim[0]

    # Anchor point
    if loc == 'lower left':
        x0 = xlim[0] + inset * xr
        y0 = ylim[0] + inset * yr
    else:  # lower right
        x0 = xlim[1] - inset * xr
        y0 = ylim[0] + inset * yr

    # Horizontal bar (time)
    if dx is not None:
        dx = float(dx)
        if loc == 'lower right':
            x1 = x0 - dx
            ax.plot([x1, x0], [y0, y0], color='black', lw=thickness, solid_capstyle='butt')
            tx = (x0 + x1) / 2
        else:
            x1 = x0 + dx
            ax.plot([x0, x1], [y0, y0], color='black', lw=thickness, solid_capstyle='butt')
            tx = (x0 + x1) / 2

        if labelx is not None:
            ax.text(tx, y0 - text_offset * yr, labelx,
                    ha='center', va='top', fontsize=8)

    # Vertical bar (voltage)
    if dy is not None:
        dy = float(dy)
        y1 = y0 + dy
        ax.plot([x0, x0], [y0, y1], color='black', lw=thickness, solid_capstyle='butt', clip_on=False)

        if labely is not None:
            ax.text(x0 + text_offset * xr, (y0 + y1) / 2,
                    labely, ha='left', va='center', rotation=90, fontsize=8)

    
def plot_nwb_sweeps(
    nwb_path,
    sweeps,                      # list: sweep numbers OR acquisition keys like 'data_00090_AD0'
    out_dir='figs/traces',
    fmt='pdf',
    dpi=300,
    width=7.0, height=2.2,
    lw=0.8,
    tlim=None,                   # applies after concatenation
    title=None,
    fonttype=42,

    # Time handling
    re_zero=True,                # force each sweep to start at 0 before concatenation
    gap_s=0.0,                   # optional gap between sweeps in concatenated x

    # TTL marker
    ttl_time=None,               # seconds after start of stim sweep
    ttl_on_sweep_index=1,        # 0-based index into sweeps list (your stim sweep typically 1 if baseline+stim)

    # Publication style
    pub_void_axes=False,
    scalebar_x=None,             # e.g. 1.0 (seconds)
    scalebar_y=None,             # e.g. 10.0 (mV) if converted
    scalebar_label_x=None,       # e.g. '1 s'
    scalebar_label_y=None,       # e.g. '10 mV'
    scalebar_loc='lower right',

    # Units convenience
    y_unit_target=None           # e.g. 'mV' to convert from V
):
    import matplotlib.pyplot as plt

    os.makedirs(out_dir, exist_ok=True)
    fmt = fmt.lower().strip()
    if fmt not in ('pdf', 'png'):
        raise ValueError("fmt must be 'pdf' or 'png'")

    plt.rcParams['pdf.fonttype'] = fonttype
    plt.rcParams['ps.fonttype']  = fonttype
    plt.rcParams['axes.linewidth'] = 0.6

    # normalize sweeps list
    if not isinstance(sweeps, (list, tuple)) or len(sweeps) < 1:
        raise ValueError('sweeps must be a non-empty list/tuple')

    with NWBHDF5IO(nwb_path, 'r', load_namespaces=True) as io:
        nwbfile = io.read()

        pieces_t = []
        pieces_y = []
        y_units = None

        # We track where each sweep begins in concatenated time
        sweep_start_offsets = []

        current_offset = 0.0

        for i, sw in enumerate(sweeps):
            resp_ts = None
            sweep_label = None

            # Case 1: acquisition key string (e.g., 'data_00090_AD0')
            if isinstance(sw, str):
                got = _get_series_by_acquisition_key(nwbfile, sw)
                if got is None:
                    raise ValueError(f"Acquisition key not found: {sw}")
                sweep_label, resp_ts = got

            # Case 2: numeric sweep number
            else:
                sw_num = int(sw)
                sw_info = _get_sweep_from_sweep_table(nwbfile, sw_num)
                if sw_info is not None:
                    used, sw_obj = sw_info
                    sweep_label = f'sweep_table:{used}'
                    if isinstance(sw_obj, dict):
                        resp_list = sw_obj.get('response', None)
                        if resp_list is not None and len(resp_list) > 0:
                            resp_ts = resp_list[0]

                if resp_ts is None:
                    acq_info = _get_series_from_acquisition_patterns(nwbfile, sw_num)
                    if acq_info is None:
                        raise ValueError(f'Could not locate sweep {sw_num} via sweep_table or acquisition patterns.')
                    sweep_label, resp_ts = acq_info

            t, y, units = _timeseries_to_xy(resp_ts)
            y, units2 = _maybe_convert_units(y, units, y_unit_target)
            if y_units is None:
                y_units = units2

            # re-zero each sweep
            if re_zero:
                t = t - t[0]

            # record offset for TTL positioning
            sweep_start_offsets.append(current_offset)

            # shift into concatenated axis
            t_shift = t + current_offset

            pieces_t.append(t_shift)
            pieces_y.append(y)

            # advance offset
            current_offset = float(t_shift[-1]) + float(gap_s)

        T = np.concatenate(pieces_t)
        Y = np.concatenate(pieces_y)

        # Apply global tlim after concat
        if tlim is not None:
            T, Y = _apply_tlim(T, Y, tlim)

        fig = plt.figure(figsize=(width, height), constrained_layout=True)
        ax = fig.add_subplot(1,1,1)
        ax.plot(T, Y, linewidth=lw)

        if ttl_time is not None:
            idx = int(ttl_on_sweep_index)
            if idx < 0 or idx >= len(sweep_start_offsets):
                raise ValueError('ttl_on_sweep_index out of range for sweeps list')
        
            x_ttl = sweep_start_offsets[idx] + float(ttl_time)
        
            ymin, ymax = ax.get_ylim()
            yr = (ymax - ymin)
            y_marker = ymin + 0.02 * yr          # 2% above bottom
            stem_h   = 0.05 * yr                 # 5% tall stem
        
            # stem first
            ax.plot([x_ttl, x_ttl],
                    [y_marker, y_marker + stem_h],
                    color='red', lw=0.8, alpha=0.8, zorder=9, clip_on=False)
        
            # triangle on top
            ax.plot([x_ttl], [y_marker],
                    marker='v',
                    color='red',
                    markersize=7,
                    markeredgewidth=0,
                    zorder=10,
                    clip_on=False)



        # labels / title
        if not pub_void_axes:
            ax.set_xlabel('Time (s)')
            ax.set_ylabel(y_units if y_units else 'Signal')
        if title is None:
            base = os.path.basename(nwb_path)
            ax.set_title(f'{base} | sweeps {sweeps}', fontsize=9)
        else:
            ax.set_title(title, fontsize=9)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # publication style
        if pub_void_axes:
            _void_axes(ax)
            _add_scalebar(
                ax,
                dx=scalebar_x,
                dy=scalebar_y,
                labelx=scalebar_label_x,
                labely=scalebar_label_y,
                loc=scalebar_loc,
                thickness=1.0
            )

        # Save
        stem = os.path.splitext(os.path.basename(nwb_path))[0]
        tag = '__'.join([str(s) for s in sweeps])
        out_path = os.path.join(out_dir, f'{stem}__{tag}.{fmt}')

        if fmt == 'pdf':
            fig.savefig(out_path, format='pdf')
        else:
            fig.savefig(out_path, format='png', dpi=dpi)

        plt.close(fig)
        return out_path
