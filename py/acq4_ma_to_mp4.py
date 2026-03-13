#!/usr/bin/env python3
"""
acq4_ma_to_mp4.py

Convert an ACQ4 MetaArray (.ma) movie to an MP4 suitable for PowerPoint.

Typical usage:
  python acq4_ma_to_mp4.py input.ma output.mp4 --fps 5

If your movie is time + z-stack:
  python acq4_ma_to_mp4.py input.ma output.mp4 --fps 5 --z 0
  python acq4_ma_to_mp4.py input.ma output.mp4 --fps 5 --z max
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def _import_metaarray():
    """
    ACQ4 commonly uses MetaArray via pyqtgraph.metaarray.MetaArray.
    Some installations have MetaArray as a top-level module.
    """
    try:
        from pyqtgraph.metaarray import MetaArray  # type: ignore
        return MetaArray
    except Exception:
        pass

    try:
        from MetaArray import MetaArray  # type: ignore
        return MetaArray
    except Exception as e:
        raise ImportError(
            "Could not import MetaArray.\n"
            "Try:\n"
            "  pip install pyqtgraph\n"
            "or ensure you're using the same Python environment as ACQ4.\n"
            f"Original error: {e}"
        )


def _load_ma(path, MetaArray):
    ma = MetaArray(file=path)
    arr = np.array(ma)  # strips metadata, keeps numeric array
    return arr, ma


def _normalize_to_uint8(frame, p_lo=1.0, p_hi=99.0):
    """
    Percentile-based contrast stretch to 8-bit for video encoding.
    """
    f = frame.astype(np.float32, copy=False)
    if not np.isfinite(f).any():
        return np.zeros_like(frame, dtype=np.uint8)

    lo = np.nanpercentile(f, p_lo)
    hi = np.nanpercentile(f, p_hi)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        # Fallback: min/max
        lo = np.nanmin(f)
        hi = np.nanmax(f)
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            return np.zeros_like(frame, dtype=np.uint8)

    f = (f - lo) / (hi - lo)
    f = np.clip(f, 0.0, 1.0)
    return (f * 255.0 + 0.5).astype(np.uint8)


def _select_frames(arr, z_mode=None):
    """
    Return frames shaped (t, y, x).
    Supported input shapes:
      (t, y, x)
      (t, z, y, x)
      (y, x, t) or (x, y, t) are not handled automatically—ACQ4 usually stores t first.
    """
    if arr.ndim == 3:
        # assume (t, y, x)
        return arr
    if arr.ndim == 4:
        # assume (t, z, y, x)
        if z_mode is None:
            raise ValueError(
                "Input appears to be (t, z, y, x). Provide --z <index> or --z max."
            )
        if z_mode == "max":
            return np.nanmax(arr, axis=1)
        else:
            z_idx = int(z_mode)
            return arr[:, z_idx, :, :]
    raise ValueError(f"Unsupported array shape {arr.shape}. Expected 3D or 4D.")


def _write_png_sequence(frames_tyxt, out_dir, p_lo, p_hi):
    """
    Write frames as PNG files: frame_000001.png, ...
    Uses imageio if available; falls back to PIL.
    """
    try:
        import imageio.v2 as imageio  # type: ignore
        writer = "imageio"
    except Exception:
        imageio = None
        writer = "pil"

    if writer == "pil":
        try:
            from PIL import Image  # type: ignore
        except Exception as e:
            raise ImportError(
                "Need either imageio or pillow to write PNG frames.\n"
                "Try:\n"
                "  pip install imageio pillow\n"
                f"Original error: {e}"
            )

    t = frames_tyxt.shape[0]
    for i in range(t):
        f8 = _normalize_to_uint8(frames_tyxt[i], p_lo=p_lo, p_hi=p_hi)
        fn = os.path.join(out_dir, f"frame_{i:06d}.png")
        if writer == "imageio":
            imageio.imwrite(fn, f8)
        else:
            Image.fromarray(f8, mode="L").save(fn)


def _run_ffmpeg(png_dir, fps, out_mp4, scale=None):
    """
    Build an MP4 using H.264 and yuv420p (PowerPoint-friendly).
    """
    if shutil.which("ffmpeg") is None:
        raise RuntimeError(
            "ffmpeg not found on PATH.\n"
            "Install it (Ubuntu):\n"
            "  sudo apt-get update && sudo apt-get install -y ffmpeg"
        )

    inp = os.path.join(png_dir, "frame_%06d.png")
    cmd = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(fps),
        "-i",
        inp,
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-movflags",
        "+faststart",
    ]

    # optional scaling (good for huge frames)
    if scale is not None:
        # scale like "1920:-2" or "1280:-2"
        cmd.extend(["-vf", f"scale={scale}"])

    cmd.append(out_mp4)

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "ffmpeg failed.\n\nSTDERR (last 60 lines):\n"
            + "\n".join(proc.stderr.splitlines()[-60:])
        )


def main():
    ap = argparse.ArgumentParser(description="Convert ACQ4 .ma (MetaArray) to MP4.")
    ap.add_argument("input_ma", help="Path to ACQ4 .ma file")
    ap.add_argument("output_mp4", help="Output MP4 path")
    ap.add_argument("--fps", type=float, default=5.0, help="Frames per second (default: 5)")
    ap.add_argument(
        "--z",
        default=None,
        help="For (t,z,y,x) data: z index (e.g. 0) or 'max' for max-projection.",
    )
    ap.add_argument(
        "--p_lo",
        type=float,
        default=1.0,
        help="Low percentile for contrast stretch (default: 1)",
    )
    ap.add_argument(
        "--p_hi",
        type=float,
        default=99.0,
        help="High percentile for contrast stretch (default: 99)",
    )
    ap.add_argument(
        "--scale",
        default=None,
        help="Optional ffmpeg scale, e.g. '1920:-2' to downscale large frames.",
    )
    args = ap.parse_args()

    MetaArray = _import_metaarray()
    arr, _ma = _load_ma(args.input_ma, MetaArray)

    frames = _select_frames(arr, z_mode=args.z)

    # Ensure (t, y, x) ordering; ACQ4 typically does this already.
    if frames.ndim != 3:
        raise ValueError(f"Internal error: frames ndim {frames.ndim}")

    with tempfile.TemporaryDirectory(prefix="acq4_ma_frames_") as td:
        _write_png_sequence(frames, td, p_lo=args.p_lo, p_hi=args.p_hi)
        _run_ffmpeg(td, fps=args.fps, out_mp4=args.output_mp4, scale=args.scale)

    print(f"Done: {args.output_mp4}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
