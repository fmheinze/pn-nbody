#!/usr/bin/env python3
"""
trajectory_anim3.py

Animate N-body trajectories from a whitespace-delimited data file.

Input format:
- First column: time t
- Remaining columns: coordinates for each body.
  2D:  t x0 y0 x1 y1 ...
  3D:  t x0 y0 z0 x1 y1 z1 ...

This script:
- Auto-detects 2D vs 3D primarily from the header (presence of z0/z1/...).
  If no header is present and the column count is ambiguous (divisible by 2 and 3),
  it defaults to 2D unless overridden via --dim.
- Skips frames automatically to keep playback time reasonable for long runs.
- Ensures each body's marker color matches its trajectory line color.
- Enforces equal axis ratios (square in 2D, cube in 3D).

Usage:
  python trajectory_anim3.py output_pos.dat
  python trajectory_anim3.py 2d.dat

Optional:
  --dim 2|3            Force dimensionality (useful if no header).
  --duration SECONDS   Target playback length (default 10).
  --max-frames N       Hard cap on animation frames (default 2000).
  --interval MS        Delay between frames in ms (default 30).
"""

import argparse
import math
from typing import Tuple, List, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def _load_data(path: str) -> Tuple[np.ndarray, Optional[List[str]]]:
    """
    Load numeric data and (optionally) header names.
    Returns (data, names) where names is None if no header.
    """
    # Try to read header line; if it starts with a letter, treat as header
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline()

    has_header = False
    names = None
    if first.strip() and any(ch.isalpha() for ch in first.strip().split()[0]):
        # e.g. "t x0 y0 z0 ..."
        has_header = True
        names = first.strip().split()

    if has_header:
        data = np.genfromtxt(path, comments="#", skip_header=1)
    else:
        data = np.genfromtxt(path, comments="#")

    if data.ndim == 1:
        # Single row -> make it 2D for consistency
        data = data[None, :]

    if data.shape[1] < 3:
        raise ValueError(f"Expected at least 3 columns (t and coords). Got {data.shape[1]}.")

    return data, names


def _infer_dim(ncols_total: int, names: Optional[List[str]], forced_dim: Optional[int]) -> int:
    """
    Infer dimensionality (2 or 3).
    - If forced_dim is provided, return it.
    - Else if names include any 'z' coordinate, return 3.
    - Else infer from coordinate column count; default to 2 if ambiguous.
    """
    if forced_dim in (2, 3):
        return forced_dim

    if names:
        # Any column name like z0, z1, ...
        for nm in names:
            if nm.startswith("z"):
                return 3
        return 2

    ncoords = ncols_total - 1
    div2 = (ncoords % 2 == 0)
    div3 = (ncoords % 3 == 0)
    if div3 and not div2:
        return 3
    if div2 and not div3:
        return 2
    # Ambiguous or neither -> default to 2 (user can override)
    return 2


def _compute_frame_indices(t: np.ndarray, interval_ms: int, target_duration_s: float, max_frames: int) -> np.ndarray:
    """
    Choose a subset of frame indices to keep the animation watchable.
    We aim for ~target_duration_s of playback at the chosen interval, but
    we also cap the number of frames by max_frames.
    """
    n = len(t)
    if n <= 1:
        return np.arange(n, dtype=int)

    fps = 1000.0 / max(interval_ms, 1)
    target_frames = int(max(2, round(target_duration_s * fps)))
    keep = min(n, max_frames, target_frames)

    if keep >= n:
        return np.arange(n, dtype=int)

    # Sample uniformly in index space (monotonic time is assumed)
    idx = np.linspace(0, n - 1, keep)
    idx = np.unique(np.round(idx).astype(int))
    return idx


def _reshape_positions(data: np.ndarray, dim: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract time and positions.
    Returns t (N,), pos (N, nbodies, dim)
    """
    t = data[:, 0]
    coords = data[:, 1:]
    ncoords = coords.shape[1]
    if ncoords % dim != 0:
        raise ValueError(
            f"Coordinate columns ({ncoords}) not divisible by dim={dim}. "
            f"Try forcing --dim 2 or --dim 3."
        )
    nb = ncoords // dim
    pos = coords.reshape((-1, nb, dim))
    return t, pos


def _set_equal_axes_2d(ax, pos: np.ndarray, pad_frac: float = 0.05):
    """
    Square plot with equal scaling in x/y, spanning all trajectories.
    """
    x = pos[:, :, 0].ravel()
    y = pos[:, :, 1].ravel()
    xmin, xmax = np.nanmin(x), np.nanmax(x)
    ymin, ymax = np.nanmin(y), np.nanmax(y)

    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    span = max(xmax - xmin, ymax - ymin)
    if span == 0:
        span = 1.0
    span *= (1.0 + pad_frac)

    half = 0.5 * span
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal", adjustable="box")


def _set_equal_axes_3d(ax, pos: np.ndarray, pad_frac: float = 0.05):
    """
    Cube plot with equal scaling in x/y/z, spanning all trajectories.
    Works across matplotlib versions (uses set_box_aspect when available).
    """
    x = pos[:, :, 0].ravel()
    y = pos[:, :, 1].ravel()
    z = pos[:, :, 2].ravel()
    xmin, xmax = np.nanmin(x), np.nanmax(x)
    ymin, ymax = np.nanmin(y), np.nanmax(y)
    zmin, zmax = np.nanmin(z), np.nanmax(z)

    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    cz = 0.5 * (zmin + zmax)
    span = max(xmax - xmin, ymax - ymin, zmax - zmin)
    if span == 0:
        span = 1.0
    span *= (1.0 + pad_frac)
    half = 0.5 * span

    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_zlim(cz - half, cz + half)

    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect((1, 1, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", help="Path to trajectory data file (whitespace-delimited).")
    ap.add_argument("--dim", type=int, choices=[2, 3], default=None, help="Force dimensionality (2 or 3).")
    ap.add_argument("--duration", type=float, default=10.0, help="Target playback length in seconds.")
    ap.add_argument("--max-frames", type=int, default=2000, help="Maximum number of frames to animate.")
    ap.add_argument("--interval", type=int, default=30, help="Delay between frames in ms.")
    ap.add_argument("--trail", type=int, default=0,
                    help="If >0, show only the last TRAIL points of each trajectory (0 = full).")
    args = ap.parse_args()

    data, names = _load_data(args.file)
    dim = _infer_dim(data.shape[1], names, args.dim)
    t, pos = _reshape_positions(data, dim)

    frame_idx = _compute_frame_indices(t, args.interval, args.duration, args.max_frames)
    t_anim = t[frame_idx]
    pos_anim = pos[frame_idx, :, :]

    nb = pos_anim.shape[1]

    # Figure / axes
    if dim == 3:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        _set_equal_axes_3d(ax, pos_anim)
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        _set_equal_axes_2d(ax, pos_anim)
        ax.set_facecolor("whitesmoke")
        ax.grid(ls="--", alpha=0.5)

    # Artists: one line + one point per body, sharing color
    lines = []
    points = []

    for i in range(nb):
        if dim == 3:
            (ln,) = ax.plot([], [], [], lw=1)
            color = ln.get_color()
            (pt,) = ax.plot([], [], [], marker="o", linestyle="None", markersize=5, color=color)
        else:
            (ln,) = ax.plot([], [], lw=1)
            color = ln.get_color()
            (pt,) = ax.plot([], [], marker="o", linestyle="None", markersize=5, color=color)

        lines.append(ln)
        points.append(pt)

    title = ax.set_title("")

    def init():
        for ln, pt in zip(lines, points):
            if dim == 3:
                ln.set_data([], [])
                ln.set_3d_properties([])
                pt.set_data([], [])
                pt.set_3d_properties([])
            else:
                ln.set_data([], [])
                pt.set_data([], [])
        title.set_text("")
        return [*lines, *points, title]

    def update(frame: int):
        p = pos_anim[frame]
        cur_t = t_anim[frame]

        for i in range(nb):
            if args.trail and args.trail > 1:
                start = max(0, frame - args.trail + 1)
                trail = pos_anim[start:frame + 1, i, :]
            else:
                trail = pos_anim[:frame + 1, i, :]

            if dim == 3:
                lines[i].set_data(trail[:, 0], trail[:, 1])
                lines[i].set_3d_properties(trail[:, 2])

                points[i].set_data([p[i, 0]], [p[i, 1]])
                points[i].set_3d_properties([p[i, 2]])
            else:
                lines[i].set_data(trail[:, 0], trail[:, 1])
                points[i].set_data([p[i, 0]], [p[i, 1]])

        title.set_text(f"t = {cur_t:g}")
        return [*lines, *points, title]

    ani = FuncAnimation(fig, update, frames=len(t_anim), init_func=init,
                        interval=args.interval, blit=False)

    plt.show()


if __name__ == "__main__":
    main()
