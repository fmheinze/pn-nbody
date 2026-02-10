#!/usr/bin/env python3
"""trajectory_plot_final.py

Plot N-body trajectories (no animation) from a whitespace-delimited data file.

Input format:
- First column: time t
- Remaining columns: coordinates for each body.
  2D:  t x0 y0 x1 y1 ...
  3D:  t x0 y0 z0 x1 y1 z1 ...

This script:
- Auto-detects 2D vs 3D primarily from the header (presence of z0/z1/...).
  If no header is present and the column count is ambiguous (divisible by 2 and 3),
  it defaults to 2D unless overridden via --dim.
- Draws full trajectories as lines.
- Places a marker at each body's FINAL location.
- Enforces equal axis ratios (square in 2D, cube in 3D).
- Optionally decimates very long trajectories for faster plotting.

Usage:
  python trajectory_plot_final.py output_pos.dat
  python trajectory_plot_final.py 2d.dat

Optional:
  --dim 2|3            Force dimensionality (useful if no header).
  --max-points N       If set, uniformly downsample to at most N time points for plotting.
                       (Markers still reflect the true final sample.)
  --pad-frac F         Axis padding fraction (default 0.05).
"""

import argparse
from typing import Tuple, List, Optional

import numpy as np
import matplotlib.pyplot as plt


def _load_data(path: str) -> Tuple[np.ndarray, Optional[List[str]]]:
    """Load numeric data and (optionally) header names."""
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
        data = data[None, :]

    if data.shape[1] < 3:
        raise ValueError(f"Expected at least 3 columns (t and coords). Got {data.shape[1]}.")

    return data, names


def _infer_dim(ncols_total: int, names: Optional[List[str]], forced_dim: Optional[int]) -> int:
    """Infer dimensionality (2 or 3)."""
    if forced_dim in (2, 3):
        return forced_dim

    if names:
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

    # Ambiguous or neither -> default to 2
    return 2


def _reshape_positions(data: np.ndarray, dim: int) -> Tuple[np.ndarray, np.ndarray]:
    """Extract time and positions: t (N,), pos (N, nbodies, dim)."""
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


def _decimate_indices(n: int, max_points: Optional[int]) -> np.ndarray:
    """Uniformly downsample indices to length <= max_points."""
    if not max_points or max_points <= 0 or n <= max_points:
        return np.arange(n, dtype=int)

    idx = np.linspace(0, n - 1, int(max_points))
    idx = np.unique(np.round(idx).astype(int))
    return idx


def _set_equal_axes_2d(ax, pos: np.ndarray, pad_frac: float = 0.05):
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
    ap.add_argument("--max-points", type=int, default=0,
                    help="If >0, uniformly downsample to at most this many time points for plotting.")
    ap.add_argument("--pad-frac", type=float, default=0.05, help="Axis padding fraction.")
    args = ap.parse_args()

    data, names = _load_data(args.file)
    dim = _infer_dim(data.shape[1], names, args.dim)
    t, pos = _reshape_positions(data, dim)

    # True final position (even if we decimate for plotting)
    final_pos = pos[-1, :, :]

    idx = _decimate_indices(len(t), args.max_points if args.max_points > 0 else None)
    pos_plot = pos[idx, :, :]

    nb = pos.shape[1]

    if dim == 3:
        fig = plt.figure(figsize=(6, 6), constrained_layout=True)
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_facecolor("whitesmoke")
        ax.grid(ls="--", alpha=0.5)

    # One line per body; marker uses same color as the line
    for i in range(nb):
        if dim == 3:
            (ln,) = ax.plot(pos_plot[:, i, 0], pos_plot[:, i, 1], pos_plot[:, i, 2], lw=1)
            c = ln.get_color()
            ax.plot([final_pos[i, 0]], [final_pos[i, 1]], [final_pos[i, 2]],
                    marker="o", linestyle="None", markersize=5, color=c)
        else:
            (ln,) = ax.plot(pos_plot[:, i, 0], pos_plot[:, i, 1], lw=1)
            c = ln.get_color()
            ax.plot([final_pos[i, 0]], [final_pos[i, 1]],
                    marker="o", linestyle="None", markersize=5, color=c)

    # Equal aspect / cube
    if dim == 3:
        _set_equal_axes_3d(ax, pos_plot, pad_frac=args.pad_frac)
    else:
        _set_equal_axes_2d(ax, pos_plot, pad_frac=args.pad_frac)

    plt.show()


if __name__ == "__main__":
    main()
