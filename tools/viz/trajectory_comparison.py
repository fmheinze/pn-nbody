#!/usr/bin/env python3
"""trajectory_plot_compare.py

Plot N-body trajectories from TWO whitespace-delimited data files in the SAME plot
for visual comparison (no animation).

Input format:
- First column: time t
- Remaining columns: coordinates for each body.
  2D:  t x0 y0 x1 y1 ...
  3D:  t x0 y0 z0 x1 y1 z1 ...

Behavior:
- Auto-detects 2D vs 3D primarily from the header (presence of z0/z1/...).
  If no header is present and the column count is ambiguous (divisible by 2 and 3),
  it defaults to 2D unless overridden via --dim.
- Plots trajectories from file1 with SOLID lines.
- Plots trajectories from file2 with DASHED lines, using the SAME colors per body
  as the corresponding body in file1.
- Marks final positions (optional): by default, file1 uses circle markers ("o"),
  file2 uses "x" markers, both in the matching body color.
- Enforces equal axis ratios (square in 2D, cube in 3D).
- Optionally decimates very long trajectories for faster plotting.

Usage:
  python trajectory_plot_compare.py runA.dat runB.dat

Optional:
  --dim 2|3            Force dimensionality (useful if no header).
  --max-points N       If >0, uniformly downsample to at most N time points per file for plotting.
                       (Markers still reflect the true final sample.)
  --pad-frac F         Axis padding fraction (default 0.05).
  --no-markers         Do not place final-position markers.
"""

import argparse
from typing import Tuple, List, Optional

import numpy as np
import matplotlib.pyplot as plt


def _load_data(path: str) -> Tuple[np.ndarray, Optional[List[str]]]:
    """Load numeric data and (optionally) header names."""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        first = f.readline()

    names: Optional[List[str]] = None
    has_header = False
    if first.strip():
        # Header if the first token starts with a letter (e.g. "t") OR if any token contains letters.
        toks = first.strip().split()
        if any(any(ch.isalpha() for ch in tok) for tok in toks):
            has_header = True
            names = toks

    if has_header:
        data = np.genfromtxt(path, comments="#", skip_header=1)
    else:
        data = np.genfromtxt(path, comments="#")

    if data.ndim == 1:
        data = data[None, :]

    if data.shape[1] < 3:
        raise ValueError(f"{path}: expected at least 3 columns (t and coords). Got {data.shape[1]}.")

    return data, names


def _infer_dim(ncols_total: int, names: Optional[List[str]], forced_dim: Optional[int]) -> int:
    """Infer dimensionality (2 or 3)."""
    if forced_dim in (2, 3):
        return forced_dim

    if names:
        # If any column name starts with 'z', treat as 3D; else 2D.
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

    # Ambiguous or neither -> default to 2D
    return 2


def _reshape_positions(data: np.ndarray, dim: int, path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Extract time and positions: t (N,), pos (N, nbodies, dim)."""
    t = data[:, 0]
    coords = data[:, 1:]
    ncoords = coords.shape[1]

    if ncoords % dim != 0:
        raise ValueError(
            f"{path}: coordinate columns ({ncoords}) not divisible by dim={dim}. "
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


def _set_equal_axes_2d(ax, pos_all: np.ndarray, pad_frac: float = 0.05):
    x = pos_all[:, :, 0].ravel()
    y = pos_all[:, :, 1].ravel()
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


def _set_equal_axes_3d(ax, pos_all: np.ndarray, pad_frac: float = 0.05):
    x = pos_all[:, :, 0].ravel()
    y = pos_all[:, :, 1].ravel()
    z = pos_all[:, :, 2].ravel()
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
    ap.add_argument("file1", help="First trajectory data file (solid lines).")
    ap.add_argument("file2", help="Second trajectory data file (dashed lines, same colors).")
    ap.add_argument("--dim", type=int, choices=[2, 3], default=None, help="Force dimensionality (2 or 3).")
    ap.add_argument("--max-points", type=int, default=0,
                    help="If >0, uniformly downsample to at most this many time points per file for plotting.")
    ap.add_argument("--pad-frac", type=float, default=0.05, help="Axis padding fraction.")
    ap.add_argument("--no-markers", action="store_true", help="Do not draw final-position markers.")
    args = ap.parse_args()

    data1, names1 = _load_data(args.file1)
    data2, names2 = _load_data(args.file2)

    dim1 = _infer_dim(data1.shape[1], names1, args.dim)
    dim2 = _infer_dim(data2.shape[1], names2, args.dim)
    if dim1 != dim2:
        raise ValueError(
            f"Dimensionality mismatch: {args.file1} inferred dim={dim1}, "
            f"{args.file2} inferred dim={dim2}. Use --dim to force."
        )
    dim = dim1

    t1, pos1 = _reshape_positions(data1, dim, args.file1)
    t2, pos2 = _reshape_positions(data2, dim, args.file2)

    nb1, nb2 = pos1.shape[1], pos2.shape[1]
    if nb1 != nb2:
        raise ValueError(
            f"Number of bodies mismatch: {args.file1} has {nb1}, {args.file2} has {nb2}. "
            "To compare, both files should contain the same bodies in the same order."
        )
    nb = nb1

    final1 = pos1[-1, :, :]
    final2 = pos2[-1, :, :]

    idx1 = _decimate_indices(len(t1), args.max_points if args.max_points > 0 else None)
    idx2 = _decimate_indices(len(t2), args.max_points if args.max_points > 0 else None)
    pos1p = pos1[idx1, :, :]
    pos2p = pos2[idx2, :, :]

    if dim == 3:
        fig = plt.figure(figsize=(6, 6), constrained_layout=True)
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
    else:
        fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_facecolor("whitesmoke")
        ax.grid(ls="--", alpha=0.5)

    # Plot file1 (solid) and capture colors per body
    colors: List[str] = []
    for i in range(nb):
        if dim == 3:
            (ln1,) = ax.plot(pos1p[:, i, 0], pos1p[:, i, 1], pos1p[:, i, 2], lw=1, linestyle='-')
        else:
            (ln1,) = ax.plot(pos1p[:, i, 0], pos1p[:, i, 1], lw=1, linestyle='-')
        c = ln1.get_color()
        colors.append(c)

        if not args.no_markers:
            if dim == 3:
                ax.plot([final1[i, 0]], [final1[i, 1]], [final1[i, 2]],
                        marker="o", linestyle="None", markersize=5, color=c)
            else:
                ax.plot([final1[i, 0]], [final1[i, 1]],
                        marker="o", linestyle="None", markersize=5, color=c)

    # Plot file2 (dashed) using same colors
    for i in range(nb):
        c = colors[i]
        if dim == 3:
            ax.plot(pos2p[:, i, 0], pos2p[:, i, 1], pos2p[:, i, 2],
                    lw=1, linestyle="--", color=c)
            if not args.no_markers:
                ax.plot([final2[i, 0]], [final2[i, 1]], [final2[i, 2]],
                        marker="x", linestyle="None", markersize=6, color=c)
        else:
            ax.plot(pos2p[:, i, 0], pos2p[:, i, 1],
                    lw=1, linestyle="--", color=c)
            if not args.no_markers:
                ax.plot([final2[i, 0]], [final2[i, 1]],
                        marker="x", linestyle="None", markersize=6, color=c)

    # Equal aspect / cube based on combined extents
    # Combine a decimated set for bounds to keep it cheap
    pos_all = np.concatenate([pos1p, pos2p], axis=0)
    if dim == 3:
        _set_equal_axes_3d(ax, pos_all, pad_frac=args.pad_frac)
    else:
        _set_equal_axes_2d(ax, pos_all, pad_frac=args.pad_frac)

    plt.show()


if __name__ == "__main__":
    main()
