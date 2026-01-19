#!/usr/bin/env python3
"""
Plot particle trajectories from two position-time data files
on top of each other in a single 3D figure, with matching colors
for the same particles across files.

Usage:
    python plot_two_trajectories.py file1.dat file2.dat
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# ---------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------
def load_data(filename):
    """
    Load position-time data.
    Expected format:
        t  x0  y0  z0  x1  y1  z1  x2  y2  z2 ...
    """
    with open(filename, "r") as f:
        header = f.readline().strip().split()

    data = np.loadtxt(filename, skiprows=1)
    if data.ndim == 1:
        data = data[None, :]

    n_steps, n_cols = data.shape
    coord_cols = n_cols - 1
    if coord_cols % 3 != 0:
        raise ValueError(
            f"Expected 3 coordinates per particle, got {coord_cols} in {filename}."
        )

    n_particles = coord_cols // 3
    t = data[:, 0]
    coords = data[:, 1:].reshape(n_steps, n_particles, 3)
    return t, coords, n_particles


# ---------------------------------------------------------------
# 3D equal aspect helper
# ---------------------------------------------------------------
def set_equal_aspect_3d(ax, X, Y, Z):
    """Set equal aspect ratio for a 3D plot."""
    x_min, x_max = np.min(X), np.max(X)
    y_min, y_max = np.min(Y), np.max(Y)
    z_min, z_max = np.min(Z), np.max(Z)

    pad = 0.05
    x_pad = pad * (x_max - x_min or 1.0)
    y_pad = pad * (y_max - y_min or 1.0)
    z_pad = pad * (z_max - z_min or 1.0)
    x_min, x_max = x_min - x_pad, x_max + x_pad
    y_min, y_max = y_min - y_pad, y_max + y_pad
    z_min, z_max = z_min - z_pad, z_max + z_pad

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    max_range = max(x_range, y_range, z_range)
    x_mid = 0.5 * (x_max + x_min)
    y_mid = 0.5 * (y_max + y_min)
    z_mid = 0.5 * (z_max + z_min)

    ax.set_xlim(x_mid - max_range / 2, x_mid + max_range / 2)
    ax.set_ylim(y_mid - max_range / 2, y_mid + max_range / 2)
    ax.set_zlim(z_mid - max_range / 2, z_mid + max_range / 2)

    try:
        ax.set_box_aspect([1, 1, 1])  # Matplotlib >= 3.3
    except AttributeError:
        pass


# ---------------------------------------------------------------
# Plot function
# ---------------------------------------------------------------
def plot_two_runs(t1, coords1, t2, coords2, n_particles, label1, label2):
    """Plot two runs' trajectories in a single 3D plot with matching colors."""
    # Combine coordinates to find overall bounds
    x_all = np.concatenate([coords1[:, :, 0], coords2[:, :, 0]], axis=None)
    y_all = np.concatenate([coords1[:, :, 1], coords2[:, :, 1]], axis=None)
    z_all = np.concatenate([coords1[:, :, 2], coords2[:, :, 2]], axis=None)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    set_equal_aspect_3d(ax, x_all, y_all, z_all)

    # Color map — consistent colors for same particle
    cmap = plt.get_cmap("tab10" if n_particles <= 10 else "hsv")
    colors = [cmap(i / n_particles) for i in range(n_particles)]

    for i in range(n_particles):
        color = colors[i]
        # Run 1: solid
        ax.plot(
            coords1[:, i, 0],
            coords1[:, i, 1],
            coords1[:, i, 2],
            lw=1.5,
            color=color,
            label=f"{label1} p{i}" if i == 0 else None,
        )
        # Run 2: dashed
        ax.plot(
            coords2[:, i, 0],
            coords2[:, i, 1],
            coords2[:, i, 2],
            lw=1.5,
            ls="--",
            color=color,
            label=f"{label2} p{i}" if i == 0 else None,
        )

    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Plot trajectories from two files with matching colors per particle."
    )
    parser.add_argument(
        "file1",
        nargs="?",
        default="output_pos_1.dat",
        help="First input data file (default: output_pos_1.dat)",
    )
    parser.add_argument(
        "file2",
        nargs="?",
        default="output_pos_2.dat",
        help="Second input data file (default: output_pos_2.dat)",
    )

    args = parser.parse_args()

    t1, coords1, n1 = load_data(args.file1)
    t2, coords2, n2 = load_data(args.file2)

    if n1 != n2:
        raise ValueError(
            f"Files have different numbers of particles: {n1} vs {n2} "
            f"({args.file1} vs {args.file2})"
        )

    print(
        f"Loaded {len(t1)} steps, {n1} particles from {args.file1!r};\n"
        f"       {len(t2)} steps, {n2} particles from {args.file2!r}"
    )

    plot_two_runs(t1, coords1, t2, coords2, n1, label1=args.file1, label2=args.file2)


if __name__ == "__main__":
    main()
