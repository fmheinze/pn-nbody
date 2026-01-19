#!/usr/bin/env python3
"""
Animate particle trajectories from position-time data.

Usage:
    python animate_particles.py output_pos.dat

If no filename is given, defaults to 'output_pos.dat'.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D projection)


def load_data(filename):
    """
    Load position-time data.
    Expected format:
        t  x0  y0  z0  x1  y1  z1  x2  y2  z2 ...
    """
    # Read header to know how many columns there are (mainly for sanity check)
    with open(filename, "r") as f:
        header = f.readline().strip().split()

    data = np.loadtxt(filename, skiprows=1)

    if data.ndim == 1:
        # Single row -> promote to 2D
        data = data[None, :]

    n_steps, n_cols = data.shape

    # First column is time, rest are coordinates
    coord_cols = n_cols - 1
    if coord_cols % 3 != 0:
        raise ValueError(
            f"Expected 3 coordinates per particle (x,y,z), "
            f"but got {coord_cols} coordinate columns."
        )

    n_particles = coord_cols // 3

    t = data[:, 0]
    coords = data[:, 1:].reshape(n_steps, n_particles, 3)  # (time, particle, xyz)

    return t, coords, n_particles


def animate_trajectories(t, coords, n_particles, save=None):
    """
    Create and show (or save) a 3D animation of particle trajectories.

    t:          (T,) array of times
    coords:     (T, N, 3) array of positions
    n_particles: number of particles
    save:       filename to save (e.g. 'traj.mp4') or None to just show
    """
    n_steps = len(t)

    # Compute overall bounds for a nice fixed axis range
    x = coords[:, :, 0]
    y = coords[:, :, 1]
    z = coords[:, :, 2]

    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    z_min, z_max = np.min(z), np.max(z)

    # Add a little padding so points aren't on the border
    pad_x = 0.05 * (x_max - x_min or 1.0)
    pad_y = 0.05 * (y_max - y_min or 1.0)
    pad_z = 0.05 * (z_max - z_min or 1.0)

    x_min, x_max = x_min - pad_x, x_max + pad_x
    y_min, y_max = y_min - pad_y, y_max + pad_y
    z_min, z_max = z_min - pad_z, z_max + pad_z

    # ---- Enforce equal aspect ratio in data space ----
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    max_range = max(x_range, y_range, z_range)

    x_mid = 0.5 * (x_max + x_min)
    y_mid = 0.5 * (y_max + y_min)
    z_mid = 0.5 * (z_max + z_min)

    x_min, x_max = x_mid - max_range / 2, x_mid + max_range / 2
    y_min, y_max = y_mid - max_range / 2, y_mid + max_range / 2
    z_min, z_max = z_mid - max_range / 2, z_mid + max_range / 2

    # Set up figure and 3D axis (only once)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    #ax.set_title("Particle trajectories")

    # If Matplotlib ≥ 3.3, this ensures a 1:1:1 box on screen too
    try:
        ax.set_box_aspect([1, 1, 1])
    except AttributeError:
        pass

    # One line + one point per particle; make marker color = line color
    lines = []
    points = []
    for i in range(n_particles):
        # Create the line first so it gets a color from the color cycle
        (line,) = ax.plot([], [], [], lw=1, label=f"Particle {i}")
        color = line.get_color()  # grab its color
        # Now create a point with the **same** color
        (point,) = ax.plot([], [], [], "o", color=color)
        lines.append(line)
        points.append(point)

    # Animation update function
    def update(frame):
        # Frame goes from 0 to n_steps-1
        for i in range(n_particles):
            x_i = coords[: frame + 1, i, 0]
            y_i = coords[: frame + 1, i, 1]
            z_i = coords[: frame + 1, i, 2]

            lines[i].set_data(x_i, y_i)
            lines[i].set_3d_properties(z_i)

            # Last point (current position)
            points[i].set_data(x_i[-1:], y_i[-1:])
            points[i].set_3d_properties(z_i[-1:])

        ax.set_title(f"Particle trajectories (t = {t[frame]:.3f})")
        return lines + points

    interval = 30
    anim = FuncAnimation(
        fig, update, frames=n_steps, interval=interval, blit=False, repeat=True
    )

    if save is not None:
        print(f"Saving animation to {save} ...")
        anim.save(save, fps=1000.0 / interval, dpi=150)
        print("Done.")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Animate particle trajectories from position-time data."
    )
    parser.add_argument(
        "filename",
        nargs="?",
        default="output_pos.dat",
        help="Input data file (default: output_pos.dat)",
    )
    parser.add_argument(
        "--save",
        metavar="FILE",
        help="Save animation to file instead of showing it (e.g. traj.mp4)",
    )

    args = parser.parse_args()

    t, coords, n_particles = load_data(args.filename)
    print(
        f"Loaded {len(t)} time steps, {n_particles} particles "
        f"from {args.filename!r}"
    )
    animate_trajectories(t, coords, n_particles, save=args.save)


if __name__ == "__main__":
    main()
