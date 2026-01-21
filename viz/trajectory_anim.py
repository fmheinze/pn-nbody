import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import axes3d, Axes3D


def read_data(filename):
    # Get the number of bodies, the dimensions and number of timesteps
    with open(filename, 'r') as file:
        file.readline()     # Ignore first line
        line = file.readline().split()
        num_bodies = int(line[-1][-1]) + 1
        dim = (len(line) - 1) // (2 * num_bodies)
        num_timesteps = sum(1 for _ in file)
    
    # Initialize arrays to store data
    times = np.empty(num_timesteps)
    pos = np.empty((num_timesteps, num_bodies, dim))
    vel = np.empty((num_timesteps, num_bodies, dim))

    # Read the data from the file and store the values in the arrays
    with open(filename, 'r') as file:
        line = file.readline().split()
        masses = np.array([float(line[3*i]) for i in range(1, num_bodies + 1)])
        file.readline()     # Ignore the second line
        for n in range(num_timesteps):
            line = np.array([float(item) for item in file.readline().split()])
            times[n] = line[0]
            pos[n, :, :] = line[1:dim*num_bodies+1].reshape(num_bodies, dim)
            vel[n, :, :] = line[1+dim*num_bodies:].reshape(num_bodies, dim)
    return times, masses, pos, vel
    

def animate_trajectories(times, masses, pos, vel, com_frame=False, com_frame_bodies='all', skip_factor=1, lim=None):
    num_timesteps = pos.shape[0]
    num_bodies = pos.shape[1]
    dim = pos.shape[2]
    if com_frame and com_frame_bodies == 'all':
        com_frame_bodies = np.arange(0, num_bodies, 1)

    # Compute center of mass
    if com_frame:
        total_mass = np.sum(masses[com_frame_bodies])
        com = np.zeros_like(pos[:, 0, :])
        for i in com_frame_bodies:
            com += masses[i] * pos[:, i, :]
        com /= total_mass

    # Setting up the plot window
    fig = plt.figure(figsize=(6, 6), constrained_layout=True)
    if dim == 2:
        ax = fig.add_subplot()
        ax.set_facecolor('ghostwhite')
        ax.grid(ls='--')
    else:
        ax = fig.add_subplot(projection='3d')
    ax.set_aspect('equal')

    if lim is None:
        if com_frame:
            lim = 0.0
            for i in range(num_bodies):
                if(1.2 * np.max(np.abs(pos[:, i, :] - com)) > lim):
                    lim = 1.2 * np.max(np.abs(pos[:, i, :] - com))
        else:
            lim = 1.2 * np.max(np.abs(pos))
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    if dim == 3:
        ax.set_zlim([-lim, lim])
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if dim == 3:
        ax.set_zlabel('z')

    # Initializing the animated objects
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
              'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan']
    markersizes = [(2 * np.log(0.001 * (masses[i] + 0.5)) + 20) for i in range(num_bodies)]
    points = [plt.plot([], [], color=colors[i % len(colors)], ls='', marker='o', markersize=markersizes[i])[0] for i in range(num_bodies)]
    lines = [plt.plot([], [], color=colors[i % len(colors)], lw=2)[0] for i in range(num_bodies)]
    artists = [*points, *lines]

    def init():
        for i in range(2 * num_bodies):
            artists[i].set_data([], [])
            if dim == 3:
                artists[i].set_3d_properties([])
        return artists
    
    def update(n):
        n *= skip_factor
        for i in range(num_bodies):
            if com_frame:
                artists[i].set_data([pos[n, i, 0] - com[n, 0]], [pos[n, i, 1] - com[n, 1]])
                artists[num_bodies + i].set_data(pos[:n, i, 0] - com[:n, 0], pos[:n, i, 1] - com[:n, 1])
            else:
                artists[i].set_data([pos[n, i, 0]], [pos[n, i, 1]])
                artists[num_bodies + i].set_data(pos[:n, i, 0], pos[:n, i, 1])
            if dim == 3:
                if com_frame:
                    artists[i].set_3d_properties([pos[n, i, 2] - com[n, 2]])
                    artists[num_bodies + i].set_3d_properties([pos[:n, i, 2] - com[:n, 2]])
                else:
                    artists[i].set_3d_properties([pos[n, i, 2]])
                    artists[num_bodies + i].set_3d_properties([pos[:n, i, 2]])
        return artists
    
    ani = FuncAnimation(fig, update, init_func=init, interval=30, frames=int(num_timesteps/skip_factor), blit=True)
    plt.show()


def plot_trajectories(masses, pos, com_frame=False, com_frame_bodies='all', com_trajectories=None, plot_only_coms=True, lim=None):
    num_timesteps = pos.shape[0]
    num_bodies = pos.shape[1]
    dim = pos.shape[2]
    if com_frame and com_frame_bodies == 'all':
        com_frame_bodies = np.arange(0, num_bodies, 1)

    # Compute center of mass
    if com_frame:
        com = 0
        for i in com_frame_bodies:
            com += masses[i] * pos[:, i, :] / np.sum(masses[com_frame_bodies])
    
    if com_trajectories is not None:
        coms = np.zeros((num_timesteps, com_trajectories.shape[0], dim))
        for i in range(com_trajectories.shape[0]):
            for j in com_trajectories[i, :]:
                coms[:, i, :] += masses[j] * pos[:, j, :] / np.sum(masses[com_trajectories[i, :]])

    # Setting up the plot window
    fig = plt.figure(figsize=(6, 6), constrained_layout=True)
    if dim == 2:
        ax = fig.add_subplot()
        ax.set_facecolor('ghostwhite')
        ax.grid(ls='--')
    else:
        ax = fig.add_subplot(projection='3d')
    ax.set_aspect('equal')

    if lim is None:
        if com_frame:
            lim = 0.0
            for i in range(num_bodies):
                if(1.2 * np.max(np.abs(pos[:, i, :] - com)) > lim):
                    lim = 1.2 * np.max(np.abs(pos[:, i, :] - com))
        else:
            lim = 1.2 * np.max(np.abs(pos))
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    if dim == 3:
        ax.set_zlim([-1.0, 0.4])
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if dim == 3:
        ax.set_zlabel('z')

    print(pos[:, 0, 0])
    
    # Plotting the trajectories
    for i in range(num_bodies):
        if com_frame:
            if dim == 2:
                ax.plot(pos[:, i, 0] - com[:, 0], pos[:, i, 1] - com[:, 1], lw=2)
            else:
                ax.plot(pos[:, i, 0] - com[:, 0], pos[:, i, 1] - com[:, 1], pos[:, i, 2] - com[:, 2], lw=2)
        else:
            if dim == 2:
                ax.plot(pos[:, i, 0], pos[:, i, 1], lw=2)
            else:
                ax.plot(pos[:, i, 0], pos[:, i, 1], pos[:, i, 2], lw=2, color='tab:orange', visible=(not(i in com_trajectories and plot_only_coms)))

    # Plotting the center of mass trajectories
    if com_trajectories is not None:
        for i in range(com_trajectories.shape[0]):
            if com_frame:
                if dim == 2:
                    ax.plot(coms[:, i, 0] - com[:, 0], coms[:, i, 1] - com[:, 1], lw=2)
                else:
                    ax.plot(coms[:, i, 0] - com[:, 0], coms[:, i, 1] - com[:, 1], coms[:, i, 2] - com[:, 2], lw=2)
            else:
                if dim == 2:
                    ax.plot(coms[:, i, 0], coms[:, i, 1], lw=2)
                    pass
                else:
                    ax.plot(coms[:, i, 0], coms[:, i, 1], coms[:, i, 2], lw=2, color='tab:purple')
                    pass
    plt.show()


times, masses, pos, vel = read_data("/ssd/zo86hec/pn-nbody/output/hierarchical4/output_pos.dat")

animate_trajectories(times, masses, pos, vel, com_frame=False, com_frame_bodies=[0, 1], skip_factor=5, lim=20)

plot_trajectories(masses, pos, com_frame=False, com_frame_bodies=[0, 1], com_trajectories=np.array([[0, 1]]),
                  plot_only_coms=True, lim=20)
