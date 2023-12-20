import MDAnalysis as mda
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.animation as animation

def get_pbc(water):
    """
    Calculate the size of the periodic boundary box (PBC) in the x and y directions based on the water molecule coordinates.

    Args:
        water (variable): A variable representing the water molecule coordinates.

    Returns:
        tuple: A tuple containing the size of the periodic boundary box in the x direction and the size of the periodic boundary box in the y direction.
    """
    box_dims = water.bbox()
    xmin, ymin = box_dims[0,0:2]
    xmax, ymax = box_dims[1,0:2]
    xsize = xmax-xmin
    ysize = ymax-ymin

    return xsize, ysize

def shift_and_wrap(positions, maxval):
    """
    Shifts the positions so that the minimum value becomes zero, and then wraps the shifted positions within the range of zero to the maximum value.

    Args:
        positions (array-like): An array of positions.
        maxval (int): The maximum value to wrap the positions within.

    Returns:
        array-like: An array of positions wrapped within the range of zero to maxval.
    """
    minval = np.min(positions)
    shifted = positions - minval
    wrapped = shifted % maxval
    return wrapped

def get_counts(membrane, water, leaflet, width):
    """
    Calculates the counts of particles in a 2D grid based on their positions in a simulation.

    Args:
        membrane (object): A variable representing the membrane in the simulation.
        water (object): A variable representing the water molecules in the simulation.
        leaflet (object): A variable representing the leaflet in the simulation.
        width (float): The width of each grid cell.

    Returns:
        numpy.ndarray: An array of shape (simlength, nxs, nys) containing the counts of particles in each grid cell for each frame in the simulation.
    """

    simlength = membrane.trajectory.n_frames-1
    xsize, ysize = get_pbc(water)
    nxs = int(xsize/width)
    nys = nxs
    all_counts = np.zeros((simlength, nxs, nys))
    for i in range(simlength):
        membrane.trajectory.next()
        xsize, ysize = get_pbc(water)
        pos = leaflet.positions
        initxs = pos[:,0]
        pos[:,0] = shift_and_wrap(initxs, xsize)
        pos[:,1] = shift_and_wrap(pos[:,1], ysize)
        cts, xs, ys = np.histogram2d(pos[:,0], pos[:,1], [nxs, nys])
        all_counts[i] = cts

    return all_counts

def animate_counts(all_counts, fname):
    """
    Create an animation of particle counts in a 2D grid based on the input array `all_counts`.

    Args:
        all_counts (numpy.ndarray): An array of shape (simlength, nxs, nys) containing the counts of particles in each grid cell for each frame in the simulation.
        fname (str): The filename of the output animation file.

    Returns:
        None. The function saves the animation as a video file.
    """
    fig, ax = plt.subplots()

    ims = []
    for counts in all_counts:
        im = ax.imshow(counts, animated=True)
        ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=1000)

    ani.save(fname)