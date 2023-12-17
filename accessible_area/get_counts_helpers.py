import MDAnalysis as mda
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import matplotlib.animation as animation
import nglview as nv

def get_pbc(water):
    box_dims = water.bbox()
    xmin, ymin = box_dims[0,0:2]
    xmax, ymax = box_dims[1,0:2]
    xsize = xmax-xmin
    ysize = ymax-ymin

    return xsize, ysize

def shift_and_wrap(positions, max):
    min = np.min(positions)
    shifted = positions-min
    wrapped = shifted%max
    return wrapped

def get_counts(membrane, width):
    simlength = membrane.trajectory.n_frames-1
    xsize, ysize = get_pbc(water)
    nxs = int(xsize/width)
    nys = nxs
    all_counts = np.zeros((simlength, nxs, nys))
    for i in range(simlength):
        membrane.trajectory.next()
        xsize, ysize = get_pbc(water)
        pos = upper.positions
        initxs = pos[:,0]
        pos[:,0] = shift_and_wrap(initxs, xsize)
        pos[:,1] = shift_and_wrap(pos[:,1], ysize)
        cts, xs, ys = np.histogram2d(pos[:,0], pos[:,1], [nxs, nys])
        all_counts[i] = cts

    return all_counts

def animate_counts(all_counts, fname):
    fig, ax = plt.subplots()

    ims = []
    for counts in all_counts:
        im = ax.imshow(counts, animated=True)
        ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=1000)

    ani.save(fname)