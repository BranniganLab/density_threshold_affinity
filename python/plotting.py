#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:41:24 2024.

@author: js2746
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from matplotlib.colors import ListedColormap, Normalize
from Site import Site
# from Symmetric_Site import Symmetric_Site
from utils import calculate_hist_mode


def make_custom_colormap():
    """
    Make a custom colormap for plotting. The colormap starts at 0.35 to ensure \
    that there is a visual difference between 0 (no signal) and depleted (near-\
    zero enrichment).

    Returns
    -------
    my_cmap : matplotlib colormap
        The colormap to use when plotting.

    """
    depleted = plt.colormaps['RdGy_r']
    enriched = plt.colormaps['bwr']
    newcolors = np.concatenate([depleted(np.linspace(0.35, 0.5, 128)), enriched(np.linspace(0.5, 1, 128))])
    my_cmap = ListedColormap(newcolors)
    return my_cmap


def outline_site_new(ax, site, grid_dims):
    """
    Draw an outline around each bin in this Site.

    Parameters
    ----------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.
    site : Site or Symmetric_Site
        The site you want to outline.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.

    Returns
    -------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.

    """
    if isinstance(site, Site):
        for each_bin in site.bin_coords:
            ax = outline_bin(ax, each_bin, grid_dims)
    # elif isinstance(site, Symmetric_Site):
    #     for each_site in site.site_list:
    #         for each_bin in each_site.bin_coords:
    #             ax = outline_bin(ax, each_bin, grid_dims)
    else:
        Exception("site must be a Site or Symmetric_Site.")
    return ax


def outline_bin(ax, bin_coords, grid_dims):
    """
    Draw an outline around this bin.

    Parameters
    ----------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.
    bin_coords : tuple
        The tuple of bin coordinates stored in (r_bin, theta_bin) format.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.

    Returns
    -------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.

    """
    assert isinstance(bin_coords, tuple), "bin_coords must be a tuple of bin coordinates."
    dr, _, dtheta, _, _ = grid_dims
    start_theta = bin_coords[1] * dtheta
    end_theta = start_theta + dtheta
    inner_r = bin_coords[0] * dr
    outer_r = inner_r + dr
    theta_range = np.linspace(start_theta, end_theta, 100)
    ax.fill_between(theta_range, inner_r, outer_r, facecolor=(0, 0, 0, 0), edgecolor='k')
    return ax


def create_heatmap_figure_and_axes(lipids, cmap, v_vals, figwidth, figheight, helices):
    """
    Create the heatmap figure and enough axes to accommodate all the lipids and\
    leaflets.

    Parameters
    ----------
    lipids : list
        The names of the lipids you intend to plot.
    cmap : matplotlib ListedColormap
        A custom colormap.
    v_vals : tuple
        The (vmin, vmid, and vmax)
    figwidth : int, optional
        Figure width.
    figheight : int, optional
        Figure height.
    helices : list
        The outer and inner helix coordinates, in that order.

    Returns
    -------
    fig : matplotlib Fig object
        The figure you just created.
    matplotlib Axes objects
        The polar projection axes that were created.

    """
    assert isinstance(lipids, list), "lipids must be a list of strings"
    assert len(lipids) > 0, "lipids cannot be an empty list"
    assert isinstance(helices, list), "helices must be a list"
    assert len(helices) == 2, "helices must contain [outer helices, inner helices]"
    numlipids = len(lipids)
    vmin, vmid, vmax = v_vals
    fig = plt.figure(figsize=(figwidth, figheight))
    gs = gridspec.GridSpec(numlipids, 2, figure=fig, wspace=0.15, hspace=0.15)
    for gridbox in range(numlipids * 2):
        ax = plt.subplot(gs[gridbox], projection='polar')
        if gridbox % (numlipids * 2) == 0:
            # put the lipid name to the left of the axes object
            trans = mtransforms.ScaledTranslation(-40 / 72, -1.5, fig.dpi_scale_trans)
            print(gridbox)
            ax.text(0.0, 1.0, lipids[gridbox // (numlipids * 2)], transform=ax.transAxes + trans, fontsize='medium', va='bottom', fontfamily='serif')
        if gridbox == 0:
            ax.set_title("Outer")
        elif gridbox == 1:
            ax.set_title("Inner")
        if gridbox % 2 == 0:
            ax = plot_helices(helices[0], False, ax, 50)
        else:
            ax = plot_helices(helices[1], False, ax, 50)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.21, 1, 0.5, 0.02])
    sm = mpl.cm.ScalarMappable(cmap=cmap)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks(np.linspace(0, 1, 5))
    cbar.ax.set_xticklabels([vmin, (vmin + vmid) / 2, vmid, (vmid + vmax) / 2, vmax])
    return fig, fig.axes


def bin_prep(bin_info):
    """
    Configure the arrays needed for plotting polar heatmaps.

    Parameters
    ----------
    bin_info : namedtuple
        Contains Nr, Ntheta, dr, and dtheta information.

    Returns
    -------
    list
        The two numpy ndarrays needed for plotting a heatmap.

    """
    r_vals = np.linspace(0, bin_info.Nr * bin_info.dr, bin_info.Nr + 1)
    theta_vals = np.linspace(0, 2 * np.pi, bin_info.Ntheta + 1)
    r_vals, theta_vals = np.meshgrid(r_vals, theta_vals, indexing='ij')
    return [r_vals, theta_vals]


def plot_heatmap(ax, data, grid, cmap, v_vals):
    """
    Plot a heatmap on a pre-existing axes object.

    Parameters
    ----------
    ax : matplotlib.pyplot axes object
        The pre-existing axes object you wish to plot a heatmap on.
    data : ndarray
        The heatmap heat values (probably density enrichment).
    grid : 2-tuple of ndarrays
        Contains the radius and theta values for plotting the polar projection.
    cmap : colorbar object
        Custom colorbar.
    v_vals : 3-tuple
        (colorbar vmin, vmid, and vmax).

    Returns
    -------
    ax : matplotlib.pyplot axes object
        The axes object, which now contains your heatmap.

    """
    vmin, vmid, vmax = v_vals
    norm = MidpointNormalize(midpoint=vmid, vmin=vmin, vmax=vmax)
    ax.grid(False)
    plt.axis('off')
    radius, theta = grid
    ax.pcolormesh(theta, radius, data, cmap=cmap, norm=norm, zorder=0, edgecolors='face', linewidth=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    return ax


def plot_histogram(ax, data, area, bulk_mode="NULL", plot_probability=False):
    """
    Plot a histogram.

    Parameters
    ----------
    ax : matplotlib.pyplot axes object
        The pre-existing axes object you wish to plot a histogram on.
    data : ndarray or list
        The histogram counts.
    area : float
        The accessible area of the site.
    bulk_mode : float or "NULL", optional
        If not "NULL", add a dashed red line showing n_peak. The default is "NULL".
    plot_probability : boolean, optional
        If True, turn the y axis into probability percentages, rather than \
        raw counts. The default is False.

    Returns
    -------
    ax : matplotlib.pyplot axes object
        The axes object, now with your histogram plotted on it.

    """
    if plot_probability:
        data = data / np.sum(data)
    ax.plot(range(len(data)), data)
    ax.set_ylabel("Probability")
    ax.set_xlabel(f"Number of beads in an area about {area} " + r"$\AA^2$")
    mode = calculate_hist_mode(data)
    ax.vlines([mode], 0, np.max(data), color='black', linestyles='dashed', label=f"mode={mode}")
    if bulk_mode != "NULL":
        ax.vlines([bulk_mode], 0, np.max(data), color='red', linestyles='dashed', label=f"bulk mode={bulk_mode}")
    ax.legend()
    return ax


# This class comes from Liam Sharp and could potentially be rewritten to be
# more clear.
class MidpointNormalize(Normalize):
    """
    Normalise the colorbar so that diverging bars work their way either side from a prescribed midpoint value).

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def plot_helices(helices, colorbychain, ax, markersize=3, sub=["tab:blue", "tab:cyan", "tab:green", "tab:purple", "tab:brown", "tab:olive"]):
    """
    Plot helices on a polar plot.

    Parameters
    ----------
    - helices (list or array): The helix data to be plotted. Each element in the list represents a set of helices, where each set is a list of angles and radii.
    - colorbychain (bool): A flag to determine if the helices should be colored by chain.
    - ax (matplotlib.axes.Axes): The polar subplot axis on which to plot the helices.
    - markersize (int, optional): The size of the scatter markers. Defaults to 3.
    - sub (list, optional): The list of colors to use for coloring the helices. Defaults to a predefined list of colors.

    Returns
    -------
    - ax (matplotlib.axes.Axes): The modified polar subplot axis with the helices plotted.

    """
    if len(np.shape(helices)) == 1:
        helices = np.reshape(helices, (1, len(helices)))
    for i, pro in enumerate(helices[:]):
        if colorbychain:
            colors = sub[i]
        else:
            colors = sub[:len(pro[::2])]
        ax.scatter(np.deg2rad(pro[1::2]), pro[::2], color=colors, linewidth=None,
                   zorder=1, s=markersize)
    return ax
