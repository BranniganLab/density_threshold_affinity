#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:41:24 2024.

@author: js2746
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, Normalize
from scipy import constants
from DTA.utils import calculate_hist_mode
from DTA.Site import Site
from DTA.SymmetricSite import SymmetricSite


class HeatmapSettings:
    """
    Class for heatmap Figure settings in DTA.

    Attributes
    ----------
    row_names : list
        A list of the names you wish to appear to the left of each row of \
        heatmaps. Frequently will be the name of the lipid species.
    col_names : list
        A list of the names you wish to appear above each column of heatmaps. \
        Frequently will be "outer leaflet" and "inner leaflet".
    fig_height : float
        Figure height, in inches.
    fig_width : float
        Figure width, in inches.
    colormap : matplotlib cmap object
        The colormap to use.
    max_enrichment : float or int
        How high you want your colorbar to go. The minimum will scale proportionally.
    r_vals, theta_vals : numpy ndarrays
        The numpy meshgrids needed to plot a heatmap using polar coordinates.
    """

    def __init__(self, row_names, col_names, fig_dims, colormap, max_enrichment, grid_dims=None):
        """
        Create a HeatmapSettings object.

        Parameters
        ----------
        row_names : list
            A list of the names you wish to appear to the left of each row of \
            heatmaps. Frequently will be the name of the lipid species.
        col_names : list
            A list of the names you wish to appear above each column of heatmaps. \
            Frequently will be "outer leaflet" and "inner leaflet".
        fig_dims : 2-tuple or 2-list
            The height and width, in inches, of your Figure.
        colormap : matplotlib cmap object
            The colormap to use.
        max_enrichment : float or int
            How high you want your colorbar to go. The minimum will scale proportionally.
        grid_dims : namedtuple, optional
            Contains Nr, Ntheta, dr, and dtheta information.. The default is None.

        Raises
        ------
        TypeError
            Makes sure that row_names and col_names are lists.

        Returns
        -------
        None.

        """
        if not isinstance(col_names, list):
            raise TypeError(f"{col_names} must be a list instead of a {type(col_names)}.")
        if not isinstance(row_names, list):
            raise TypeError(f"{row_names} must be a list instead of a {type(row_names)}.")
        self.row_names = row_names
        self.col_names = col_names
        self.fig_height = fig_dims[0]
        self.fig_width = fig_dims[1]
        self.colormap = colormap
        self.colorbar_range = (1 / max_enrichment, 1, max_enrichment)
        if grid_dims is not None:
            self.r_vals, self.theta_vals = bin_prep(grid_dims)
        else:
            self.r_vals, self.theta_vals = None, None

    def add_grid_dims(self, grid_dims):
        """
        Calculate the meshgrids needed to plot a heatmap in polar coordinates.

        Parameters
        ----------
        grid_dims : namedtuple
            Contains Nr, Ntheta, dr, and dtheta information.. The default is None.

        Returns
        -------
        None.

        """
        self.r_vals, self.theta_vals = bin_prep(grid_dims)


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


def outline_site(ax, site, grid_dims, linewidth=1, color='black'):
    """
    Draw an outline around a Site or around each site in a SymmetricSite.

    Parameters
    ----------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.
    site : Site or Symmetric_Site
        The site you want to outline.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.
    linewidth : float, optional
        The width of the outline. Default is 1.
    color : str, optional
        The color you want to outline the site in. Default is 'black'.

    Returns
    -------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.

    """
    assert isinstance(site, (Site, SymmetricSite)), "site must be a Site or SymmetricSite."
    assert site.bin_coords is not None, "Site must be fully defined first. Please add bin_coords."
    edges = isolate_unique_site_edges(site.bin_coords, grid_dims)
    for edge_tuple in edges:
        r, theta = edge_tuple
        ax.plot(theta, r, color=color, linewidth=linewidth, marker=None)
    return ax


def isolate_unique_site_edges(bin_coords_list, grid_dims):
    """
    List all of the edges that need to be drawn in order to enclose the site.\
    In this lattice space, if edges are repeated it means they are internal \
    edges and should not be drawn at all.

    Parameters
    ----------
    bin_coords_list : list or tuples
        The bins that belong to this site in (r, theta) format. e.g. \
        [(2, 10), (2, 11), (2, 12)] would correspond to the 11th, 12th, and \
        13th theta bins in the 3rd radial bin from the origin. Bin coordinates \
        are zero-indexed by convention.
        The list of bin.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.

    Returns
    -------
    line_list : list of tuples
        The list of exterior bin edges needed in order to draw the site outline.\
        Each bin edge is a tuple, with each tuple value being an ndarray.

    """
    line_list = []
    for coord_pair in bin_coords_list:
        edges = compile_bin_edges(coord_pair, grid_dims)
        for edge in edges:
            assert isinstance(edge, tuple), "edge must be a tuple"
            index = _find_edge_in_list(edge, line_list)
            if index != -1:
                line_list.pop(index)
            else:
                line_list.append(edge)
    return line_list


def _find_edge_in_list(edge, line_list):
    """
    Iterate through list and compare elements. Need to do this manually rather\
    than with 'in' because list contents are numpy ndarrays.

    Parameters
    ----------
    edge : tuple of ndarrays
        The line you are looking for in line_list.
    line_list : list of tuples of ndarrays
        The list of all lines to be potentially drawn.

    Returns
    -------
    int
        The index in line_list at which edge is found.

    """
    for index, line in enumerate(line_list):
        if np.allclose(edge, line):
            return index
    return -1


def compile_bin_edges(bin_coords, grid_dims):
    """
    Determine the 4 lines that outline this bin.

    Parameters
    ----------
    bin_coords : tuple
        The tuple of bin coordinates stored in (r_bin, theta_bin) format.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.

    Returns
    -------
    tuple of tuples
        The four lines corresponding to the bin edges. Each tuple contains two\
        ndarrays of equal lengths.

    """
    assert isinstance(bin_coords, tuple), "bin_coords must be a tuple of bin coordinates."
    start_theta = bin_coords[1] * grid_dims.dtheta
    end_theta = start_theta + grid_dims.dtheta
    inner_r = bin_coords[0] * grid_dims.dr
    outer_r = inner_r + grid_dims.dr
    theta_range = np.linspace(start_theta, end_theta, 100)
    r_range = np.linspace(inner_r, outer_r, 100)
    line1 = (np.linspace(inner_r, inner_r, 100), theta_range)
    line2 = (np.linspace(outer_r, outer_r, 100), theta_range)
    if np.allclose(start_theta, 2 * np.pi):
        start_theta = 0
    if np.allclose(end_theta, 2 * np.pi):
        end_theta = 0
    line3 = (r_range, np.linspace(start_theta, start_theta, 100))
    line4 = (r_range, np.linspace(end_theta, end_theta, 100))
    return (line1, line2, line3, line4)


def create_heatmap_figure_and_axes(heatmap_settings):
    """
    Create the heatmap figure and enough axes to accommodate all the lipids and\
    leaflets.

    Parameters
    ----------
    heatmap_settings : HeatmapSettings object
        All of the auxiliary (non-data) information needed to plot a heatmap.

    Returns
    -------
    fig : matplotlib Fig object
        The figure you just created.

    """
    num_rows = len(heatmap_settings.row_names)
    num_cols = len(heatmap_settings.col_names)

    fig = plt.figure(figsize=(heatmap_settings.fig_width, heatmap_settings.fig_height))
    gs = gridspec.GridSpec(num_rows, num_cols, figure=fig, wspace=0.15, hspace=0.15)
    for gridbox in range(num_rows * num_cols):
        ax = plt.subplot(gs[gridbox], projection='polar')
        if gridbox < num_cols:
            ax.set_title(heatmap_settings.col_names[gridbox], fontsize='medium')
        if gridbox % num_cols == 0:
            # put the row name to the left of the axes object
            ax.text(-0.5, 0.5, heatmap_settings.row_names[gridbox // num_cols], transform=ax.transAxes, fontsize='medium', va='center', fontfamily='serif')
    return fig


def plot_helices_on_panels(fig, helices):
    """
    Plot helix locations on each panel present in the figure.

    Parameters
    ----------
    fig : matplotlib Figure object
        The Figure containing your heatmap panels.
    helices : list of numpy ndarrays
        List containing one set of helix coordinates per panel.

    Returns
    -------
    fig : matplotlib Figure
        The Figure containing your heatmap panels, now with helix locations.

    """
    if not isinstance(helices, list):
        raise TypeError(f"{helices} must be a list instead of a {type(helices)}.")
    if not all(isinstance(item, np.ndarray) for item in helices):
        raise TypeError("helices must be a list of ndarrays")
    assert len(helices) == np.ravel(fig.axes).shape[0]
    for ax, helix_set in zip(np.ravel(fig.axes), helices):
        ax = plot_helices(helix_set, False, ax, 50)
    return fig


def make_colorbar(fig, heatmap_settings):
    """
    Generate the colorbar. Must be done after plot_heatmap.

    Parameters
    ----------
    fig : fig
        The figure object from matplotlib.
    heatmap_settings : HeatmapSettings object
        All the auxiliary information needed to plot a heatmap.

    Returns
    -------
    fig

    """
    vmin, vmid, vmax = heatmap_settings.colorbar_range
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.21, 1, 0.5, 0.02])
    sm = mpl.cm.ScalarMappable(cmap=heatmap_settings.colormap)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks(np.linspace(0, 1, 3))
    cbar.ax.set_xticklabels([round(vmin, 2), vmid, vmax])
    return fig


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


def plot_heatmap(ax, data, heatmap_settings):
    """
    Plot a heatmap on a pre-existing axes object.

    Parameters
    ----------
    ax : matplotlib.pyplot axes object
        The pre-existing axes object you wish to plot a heatmap on.
    data : ndarray
        The heatmap heat values (probably density enrichment).
    heatmap_settings : HeatmapSettings object
        All the auxiliary information needed to plot a heatmap.

    Returns
    -------
    ax : matplotlib.pyplot axes object
        The axes object, which now contains your heatmap.

    """
    vmin, vmid, vmax = heatmap_settings.colorbar_range
    norm = MidpointNormalize(midpoint=vmid, vmin=vmin, vmax=vmax)
    ax.grid(False)
    ax.pcolormesh(heatmap_settings.theta_vals, heatmap_settings.r_vals, data, cmap=heatmap_settings.colormap, norm=norm, zorder=0, edgecolors='face', linewidth=0)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    return ax


def plot_histogram(ax, data, area, bulk_mode="NULL", plot_probability=False, show_target=False):
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
    show_target : boolean, optional
        If bulk_mode is True and show_target is True, label the red line "target\
        mode". If bulk_mode is True and show_target is False, label the red line\
        "bulk mode".

    Returns
    -------
    ax : matplotlib.pyplot axes object
        The axes object, now with your histogram plotted on it.

    """
    if plot_probability:
        data = data / np.sum(data)
        ax.set_ylabel("Probability")
    else:
        ax.set_ylabel("Counts")
    ax.plot(range(len(data)), data)
    ax.set_xlabel(f"Number of beads in an area about {area} " + r"$\AA^2$")
    mode = calculate_hist_mode(data)
    ax.vlines([mode], 0, np.max(data), color='black', linestyles='dashed', label=f"mode={mode}")
    if bulk_mode != "NULL":
        if show_target:
            ax.vlines([bulk_mode], 0, np.max(data), color='red', linestyles='dotted', label=f"target mode={bulk_mode}")
        else:
            ax.vlines([bulk_mode], 0, np.max(data), color='red', linestyles='dotted', label=f"bulk mode={bulk_mode}")
    ax.legend()
    return ax


def plot_titration_curve(ax, deltaG, deltaG_std, temperature, label, plot_error=True, error_type='std', n_replicas=None, color=None, linestyle='solid'):
    """
    Plot a titration curve on an existing figure.

    Parameters
    ----------
    ax : matplotlib Axes object
        The axis you want this plot to appear on.
    deltaG : float
        The deltaG_bind, in kcal/mol.
    deltaG_std : float
        The standard deviation of the mean for your deltaG.
    temperature : float
        The temperature of your system.
    label : str
        The label to add to the legend.
    plot_error : boolean, optional
        If True, also plot the error as a shaded region around the titration \
        curve. If False, do not plot the error; error_type and n_replicas are \
        ignored with this option. The default is True.
    error_type : str, optional
        Use 'std' for standard deviation. This is the default.
        Use 'ste' for standard error. User must specify n_replicas with this \
        option.
        Use 'CI' for a 95% confidence interval. User must specify n_replicas \
        with this option.
    n_replicas : int, optional
        The number of samples/replicas in your deltaG calculation. Used to \
        calculate the standard error and/or confidence interval. The default is\
        None.
    color : str, optional
        A matplotlib-recognizeable color string. Applied to curve. Default is \
        None, which results in matplotlib choosing your color from its default\
        palette.
    linestyle : str, optional
        A matplotlib-recognizable linestyle string. Applied to curve. Default is\
        'solid'.

    Returns
    -------
    ax : matplotlib Axes object
        The axis, now with your plot drawn on it.

    """
    RT = temperature * constants.R / 4184.  # 4184 converts J to kcal
    mol_pcts = np.linspace(0, 1, 10000)
    P_occ = mol_pcts / (np.exp(deltaG / RT) + mol_pcts)
    if color is not None:
        ax.plot(mol_pcts, P_occ, label=label, linewidth=3, linestyle=linestyle, color=color)
    else:
        ax.plot(mol_pcts, P_occ, label=label, linewidth=3, linestyle=linestyle)
    if plot_error:
        assert error_type in ['std', 'ste', 'CI'], "error_type must be std, ste, or CI."
        error = deltaG_std
        if error_type in ['ste', 'CI']:
            assert n_replicas is not None, "Need n_replicas to calculate standard error or confidence interval"
            error = error / np.sqrt(n_replicas)
            if error_type == 'CI':
                error = error * 1.96
        lwr_bound = mol_pcts / (np.exp((deltaG - error) / RT) + mol_pcts)
        upr_bound = mol_pcts / (np.exp((deltaG + error) / RT) + mol_pcts)
        if color is not None:
            ax.fill_between(mol_pcts, lwr_bound, upr_bound, alpha=0.2, color=color)
        else:
            ax.fill_between(mol_pcts, lwr_bound, upr_bound, alpha=0.2)
    return ax


def calc_x_50(deltaG, error, temperature):
    """
    Calculate an x_50 and error.

    Parameters
    ----------
    deltaG : float
        The deltaG_bind, in kcal/mol.
    error : float
        The error (can be standard deviation, standard error, etc) for your deltaG.
    temperature : float
        The temperature of your system.

    Returns
    -------
    x_50 : float
        The mol % at which you can expect your site to be occupied 50% of the time.
    err_on_x_50 : float
        The error.

    """
    RT = temperature * constants.R / 4184.  # 4184 converts J to kcal
    x_50 = np.exp(deltaG / RT)
    err_on_x_50 = x_50 - np.exp((deltaG - error) / RT)
    return x_50, err_on_x_50


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


def plot_helices(helices, colorbychain, ax, markersize=3, colorlist=None):
    """
    Plot helices on a polar plot.

    Parameters
    ----------
    - helices (list or array): The helix data to be plotted. Each element in the list represents a set of helices, where each set is a list of angles and radii.
    - colorbychain (bool): A flag to determine if the helices should be colored by chain.
    - ax (matplotlib.axes.Axes): The polar subplot axis on which to plot the helices.
    - markersize (int, optional): The size of the scatter markers. Defaults to 3.
    - colorlist (list, optional): The list of colors to use for coloring the helices. Defaults to a predefined list of colors.

    Returns
    -------
    - ax (matplotlib.axes.Axes): The modified polar subplot axis with the helices plotted.

    """
    if colorlist is None:
        colorlist = ["tab:blue", "tab:cyan", "tab:green", "tab:purple", "tab:brown", "tab:olive", "tab:orange"]
    if len(np.shape(helices)) == 1:
        helices = np.reshape(helices, (1, len(helices)))
    for i, pro in enumerate(helices[:]):
        if colorbychain:
            colors = colorlist[i]
        else:
            colors = colorlist[:len(pro[::2])]
        ax.scatter(np.deg2rad(pro[1::2]), pro[::2], color=colors, linewidth=None,
                   zorder=1, s=markersize)
    return ax


def make_density_enrichment_heatmap(enrichments_list, helices, heatmap_settings):
    """
    Make a figure and axes objects. Plot heatmaps of density enrichment for each\
    system on each axes object. Return the figure and axes.

    Parameters
    ----------
    enrichments_list : list
        A list of 2d ndarrays containing enrichment values for each bin in the\
        lattice. One list item per heatmap.
    helices : list of ndarrays
        Each ndarray in the list contains helix coordinates. There should be one\
        ndarray per heatmap.
    heatmap_settings : HeatmapSettings object
        Contains all of the auxiliary information needed to plot heatmaps.

    Returns
    -------
    fig : Figure object
        The matplotlib Figure object containing your plots.
    axes : Axes object or list of Axes objects
        The matplotlib Axes object(s) containing your plot(s).

    """
    if not isinstance(enrichments_list, list):
        raise TypeError(f"{enrichments_list} must be a list instead of a {type(enrichments_list)}.")
    if not all(isinstance(item, np.ndarray) for item in enrichments_list):
        raise TypeError("enrichments_list must be a list of ndarrays")
    if not all(len(arr.shape) == 2 for arr in enrichments_list):
        raise TypeError("enrichments_list must contain 2d ndarrays")

    fig = create_heatmap_figure_and_axes(heatmap_settings)
    axes = np.ravel(fig.axes)
    if len(enrichments_list) != axes.shape[0]:
        raise IndexError(f"Number of enrichments_list items ({len(enrichments_list)}) does not match number of figure panels ({axes.shape[0]}).")

    fig = plot_helices_on_panels(fig, helices)
    for index, ax in enumerate(axes):
        ax = plot_heatmap(ax, enrichments_list[index], heatmap_settings)
    fig = make_colorbar(fig, heatmap_settings)

    return fig, axes
