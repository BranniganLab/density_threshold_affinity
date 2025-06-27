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
from matplotlib.colors import ListedColormap, Normalize
from scipy import constants
from DTA.utils import calculate_hist_mode, load_inclusion_helices
from DTA.Site import Site
from DTA.SymmetricSite import SymmetricSite
from DTA.density import parse_tcl_dat_file, calculate_density_enrichment, calculate_density
from pathlib import Path


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


def outline_site(axes, site, grid_dims, index, linewidth=1, color='black'):
    """
    Draw an outline around a Site or around each site in a SymmetricSite.

    Parameters
    ----------
    axes : list of matplotlib.pyplot Axes objects
        The list of Axes objects containing plots.
    site : Site or Symmetric_Site
        The site you want to outline.
    grid_dims : namedtuple
        Contains dr, number of r bins, dtheta, number of theta bins, and number\
        of frames contained in file.
    index : int
        The index of the axes list corresponding to the plot you want to draw \
        a site outline on.
    linewidth : float, optional
        The width of the outline. Default is 1.
    color : str, optional
        The color you want to outline the site in. Default is 'black'.

    Returns
    -------
    axes : list of matplotlib.pyplot Axes objects
        The Axes objects, now with a site outlined.

    """
    assert isinstance(site, (Site, SymmetricSite)), "site must be a Site or SymmetricSite."
    assert site.bin_coords is not None, "Site must be fully defined first. Please add bin_coords."
    assert isinstance(index, int), "index must be an int"
    assert (index < len(axes)), f"index {index} exceeds length of axes list."
    edges = isolate_unique_site_edges(site.bin_coords, grid_dims)
    for edge_tuple in edges:
        r, theta = edge_tuple
        axes[index].plot(theta, r, color=color, linewidth=linewidth, marker=None)
    return axes


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
    for index in range(len(line_list)):
        if np.allclose(edge, line_list[index]):
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
    figwidth : int
        Figure width.
    figheight : int
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
        if gridbox == 0:
            ax.set_title("Outer")
        elif gridbox == 1:
            ax.set_title("Inner")
        if gridbox % 2 == 0:
            ax = plot_helices(helices[0], False, ax, 50)
        else:
            ax = plot_helices(helices[1], False, ax, 50)
        if gridbox % 2 == 0:
            # put the lipid name to the left of the axes object
            ax.text(-0.5, 0.5, lipids[gridbox // 2], transform=ax.transAxes, fontsize='medium', va='center', fontfamily='serif')
    return fig, fig.axes


def make_colorbar(fig, v_vals, cmap):
    """
    Generate the colorbar. Must be done after plot_heatmap.

    Parameters
    ----------
    fig : fig
        The figure object from matplotlib.
    v_vals : list or tuple
        min, mid, and max values for the colorbar.
    cmap : colormap
        matplotlib colormap object.

    Returns
    -------
    fig

    """
    vmin, vmid, vmax = v_vals
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.21, 1, 0.5, 0.02])
    sm = mpl.cm.ScalarMappable(cmap=cmap)
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


def plot_heatmap(ax, data, grid_dims, cmap, v_vals):
    """
    Plot a heatmap on a pre-existing axes object.

    Parameters
    ----------
    ax : matplotlib.pyplot axes object
        The pre-existing axes object you wish to plot a heatmap on.
    data : ndarray
        The heatmap heat values (probably density enrichment).
    grid_dims : namedtuple
        Contains Nr, Ntheta, dr, and dtheta information.
    cmap : colorbar object
        Custom colorbar.
    v_vals : 3-tuple
        (colorbar vmin, vmid, and vmax).

    Returns
    -------
    ax : matplotlib.pyplot axes object
        The axes object, which now contains your heatmap.

    """
    grid = bin_prep(grid_dims)
    vmin, vmid, vmax = v_vals
    norm = MidpointNormalize(midpoint=vmid, vmin=vmin, vmax=vmax)
    ax.grid(False)
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
        ax.set_ylabel("Probability")
    else:
        ax.set_ylabel("Counts")
    ax.plot(range(len(data)), data)
    ax.set_xlabel(f"Number of beads in an area about {area} " + r"$\AA^2$")
    mode = calculate_hist_mode(data)
    ax.vlines([mode], 0, np.max(data), color='black', linestyles='dashed', label=f"mode={mode}")
    if bulk_mode != "NULL":
        ax.vlines([bulk_mode], 0, np.max(data), color='red', linestyles='dashed', label=f"bulk mode={bulk_mode}")
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


def plot_helices(helices, colorbychain, ax, markersize=3, sub=["tab:blue", "tab:cyan", "tab:green", "tab:purple", "tab:brown", "tab:olive", "tab:orange"]):
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


def make_density_enrichment_heatmap(system_names, my_cmap, max_enrichment, helix_definitions, figdims, leaflets, root_path, replicas=None):
    """
    Make a figure and axes objects. Plot heatmaps of density enrichment for each\
    system on each axes object. Return the figure, axes, and lattice grid dimensions.

    Parameters
    ----------
    system_names : list
        A list of the different system names. Must correspond to system name \
        provided to PolarDensityBin.
    my_cmap : cmap
        The colormap to use.
    max_enrichment : float or int
        How high you want your colorbar to go. The minimum will scale proportionally.
    helix_definitions : str or Path
        The directory containing your helix coordinate outputs from PolarDensityBin.
    figdims : 2-tuple
        The figure height and width, in inches.
    leaflets : list
        A list containing "outer", "inner", or both, depending on what you want\
        to plot.
    root_path : str or Path
        The path to the directory containing your PolarDensityBin outputs (if \
        replicas is None), or the directory containing your replica directories.
    replicas : None or list, optional
        If you have multiple replicas that you would like averaged together, \
        ensure that the files are contained within separate subdirectories of \
        root_path. Provide the subdirectory names as a list of strs. The default\
        is None, which turns off this feature.

    Returns
    -------
    fig1 : Figure object
        That matplotlib Figure object containing your plots.
    axes : Axes object or list of Axes objects
        The matplotlib Axes object(s) containing your plot(s).
    grid_dims : named tuple
        The lattice dimension information; used by other functions.

    """
    assert isinstance(system_names, list), "system_names must be a list, even if it only contains one item."
    assert isinstance(leaflets, list), "leaflets must be a list, even if it only contains one item."
    assert (replicas is None) or (isinstance(replicas, list)), "replicas must either be None or a list."
    assert isinstance(helix_definitions, (str, Path)), "helix_definitions must be a str or Path object."
    if isinstance(helix_definitions, str):
        helix_definitions = Path(helix_definitions)
    assert helix_definitions.exists(), "helix_definitions not found."
    assert isinstance(root_path, (str, Path)), "root_path must be a str or Path object."
    if isinstance(root_path, str):
        root_path = Path(root_path)
    assert root_path.exists(), "root_path not found."
    fig_h, fig_w = figdims
    colorbar_range = (1 / max_enrichment, 1, max_enrichment)
    helices = load_inclusion_helices(helix_definitions)
    fig1, axes = create_heatmap_figure_and_axes(system_names, my_cmap, colorbar_range, figwidth=fig_w, figheight=fig_h, helices=helices)
    index = 0
    for species in system_names:
        for leaf in leaflets:
            if replicas:
                rep_list = []
                for rep in replicas:
                    rep_path = root_path.joinpath(rep, f"{species}.{leaf}.avg.dat")
                    counts, grid_dims, system_info = parse_tcl_dat_file(rep_path, bulk=False)
                    density_enrichment = calculate_density_enrichment(calculate_density(counts, grid_dims), system_info.ExpBeadDensity)
                    rep_list.append(density_enrichment)
                avg_enrichment = np.mean(np.stack(tuple(rep_list), axis=0), axis=0)
            else:
                path = root_path.joinpath(f"{species}.{leaf}.avg.dat")
                counts, grid_dims, system_info = parse_tcl_dat_file(path, bulk=False)
                density_enrichment = calculate_density_enrichment(calculate_density(counts, grid_dims), system_info.ExpBeadDensity)
                avg_enrichment = density_enrichment
            axes[index] = plot_heatmap(axes[index], avg_enrichment, grid_dims, my_cmap, colorbar_range)
            index += 1
    fig1 = make_colorbar(fig1, colorbar_range, my_cmap)
    return fig1, axes, grid_dims
