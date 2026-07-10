#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:41:24 2024.

@author: js2746
"""
from dataclasses import dataclass, InitVar, field
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib.colors import ListedColormap, Normalize
from scipy import constants
from dta.utils import calculate_hist_mode
from dta.Site import Site
from dta.SymmetricSite import SymmetricSite
from dta.SiteAcrossReplicas import SiteAcrossReplicas
from dta.bin_logic import PolarBinGrid


@dataclass
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
    fig_dims : 2-tuple
        Figure height and width, in inches.
    grid : PolarBinGrid
        Contains lattice information.
    colormap : matplotlib cmap object
        The colormap to use. Provided after init.
    colorbar_range : 3-tuple
        The min, mid, and max values of your colorbar. Automatically calculated\
        from max_enrichment InitVar.

    InitVars
    --------
    max_enrichment : float
        Gets passed to post_init and turned into colorbar_range.
    """

    row_names: list
    col_names: list
    fig_dims: tuple
    max_enrichment: InitVar[float]
    grid: PolarBinGrid
    colormap: mpl.colors.ListedColormap = field(init=False)
    colorbar_range: tuple = field(init=False)

    def __post_init__(self, max_enrichment: float):
        """
        Calculate colorbar_range and make sure row_names and col_names are lists.

        Parameters
        ----------
        max_enrichment : float
            How high you want your colorbar to go. The minimum will scale proportionally.

        Raises
        ------
        TypeError
            Makes sure that row_names and col_names are lists.

        Returns
        -------
        None.

        """
        if not isinstance(self.col_names, list):
            raise TypeError(f"{self.col_names} must be a list instead of a {type(self.col_names)}.")
        if not isinstance(self.row_names, list):
            raise TypeError(f"{self.row_names} must be a list instead of a {type(self.row_names)}.")
        if not isinstance(self.grid, PolarBinGrid):
            raise TypeError(f"grid must be a PolarBinGrid instead dof a {type(self.grid)}.")
        if max_enrichment <= 1:
            raise ValueError("max_enrichment must be greater than 1.")
        self.colorbar_range = (1 / max_enrichment, 1, max_enrichment)


def make_custom_colormap() -> ListedColormap:
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


def outline_site(
    ax: plt.Axes,
    site: Site | SymmetricSite | SiteAcrossReplicas,
    linewidth: float = 1,
    color: str = 'black'
) -> plt.Axes:
    """
    Draw an outline around a Site or around each site in a SymmetricSite.

    Parameters
    ----------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.
    site : Site, SymmetricSite, or SiteAcrossReplicas
        The site you want to outline.
    linewidth : float, optional
        The width of the outline. Default is 1.
    color : str, optional
        The color you want to outline the site in. Default is 'black'.

    Returns
    -------
    ax : matplotlib.pyplot Axes object
        The Axes object you want to draw this on.

    """
    if not isinstance(site, (Site, SymmetricSite)):
        if isinstance(site, SiteAcrossReplicas):
            site = site.base_site
        else:
            raise TypeError("site must be a Site, SymmetricSite, or SiteAcrossReplicas.")
    if site.bin_coords is None:
        raise ValueError("Site must be fully defined first. Please add bin_coords.")
    edges = site.grid.list_all_exposed_edges(site.bin_coords)
    for edge in edges:
        theta = (edge.endpoint1[1], edge.endpoint2[1])
        r = (edge.endpoint1[0], edge.endpoint2[0])
        ax.plot(theta, r, color=color, linewidth=linewidth, marker=None)
    return ax


def create_heatmap_figure_and_axes(heatmap_settings: HeatmapSettings) -> plt.Figure:
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

    fig = plt.figure(figsize=(heatmap_settings.fig_dims[1], heatmap_settings.fig_dims[0]))
    gs = gridspec.GridSpec(num_rows, num_cols, figure=fig, wspace=0.15, hspace=0.15)
    for gridbox in range(num_rows * num_cols):
        ax = plt.subplot(gs[gridbox], projection='polar')
        if gridbox < num_cols:
            ax.set_title(heatmap_settings.col_names[gridbox], fontsize='medium')
        if gridbox % num_cols == 0:
            # put the row name to the left of the axes object
            ax.text(-0.5, 0.5, heatmap_settings.row_names[gridbox // num_cols], transform=ax.transAxes, fontsize='medium', va='center', fontfamily='serif')
    return fig


def plot_helices_on_panels(fig: plt.Figure, helices: list[np.ndarray]):
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
    assert len(helices) == np.ravel(fig.axes).shape[0]
    for ax, helix_set in zip(np.ravel(fig.axes), helices):
        ax = plot_helices(helix_set, False, ax, 50)
    return fig


def make_colorbar(fig: plt.Figure, heatmap_settings: HeatmapSettings) -> plt.Figure:
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


def plot_heatmap(ax: plt.Axes, data: np.ndarray, heatmap_settings: HeatmapSettings) -> plt.Axes:
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
    r_vals = heatmap_settings.grid.r_grid
    theta_vals = heatmap_settings.grid.theta_grid
    ax.pcolormesh(theta_vals, r_vals, data, cmap=heatmap_settings.colormap, norm=norm, zorder=0, edgecolors='face', linewidth=0)
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
        If True, turn the y axis into probability percentages, rather than
        raw counts. The default is False.
    show_target : boolean, optional
        If bulk_mode is True and show_target is True, label the red line "target
        mode". If bulk_mode is True and show_target is False, label the red line
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
        If True, also plot the error as a shaded region around the titration
        curve. If False, do not plot the error; error_type and n_replicas are
        ignored with this option. The default is True.
    error_type : str, optional
        Use 'std' for standard deviation. This is the default.
        Use 'ste' for standard error. User must specify n_replicas with this
        option.
        Use 'CI' for a 95% confidence interval. User must specify n_replicas
        with this option.
    n_replicas : int, optional
        The number of samples/replicas in your deltaG calculation. Used to
        calculate the standard error and/or confidence interval. The default is
        None.
    color : str, optional
        A matplotlib-recognizeable color string. Applied to curve. Default is
        None, which results in matplotlib choosing your color from its default
        palette.
    linestyle : str, optional
        A matplotlib-recognizable linestyle string. Applied to curve. Default is
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
    if isinstance(helices, list):
        helices = np.array(helices, dtype=np.float64)
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
        A list of 2d ndarrays containing enrichment values for each bin in the
        lattice. One list item per heatmap.
    helices : list of ndarrays
        Each ndarray in the list contains helix coordinates. There should be one
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
