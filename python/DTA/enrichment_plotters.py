import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
plt.rcParams['axes.grid'] = False 

from . import polarDensity_helper as pc

#    ax = plt.subplot(111,projection="polar")
def plot_single(toplot, 
                theta,
                radius,
                ax,
                vmid=1,
                vmin=0,
                vmax=2,
                norm=None,
                yticks=[],
                xticks=[],
                cmap='RdBu'):
    """
    Create a polar plot with a color mesh representing the density of a given data set.
    
    Args:
        toplot (array-like): The density data to be plotted.
        theta (array-like): The theta values for the polar plot.
        radius (array-like): The radius values for the polar plot.
        ax (matplotlib.axes.Axes): The polar subplot on which to plot the data.
        vmid (float, optional): The midpoint value for the color mapping. Defaults to 1.
        vmin (float, optional): The minimum value for the color mapping. Defaults to 0.
        vmax (float, optional): The maximum value for the color mapping. Defaults to 2.
        norm (matplotlib.colors.Normalize, optional): The normalization object for the color mapping. Defaults to None.
        yticks (list, optional): The y-axis tick labels. Defaults to an empty list.
        xticks (list, optional): The x-axis tick labels. Defaults to an empty list.
        cmap (str or colormap, optional): The colormap for the color mapping. Defaults to 'RdBu'.
    
    Returns:
        matplotlib.axes.Axes: The modified polar subplot with the density data plotted.
    """
    if norm is None:
        norm = pc.MidpointNormalize(midpoint=vmid,vmin=vmin,vmax=vmax)
    ax.grid(False)
    plt.axis('off')
    s = ax.pcolormesh(theta, 
                        radius, 
                        toplot,
                        cmap=cmap,
                        norm=norm,
                        zorder=0,
                        edgecolors='face',
                        linewidth=0,
                        )
    s.set_edgecolor('face')

    ax.set_xticklabels(xticks)
    ax.set_yticklabels(yticks)
    
    return ax

#  orange   light blue   green       purple      amber      blue       red 
#"#E69F00"  "#56B4E9"  "#009E73"   "#CC79A7"   "#F5C710"  "#0072B2"  "#D55E00"  
#sub = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F5C710", "#0072B2", "#D55E00"]

def plot_helices(helices, colorbychain, ax, markersize=3, sub=["tab:blue", "tab:cyan", "tab:green", "tab:purple", "tab:brown", "tab:olive"]):
    """
    Plot helices on a polar plot.

    Parameters:
    - helices (list or array): The helix data to be plotted. Each element in the list represents a set of helices, where each set is a list of angles and radii.
    - colorbychain (bool): A flag to determine if the helices should be colored by chain.
    - ax (matplotlib.axes.Axes): The polar subplot axis on which to plot the helices.
    - markersize (int, optional): The size of the scatter markers. Defaults to 3.
    - sub (list, optional): The list of colors to use for coloring the helices. Defaults to a predefined list of colors.

    Returns:
    - ax (matplotlib.axes.Axes): The modified polar subplot axis with the helices plotted.
    """

    if len(np.shape(helices))==1:
        helices = np.reshape(helices, (1,len(helices)))
    for i,pro in enumerate(helices[:]):
        if colorbychain:
            colors = sub[i]
        else:
            colors = sub[:len(pro[::2])]        
        ax.scatter(np.deg2rad(pro[1::2]),
                    pro[::2],
                    color=colors,
                    linewidth=None,
                    zorder=1, 
                    s=markersize,
                    )
    return ax
def polar_plot(data_in, 
               theta, 
               radius, 
               chains_groups, 
               helices_lwr=None, 
               helices_upr=None, 
               vmax=2, 
               vmid=1, 
               vmin=0, 
               colorbychain=True,
               figwidth=20,
               figheight=20,
               xticks=[],
               yticks=[],
               cmap=plt.cm.bwr_r):
    """
    Creates a polar plot with color mesh representing the density of a given data set. 
    It also plots helices on the polar plot.

    Parameters:
    - data_in (array/list): Array or list of density data.
    - theta (array): Array of theta values for the polar plot.
    - radius (array): Array of radius values for the polar plot.
    - chains_groups (list): List of lipids to plot.
    - helices_lwr (optional): Lower helices data.
    - helices_upr (optional): Upper helices data.
    - vmax (float): Maximum value for color mapping.
    - vmid (float): Midpoint value for color mapping.
    - vmin (float): Minimum value for color mapping.
    - colorbychain (bool): Flag to determine if helices should be colored by chain.
    - figwidth (int): Width of the figure.
    - figheight (int): Height of the figure.
    - xticks (list, optional): X-axis tick labels.
    - yticks (list, optional): Y-axis tick labels.
    - cmap (colormap, optional): Colormap for color mapping.

    Returns:
    - fig (figure): The created figure.
    - axes (axes): The axes of the subplots.
    """
    data_in = pc.sum_reps(data_in)
    fig = plt.figure(figsize=(figwidth,figheight))
    gs1=gridspec.GridSpec(len(chains_groups),2,wspace=.15, hspace=0.15)
    plt.rcParams.update({'font.size': 10})
    norm1 = pc.MidpointNormalize(midpoint=vmid,vmin=vmin,vmax=vmax)
    grid = 0    
    for cg in chains_groups:
        for leaf in data_in.columns:
            toplot = data_in.at[cg,leaf]
            ax = plt.subplot(gs1[grid],projection="polar")
            ax = plot_single(toplot, theta, radius, ax, norm=norm1, xticks=xticks, yticks=yticks, cmap=cmap)

            if leaf=="Outer":
                helices = helices_upr
            else:
                helices = helices_lwr

            ax = plot_helices(helices, colorbychain, ax, markersize=50)
            if grid%2==0:
                trans = mtransforms.ScaledTranslation(-40/72, -1.5, fig.dpi_scale_trans)
                ax.text(0.0, 1.0, cg, transform=ax.transAxes + trans,
                        fontsize='medium', va='bottom', fontfamily='serif')
            if grid < 2:
                ax.set_title(leaf)
            grid = grid + 1

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.21, 1, 0.5, 0.02])
    sm = plt.cm.ScalarMappable(cmap=cmap)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks(np.linspace(0,1,5))
    cbar.ax.set_xticklabels([vmin, (vmin+vmid)/2, vmid, (vmid+vmax)/2, vmax])

    return fig, fig.axes

#    ax = plt.subplot(111,projection="polar")
def plot_single(toplot, 
                theta,
                radius,
                ax,
                vmid=1,
                vmin=0,
                vmax=2,
                norm=None,
                yticks=[],
                xticks=[],
                cmap='RdBu'):
    """
    Create a polar plot with a color mesh representing the density of a given lipid in a simulation.

    Parameters:
    - toplot (array-like): The density data to be plotted.
    - theta (array-like): The theta values for the polar plot.
    - radius (array-like): The radius values for the polar plot.
    - ax (matplotlib.axes.Axes): The polar subplot on which to plot the data.
    - vmid (float, optional): The midpoint value for the color mapping. Defaults to 1.
    - vmin (float, optional): The minimum value for the color mapping. Defaults to 0.
    - vmax (float, optional): The maximum value for the color mapping. Defaults to 2.
    - norm (matplotlib.colors.Normalize, optional): The normalization object for the color mapping. Defaults to None.
    - yticks (list, optional): The y-axis tick labels. Defaults to an empty list.
    - xticks (list, optional): The x-axis tick labels. Defaults to an empty list.
    - cmap (str or colormap, optional): The colormap for the color mapping. Defaults to 'RdBu'.

    Returns:
    - matplotlib.axes.Axes: The modified polar subplot with the density data plotted.
    """
    if norm is None:
        norm = pc.MidpointNormalize(midpoint=vmid,vmin=vmin,vmax=vmax)
    ax.grid(False)
    plt.axis('off')
    s = ax.pcolormesh(theta, 
                        radius, 
                        toplot,
                        cmap=cmap,
                        norm=norm,
                        zorder=0,
                        edgecolors='face',
                        linewidth=0,
                        )
    s.set_edgecolor('face')

    ax.set_xticklabels(xticks)
    ax.set_yticklabels(yticks)
    
    return ax


def get_helices(the_path):
    """
    Get helix locations from the given path.

    Parameters:
    - the_path (str): The path to the directory containing the helix coordinate files.

    Returns:
    - helices_upr (numpy.ndarray or None): The upper helix coordinates, or None if the file is not found.
    - helices_lwr (numpy.ndarray or None): The lower helix coordinates, or None if the file is not found.
    """
    try:
        helices_lwr = np.loadtxt(the_path.joinpath("Protein_coords_lwr.dat"))
        helices_upr = np.loadtxt(the_path.joinpath("Protein_coords_upr.dat"))
    except FileNotFoundError:
        helices_lwr = None
        helices_upr = None
        print("Protein coordinates not found")

    return helices_upr, helices_lwr

def read_rep(file_list, chains_groups, leaflets=['low','upp'], enrich=True):
    """
    Reads a list of files and performs analysis on each file to generate counts and enrichments data.

    Parameters:
    - file_list (list): A list of file paths to be read and analyzed.
    - chains_groups (list): A list of chain groups.
    - leaflets (list, optional): A list of leaflets. Defaults to ['low', 'upp'].
    - enrich (bool, optional): A flag to determine if enrichment analysis should be performed. Defaults to True.

    Returns:
    - counts (pandas.DataFrame): A DataFrame containing the counts data.
    - enrichments (pandas.DataFrame): A DataFrame containing the enrichments data.
    """
    enrichments = pd.DataFrame(index=chains_groups, columns=leaflets)
    counts = pd.DataFrame(index=chains_groups, columns=leaflets)

    idx = 0
    for fl in file_list:
        filename = fl.name
        tmp_chain = filename.split('.')[0]
        tmp_nm = filename.split('.')[1]

        if idx == 0:
            rad, dr, dth, theta, radius, frames, Ntheta = pc.Coord_Get(fl)
            A = (radius * dr * dth)
            thetas = np.unique(theta)
        idx+=1
        toadd = np.loadtxt(fl, skiprows=1)
        toadd = toadd[:,3:-1]
        data, num_mol,avg_A,beads,exrho,avg_chain = pc.get_polar_data(fl)
        counts.at[tmp_chain,tmp_nm] = toadd
        enrichments.at[tmp_chain,tmp_nm] = pc._analysis_call_(fl, radius, dr, dth, frames, enrich=enrich)

    return counts, enrichments

def get_file_list(root, lipids):
    """
    Get a list of file paths for a given root directory and a list of lipids.

    Parameters:
    - root (str or pathlib.Path): The root directory to search for files.
    - lipids (list): A list of lipids to search for.

    Returns:
    - file_list (list): A list of file paths matching the given lipids in the root directory.
    """
    file_list = []
    for lip in lipids:
        toadd = list(root.glob(f"{lip}*avg.dat") )
        file_list = np.append(file_list,toadd)
    return file_list

