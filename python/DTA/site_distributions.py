from dataclasses import dataclass
import numpy as np
from string import ascii_uppercase
@dataclass
class Site:
    """
    Represents a binding site with various properties such as inner and outer radii, number of theta bins, and title.
    Calculates and stores additional attributes such as densities, Npeak, and mean.

    Attributes:
        inner_r (float): The inner radius of the binding site.
        outer_r (float): The outer radius of the binding site.
        thetabins (int): The number of theta bins.
        title (str): The title of the binding site.
        counts (list): A list of counts.
        area (float): The calculated area of the site.
        densities (list): A list of densities.
        Npeak (float): The calculated Npeak value.
        mean (float): The calculated mean value.
    """
    inner_r: float=0
    outer_r: float=0
    thetabins: int=0
    title: str=""
    counts: list=None
    area: float=0
    densities: list=None
    Npeak: float=0
    peak: int=0
    mean: float=0
    expPunocc: float=0
    theta_vals: list=None
    dr: float=0

def make_simple_site(the_data, inner_r=0, outer_r=0, nth=1, Ntheta=1, dr=1, dth=1, exrho=0, frames=1, the_thetas=None, title="", Npeak=None, accessible_area=None):
    """
    Create a Site object with specified inner and outer radii, number of theta bins, and title.
    Calculate the shell data based on the given parameters using the get_shell function.
    Flatten the shell data and store it in the counts attribute of the Site object.
    Calculate the area of the site using the get_area function and store it in the area attribute.
    Calculate additional attributes such as densities, Npeak, and mean using the get_site_stats function.
    
    Args:
        the_data (DataFrame): The data (bin counts) for the whole system.
        inner_r (float): The inner radius of the binding site.
        outer_r (float): The outer radius of the binding site.
        nth (int): The number of theta bins.
        Ntheta (int): The number of theta bins in the shell data.
        dr (float): The width of each shell.
        dth (float): The width of each theta bin.
        exrho (float): A scaling factor.
        frames (int): The number of frames in the trajectory.
        the_thetas (list): A list of indices representing the desired theta bins.
        title (str, optional): The title of the binding site.
        
    Returns:
        Site: A Site object with the calculated attributes: counts, area, densities, Npeak, and mean.
    """
    the_site = Site(inner_r, outer_r, nth, title)

    if len(the_data)>0:
        tmp_shell = get_shell(the_data, the_site, frames, Ntheta, dr)
        the_site.counts = tmp_shell[:,the_thetas]
        if nth > 1:
            the_site.counts = np.sum(the_site.counts, axis=-1)
    else:
        the_site.counts = np.ndarray([])

    if accessible_area is not None:
        the_site.area = accessible_area
    else:
        the_site.area = get_area(the_site, dth)
    the_site = get_site_stats(the_site)

    the_site.theta_vals = (the_thetas*dth)%(2*np.pi)
    the_site.dr = dr
    if Npeak is not None:
        the_site.Npeak = Npeak
    else:
        the_site.Npeak = the_site.peak
    return the_site

def combine_sites(list_of_sites,
                  exrho,
                  newtitle="composite site",
                  custom_area=None,
                  symmetric=False):
    """
    Combines the properties of the input sites to create a new composite site.

    Args:
        list_of_sites (list): A list of `Site` objects representing the individual sites to be combined.
        exrho (float): The expected density of the bin from the bulk continuum approximation.
        newtitle (str, optional): The title to be assigned to the new composite site. Default is "composite site".
        custom_area (float, optional): A custom area definition
        symmetric (bool, optional): Whether or not the combined sites are symmetric copies of each other. Default False.
    Returns:
        Site: A new `Site` object representing the composite site with combined properties and updated attributes.
    """
    new_site = Site()
    new_site.counts = np.zeros_like(list_of_sites[0].counts)
    for site in list_of_sites:
        new_site.inner_r = np.min([new_site.inner_r, site.inner_r])
        new_site.outer_r = np.max([new_site.outer_r, site.outer_r])
        if symmetric:
            new_site.counts = np.concatenate([new_site.counts, site.counts])
        else:
            new_site.counts = new_site.counts + site.counts
        if not symmetric and custom_area is None:
            new_site.area = new_site.area + site.area

    if symmetric:
        new_site.area = list_of_sites[0].area
    if custom_area is not None:
        print("Warning: custom area overrides other area calculations and assignments")
        new_site.area = custom_area

    new_site = get_site_stats(new_site)
    new_site.title = newtitle

    return new_site

def get_shell(some_data, site:Site, frames, Ntheta, dr):
    """
    Returns the flattened shell data based on the given parameters.

    Parameters:
    some_data (numpy.ndarray): An array containing the data (either counts or entrichments).
    site (Site): A Site object representing the inner and outer radii of the shell.
    frames (int): The number of frames in the data.
    Ntheta (int): The number of theta bins.
    dr (float): The width of each shell.

    Returns:
    numpy.ndarray: The flattened shell data with shape (frames, Ntheta).
    """

    mask = np.logical_and(some_data[:,0] >= site.inner_r, some_data[:,1] <= site.outer_r)
    meta_shell = some_data[mask, 3:-1]
    nshells = (site.outer_r-site.inner_r)/dr
    assert nshells == int(nshells), "Error: non-integer number of shells"
    reshaped = np.reshape(meta_shell, (int(nshells),frames, Ntheta))
    flattened = np.sum(reshaped, axis=0)
    return flattened

def plot_density(site: Site, ax):
    """
    Plot the probability density of a site object.

    Args:
        site (Site): A Site object that contains information about the site.
        ax: An ax object from the matplotlib.pyplot library that represents the plot axes.

    Returns:
        ax: The modified ax object with the plot.
    """
    dG = get_dg(site)

    ax.vlines(site.Npeak, 0, np.max(site.densities), color='red', linestyles='solid', label=r"$N_\mathrm{peak}$")
    #ax.vlines(site.mean, 0, np.max(site.densities), color='red', linestyles='solid', label="actual mean")
    ax.plot(np.arange(len(site.densities)), site.densities, label=r"$\Delta G =$"+f"{np.round(dG,2)}kcal/mol")

    ax.legend()
    ax.set_xlabel(r"$N_\mathrm{beads}$")
    ax.set_ylabel("Probability Density")
    ax.set_title(site.title)

    return ax

def get_site_stats(site):
    """
    Calculate the statistics of a given site object and return an updated new_site object with additional attributes.

    Args:
        site (Site): A Site object representing a site with certain properties.
        exrho (float): A scaling factor.

    Returns:
        Site: An updated Site object with additional attributes densities, Npeak, and mean.
    """
    new_site = site
    frequencies = np.bincount(site.counts.astype(int).flatten())
    new_site.densities = frequencies/np.sum(frequencies)
    
    new_site.peak = np.argmax(new_site.densities)
    new_site.mean = np.mean(site.counts)
    return new_site

def get_area(site: Site, dth: float) -> float:
    """
    Calculate the area of a site object based on its inner and outer radii, the number of theta bins, and the width of each bin.

    Args:
        site (Site): A Site object representing the inner and outer radii of the shell.
        dth (float): The width of each theta bin.

    Returns:
        float: The calculated area of the site.
    """
    return np.mean([site.inner_r, site.outer_r]) * (site.outer_r - site.inner_r) * dth * site.thetabins

def get_dg(site, RT=0.6427483):
    """
    Calculates the difference in free energy (dG) between two probability distributions based on the densities of a given site object.

    Args:
        site (Site): A Site object representing a binding site
        RT (float, optional): The gas constant multiplied by the temperature. Defaults to 0.6427483, which gives kcal/mol at 323K, default for MARTINI 

    Returns:
        float: The difference in free energy between the two probability distributions, p_< and p_>
    """
    if site.Npeak == 0:
        dG = -RT * np.log(site.expPunocc/site.densities[0])

    else:
        beads = np.arange(len(site.densities))
        mask_lt = beads<=site.Npeak
        mask_gt = beads>site.Npeak
        p_lessthan = np.sum(site.densities[mask_lt])
        p_greaterthan = np.sum(site.densities[mask_gt])
        dG = -RT * np.log(p_greaterthan / p_lessthan)
    return dG


def plot_site(ax, thetas, mintheta, maxtheta, rmin, rmax):
    """
    Plot a filled region on a polar plot, representing a binding site within a specified range of angles and radii.

    Parameters:
    ax (matplotlib.axes.Axes): The matplotlib axes object representing the polar plot.
    thetas (numpy.ndarray): An array of angles defining the polar plot.
    mintheta (int): The index of the starting angle for the binding site.
    maxtheta (int): The index of the ending angle for the binding site.
    rmin (float): The minimum radius of the binding site.
    rmax (float): The maximum radius of the binding site.

    Returns:
    matplotlib.axes.Axes: The modified matplotlib axes object with the filled region representing the binding site.
    """
    minangle = thetas[mintheta]+thetas[1]/2
    maxangle = thetas[maxtheta%len(thetas)]+thetas[1]/2
    if maxangle<minangle:
        maxangle = minangle+maxangle+thetas[1]*1.5
    angles = np.linspace(minangle, maxangle,50)
    ax.fill_between(angles, rmin, rmax, edgecolor='black', facecolor='none')
    return ax

def make_symmetric_sites(the_data, theta_start, width, inner_r, outer_r, Ntheta, dr, dth, exrho, frames, binspersubunit=10, subunits=5, sitename="site"):
    """
    Create a list of symmetric binding sites based on the given parameters.

    Args:
        the_data (array-like): The data (bin counts) for the whole system.
        theta_start (int): The starting theta bin index.
        width (int): The width of each theta bin.
        rmin (float): The minimum radius of the binding site.
        rmax (float): The maximum radius of the binding site.
        binspersubunit (int, optional): The number of bins per subunit. Defaults to 10.
        subunits (int, optional): The number of subunits. Defaults to 5.
        sitename (str, optional): The base name for the binding sites. Defaults to "site".

    Returns:
        list: A list of symmetric binding sites.
    """
    sites = []
    n_radial_bins = binspersubunit*subunits
    ts = np.arange(theta_start,n_radial_bins,binspersubunit)
    units = ascii_uppercase[:subunits]
    for unit, tmin in zip(units, ts):
        the_thetas = np.arange(tmin,tmin+width)%n_radial_bins
        toadd = make_simple_site(the_data, inner_r, outer_r, width, Ntheta, dr, dth, exrho, frames, the_thetas, sitename+unit)
        sites.append(toadd)
    return sites

def outline_site(ax, site):
    """
    Outline a binding site on a plot with a black border.

    Parameters:
    - ax (matplotlib.axes.Axes): The matplotlib axes object representing the plot.
    - site (Site): A Site object representing the binding site.

    Returns:
    - matplotlib.axes.Axes: The modified matplotlib axes object with the binding site outlined.
    """

    sorted_thetas = np.sort(site.theta_vals)
    dtheta = sorted_thetas[1]-sorted_thetas[0]
    the_thetas = site.theta_vals-dtheta/2
    the_thetas = np.append(the_thetas, the_thetas[-1]+dtheta)
    ax.fill_between(the_thetas, site.inner_r, site.outer_r, facecolor=(0,0,0,0), edgecolor='k')
    return ax
