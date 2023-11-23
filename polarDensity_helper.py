#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 17:38:26 2019

@author: liam
"""

import numpy as np
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pdq
import os
from matplotlib.colorbar import Colorbar


class MidpointNormalize(mcol.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		mcol.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def spec_binding(ddG_in,leaf,ind_lip,bind_sight):
    
    ddG_out = np.nan_to_num(ddG_in[leaf][ind_lip] * bind_sight,nan=0.0)
    ddG_out [ddG_out > 10**(100)] = 0.0
    return ddG_out
   
def load_prob0():
    return pd.read_csv("../Data/polar/Probability_0_Synapse.csv",header=0, index_col=0)

# parses out names and replicas from file names
# a bit overboard
def get_col_names(files):
    col = []
    for fl in files:
        if len(fl) == 30:
            col.append(fl[18:26])
        elif len(fl) == 28:
            col.append(fl[16:24]) 
        elif len(fl) == 29:
            tmp_nm = fl[17:25]
            if tmp_nm.startswith("a") ==  False:
                col.append("a"+tmp_nm) 
        elif len(fl) > 30:
             col.append(fl[21:-4])
    return np.unique(col)

# confirms file exists
def file_check(files):
    flj = 0
    for fli, fl in enumerate(files):
        if os.path.isfile("../Data/polar/%s"%fl) == True:
            continue
        
        elif flj >= len(files) - 2:
            return False
        
        else:
            print("%s -- Does not exist"%fl)
            flj =+ 1
    return True

# carteasian to polar
def cart2pol(x, y):
    thta = np.arctan2(y, x)
    rho = np.hypot(x, y)
    return rho, thta


# selects lipids based on membrane
# no longer used
def Syn_ooc_comp():
    oo = ["POPE", "POPC", "PIPG", "PIPC", "PAPS", "PAPI", "OUPE", "IUPE", 
          "IEPE", "IAPC", "DXSM", "DPSM", "DOPS", "DLSM", "CHOL"]
    
    sn = ["CHOL","DOPC", "DPPC", "DPPS", "DPSM", "OAPE", "OIPC",
                     "OIPE", "OUPC", "OUPE", "OUPS", "PAP1", "PAP2", "PAP3", 
                     "PAPA", "PAPC", "PAPE", "PAPI", "PAPS", "PBSM", "PFPC", 
                     "PIPI", "PNSM", "POP1", "POP2", "POP3", "POPA", "POPC", 
                     "POPE", "POPI", "POPS", "POSM", "PUPC", "PUPE", "PUPI",
                     "PUPS"]
    #["PUPS", "PUPG", "PUPE", "PUPC", "POSM", "POPS", "POPI", "POPE",
         # "POPC", "PIPE", "PIPC", "PGPE", "PAPI", "PAPG", "PAPE", "PAPC", 
         # "OUPE", "OUPC", "ONPS", "OGPE", "OAPE", "DUPG", "DPPC", "DOPE", "CHOL"]
    sn_oo = set(sn).intersection(oo)
    return oo,sn,sn_oo

# calculates lipid ratios in a membrane based on a file of membrane composition
# Used by Liam from build files
def get_rat_val(satfl):
    rat_tmp = pd.read_csv(satfl, index_col=0).T
    rat_up = rat_tmp[["m1_a%i_up"%rep for rep in range(1,11)]].mean(axis=1)
    rat_low = rat_tmp[["m1_a%i_low"%rep for rep in range(1,11)]].mean(axis=1)
    rat = pd.concat([rat_up, rat_low], axis=1,keys=["Outer","Inner"])
    return rat

# takes log of density map
def logify(in_arr):
    in_arr = np.log(in_arr)
    in_arr[in_arr==-np.inf] = np.nan
    return in_arr    
 
# "Computs" center of "mass" for radial bins
def calc_rad(dat):
    # TODO Something is wrong here
    # starts at 0, didn't it start at 4?
    rad1,counts = np.unique(dat[:,0], return_counts=True)
    rad2 = np.unique(dat[:,1])
    rad = rad2+(rad2-rad1)/2.0
    counts = counts[0]
    return rad,counts
    
    
   
def Coord_Get(fl_in):
    '''
    Function to determine the radius and theta matrix
    plt_size is dependent on the size of the membrane plotted
    *** plt_size is no longer used remove ***
    
    dr, dth, theta, radius = Coord_Get(50)
    '''
    dat = np.loadtxt("../Data/polar/%s"%fl_in,skiprows=1)
    frames = 1
    #rad = dat[:,1]+(dat[:,1]-dat[:,0])/2.0    
    rad, frames = calc_rad(dat)
    dr = rad[1]-rad[0]
    Ntheta = 50
    the = np.linspace(0,2*np.pi,Ntheta +1)
    dth = the[1]-the[0]
    theta,radius=np.meshgrid(the,rad)

    return rad, dr, dth, theta, radius, frames

def get_header_info(sat_files):
    '''
    Parameters
    ----------
    sat_files : str
        path to a data file. extracts extra information
        used for analysis

    Returns
    -------
    list
        Number of molecules, number of beads per mol, and expected
        density.
    '''
    line = None
    if (sat_files.find("Neutral")!=-1 or sat_files.find("Anionic")!=-1):
        sat_files = sat_files.replace(".1.",".header.")
        line = open("../Data/polar/%s"%sat_files,'r').readline()
    else:
        line = open("../Data/polar/%s"%sat_files,'r').readline()
    num_mol = float(line.split(",")[0].split(":")[1].split(" ")[1])
    beads = float(line.split(",")[1].split(":")[1].split(" ")[1])
    avg_A = float(line.split(',')[2].split(":")[1].split(" ")[1])
    exrho = float(line.split(",")[3].split(":")[1].split("/")[0])
    avg_chain = float(line.split(",")[4].split(":")[1].split(" ")[1])
    return num_mol,avg_A,beads,exrho,avg_chain

def get_polar_data(sat_files):
    """
    _get_polar_data_("input")
    
    Function to get the data and informaiton in the header. The
    header information is used for density and for enrichment    
    """
    num_mol,avg_A,beads,exrho,avg_chain = get_header_info(sat_files)
    try:
        dat = np.loadtxt("../Data/polar/%s"%sat_files,skiprows=1)#/avg_chain
    except:
        return 0
  
    dat_set = dat[:,3:] 
    return dat_set,num_mol,avg_A,beads,exrho,avg_chain
    
def _analysis_call_(sat_file, rad, dr, dth, frames, enrich=False,):#cbs=None,w3bs=None, 
                    #ddg=False, rat=None):
    '''

    Parameters
    ----------
    sat_files : str
        Path to file.
    rad : np.array
        Radius array (1D).
    dr : float
        Change in radius between r1 and r2.
    dth : float
        Change in angle (radians).

    Returns
    -------
    rho : np.array 
        2D array of density data.
    '''
    data, num_mol,avg_A,beads,exrho,avg_chain = get_polar_data(sat_file)
    if (frames > 1):
        assert np.shape(data[::frames,:]) == np.shape(rad), "Radius matrix different shape form data"
    else:
        assert np.shape(data) == np.shape(rad), "Radius matrix different shape form data"
    
    # default... don't have a way at the moment
    # to skip this step....
    A = (rad * dr * dth)
    rho = data.copy()
    if (frames > 1):
        for f in range(frames):
            assert np.shape(data[f::frames]) == np.shape(A), "Data array wrong shape"
            # row 1 is 0 for data and A -> nan
            rho[f::frames] = data[f::frames] / A
    else:
        rho = data / A #(rad * dr * dth)
    if enrich == True:
        if exrho > 0:
            rho = rho / exrho
        else:
            rho = np.zeros_like(rho)#rho / 0.0000312740
        
    return rho

def sum_reps(rho_data):
	# averages data over replicas (#TODO Name change)
    upper = [labl for labl in rho_data.columns if (labl.endswith("upp") and not labl.startswith(".1.upp"))]
    lower = [labl for labl in rho_data.columns if (labl.endswith("low") and not labl.startswith(".1.low"))]
    rho_up = rho_data[upper].sum(axis=1)/len(upper)
    rho_lo = rho_data[lower].sum(axis=1)/len(lower)
    return pd.concat([rho_up,rho_lo],axis=1,keys=["Outer","Inner"])
 
# opens protein alpha helecies center of mass files
def prot_coord():
# Main plotting routine
    protein = np.loadtxt("../Data/Protein_coords_upr.dat")
    Au = protein[0,:]
    Bu = protein[1,:]
    Cu = protein[2,:]
    Du = protein[3,:]
    Eu = protein[4,:]
    chain_up = [Au,Bu,Cu,Du,Eu]
    
    protein = np.loadtxt("../Data/Protein_coords_lwr.dat")
    Al = protein[0,:]
    Bl = protein[1,:]
    Cl = protein[2,:]
    Dl = protein[3,:]
    El = protein[4,:]
    chain_lw = [Al,Bl,Cl,Dl,El]
    return chain_up, chain_lw

def cholesterol_binding_selection(protein,radian,theta):
    '''
    cholesterol_binding_selection 
    Loops across protein subunis, finds specified values for r and theta
    and returns a matrix of boolean values True for cholesterol binding
    regions, false otherwise.
    '''    
    protein.append(protein[0]) #assumes protien is a list
    rad = radian.copy()
    the = theta.copy()
    dtheta = theta[0,1] - theta[0,0]
    mod = .5
    out = []
    for sub_ind in range(len(protein)-1):
        flag = None
        p_theta = [np.deg2rad(protein[sub_ind+1][5]),np.deg2rad(protein[sub_ind][1])]
        if p_theta[0] > p_theta[1]:
            flag = True
            p_theta = [[0,p_theta[1]],[p_theta[0],2*np.pi]]
        # because alpha_beta is after 0, there was
        # an issue circling back to 0
        if p_theta[0] > p_theta[1]:
            p_theta = [p_theta[1],-p_theta[0]]
        if sub_ind == 0:
            the = (theta >= np.min(p_theta-dtheta*mod)) * (theta <= np.max(p_theta-dtheta*mod))
            rad = (radian >= 10) * (radian <= 32)
        elif sub_ind > 0:
            if flag == True:
                the = ((theta >= np.min(p_theta[0]-dtheta*mod)) * (theta <= np.max(p_theta[0]-dtheta*mod))) + ((theta >= np.min(p_theta[1] - dtheta*mod)) * (theta <= np.max(p_theta[1])))
                rad = (radian >= 10) * (radian <= 32)
            else:
                the = (theta >= np.min(p_theta-dtheta*mod)) * (theta <= np.max(p_theta-dtheta*mod))
                rad = (radian >= 10) * (radian <= 32)
        out.append(rad*the)
    return np.array(out)#rad*the

def omega3_binding_selection(protein,radian,theta):
    '''
    omega3_binding_selection 
    
    Loops across protein subunis, finds specified values for r and theta
    and returns a matrix of boolean values True for w3 binding
    regions, false otherwise.
    '''    
    rad = radian.copy()
    the = theta.copy()
    dtheta = theta[0,1] - theta[0,0]
    mod = .5

    out = []
    for sub_ind in range(len(protein)):
        flag = None
        p_theta = [np.deg2rad(protein[sub_ind][1]),np.deg2rad(protein[sub_ind][5])]
        if p_theta[0] > p_theta[1]:
            flag = True
            p_theta = [[0,p_theta[1]],[p_theta[0],2*np.pi]]
        if sub_ind == 0:
            the = (theta >= np.min(p_theta-dtheta*mod)) * (theta <= np.max(p_theta-dtheta*mod))
            rad = (radian >= 10) * (radian <= 44)
        elif sub_ind > 0:
            if flag == True:
                the = ((theta >= np.min(p_theta[0]-dtheta*mod)) * (theta <= np.max(p_theta[0]-dtheta*mod))) + ((theta >= np.min(p_theta[1]-dtheta*mod)) * (theta <= np.max(p_theta[1]-dtheta*mod)))
                rad = (radian >= 10) * (radian <= 44)
            else:
                the = (theta >= np.min(p_theta-dtheta*mod)) * (theta <= np.max(p_theta-dtheta*mod))
                rad = (radian >= 10) * (radian <= 44)
        out.append(rad*the)
    return np.array(out)

def _Polar_Plot_(data_in, theta, radius, chains_groups,memb):
	# plots densities
	# data_in = array/list of density data
	# theta, radius = arrays of position bins
	# chains_groups = old name, really lipids to plot
	# memb = redundant
	
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import LogFormatterExponent
    data_in = sum_reps(data_in)
    fig = plt.figure(figsize=(5,22.5))
    gs1=gridspec.GridSpec(len(chains_groups)+2,2,wspace=.15, hspace=0.15)
    max_norm = 0.003#
    plt.rcParams.update({'font.size': 10})
    norm1 = MidpointNormalize(midpoint=max_norm/2,vmin=0,vmax=max_norm)
    norm1 = MidpointNormalize(midpoint=1,vmin=0,vmax=2)
    cmap = plt.cm.RdBu_r#PuOr
    cmap.set_bad(color='black')
    grid = 0
    #chains_up, chains_lo = prot_coord()
    sub = ["g",'m','grey','green','cyan']

    for cg in chains_groups:
        
        ax = plt.subplot(gs1[grid],projection="polar")
        try:
            s = ax.pcolormesh(theta, radius, ((data_in.at[cg,"Outer"]  )),
                              cmap=cmap,norm=norm1,zorder=0)
        except:
            grid = grid + 1
            continue
        #for i,pro in enumerate(chains_up[:5]):
        #    ax.scatter(np.deg2rad(pro[1::2]),pro[::2],edgecolor=sub[i],
        #               facecolors=sub[i],linewidth=1,zorder=2,s=np.shape(data_in)[0]*10)
        
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        grid = grid + 1
        
        ax = plt.subplot(gs1[grid],projection="polar")
        try:
            ax.pcolormesh(theta, radius, (data_in.at[cg,"Inner"]) ,
                              cmap=cmap,norm=norm1,zorder=0)
        except:
            grid = grid + 1
            continue
        for i,pro in enumerate(chains_lo[:5]):
            ax.scatter(np.deg2rad(pro[1::2]),pro[::2],edgecolor=sub[i],
                       facecolors=sub[i],linewidth=1,zorder=2,s=np.shape(data_in)[0]*10)
        
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        grid = grid + 1

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.21, .89, 0.5, 0.008])
    cb=fig.colorbar(s, cax=cbar_ax,ticks=[0,0.5,1,1.5,2],orientation="horizontal")

    plt.tight_layout()
    plt.show()

