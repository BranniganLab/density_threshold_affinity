#!/usr/bin/env python
# coding: utf-8

# In[86]:


import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib import colormaps as cmaps
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterExponent
from DTA.enrichment_plotters import polar_plot, get_file_list, get_helices, read_rep
from DTA.polarDensity_helper import Coord_Get, get_header_info,_analysis_call_
from DTA.site_distributions import outline_site, make_symmetric_sites, combine_sites, plot_density,make_simple_site, Site, plot_bulk_counts
plt.rcParams['axes.grid'] = False 


# # Important User Settings:

# In[87]:


# root = Path("/home/ems363/Projects/ELIC_PE_PG/trajectory_version")
# chains_groups = ["POPE", "POPG"]
#root = Path("/home/ems363/Projects/ELIC_PE_PG/aggregated/")
root = Path("/home/liam/Documents/corringer/AA/5/ACA/")
#root = Path("/home/ems363/Projects/ELIC_PE_PG/liam/trajectory_version/")
chains_groups = ["ACA"]#,"POPG", "POPE"]
lipids = chains_groups


# In[88]:


enrich = True

# get files to use

file_list = []
for lip in lipids:
    #toadd = list(root.glob(f"{lip}.dat*avg.dat") )
    toadd = list(root.glob(f"{lip}*avg.dat") )
    file_list = np.append(file_list,toadd)

leaflets = ['low', 'upp']
print(file_list)


# In[77]:


enrichments = pd.DataFrame(index=chains_groups, columns=leaflets)
counts = pd.DataFrame(index=chains_groups, columns=leaflets)

idx = 0
for fl in file_list:
    if idx == 0:
        rad, dr, dth, theta, radius, frames, Ntheta = Coord_Get(fl)

    filename = fl.name

    tmp_chain = filename.split('.')[0]
    ## 
    tmp_nm = filename.split('.')[1]

    # This is a hack. The above part does not have a "flexible"
    # method to consider sim type (a, b ...)
    idx+=1
    toadd = np.loadtxt(fl, skiprows=1)
    toadd = toadd[:,3:-1]
    counts.at[tmp_chain,tmp_nm] = toadd
    enrichments.at[tmp_chain,tmp_nm] = _analysis_call_(fl, radius, dr, dth, frames, enrich=enrich)

print(enrichments)

# Optional helix locations
try:
    helices_lwr = np.loadtxt(root.joinpath("Protein_coords_lwr.dat"))
    helices_upr = np.loadtxt(root.joinpath("Protein_coords_upr.dat"))
except FileNotFoundError:
    helices_lwr = None
    helices_upr = None
    print("Protein coordinates not found")


# In[78]:


# theta is not defined
thetas = np.unique(theta.flatten())


# In[79]:


bins = {1:(range(4,30,6), np.repeat(27.5,5)), 
        2:(range(5,30,6), np.repeat(27.5,5)),
        3:(range(4,30,6), np.repeat(32.5,5)),
        4:(range(5,30,6), np.repeat(32.5,5)),
        5:(range(4,30,6), np.repeat(22.5,5))}


# In[80]:


from matplotlib.colors import ListedColormap
depleted = cmaps.get_cmap('RdGy_r')
enriched = cmaps.get_cmap('bwr_r')
newcolors = depleted(np.linspace(0, 1, 256))
newcolors[128:] = enriched(np.linspace(0.5,1,128))
my_cmap = ListedColormap(newcolors)


# In[81]:


from matplotlib.colors import ListedColormap
depleted = cmaps.get_cmap('PiYG')
enriched = cmaps.get_cmap('bwr_r')
newcolors = np.concatenate([depleted(np.linspace(0.4, 0.5, 128)), enriched(np.linspace(0.5,1,128))])
my_cmap = ListedColormap(newcolors)


# In[82]:


enrichments.at['ACA', 'upp'][0:1,:]=0
fig, axes = polar_plot(enrichments, 
                       theta, 
                       radius, 
                       chains_groups, 
                       helices_lwr, 
                       helices_upr, 
                       colorbychain=False, 
                       vmin=0.75, 
                       vmax=1.5, 
                       vmid=1,
                       figheight=8,
                       figwidth=8,
                       cmap=my_cmap)
fig.tight_layout()
plt.savefig(root.joinpath("ELIC_enrichments.pdf"))
plt.show()


# In[83]:

enrichments.at['ACA', 'upp']
enrichments.at['ACA', 'low']


# In[ ]:



