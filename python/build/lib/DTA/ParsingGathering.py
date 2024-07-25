#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 22:17:41 2020

@author: liam
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from . import polarRadial_Complex_helper as pc

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


def parse_expected(leaf,A_acc,lip):
  
  
    # Runs for the non-ideal system 
    # memb is an old var and should be deleted
    list_expect = []
    lipids=None
    frames=None
    # Files contain all lipids, name is the area size
    # should first make sure the replicas are averaged together
    #A_acc = areas
    for ai, a in enumerate(A_acc):#range(1,11):
        # At current, exp_1... is the previous area data 
        # exp_2 is the current area data
        tmp_exp = pd.read_csv("../Data/polar/exp_3_synpR_%s_area_expected.%s.dat"%(str(a),leaf),sep=' ', index_col=0, header=0)
        #tmp_exp = pd.read_csv("../Data/polar/exp_1_nach_de4.2_%s_area_expected.%s.dat"%(str(a),leaf),sep=' ', index_col=0, header=None)
        if ai == 0:
            #lipids = tmp_exp.index
            frames = tmp_exp.columns
        # TODO BUG IS HERE. ONLY OUTPUTS THE FIRST ROW ## Don't think this is true... prety sure I got this.
        # Need a data frame for each area
        list_expect.append(tmp_exp.T[lip])
    #exp_val = pd.DataFrame(np.mean(list_expect,axis=0), index=lipids,columns=frames)
    exp_val = pd.DataFrame(list_expect, index=A_acc, columns=frames)
    return exp_val

def set_expected(memb,lip,areas):
    exp_upp = parse_expected("upp",areas,lip)
    exp_low = parse_expected("low",areas,lip)
    #TODO terat this as analyzing series THEN make a dataframe
    
    upp_hist,upp_edg = [],[]#pd.DataFrame(index=exp_upp.index,columns=np.arange(-1,54,1)),pd.DataFrame(index=exp_upp.index,columns=np.arange(-1,55,1))
    low_hist,low_edg = [],[]#pd.DataFrame(index=exp_low.index,columns=np.arange(-1,54,1)),pd.DataFrame(index=exp_low.index,columns=np.arange(-1,55,1))
    
    for key_id, key in enumerate(exp_upp.index):
        low_hist.append(expected_hist(exp_low.iloc[key_id])) #np.histogram(exp_low.iloc[key_id], np.arange(-1,55,1))[0])
        upp_hist.append(expected_hist(exp_upp.iloc[key_id])) #np.histogram(exp_upp.iloc[key_id], np.arange(-1,55,1))[0])
    # no 1/10 because we already takethe mean
    
    return pd.DataFrame(low_hist,index=exp_low.index), pd.DataFrame(upp_hist,
                        index=exp_upp.index), len(exp_upp.columns)    

def parse_expected_A(memb,leaf,areas):
    #mods to run DPPC only membrane with multip areas
    list_expect = []
    lipids=None
    frames=None
    #A_exp = [645,516,387,258,129]
    A_exp = areas
    for ai, a in enumerate(A_exp):#range(1,11):
        tmp_exp = pd.read_csv("../Data/polar/exp_3_synpR_a%i_area_expected.%s.dat"%(a,leaf),sep=' ', index_col=0, header=None)
        #tmp_exp = pd.read_csv("../Data/polar/exp_1_nach_de4.2_%s_area_expected.%s.dat"%(str(a),leaf),sep=' ', index_col=0, header=None)
        if ai == 0:
            lipids = tmp_exp.index
            frames = tmp_exp.columns
        list_expect.append(tmp_exp.values[0])
    #exp_val = pd.DataFrame(np.mean(list_expect,axis=0), index=lipids,columns=frames)
    exp_val = pd.DataFrame(list_expect, index=A_exp, columns=frames)
    return exp_val

# parse_expected_A("", "upp")

# moded for single lipids
def expected_hist(inval):
    return np.histogram(inval, np.arange(0,200,1))[0]
    
def calc_expected_prob(inval, norm):
    return inval/norm

def set_expected_prob(histUp, histLow,normUp, normLow):
    
    upp_prob = calc_expected_prob(np.asarray(histUp), normUp)
    low_prob = calc_expected_prob(np.asarray(histLow), normLow)
    
    return pd.DataFrame(low_prob,index=histLow.index), pd.DataFrame(upp_prob,
                        index=histUp.index)    

def get_area_accessable(Area, counts, ext_up, ext_low, ar_val):
    for j,i in enumerate(["alpha-beta","beta-delta",
        "delta-alpha","alpha-gamma","gamma-alpha","alpha_b","beta",
        "delta","alpha_g","gamma", ]):     
        Aup = [n for n in ar_val if str(n)[:3]==ext_up[j].name][0]
        Alo = [n for n in ar_val if str(n)[:3]==ext_low[j].name][0]
        Area["Outer_%s"%i]["A"] = Aup
        Area["Inner_%s"%i]["A"] = Alo
        Area["Outer_%s"%i]["acc"] = np.argmax(counts["Outer_%s"%i]) * Aup / np.argmax(ext_up[j])
        Area["Inner_%s"%i]["acc"] = np.argmax(counts["Inner_%s"%i]) * Alo / np.argmax(ext_low[j])
        Area["Outer_%s"%i]["acc_root"] = np.sqrt(np.argmax(counts["Outer_%s"%i]) * Aup / np.argmax(ext_up[j]))
        Area["Inner_%s"%i]["acc_root"] = np.sqrt(np.argmax(counts["Inner_%s"%i]) * Alo / np.argmax(ext_low[j]))
        Area["Outer_%s"%i]["res"] = Area["Outer_%s"%i]["A"] - Area["Outer_%s"%i]["acc"]   
        Area["Inner_%s"%i]["res"] = Area["Inner_%s"%i]["A"] - Area["Inner_%s"%i]["acc"]     
    return Area

def set_expected_A(memb,areas):
    exp_upp = parse_expected_A(memb,"upp",areas)
    exp_low = parse_expected_A(memb,"low",areas)
    #TODO terat this as analyzing series THEN make a dataframe
    
    upp_hist,upp_edg = [],[]#pd.DataFrame(index=exp_upp.index,columns=np.arange(-1,54,1)),pd.DataFrame(index=exp_upp.index,columns=np.arange(-1,55,1))
    low_hist,low_edg = [],[]#pd.DataFrame(index=exp_low.index,columns=np.arange(-1,54,1)),pd.DataFrame(index=exp_low.index,columns=np.arange(-1,55,1))
    
    for key_id, key in enumerate(exp_upp.index):
        low_hist.append(expected_hist(exp_low.iloc[key_id])) #np.histogram(exp_low.iloc[key_id], np.arange(-1,55,1))[0])
        upp_hist.append(expected_hist(exp_upp.iloc[key_id])) #np.histogram(exp_upp.iloc[key_id], np.arange(-1,55,1))[0])
    # no 1/10 because we already takethe mean
    
    return pd.DataFrame(low_hist,index=exp_low.index), pd.DataFrame(upp_hist,
                        index=exp_upp.index), len(exp_upp.columns) 

def calc_dG(prob_0, prob_n):
    dP = np.array(prob_n.values.tolist()) / np.array(prob_0.values.tolist())
    RT = 0.642748322147651006711409395973154362416107382550335570469
    dG = pd.DataFrame(-RT*np.log(dP)
        ,index=prob_0.index, columns=prob_0.columns)
    return dG

def calc_dG_exp(prob_0, prob_n):
    dP = np.array(prob_n.values.tolist()) / np.array(prob_0.values.tolist())
    RT = 0.642748322147651006711409395973154362416107382550335570469
    dG = pd.DataFrame(-RT*np.log(dP)
        ,index=prob_0.index,)
    return dG

def plot_dG(var_in, flnm, vmin=-4.5,vmax=2.5):

    cmap = plt.cm.RdBu#PuOr
    cmap.set_bad(color='black')
    norm = MidpointNormalize(midpoint=0,vmin=-4.5,vmax=2.5)
    c = plt.pcolormesh(var_in,norm=norm,cmap=cmap)
    
    plt.yticks(np.linspace(0.5,35.5,36),labels=var_in.index.tolist(),fontsize=7)
    plt.xticks(np.linspace(0.5,len(var_in.columns)-1.5,len(var_in.columns)-1),
               labels=var_in.columns.tolist(),rotation=90, fontsize=7)
    plt.colorbar(c)
    plt.tight_layout()
    plt.show()
    #plt.savefig("%s.pdf"%flnm)
    plt.close()