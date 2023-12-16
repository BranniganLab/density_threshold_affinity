#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:24:11 2020

@author: liam
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from . import ParsingGathering as pg
from . import polarRadial_Complex_helper as pc
#import Calc_Error_dG as ceg


######################## ERROR CALC ###########################################
def get_prob_error(chains_groups, rep):
        
    inter_sites = get_intersites()
    counts = None
    for cgid, cg in enumerate(chains_groups):
        #prob_for_sem = None
        for a in range(1,rep+1):
            #files = ["polar1_nach_%s_de%i.4.upp.dat"%(cg,4),"polar1_nach_%s_de%i.4.low.dat"%(cg,4)]
            files = [ "polar1_synpR_%s_a%i.4.upp.dat"%(cg,a),"polar1_synpR_%s_a%i.4.low.dat"%(cg,a),"polar1_synpR_%s_a%i.4.upp.dat"%("SM",a),"polar1_synpR_%s_a%i.4.low.dat"%("SM",a)  ]
            den_tmp = pd.DataFrame(index=chains_groups, columns=["Outer","Inner"])
            rad, dr, dth, theta, radius, frames = pc.Coord_Get(files[0])
            if (a == 1):
                counts = pd.DataFrame(index=(np.linspace(0, 199,200)), columns=(inter_sites)).fillna(0.0)
            
            chains_up = pc.prot_coord()
            bindsites = pc.cholesterol_binding_selection(chains_up[0],radius,theta), pc.cholesterol_binding_selection(chains_up[1],radius,theta),pc.omega3_binding_selection(chains_up[0],radius,theta),pc.omega3_binding_selection(chains_up[1],radius,theta)
            print(cg,a)
            den_tmp["Outer"][cg] = pc._analysis_call_(files[0], radius, dr, dth, frames, enrich=False) # not sure how the output will look
            den_tmp["Outer"][cg] = den_tmp["Outer"][cg] + pc._analysis_call_(files[2], radius, dr, dth, frames, enrich=False)
            den_tmp["Inner"][cg] = pc._analysis_call_(files[1], radius, dr, dth, frames, enrich=False) 
            den_tmp["Inner"][cg] = den_tmp["Inner"][cg] + pc._analysis_call_(files[3], radius, dr, dth, frames, enrich=False) # runn that line!
            
            #This is only going to work short term....
            den_tmp = densit_frames(den_tmp, cg, frames)
            n_avg = nAvg(den_tmp, radius, dr, dth, frames, bindsites)
            counts_temp,edge = histogram_sites(n_avg)
            prob_for_sem = counts_temp / len(n_avg.iloc[:-1,0])
            for s in prob_for_sem:
                prob_for_sem[s].to_csv("../Data/deltaG/stats/%s_%s_m1_a%i.csv"%(cg,s,a))               

def calc_error(chains_groups,rep):
    inter_sites = get_intersites()
    outError = pd.DataFrame(index=chains_groups, columns=inter_sites).fillna(0.0)
    for cgid, cg in enumerate(chains_groups):
        semDist = pd.DataFrame(index=(np.linspace(0, 199,200)), columns=inter_sites).fillna(0.0)
        for s in inter_sites:
            standError = pd.DataFrame(index=(np.linspace(0, 199,200)), columns=["a%i"%a for a in range(1,rep+1)]).fillna(0.0)
            for a in range(1,rep+1):
                standError["a%i"%a] = pd.read_csv("../Data/deltaG/stats/%s_%s_m1_a%i.csv"%(cg,s,a),index_col=0,header=0)
            semDist[s] = standError.sem(axis=1)
        semDist.to_csv("../Data/deltaG/stats/sem_%s_m1_a%i.csv"%(cg,a))
        outError.T[cg] = semDist.mean()
    outError.to_csv("../Data/deltaG/stats/mean_sem.csv")

def get_area_error(areas):
    for a,s in zip(areas,get_intersites()):
        for leaf in ["upp","low"]:
            for rep in range(1,11):
                exp_bulk = pd.read_csv("../Data/polar/exp_1_synpR_a%i.2_%s_area_expected.%s.dat"%(rep,a,leaf),
                                      sep=" ",index_col=0,header=None)
                for exb in exp_bulk.T:
                    (pd.Series(np.histogram(exp_bulk.T[exb],bins=200,range=(-1,
                        199))[0])/len(exp_bulk.columns)).to_csv("../Data/deltaG/stats/%s_%s_m1_a%i_exp.csv"%(exb,s,rep))   
                    
def calc_area_error(chains_groups,rep):
    inter_sites = get_intersites()
    outError = pd.DataFrame(index=chains_groups, columns=inter_sites).fillna(0.0)
    for cgid, cg in enumerate(chains_groups):
        semDist = pd.DataFrame(index=(np.linspace(0, 199,200)), columns=inter_sites).fillna(0.0)
        for s in inter_sites:
            standError = pd.DataFrame(index=(np.linspace(0, 199,200)), columns=["a%i"%a for a in range(1,rep+1)]).fillna(0.0)
            for a in range(1,rep+1):
                standError["a%i"%a] = pd.read_csv("../Data/deltaG/stats/%s_%s_m1_a%i_exp.csv"%(cg,s,a),index_col=0,header=0)
            semDist[s] = standError.sem(axis=1)
        semDist.to_csv("../Data/deltaG/stats/sem_%s_m1_a%i_exp.csv"%(cg,a))
        outError.T[cg] = semDist.mean()
    outError.to_csv("../Data/deltaG/stats/mean_sem_exp.csv")
    
def get_error_ratio(lipids,rep):
    intersites = get_intersites()
    RT = 0.642748322147651006711409395973154362416107382550335570469
    # PL_err = pd.DataFrame(index=lipids, columns=intersites)
    # PG_err = pd.DataFrame(index=lipids, columns=intersites)
    dG_err = pd.DataFrame(index=lipids, columns=intersites)
    dG_avg = pd.DataFrame(index=lipids, columns=intersites)
    dG = pd.DataFrame(index=lipids, columns=intersites)
    for cg in lipids:
        for s in intersites:
            temp_PL = []
            temp_PG = []
            for a in range(1,rep+1):
                prob = pd.read_csv("../Data/deltaG/stats/%s_%s_m1_a%i.csv"%(cg,s,a),index_col=0,header=0)
                exptP = pd.read_csv("../Data/deltaG/stats/%s_%s_m1_a%i_exp.csv"%(cg,s,a),index_col=0,header=0)
                #ex_max = get_std(exptP)
                ex_max = exptP.values.argmax()      
                temp_PL.append(prob.values[:ex_max].sum())
                temp_PG.append(prob.values[ex_max+1:].sum())
            # PL_err[s][cg] = pd.Series(temp_PL).sem()
            # PG_err[s][cg] = pd.Series(temp_PG).sem()
            dP_err = pd.Series(np.asarray(temp_PG)/np.asarray(temp_PL))
            dG[s][cg] = (-RT * np.log(pd.Series(np.asarray(temp_PG)/np.asarray(temp_PL))))
            dG_avg[s][cg] = (-RT * np.log(pd.Series(np.asarray(temp_PG)/np.asarray(temp_PL)))).mean()
            dG_err[s][cg] = (-RT * np.log(pd.Series(np.asarray(temp_PG)/np.asarray(temp_PL))).replace(np.inf,0).replace(-np.inf,0)).sem()
    dG_err.to_csv("../Data/deltaG/stats/SEM_dG.csv")
    print(dG_err)
    #return dG, dG_avg, dG_err#PG_err, PL_err

def get_dG_error(PG_err, PL_err):
    RT = 0.642748322147651006711409395973154362416107382550335570469
    # TODO PG/PL has issues with 0/0?!?!?!
    dP_err = (PG_err / PL_err)
    #dP_err[dP_err == 0] = np.nan
    dG_calc_err = -1.0 * RT * np.log(dP_err.replace(0,np.nan))
    dG_err = pd.DataFrame(dG_calc_err,index=PG_err.index,columns=PG_err.columns)
    return
            #     for icol,col in enumerate(prob.columns):
            #         PL[col][lip] = prob[col][:ex_max[icol]].sum()
            #         PG[col][lip] = prob[col][ex_max[icol]+1:].sum()
            # dP = calc_prob_ratio(PL,PG)
    # return PL, PG, dP    
###############################################################################

############################ SUBROUTINES FOR DISTRIBUTIONS ####################
# TODO the averages don't look like they match the probability
def densit_frames(density_df, cg, frames): 
    den_frms = pd.DataFrame(index=[f for f in range(frames)], columns=["Outer","Inner"])
    for f in range(frames):
        den_frms["Outer"][f] = density_df["Outer"][cg][f::frames,:]
        den_frms["Inner"][f] = density_df["Inner"][cg][f::frames,:]  
    return den_frms

def nAvg_bin(density_df, bind_site, A):
    density_df = density_df  * bind_site * A 
    # Lower data seems to be added to each map...
    return  density_df 

def iterate_sites(n_avg_frm, density_frm, bind_site, A):
    bind_ind = 0
    bind_count = 0
    for  key in n_avg_frm.index:
        k = None
        if key.startswith("Outer")==True:
            k = ["Outer",0]
            if bind_count == 1:
                k = ["Outer",2]
        else:
            k = ["Inner",1]
            if bind_count == 3:
                k = ["Inner",3]
                
        n_avg_frm[key] = nAvg_bin(density_frm[k[0]], bind_site[k[1]][bind_ind], A).sum() 
        bind_ind = bind_ind + 1
        # sloppy
        if (bind_ind >= len(n_avg_frm.index)/4): 
            bind_ind = 0
            bind_count = bind_count + 1
    
    return n_avg_frm

def nAvg(density_df, radius, dr, dth, frames, bind_site):
    n_avg = pd.DataFrame(index=[f for f in range(frames)], columns=["%s_%s"%(i,j) 
        for i in ["Outer","Inner"] for j in ["alpha-beta","beta-delta",
        "delta-alpha","alpha-gamma","gamma-alpha","alpha_b","beta",
        "delta","alpha_g","gamma"] ])
    for  frm in range(frames): 
        n_avg.iloc[frm] = iterate_sites(n_avg.iloc[frm], density_df.iloc[frm], bind_site, radius * dr * dth) 
    return n_avg
    

def histogram_sites(n_avg_sites):
    hist, edge = pd.DataFrame(columns=(n_avg_sites.columns)),pd.DataFrame(columns=(n_avg_sites.columns))
    for key in n_avg_sites:
        hist[key], edge[key] = np.histogram(n_avg_sites[key][:-1], np.arange(0,200,1))
    return hist, edge

def max_count(df_in,lim=1):
    out = pd.Series(index=df_in.columns,dtype=object)
    for key in df_in.columns:
        out[key] = np.where(df_in[key][lim:] == np.max(df_in[key][lim:]))[0][0]+lim
    return out

def mean_count(df_in, lim=1):
    return df_in.iloc[lim:, :].mean()


def get_intersites():
        inter_sites = ["%s_%s"%(i,j) 
        for i in ["Outer","Inner"] for j in ["alpha-beta","beta-delta",
        "delta-alpha","alpha-gamma","gamma-alpha","alpha_b","beta",
        "delta","alpha_g","gamma"] ]
        return inter_sites

def get_prob(chains_groups, rep):
        
    inter_sites = get_intersites()
    counts = None 
    
    for cgid, cg in enumerate(chains_groups):
        #prob_for_sem = None
        for a in range(1,rep+1):
            #files = ["polar1_nach_%s_de%i.4.upp.dat"%(cg,4),"polar1_nach_%s_de%i.4.low.dat"%(cg,4)]
            files = [ "polar1_synpR_%s_a%i.4.upp.dat"%(cg,a),"polar1_synpR_%s_a%i.4.low.dat"%(cg,a) ]
            den_tmp = pd.DataFrame(index=chains_groups, columns=["Outer","Inner"])
            rad, dr, dth, theta, radius, frames = pc.Coord_Get(files[0])
            # exp_low, exp_upp, norm = pg.set_expected("",)  
            # exp_lowP, exp_uppP=pg.set_expected_prob(exp_upp, exp_low, norm, norm) 
            # expA_low, expA_upp, norm=pg.set_expected_A("")
            # expA_lowP, expA_uppP=pg.set_expected_prob(expA_upp, expA_low, norm, norm)
            if (a == 1):
                counts = pd.DataFrame(index=(np.linspace(0, 199,200)), columns=(inter_sites)).fillna(0.0)
            
            chains_up = pc.prot_coord()
            bindsites = pc.cholesterol_binding_selection(chains_up[0],radius,theta), pc.cholesterol_binding_selection(chains_up[1],radius,theta),pc.omega3_binding_selection(chains_up[0],radius,theta),pc.omega3_binding_selection(chains_up[1],radius,theta)
            print(cg,a)
            den_tmp["Outer"][cg] = pc._analysis_call_(files[0], radius, dr, dth, frames, enrich=False) # not sure how the output will look
            den_tmp["Inner"][cg] = pc._analysis_call_(files[1], radius, dr, dth, frames, enrich=False) # runn that line!
            
            #This is only going to work short term....
            den_tmp = densit_frames(den_tmp, cg, frames)
            n_avg = nAvg(den_tmp, radius, dr, dth, frames, bindsites)
            #counts = pd.DataFrame(columns=inter_sites)
            counts_temp,edge = histogram_sites(n_avg)
            #prob_for_sem = counts / len(n_avg.iloc[:-1,0])
            #prob_for_sem.to_csv("../Data/boundary/stats/%s_m1_a%i.csv"%(cg,a))
            counts = counts.add(counts_temp,axis=1,fill_value=0)
                    
        prob = counts / (rep * len(n_avg.iloc[:-1,0]))
        # if prob.sum().sum() < 19.5 or prob.sum().sum() > 20.5:
        #     print("\nSum of Probabilities do not diverge from expected value:")
        #     print("(%f / 20)"%prob.sum().sum())
        #     print("Will not write output files.")
        #     print("Endding\n")
        #     return -1
        counts.to_csv("../Data/deltaG/%s_counts.csv"%cg)
        prob.to_csv("../Data/deltaG/%s_prob.csv"%cg)

def area_info_lim_ex_size(in_list):
    size = (199,)
    for j,i in enumerate(in_list):
        if np.shape(i) != size:
            in_list[j] = i.iloc[:,0]
    return in_list
    

def get_area_info_exp(lipid,areas,ar_vals):
    inter_sites = get_intersites()
    counts = pd.read_csv("../Data/deltaG/DPPC_counts.csv",index_col=0,header=0)

        #files = [ "polar1_nach_%s_de%i.2.upp.dat"%(lip,1),"polar1_snach_%s_de%i.2.low.dat"%(lip,1) ]
    files = [ "polar1_nach_%s_de%i.4.upp.dat"%("DPPC",4),"polar1_nach_%s_de%i.4.low.dat"%("DPPC",4) ]
    
    rad, dr, dth, theta, radius, frames = pc.Coord_Get(files[0])
    for lip in lipid:
        exp_low, exp_upp, norm = pg.set_expected_A("",areas)
        exp_lowP, exp_uppP=pg.set_expected_prob(exp_upp, exp_low, norm, norm) 
            #expA_low, expA_upp, norm=pg.set_expected_A("")
            #expA_lowP, expA_uppP=pg.set_expected_prob(expA_upp, expA_low, norm, norm)

        ex_up = [exp_uppP.T[areas[0]],exp_uppP.T[areas[1]],exp_uppP.T[areas[2]],
                 exp_uppP.T[areas[3]], exp_uppP.T[areas[4]], exp_uppP.T[areas[5]],
                 exp_uppP.T[areas[6]], exp_uppP.T[areas[7]], exp_uppP.T[areas[8]],
                 exp_uppP.T[areas[9]]]
        ex_up = area_info_lim_ex_size(ex_up)
        ex_low = [exp_lowP.T[areas[10]],exp_lowP.T[areas[11]], exp_lowP.T[areas[12]],
                  exp_lowP.T[areas[13]], exp_lowP.T[areas[14]],exp_lowP.T[areas[15]],
                  exp_lowP.T[areas[16]],exp_lowP.T[areas[17]],exp_lowP.T[areas[18]],
                  exp_lowP.T[areas[19]]]
        ex_low = area_info_lim_ex_size(ex_low)
        
        
        Area = pd.DataFrame(columns=(inter_sites), index=["A","acc", "acc_root","res"])
        Area = pg.get_area_accessable(Area, counts, ex_up, ex_low,ar_vals)
        # ex_up,ex_low = pd.DataFrame(ex_up,index=inter_sites[:10]).T, pd.DataFrame(ex_low,index=inter_sites[10:]).T
        # exctP = pd.concat([ex_up,ex_low],axis=1)
        Area.to_csv("../Data/deltaG/%s_Area_curent.csv"%lip)
        # exctP.to_csv("../Data/deltaG/%s_Area_distribution.csv"%lip)

def get_prob_exp(chains_groups="DPPC"):
        
    inter_sites = get_intersites()
    counts = None 
    files = ["polar1_nach_%s_de%i.4.upp.dat"%("CHOL",4),"polar1_nach_%s_de%i.4.low.dat"%("CHOL",4)]
            #files = [ "polar1_synpR_%s_a%i.2.upp.dat"%(cg,a),"polar1_synpR_%s_a%i.2.low.dat"%(cg,a) ]
    den_tmp = pd.DataFrame(index=[chains_groups], columns=["Outer","Inner"])
    rad, dr, dth, theta, radius, frames = pc.Coord_Get(files[0])

    #counts = pd.DataFrame(index=(np.linspace(0, 199,200)), columns=(inter_sites)).fillna(0.0)
            
    chains_up = pc.prot_coord()
    bindsites = pc.cholesterol_binding_selection(chains_up[0],radius,theta), pc.cholesterol_binding_selection(chains_up[1],radius,theta),pc.omega3_binding_selection(chains_up[0],radius,theta),pc.omega3_binding_selection(chains_up[1],radius,theta)
    den_tmp["Outer"][chains_groups] = pc._analysis_call_(files[0], radius, dr, dth, frames, enrich=False) # not sure how the output will look
    den_tmp["Inner"][chains_groups] = pc._analysis_call_(files[1], radius, dr, dth, frames, enrich=False) # runn that line!
            
    den_tmp = densit_frames(den_tmp, chains_groups, frames)
    n_avg = nAvg(den_tmp, radius, dr, dth, frames, bindsites)
    counts,edge = histogram_sites(n_avg)

    #counts = counts.add(counts_temp,axis=1,fill_value=0)
                    
    prob = counts / (1 * len(n_avg.iloc[:-1,0]))
    counts.to_csv("../Data/deltaG/%s_counts.csv"%chains_groups)
    prob.to_csv("../Data/deltaG/%s_prob.csv"%chains_groups)

def get_area_info(lipid,areas):
    inter_sites = get_intersites()
    #counts = None
    # for lip_id, lip in enumerate(lipid):
    #     count = pd.read_csv("../Data/deltaG/%s_counts.csv"%lip,index_col=0,header=0)
    #     if lip_id == 0:
    #         counts = count.copy()
    #     else:
    #         counts = counts + count
    #files = [ "polar1_nach_%s_de%i.4.upp.dat"%("DPPC",4),"polar1_nach_%s_de%i.4.low.dat"%("DPPC",4) ]
    files = [ "polar1_synpR_%s_a%i.4.upp.dat"%("CHOL",1),"polar1_synpR_%s_a%i.4.low.dat"%("CHOL",1) ]
    
    rad, dr, dth, theta, radius, frames = pc.Coord_Get(files[0])
    for lip in lipid:
        exp_low, exp_upp, norm = pg.set_expected("",lip,areas)
        exp_lowP, exp_uppP=pg.set_expected_prob(exp_upp, exp_low, norm, norm) 
            #expA_low, expA_upp, norm=pg.set_expected_A("")
            #expA_lowP, expA_uppP=pg.set_expected_prob(expA_upp, expA_low, norm, norm)

        ex_up = [exp_uppP.T[areas[0]],exp_uppP.T[areas[1]],exp_uppP.T[areas[2]],
                 exp_uppP.T[areas[3]], exp_uppP.T[areas[4]], exp_uppP.T[areas[5]],
                 exp_uppP.T[areas[6]], exp_uppP.T[areas[7]], exp_uppP.T[areas[8]],
                 exp_uppP.T[areas[9]]]
        # area_info_lim is NOT removing duplicated values
        ex_up = area_info_lim_ex_size(ex_up)
        ex_low = [exp_lowP.T[areas[10]],exp_lowP.T[areas[11]], exp_lowP.T[areas[12]],
                  exp_lowP.T[areas[13]], exp_lowP.T[areas[14]],exp_lowP.T[areas[15]],
                  exp_lowP.T[areas[16]],exp_lowP.T[areas[17]],exp_lowP.T[areas[18]],
                  exp_lowP.T[areas[19]]]
        ex_low = area_info_lim_ex_size(ex_low)
        
        
        #Area = pd.DataFrame(columns=(inter_sites), index=["A","acc", "acc_root","res"])
        #Area = pg.get_area_accessable(Area, counts, ex_up, ex_low)
        ex_up,ex_low = pd.DataFrame(ex_up,index=inter_sites[:10]).T, pd.DataFrame(ex_low,index=inter_sites[10:]).T
        exctP = pd.concat([ex_up,ex_low],axis=1)
        #Area.to_csv("../Data/deltaG/%s_Area_curent.csv"%lip)
        exctP.to_csv("../Data/deltaG/%s_Area_distribution3.csv"%lip)

def average_area_files(areas):
    avg = None
    for a in areas:
        avg = None
        for leaf in ["upp", "low"]:
            for rep in range(1,11):
                if rep == 1:
                    avg = pd.read_csv("../Data/polar/exp_1_synpR_a%i.2_%s_area_expected.%s.dat"%(rep,a,leaf),
                                      sep=" ",index_col=0,header=None)
                else:
                    avg = avg + pd.read_csv("../Data/polar/exp_1_synpR_a%i.2_%s_area_expected.%s.dat"%(rep,a,leaf),
                                      sep=" ",index_col=0,header=None)
            avg = avg / rep
            avg.to_csv("../Data/polar/exp_3_synpR_%s_area_expected.%s.dat"%(a,leaf),sep=" ",index=True)
    
# # Compress does not appear to be nessisary
# def average_compresss(areas):
#     for a in areas:
#         for leaf in ["upp", "low"]:
#             inpt = pd.read_csv("../Data/polar/exp_1_synpR_%s_area_expected.%s.dat"%(a,leaf),sep=" ",index_col=0)
#             avg = inpt.sum()
#             avg.to_csv("../Data/polar/exp_1_synpR_%s_area_compressed.%s.dat"%(a,leaf),sep=" ",index=True)
    
    


def plot_prob_dist(lipids):
    np.random.seed(0)
    intersites = get_intersites()
    # clrs = ["darkgray","darkslategray","darkolivegreen","forestgreen","darkred","olive"
    #         ,"darkslateblue","steelblue","chocolate","yellowgreen","indianred",
    #         "darkblue","limegreen","darkseagreen","purple", "orangered",
    #         "darkturquoise","orange","gold","mediumvioletred","lawngreen",
    #         "mediumspringgreen","blueviolet","crimson","blue","fuchsia","dodgerblue",
    #         "khaki","plum","mediumslateblue","lightsalmon","violet","lightskyblue",
    #         "aquamarine","bisque","lightpink"]
    fig,ax = plt.subplots(2,int(len(intersites)/2), sharex=True, sharey=True,
                          figsize=(10,2.5))

    for lipI,lip in enumerate(lipids):
        pro_dist = pd.read_csv("../Data/deltaG/%s_prob.csv"%lip,index_col=0)
        print(pro_dist.sum().sum())
        area_dist = pd.read_csv("../Data/deltaG/%s_Area_distribution3.csv"%lip,index_col=0)
        #print(area_dist.sum())
        insJ, insI = 0,0
        for ins in intersites:          
            if insI>=10:
                insI=0
                insJ=1
            ax[insJ,insI].plot(pro_dist[ins],label=lip)#,c=clrs[lipI])
            ax[insJ,insI].plot(area_dist[ins],'--',)#c=clrs[lipI])
            ax[insJ,insI].set_xlim(-5,65)
            if insJ == 0:
                ax[insJ,insI].set_title(ins[6:]) 
            insI = insI+1 
    ax[0,-1].legend(lipids,bbox_to_anchor=(.14, 1.25),ncol=12)
    #plt.tight_layout()
    # plt.savefig("prob_CHOL.pdf")
    # plt.close()
    
def plot_mean_dist(lipids):
    np.random.seed(0)
    intersites = get_intersites()
    mean_sites = ["Out_inter","Out_M4", "In_inter","In_M4"]
    color_list = ["#1f77b4","#ff7f0e","#2ca02c","#d62728", "#9467bd", "#8c564b", 
                  "#e377c2","#7f7f7f"]    
    fig,ax = plt.subplots(len(lipids),4, sharex=True, sharey=True,
                          figsize=(10,10))
    for lipI,lip in enumerate(lipids):
        pro_dist = pd.read_csv("../Data/deltaG/%s_prob.csv"%lip,index_col=0)
        #area_dist = pd.read_csv("../Data/deltaG/%s_Area_distribution3.csv"%lip,index_col=0)
        area_dist = pd.read_csv("../Data/deltaG/%s_counts.csv"%lip,index_col=0)
        insJ, insI = 0,0
        for ins in range(0,4):          
            if insI%4==0:
                insI= 0 
            print(ins*5,(ins+1)*5)
            #ax[lipI,insI].plot(pro_dist[intersites[ins*5:(ins+1)*5]].mean(axis=1),c=color_list[lipI])
            ax[lipI,insI].plot(area_dist[intersites[ins*5:(ins+1)*5]].mean(axis=1),"--",c=color_list[lipI])
            ax[lipI,insI].set_xlim(0,30)
            ax[lipI,insI].set_ylim(0,35)
            #ax[insI,lipI].set_title(mean_sites[ins]) 
            insI = insI + 1
    #ax[2].legend(bbox_to_anchor=(.14, 1.25),ncol=12)
    plt.tight_layout()
    plt.savefig("Anionic_Counts.pdf")#"Hypoth3.pdf")
    plt.close()

def get_std(probs):
    std = []
    for prob in probs:
        std.append(probs[prob].std())
    print(std)
    return std
###############################################################################

##################### Delta G Routines ########################################
def get_prob_ratio(lipids,areas):
    intersites = get_intersites()
    PL = pd.DataFrame(index=lipids, columns=intersites)
    PG = pd.DataFrame(index=lipids, columns=intersites)
    for lip in lipids:
        prob = pd.read_csv("../Data/deltaG/%s_prob.csv"%lip,index_col=0)
        exptP = pd.read_csv("../Data/deltaG/%s_Area_distribution3.csv"%lip,index_col=0)
        #ex_max = get_std(exptP)
        ex_max = [exptP[exp].argmax() for exp in exptP]        
        for icol,col in enumerate(prob.columns):
            PL[col][lip] = prob[col][:ex_max[icol]].sum()
            PG[col][lip] = prob[col][ex_max[icol]+1:].sum()
    dP = calc_prob_ratio(PL,PG)
    return PL, PG, dP        
        
def calc_prob_ratio(PL,PG):
    return pd.DataFrame(np.array(PG.values.tolist()) / np.array(PL.values.tolist()),
                                  index=PL.index, columns=PL.columns)


def write_prob_ratios(lipids,areas):
    PL, PG, dP = get_prob_ratio(lipids,areas)
    PL.to_csv("../Data/deltaG/PL_stand.csv")
    PG.to_csv("../Data/deltaG/PG_stand.csv")
    dP.to_csv("../Data/deltaG/dP_stand.csv")

def write_dG(lipids):
    intersites = get_intersites()
    dP = pd.read_csv("../Data/deltaG/dP_stand.csv",index_col=0)
    RT = 0.642748322147651006711409395973154362416107382550335570469
    dG = pd.DataFrame(-RT*np.log(dP) ,index=lipids, columns=intersites)
    dG.to_csv("../Data/deltaG/dG_stand.csv")

def plot_dG_map(lipids,lab,title):
    fig,ax = plt.subplots(1,1,figsize=(15,10))
    dG = None
    if lab == "standard":
        dG = pd.read_csv("../Data/deltaG/dG_test.csv",index_col=0)
        minV = -1.5
        maxV = 5
    if lab == "order_sites":
        dG = pd.read_csv("../Data/deltaG/Average_dG.csv",index_col=0)
        minV = -1
        maxV = 2.5
    dG = dG.T[lipids].T.sort_values(by=dG.columns[0],ascending=False)
    cmap = plt.cm.RdBu#PuOr
    minV = -1 - 0.5
    maxV = 3+ 0.5
    cmap.set_bad(color='black')
    norm = pg.MidpointNormalize(midpoint=0,vmin=minV,vmax=maxV)    
    c = plt.pcolormesh(dG,norm=norm,cmap=cmap)
    for (i, j), z in np.ndenumerate(dG):
        plt.text(j+.5, i+.5, '{:0.2f}'.format(z), ha='center', va='center')
    plt.yticks(np.linspace(0.5,len(dG.index)-.5,len(dG.index)),labels=dG.index.tolist(),fontsize=12)
    plt.xticks(np.linspace(0.5,len(dG.columns)-.5,len(dG.columns)),
               labels=dG.columns.tolist(),rotation=90, fontsize=12)
    plt.colorbar(c)
    #plt.gca().set_aspect('equal', adjustable='box')
    #plt.tight_layout()
    plt.savefig("dG_%s.pdf"%title)
    plt.close()
    
def get_dg_histogram(lipids, sort_key):
    
    # 1. how to split up head groups (ends with Ne/Zw)
    # 2. how to split up by chains (starts with n0,....)
    # 3. how to split up by bins (use the col names)
    
    dg = pd.read_csv("../Data/deltaG/Average_dG.csv",index_col=0,header=0)
    dg_err = pd.read_csv("../Data/deltaG/stats/SEM_dG_Mean.csv",index_col=0,header=0)
    lip = None
    if sort_key == 1:
        fig,ax = plt.subplots(2,1, sharex=True,sharey=True)
        dg = dg.T
        lip = [["n0Zw","n9Zw","n6Zw","n3Zw"],["n0Ne","n9Ne","n6Ne","n3Ne"]]
        i,j =0,0
        for li, l in enumerate(lip):
            dg[l].plot(ax=ax[li], kind='bar',legend=False)
            dg[l].plot(ax=ax[li], kind='bar',legend=False)
        ax[0].legend(["Sat","Mono","n-6","n-3"],ncol=len(l))
        plt.savefig("HG_sort.pdf")
        plt.close()
    if sort_key == 2:
        fig,ax = plt.subplots(2,2, sharex=True,sharey=True)
        dg = dg.T
        lip = [["n0Zw","n0Ne"],["n9Zw","n9Ne"],["n6Zw","n6Ne"],["n3Zw","n3Ne",]]
        i,j =0,0
        for li, l in enumerate(lip):
            if li == 2:
                i = 0
                j = 1
            dg[l].plot(ax=ax[i,j], kind='bar')
            i = i+1
        plt.savefig("Acyl_sort.pdf")
        plt.close()
    if sort_key == 3:
        clrs = ['#1f77b4','#ff7f0e']
        lip = ["n0Zw","n9Zw","n6Zw","n3Zw","n0Ne","n9Ne","n6Ne","n3Ne"]
        fig,ax = plt.subplots(2,2, sharex=True,sharey=True)
        #lip = [["n0Zw","n0Ne"],["n9Zw","n9Ne"],["n6Zw","n6Ne"],["n3Zw","n3Ne",]]
        i,j =0,0
        for li, l in enumerate(dg):
            if li == 2:
                i = 0
                j = 1
            dg[l][lip].plot(ax=ax[i,j], kind='bar',color=clrs[i])
            i = i+1
        plt.savefig("Occ_sort.pdf")
        plt.close()
        
    if sort_key == 4:
        color_list = ["#1f77b4","#2ca02c","#ff7f0e","#d62728"]
        lip = ["CHOL","n0", "n9", "n6", "n3"]
        fig,ax = plt.subplots(2,1, sharex=True,sharey=True, figsize=(10,5))
        dg_neut = pd.DataFrame(index=["CHOL","n0Zw","n9Zw","n6Zw","n3Zw"],
                              columns=(dg.columns))
        dg_anio = pd.DataFrame(index=["CHOL", "n0Ne","n9Ne","n6Ne","n3Ne"],
                              columns=(dg.columns))
        # er_neut = pd.DataFrame(index=["n0Zw","n9Zw","n6Zw","n3Zw","CHOL"],
        #                       columns=(dg.columns))
        # er_anio = pd.DataFrame(index=["n0Ne","n9Ne","n6Ne","n3Ne","CHOL"],
        #                       columns=(dg.columns))
        for si, site in enumerate(dg):
            for li, l in enumerate(lip):
                if l == "CHOL":
                    dg_neut[site]["%s"%l] = dg[site]["%s"%l]
                    dg_anio[site]["%s"%l] = 0.0
                    # er_neut[site]["%s"%l] = dg_err[site]["%s"%l]
                    # er_anio[site]["%s"%l] = 0.0
                else:
                    dg_neut[site]["%sZw"%l] = dg[site]["%sZw"%l]
                    if site.startswith("Outer"):
                        dg_anio[site]["%sNe"%l] = np.nan
                        dg_err[site]["%sNe"%l] = np.nan
                    else:
                        dg_anio[site]["%sNe"%l] = dg[site]["%sNe"%l]
                    # er_neut[site]["%s"%l] = dg_err[site]["%sZw"%l]
                    # er_anio[site]["%s"%l] = dg_err[site]["%sNe"%l]
        #er_neut = pd.DataFrame(columns=dg.columns,index=["n0Zw","n9Zw","n6Zw","n3Zw","CHOL"])
        dg_neut.plot(ax=ax[0],kind="bar",color=color_list,legend=False, use_index=True, rot=0,
                     yerr=dg_err.T[["CHOL", "n0Zw","n9Zw","n6Zw","n3Zw"]].T, width=.75)
        
        dg_anio.plot(ax=ax[1],kind="bar",legend=False,color=color_list, use_index=True, rot=0,
                     yerr=dg_err.T[["n0Ne","n9Ne","n6Ne","n3Ne"]].T,width=.75)
        ax[0].set_xlabel("")
        #ax[0,0].legend(dg.columns,ncol=4)
        plt.savefig("Key_sort4_2.pdf")
        plt.close()

    if sort_key == 5:
        lip = ["n0", "n9", "n6", "n3"]
        color_list = ["#d62728","#1f77b4","#fd7082","#ff7f0e","#ebd8c5"]
        fig,ax = plt.subplots(2,2, sharex=True,sharey=True, figsize=(10,5))
        for si, s in enumerate(["Inter","M4"]):
            dg_tmp = pd.DataFrame(columns=["Neutral","Negative"],
                                  index=["CHOL", "n0","n9","n6","n3"])
            
            er_tmp = pd.DataFrame(columns=["Neutral","Negative"],
                                  index=["CHOL", "n0","n9","n6","n3"])
            
            dg_tmp["Neutral"][["CHOL","n0","n9","n6","n3"]] = dg["Outer-%s"%s][["CHOL","n0Zw","n9Zw","n6Zw","n3Zw"]]
            #dg_tmp["Negative"][["n0","n9","n6","n3"]] = dg["Outer-%s"%s][["n0Ne","n9Ne","n6Ne","n3Ne"]]
            
            er_tmp["Neutral"][["CHOL","n0","n9","n6","n3"]] = dg_err["Outer-%s"%s][["CHOL","n0Zw","n9Zw","n6Zw","n3Zw"]]
            #er_tmp["Negative"][["n0","n9","n6","n3"]] = dg_err["Outer-%s"%s][["n0Ne","n9Ne","n6Ne","n3Ne"]]
            
        
            
            #PLOT
            dg_tmp.T.plot(ax=ax[0,si],kind="bar",legend=False,rot=0,color=color_list, yerr=er_tmp.T,width=1)
            ax[0,0].set_xlabel("")
            
            
            
            
            
            dg_tmp = pd.DataFrame(columns=["Neutral","Negative"],
                                  index=["CHOL","n0","n9","n6","n3"])
            er_tmp = pd.DataFrame(columns=["Neutral","Negative"],
                                  index=["CHOL","n0","n9","n6","n3"])
            
            dg_tmp["Neutral"][["CHOL","n0","n9","n6","n3"]] = dg["Inner-%s"%s][["CHOL","n0Zw","n9Zw","n6Zw","n3Zw"]]
            dg_tmp["Negative"][["n0","n9","n6","n3"]] = dg["Inner-%s"%s][["n0Ne","n9Ne","n6Ne","n3Ne"]]
            er_tmp["Neutral"][["CHOL","n0","n9","n6","n3"]] = dg_err["Inner-%s"%s][["CHOL","n0Zw","n9Zw","n6Zw","n3Zw"]]
            er_tmp["Negative"][["n0","n9","n6","n3"]] = dg_err["Inner-%s"%s][["n0Ne","n9Ne","n6Ne","n3Ne"]]
            #PLOT
            dg_tmp.T.plot(ax=ax[1,si],kind="bar",legend=False, rot =0,color=color_list,
                           yerr=er_tmp.T,width=1)

        plt.savefig("Key_sort5.pdf")
        plt.close()

        
def write_dG_tables_latex():

    dg = pd.read_csv("../Data/deltaG/dG_stand.csv",index_col=0,header=0)
    dg_err = pd.read_csv("../Data/deltaG/stats/SEM_dG.csv",index_col=0,header=0)
    
    dg_outer_inter, err_outer_inter = dg.iloc[:,0:5], dg_err.iloc[:,0:5]
    dg_outer_m4, err_outer_m4 = dg.iloc[:,5:10], dg_err.iloc[:,0:5]
    dg_inner_inter, err_inner_inter = dg.iloc[:,10:15], dg_err.iloc[:,0:5]
    dg_inner_m4, err_inner_m4 = dg.iloc[:,15:], dg_err.iloc[:,0:5]
    
    dg_list = [dg_outer_inter, dg_outer_m4, dg_inner_inter, dg_inner_m4]
    err_list = [err_outer_inter, err_outer_m4, err_inner_inter, err_inner_m4]
    for i,j,k in zip(range(len(dg_list)),["outer","outer","inner","inner"],["inter","m4","inter","m4"]):
        dg_list[i]["Avg"] = dg_list[i].mean(axis=1)
        dg_list[i]["STD"] = (err_list[i]*np.sqrt(5)).mean(axis=1)
        dg_list[i]["SEM"] = err_list[i].mean(axis=1)
        dg_list[i].to_csv("./%s_%s.csv"%(j,k))
        dg_list[i] = dg_list[i].round(1).to_latex()
        fl = open ("./%s_%s.tex"%(j,k), "w")
        fl.write(dg_list[i])
        fl.close()
###############################################################################
############################ Runners ##########################################
# get_dg_histogram(None,3)
def build_data_set(chains_groups,area,set_area=False,ar_vals=None):
    if set_area == True:
        get_prob_exp()
        get_area_info_exp(chains_groups,area,ar_vals)
        return 0
    get_area_info(chains_groups,area)
    get_prob(chains_groups, 1)
    return 0

def build_probability_dg(chains_groups, area):
    write_prob_ratios(chains_groups, area)
    write_dG(chains_groups)
    return 0
###############################################################################

if __name__=="__main__":

    ######################## "GLOBAL" Values ######################################

    chains_groups = ["CHOL","DOPC", "DPPC", "DPPS", "DPSM", "OAPE", "OIPC",
                    "OIPE", "OUPC", "OUPE", "OUPS", "PAP1", "PAP2", "PAP3", 
                    "PAPA", "PAPC", "PAPE", "PAPI", "PAPS", "PBSM", "PFPC", 
                    "PIPI", "PNSM", "POP1", "POP2", "POP3", "POPA", "POPC", 
                    "POPE", "POPI", "POPS", "POSM", "PUPC", "PUPE", "PUPI",
                    "PUPS","PC","PE","SM","PS","PA","PI","P1a","P2a","P3a",
                    "n0Zw","n9Zw","n6Zw","n3Zw","n0Ne","n9Ne","n6Ne","n3Ne",
                    "n0","n9","n6","n3","Neutral","Anionic"]
    ar = ["104","50.","63.","81.","34.","173","117","184","195","161","68.","74.","68.","63.","40.","211","273","205","302","273"]

    #get_error_ratio(chains_groups,10)

    #get_dG_error(PG, PL)

    # chains_groups = ["CHOL", "PUPS", "DPPS","POPI"]
    #ceg.area_error(ar)
    #ceg.calc_area_error(chains_groups,10)
    #get_prob_error(chains_groups, 10)
    #calc_error(chains_groups,10)
    #get_error_ratio(chains_groups, 10)
    # get_dg_histogram(["CHOL","n0", "n9","n6","n3"],4)
    # get_dg_histogram(["CHOL","n0", "n9","n6","n3"],5)
    # build_data_set(chains_groups,ar)#,True,ar_vals)
    # build_probability_dg(chains_groups,ar)
    plot_mean_dist(["PS","PA","PI","P1a","P2a","P3a"])
