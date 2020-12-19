#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 17:56:25 2020

@author: liam
"""

import numpy as np
import pandas as pd
import polarDensity_helper as pc



def Density(memb, ddg=False, lipids=False, enrich=False):
    '''
    memb (str) should be synps or oocyt
    ddg (bool) determines if ddg analysis called
    lipids (str) points to a file to check lipids
    enrich (bool) tells script to check enrichment
    '''
    
    # initialization
    
    global Ntheta
    Ntheta = 30
    
    file_title1 = "polar1_"
    density = None
    chol_bind = None
    w3_bind = None
    rat = None
    chains_groups = None 
    col = []
    chains_up, chains_lo = pc.prot_coord()
    
    # call sanity check
    if memb=="DPPC":
        import PolarSanity
        return PolarSanity.Sanity()
        
    # Catches your using the right args
    assert memb=="synps" or memb=="oocyt", "memb var should be synps or oocyt"
    
    #sets up initial strings to collect data files
    if memb == "synps":
        file_title1 = file_title1 + "synpR"
        #rat = pc.get_rat_val("../Data/ratio/Ratio_Synapse.csv")
        if lipids == False:    
            chains_groups = pc.Syn_ooc_comp()[1]
        else:
            chains_groups = lipids
            
    if memb == "oocyt":
        file_title1 = file_title1 + "oocyR"
        #rat = pc.get_rat_val("../Data/ratio/Ratio_Oocyte.csv")
        if lipids == False:    
            chains_groups = pc.Syn_ooc_comp()[0]
        else:
            chains_groups = lipids
    
    # get files to use
    files = ["%s_%s_a%i.1.%s.dat"%(file_title1,cg,ai,leaf) for ai in range(1,11) 
        for cg in chains_groups for leaf in ["upp", "low"]]
    
    # I want to confirm all the files are being used for best
    # results; ie file check
    if (pc.file_check(files)==False): 
        print("Error: You are missing files. Please check and re-run")
        return 0
    
    # this get column names for pandas
    col = pc.get_col_names(files)
    den_tmp = pd.DataFrame(index=chains_groups, columns=col)
    
    for fli, fl in enumerate(files):
        if fli == 0:
            rad, dr, dth, theta, radius, frames = pc.Coord_Get(fl)
                #if args.ddg:
                #    rat = pd.read_csv("Ratio_%s.csv"%membrn, index_col=0)
            #if args.cbs:
            chol_bind = pc.cholesterol_binding_selection(chains_up,radius,theta)
            #if args.w3bs:
            w3_bind =  pc.omega3_binding_selection(chains_up,radius,theta)
        # brute force to parse out specific lipid times or groups (PIPs and Neutral/Anionic)
        if len(fl) == 30:
            tmp_chain = fl[13:17]
            tmp_nm = fl[18:26] 
        elif len(fl) == 28:
            tmp_chain = fl[13:15]
            tmp_nm = fl[16:24] 
        elif len(fl) == 29:
            tmp_chain = fl[13:16]
            tmp_nm = fl[17:25]
            if tmp_nm.startswith("a") ==  False:
                tmp_nm = "a"+tmp_nm 
            if  (tmp_chain == "P1") or (tmp_chain == "P2") or (tmp_chain == "P3"):
                tmp_chain =  tmp_chain + "a"
            if (tmp_chain.endswith("_")):
                tmp_chain = tmp_chain[:-1]
        elif len(fl) > 30:
            tmp_nm = fl[21:-4]
            if "Neutral" in fl:
                tmp_chain = "Neutral"
            elif "Anionic" in fl:
                tmp_chain = "Anionic"
        
        # This is a hack. The above part does not have a "flexible"
        # method to consider sim type (a, b ...)

        print(tmp_chain,tmp_nm)
        den_tmp.at[tmp_chain,tmp_nm] = pc._analysis_call_(fl, radius, dr, 
                dth, frames, enrich=enrich)
    if ddg == False:
        #plot density or enrichment density
        pc._Polar_Plot_(den_tmp, theta, radius, chains_groups,memb)
    
   
# enrich and lipids have not been implemented  

chains_groups = ["CHOL","DOPC", "DPPC", "DPPS", "DPSM", "OAPE", "OIPC",
                 "OIPE", "OUPC", "OUPE", "OUPS", "PAP1", "PAP2", "PAP3", 
                 "PAPA", "PAPC", "PAPE", "PAPI", "PAPS", "PBSM", "PFPC", 
                 "PIPI", "PNSM", "POP1", "POP2", "POP3", "POPA", "POPC", 
                 "POPE", "POPI", "POPS", "POSM", "PUPC", "PUPE", "PUPI",
                 "PUPS"]     
HG = ["PC","PE","SM","PS","PA","PI","P1a","P2a","P3a"]
acyl = ["n0","n9","n6","n3"]
ddG_sum= Density("synps", ddg=False, lipids=["Neutral","Anionic"], enrich=True) #lipids=["n0","n9","n6","n3","CHOL"]
