#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:44:55 2024

@author: js2746
"""

import numpy as np

class Site:
    """
    The basic class for a binding site on/in a protein/inclusion. User defines \
    Site with bin coordinates. Multiple symmetric Sites can be combined with \
    Symmetric_Site class. One Site (or Symmetric_Site) across multiple replicas\
    can be combined in Site_Across_Replicas class.
    
    Attributes
    ----------
    
    """