#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 11:28:04 2025.

@author: js2746
"""

import numpy as np
from density import SysInfo, Dimensions


class SpeciesCounts:
    
    def __init__(self):
        self.outer_counts = None
        self.inner_counts = None
        self.outer_sysinfo = None
        self.inner_sysinfo = None
        self.grid_dims = None

    def add_counts(self, leaflet_id, counts, sysinfo, grid_dims):
        assert leaflet_id in ['upp', 'low'], "leaflet_id must be 'upp' for outer leaflet or 'low' for inner leaflet."
        assert isinstance(counts, np.ndarray), "counts must be a numpy ndarray."
        assert len(counts.shape) == 2, "counts must be a 2D array."
        assert isinstance(sysinfo, tuple)
        assert isinstance(grid_dims, tuple)
        if leaflet_id == 1:
            self.outer_counts = counts
            self.outer_sysinfo = sysinfo
        elif leaflet_id == 2:
            self.inner_counts = counts
            self.inner_sysinfo = sysinfo
        if self.grid_dims is not None:
            assert grid_dims == self.grid_dims, "grid_dims must match"
        else:
            self.grid_dims = grid_dims


class AggregateCounts:
    """
    Aggregate SpeciesCounts objects from SAME REPLICA together into one class.

    This object CANNOT be used to aggregate objects from multiple replicas \
    together.

    """

    def __init__(self):
        self.grid_dims = None
        self.contents = []

    def add_counts(self, counts_object):
        assert isinstance(counts_object, SpeciesCounts), "counts_object must be a SpeciesCounts object."
        assert (counts_object.outer_counts is not None) or (counts_object.inner_counts is not None), "counts_object must have counts in at least one leaflet."
        if self.grid_dims is not None:
            assert counts_object.grid_dims == self.grid_dims, "grid_dims must match"
        else:
            self.grid_dims = counts_object.grid_dims
        self.contents.append(counts_object)

    def get_counts(self, leaflet_id):
        assert (self.grid_dims is not None) and (len(self.contents) != 0), "This AggregateCounts object is empty"
        assert leaflet_id in ['upp', 'low'], "leaflet_id must be 'upp' for outer leaflet or 'low' for inner leaflet."
        counts_array = np.zeros((self.grid_dims.Nr, self.grid_dims.Ntheta))
        for obj in self.contents:
            if obj.outer_counts is None:
                raise Exception(f'SpeciesCounts object is empty')
            else:
                if leaflet_id == 'upp':
                    counts_array += obj.outer_counts
                elif leaflet_id == 'low':
                    counts_array += obj.inner_counts
        return counts_array

    def get_exp_density(self, leaflet_id):
        assert (self.grid_dims is not None) and (len(self.contents) != 0), "This AggregateCounts object is empty"
        assert leaflet_id in ['upp', 'low'], "leaflet_id must be 1 for outer leaflet or 2 for inner leaflet."
        exp_density = 0
        for obj in self.contents:
            if obj.outer_sysinfo is None:
                raise Exception(f'SpeciesCounts object is empty')
            else:
                if leaflet_id == 'upp':
                    exp_density += obj.outer_sysinfo.ExpBeadDensity
                elif leaflet_id == 'low':
                    exp_density += obj.inner_sysinfo.ExpBeadDensity
        return exp_density
