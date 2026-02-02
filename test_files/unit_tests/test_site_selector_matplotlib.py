#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 16:02:37 2026

@author: js2746
"""

import numpy as np
from matplotlib.testing.decorators import check_figures_equal

from DTA.core import PolarBinGrid, BinEdge
from DTA.gui import SiteSelector
from DTA.gui.matplotlib import PolarBinRenderer


@check_figures_equal(extensions=("png",))
def test_renderer_draw_edges_matches_manual(fig_test, fig_ref):
    # test figure: uses renderer
    ax_t = fig_test.add_subplot(111, projection="polar")
    renderer = PolarBinRenderer(ax_t, plot_kwargs={"color": "k", "lw": 1, "zorder": 10})
    edges = [
        BinEdge((0.0, 1.0), (0.0, 0.0)),
        BinEdge((1.0, 1.0), (0.0, np.pi / 2)),
    ]
    renderer.draw_edges(edges)

    # reference figure: manual plotting of identical segments
    ax_r = fig_ref.add_subplot(111, projection="polar")
    for e in edges:
        ax_r.plot(e.theta_endpoints, e.r_endpoints, color="k", lw=1, zorder=10)


@check_figures_equal(extensions=("png",))
def test_selector_draw_committed_matches_expected(fig_test, fig_ref):
    theta_edges = np.linspace(0, 2 * np.pi, 9)
    r_edges = np.linspace(0, 1, 5)
    bins = {(0, 0), (0, 1), (1, 1)}

    # test figure: selector draws committed selection
    ax_t = fig_test.add_subplot(111, projection="polar")
    sel = SiteSelector(ax_t, theta_edges, r_edges, plot_kwargs={"color": "r", "lw": 2, "zorder": 20})
    sel.model.set(bins)
    sel._draw_committed()

    # reference figure: compute edges and draw manually
    ax_r = fig_ref.add_subplot(111, projection="polar")
    grid = PolarBinGrid(theta_edges, r_edges)
    edges = grid.exposed_edges(bins)
    for e in edges:
        ax_r.plot(e.theta_endpoints, e.r_endpoints, color="r", lw=2, zorder=20)
