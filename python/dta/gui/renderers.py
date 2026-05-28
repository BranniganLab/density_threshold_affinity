#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Matplotlib renderers for polar bin selection.

This module contains view-layer classes responsible for drawing selection
outlines using Matplotlib.

All event handling and gesture interpretation is handled by
Matplotlib controllers in selectors.py.
"""
from collections.abc import Iterable
import matplotlib
from dta.bin_logic.utils import BinEdge


class SelectionRenderer:
    """Renderer for drawing polar objects on a matplotlib Axes."""

    def __init__(
        self,
        ax: matplotlib.axes.Axes,
        plot_kwargs: dict = None
    ) -> None:
        """
        Create a SelectionRenderer object and tie it to an Axes instance.

        Default plot kwargs are
        {
            "color": "red",
            "lw": 2.0,
            "zorder": 20,
        }
        but can be overridden or added to with the plot_kwargs argument.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes used for drawing.
        plot_kwargs : dict
            Dictionary of matplotlib.pyplot keywords for drawing bin edges.
        """
        self.ax = ax
        default_plot_kwargs = {
            "color": "red",
            "lw": 2.0,
            "zorder": 20
        }
        if plot_kwargs is None:
            self.plot_kwargs = default_plot_kwargs
        else:
            self.plot_kwargs = default_plot_kwargs | plot_kwargs
        self.selected_artists = []
        self.hover_artists = []

    def draw_edges(
        self,
        edges: Iterable[BinEdge],
        plot_kwargs: dict = None
    ) -> list[matplotlib.artist.Artist]:
        """
        Draw a collection of bin edges.

        Parameters
        ----------
        edges : iterable of BinEdge
            Edges to draw.
        plot_kwargs : dict
            Keyword arguments passed to ``Axes.plot``.

        Returns
        -------
        list
            Matplotlib artist objects created.
        """
        artists = []

        # check if default kwargs are present; update/override if possible.
        kwargs = dict(self.plot_kwargs)
        if plot_kwargs is not None:
            kwargs.update(plot_kwargs)

        # plot each BinEdge individually
        for edge in edges:
            theta_endpoints = (edge.endpoint1[1], edge.endpoint2[1])
            r_endpoints = (edge.endpoint1[0], edge.endpoint2[0])
            artists.append(self.ax.plot(theta_endpoints, r_endpoints, **kwargs)[0])

        # return list of artists so they can be un-plotted in the future if desired.
        return artists

    def shade_interior_region(self) -> None:
        """Hold space for future method."""
