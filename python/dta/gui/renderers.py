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

    def __init__(self, ax: matplotlib.axes.Axes, plot_kwargs: dict = None) -> None:
        """
        Create a SelectionRenderer object and tie it to an Axes instance.

        Default plotting colors are set, but can be overridden or added to when
        calling drawing methods.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes used for drawing.
        plot_kwargs : dict
            Dictionary of matplotlib.pyplot keywords for drawing bin edges.
        """
        self.ax = ax
        self.preview_artists = []
        self.selection_artists = []

        self.default_preview_kwargs = {
            "color": "orange",
            "lw": 2.0,
            "zorder": 1000
        }
        default_selection_kwargs = {
            "color": "red",
            "lw": 2.0,
            "zorder": 20
        }
        if plot_kwargs is not None:
            default_selection_kwargs.update(plot_kwargs)
        self.selection_kwargs = default_selection_kwargs

    def draw_bin_edges(
        self,
        edges: list[BinEdge],
        preview: bool = True
    ) -> None:
        """
        Draw a collection of bin edges.

        Parameters
        ----------
        edges : set[BinEdge]
            Edges to draw.
        preview : bool, optional
            If True, draw an orange selection preview. If False, draw actual
            selection.

        Side Effects
        ------------
        - Removes previous artists
        - Saves drawn artists to the correct SelectionRenderer attribute.
        """
        self.clear_artists(clear_preview=True)
        self.clear_artists(clear_preview=False)
        if preview:
            kwargs = self.default_preview_kwargs
            artists = self.preview_artists
        else:
            kwargs = self.selection_kwargs
            artists = self.selection_artists

        for edge in edges:
            artists.append(self._draw_edge(edge, kwargs))

    def _draw_edge(self, edge: BinEdge, kwargs: dict) -> matplotlib.artist.Artist:
        """
        Draw a BinEdge on the polar Axes and return the artist generated.

        Parameters
        ----------
        edge : BinEdge
            The BinEdge you wish to plot.
        kwargs : dict
            The matplotlib plotting kwargs to use when plotting.

        Side Effects
        ------------
        Plots the BinEdge on the Axes.
        """
        theta_endpoints = (edge.endpoint1[1], edge.endpoint2[1])
        r_endpoints = (edge.endpoint1[0], edge.endpoint2[0])
        return self.ax.plot(theta_endpoints, r_endpoints, **kwargs)[0]

    def clear_artists(self, clear_preview: bool = True) -> None:
        """
        Remove Matplotlib artists from the Axes.

        Parameters
        ----------
        clear_preview : bool, optional
            If True, clear the preview from the Axes. If False, clear the
            selection from the Axes.

        Side Effects
        ------------
        - Calls ``artist.remove()`` on each artist.
        - Empties the preview_artists or selection_artists attribute in-place.
        """
        if clear_preview:
            artists = self.preview_artists
        else:
            artists = self.selection_artists

        for artist in artists:
            artist.remove()

        artists.clear()

    def shade_interior_region(self) -> None:
        """Hold space for future method."""
