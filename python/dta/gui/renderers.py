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

    def __init__(self, ax: matplotlib.axes.Axes) -> None:
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
        self.default_selection_kwargs = {
            "color": "red",
            "lw": 2.0,
            "zorder": 20
        }
        self.default_preview_kwargs = {
            "color": "orange",
            "lw": 2.0,
            "zorder": 1000
        }
        self.preview_artists = []
        self.selection_artists = []

    def draw_edges(
        self,
        edges: Iterable[BinEdge],
        draw_preview: bool = True,
        plot_kwargs: dict = None
    ) -> None:
        """
        Draw a collection of bin edges.

        Parameters
        ----------
        edges : iterable of BinEdge
            Edges to draw.
        draw_preview : bool, optional
            If True, draw an orange selection preview. If False, use plot_kwargs
            and draw actual selection.
        plot_kwargs : dict
            Keyword arguments passed to ``Axes.plot``.

        Side Effects
        ------------
        Saves drawn artists to the correct SelectionRenderer attribute.
        """
        if draw_preview:
            kwargs = self.default_preview_kwargs
            artists = self.preview_artists
        else:
            kwargs = self.default_selection_kwargs
            artists = self.selection_artists

        # update/override default if possible.
        if plot_kwargs is not None:
            kwargs.update(plot_kwargs)

        # plot each BinEdge individually
        for edge in edges:
            theta_endpoints = (edge.endpoint1[1], edge.endpoint2[1])
            r_endpoints = (edge.endpoint1[0], edge.endpoint2[0])
            artists.append(self.ax.plot(theta_endpoints, r_endpoints, **kwargs)[0])

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
        - Empties the hover_artists or selection_artists attribute in-place.
        """
        if clear_preview:
            artists = self.hover_artists
        else:
            artists = self.selection_artists

        for artist in artists:
            artist.remove()

        artists.clear()

    def shade_interior_region(self) -> None:
        """Hold space for future method."""
