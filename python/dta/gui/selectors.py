#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interactive polar-bin selection controllers.

This module provides a per-Axes selector (`SiteSelector`) and a figure-level
event router (`SiteSelectorManager`) for interactive selection of polar bins
drawn on Matplotlib polar Axes.

The intended usage pattern is:

1. Create one `SiteSelector` per Axes.
2. Register all selectors with a single `SiteSelectorManager` tied to the Figure.
3. The manager routes Matplotlib mouse events to the correct selector and
   ensures drag gestures are consistently owned by the selector where the drag
   began.

The implementation is designed to work in interactive backends (including
Jupyter widget backends) where modifier keys may not be reported consistently.
"""

from dta.bin_logic import PolarBinGrid, BinSelection
from dta.bin_logic.utils import unwrap_theta, Coordinate
from .selector_state import SelectionOperation, SelectorDragState
from .renderers import SelectionRenderer


class SiteSelectorManager:
    """
    Figure-level event router for multiple `SiteSelector` instances.

    The manager connects to a Matplotlib Figure's mouse events and forwards them
    to the correct selector. During a drag gesture, the selector where the drag
    began is treated as the "drag owner" and continues to receive motion/release
    events even if the mouse leaves its Axes.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure whose canvas events will be monitored and routed.
    """

    def __init__(self, fig):
        """
        Construct a manager and connect to Figure canvas events.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure whose canvas events will be monitored and routed.

        Side Effects
        ------------
        Registers Matplotlib callbacks for press, motion, and release events.
        """
        self.fig = fig
        self._selectors = {}
        self._active = {}
        self._drag_owner = None

        self._cids = [
            fig.canvas.mpl_connect("button_press_event", self._on_press_event),
            fig.canvas.mpl_connect("motion_notify_event", self._on_motion_event),
            fig.canvas.mpl_connect("button_release_event", self._on_release_event),
        ]

    def register(self, selector, *, active=False):
        """
        Register a selector with the manager.

        Parameters
        ----------
        selector : SiteSelector
            Selector instance to register.
        active : bool
            If True, activates this selector for its Axes immediately. If False,
            the selector may still become active if no selector is currently
            active for that Axes.

        Notes
        -----
        Multiple selectors may be registered for the same Axes, but only one
        selector per Axes is considered active at a time.
        """
        ax = selector.renderer.ax
        self._selectors.setdefault(ax, []).append(selector)
        if active or ax not in self._active:
            self.set_active(selector)

    def set_active(self, selector):
        """
        Set the active selector for an Axes.

        Parameters
        ----------
        selector : SiteSelector
            Selector to activate.

        Side Effects
        ------------
        - Deactivates the previously active selector for the same Axes (if any).
        - Activates the new selector by calling its :meth:`SiteSelector.on_activate`.
        """
        ax = selector.renderer.ax
        current = self._active.get(ax)
        if current is selector:
            return
        if current:
            current.on_deactivate()
        self._active[ax] = selector
        selector.on_activate()

    # ------------------------------------------------------------------
    # Matplotlib callbacks (one per event type)
    # ------------------------------------------------------------------

    def _on_press_event(self, event):
        """Matplotlib callback: mouse press."""
        ax = getattr(event, "inaxes", None)
        selector = self._active.get(ax)
        if selector is None:
            return

        # The selector that receives the press becomes the drag owner until release.
        self._drag_owner = selector

        # Latch modifiers for the duration of the gesture.
        selector.drag_tracker.mods = frozenset(self._mods_from_mouse_event(event))

        updated = bool(selector.on_press(event))
        if updated:
            self.fig.canvas.draw_idle()

    def _on_motion_event(self, event):
        """Matplotlib callback: mouse motion."""
        if self._drag_owner is None:
            return
        updated = bool(self._drag_owner.on_motion(event))
        if updated:
            self.fig.canvas.draw_idle()

    def _on_release_event(self, event):
        """Matplotlib callback: mouse release."""
        if self._drag_owner is None:
            return

        selector = self._drag_owner
        updated = bool(selector.on_release(event))
        if updated:
            self.fig.canvas.draw_idle()

        self._drag_owner = None

    def _mods_from_mouse_event(self, event) -> set[str]:
        """
        Query the MouseEvent to see if any keyboard buttons have been pressed.

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse event potentially carrying modifier state.

        Returns
        -------
        set[str]
            Modifier key names present at the time of the event. Keys include
            ``"shift"`` and ``"control"``.

        Notes
        -----
        On some interactive backends (e.g., Jupyter widget backends),
        ``event.guiEvent`` may be a dict-like object that does not always include
        modifier flags. This method attempts to read modifier flags from
        ``event.guiEvent`` first and falls back to parsing ``event.key``.
        """
        mods: set[str] = set()

        ge = getattr(event, "guiEvent", None)

        # Fast path: only attempt guiEvent if it looks like it actually has modifier fields.
        if isinstance(ge, dict) and (
            "shiftKey" in ge or "ctrlKey" in ge or "metaKey" in ge
        ):
            if ge.get("shiftKey", False):
                mods.add("shift")
            if ge.get("ctrlKey", False) or ge.get("metaKey", False):
                mods.add("control")
            if mods:
                return mods

        if ge is not None and not isinstance(ge, dict):
            # Object-style guiEvent (rare in some backends)
            try:
                if getattr(ge, "shiftKey", False):
                    mods.add("shift")
                if getattr(ge, "ctrlKey", False) or getattr(ge, "metaKey", False):
                    mods.add("control")
                if mods:
                    return mods
            except (AttributeError, TypeError):
                pass

        # Fallback: parse string representation from Matplotlib
        k = (getattr(event, "key", None) or "").lower()
        if "shift" in k:
            mods.add("shift")
        if "control" in k or "ctrl" in k or "meta" in k:
            mods.add("control")
        return mods


class SiteSelector:
    """
    Interactive controller for selecting polar bins on a single Axes.

    This class converts mouse gestures (press / drag / release) into semantic
    selection operations, updates an internal selection model, and draws both
    hover previews and committed selection outlines.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Polar axes used for interaction and drawing.
    theta_edges, r_edges : array-like
        Bin edge definitions for the polar grid.
    plot_kwargs : dict, optional
        Default Matplotlib line properties forwarded to the renderer. These
        apply to committed selection edges. Hover edges are drawn with a
        derived style.

    Notes
    -----
    The active selection operation (replace/add/subtract) is *latched* at the
    beginning of a gesture based on modifier keys and does not change until
    the gesture ends.
    """

    def __init__(self, ax, theta_edges, r_edges, plot_kwargs=None):
        """
        Construct a selector bound to a single Axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes used for interaction and drawing.
        theta_edges, r_edges : array-like
            Bin edge definitions for the polar grid.
        plot_kwargs : dict, optional
            Default Matplotlib plotting keywords for drawing committed edges.
        """
        self.grid = PolarBinGrid(theta_edges, r_edges)
        self.renderer = SelectionRenderer(ax, plot_kwargs)
        self.selection = BinSelection()
        self.drag_tracker = SelectorDragState()
        self.operation = SelectionOperation.REPLACE
        self.current_preview_bins = None

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def on_activate(self):
        """
        Activate the selector for interaction.

        This method resets transient drag/preview state, but does not modify the
        committed selection.
        """
        self.drag_tracker.drag_start = None
        self.drag_tracker.last_theta = None
        self.drag_tracker.mods = frozenset()
        self.current_preview_bins = None

    def on_deactivate(self):
        """
        Deactivate the selector and clear hover preview visuals.

        This method removes hover artists from the Axes and resets transient
        drag/preview state. It does not modify the committed selection.
        """
        self._clear_artists(self.renderer.hover_artists)
        self.drag_tracker.drag_start = None
        self.drag_tracker.last_theta = None
        self.drag_tracker.mods = frozenset()
        self.current_preview_bins = None

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------

    def on_press(self, event) -> bool:
        """
        Begin a selection gesture.

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse press event. Only presses occurring inside this selector's
            Axes are processed.

        Side Effects
        ------------
        - Initializes drag state (start location and theta history).
        - Latches the selection operation (replace/add/subtract) from the
          gesture's modifier state stored in ``drag_tracker.mods``.
        - If data coordinates are available, computes an initial preview and
          draws a hover outline.
        """
        if event.inaxes is not self.renderer.ax:
            return False

        if event.xdata is None or event.ydata is None:
            return False

        # Store drag start as (r, theta).
        click_coordinate = Coordinate(event.ydata, event.xdata)
        self.drag_tracker.drag_start = click_coordinate
        self.drag_tracker.last_theta = click_coordinate.theta_coord

        # Latch key-press modifiers.
        mods = self.drag_tracker.mods
        if "shift" in mods:
            self.operation = SelectionOperation.ADD
        elif "control" in mods:
            self.operation = SelectionOperation.SUBTRACT
        else:
            self.operation = SelectionOperation.REPLACE

        # Establish an initial preview.
        clicked_bin = self.grid.bins_in_region(
            start=click_coordinate,
            end=click_coordinate
        )
        updated_preview_bins = self._calculate_preview_bins(clicked_bin)
        self.current_preview_bins = updated_preview_bins
        self._draw_preview(updated_preview_bins)
        return True

    def on_motion(self, event) -> bool:
        """
        Update hover preview during a drag gesture.

        Parameters
        ----------
        event : matplotlib.backend_bases.MouseEvent
            Mouse motion event.

        Behavior
        --------
        - If no drag is active, this method does nothing.
        - If the cursor is outside this selector's Axes, this method does
          nothing (hover preview is frozen).
        - If data coordinates are valid, updates the preview selection based on
          the drag rectangle/sector implied by the drag start and current cursor
          location, then redraws the hover outline.
        """
        if self.drag_tracker.drag_start is None:
            return False

        # Freeze hover updates unless the cursor is inside this Axes.
        if event.inaxes is not self.renderer.ax:
            return False
        if event.xdata is None or event.ydata is None:
            return False

        current_location = Coordinate(
            event.ydata,
            unwrap_theta(self.drag_tracker.last_theta, event.xdata)
        )
        self.drag_tracker.last_theta = current_location.theta_coord

        bins = self.grid.bins_in_region(
            start=self.drag_tracker.drag_start,
            end=current_location
        )

        updated_preview_bins = self._calculate_preview_bins(bins)
        self.current_preview_bins = updated_preview_bins

        self._draw_preview(updated_preview_bins)
        return True

    def on_release(self, _event) -> bool:
        """
        Finalize a selection gesture and commit the result.

        Parameters
        ----------
        _event : matplotlib.backend_bases.MouseEvent
            Mouse release event. The event is accepted for signature
            compatibility with Matplotlib callbacks, but is not used by this
            method.

        Side Effects
        ------------
        - Commits the selection.
        - Calls :meth:`on_selection_committed` with pre/post snapshots.
        - Redraws committed selection edges and clears hover edges.
        - Resets gesture state (including latched modifiers).
        """
        if self.drag_tracker.drag_start is None:
            return False

        last_bins = frozenset(self.selection.get_bins())
        self.save_to_selection_history(last_bins)

        self.selection.set_bins(self.current_preview_bins)

        self._draw_selection()
        self._clear_artists(self.renderer.hover_artists)

        self.drag_tracker.drag_start = None
        self.drag_tracker.last_theta = None
        self.current_preview_bins = None
        self.drag_tracker.mods = frozenset()
        return True

    # ------------------------------------------------------------------
    # Selection logic
    # ------------------------------------------------------------------

    def _calculate_preview_bins(self, bins):
        """
        Compute the preview selection for the current gesture without mutation.

        Parameters
        ----------
        bins : set[tuple[int, int]]
            The bin set implied by the current drag region.

        Returns
        -------
        set[tuple[int, int]]
            The preview selection that would result if the gesture ended now.

        Notes
        -----
        The preview depends on the latched operation:

        - ``REPLACE``: preview is ``bins``.
        - ``ADD``: preview is ``current ∪ bins``.
        - ``SUBTRACT``: preview is ``current \\ bins``.
        """
        current = self.selection.get_bins()

        if self.operation is SelectionOperation.REPLACE:
            return bins
        if self.operation is SelectionOperation.ADD:
            return current | bins
        if self.operation is SelectionOperation.SUBTRACT:
            return current - bins

        return current

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def _draw_preview(self, bins):
        """
        Draw the hover preview outline for a bin set.

        Parameters
        ----------
        bins : set[tuple[int, int]]
            Bins whose boundary should be rendered as a hover preview.

        Side Effects
        ------------
        - Removes any previous hover artists.
        - Draws new hover artists on this selector's Axes.
        """
        self._clear_artists(self.renderer.hover_artists)
        edges = self.grid.exposed_edges(bins)

        hover_kwargs = {
            "color": "orange",
            "lw": 1.5,
            "zorder": self.renderer.plot_kwargs["zorder"] + 1,
        }

        self.renderer.hover_artists.extend(
            self.renderer.draw_edges(edges, hover_kwargs)
        )

    def _draw_selection(self):
        """
        Draw the committed selection outline.

        Side Effects
        ------------
        - Removes any previous committed selection artists.
        - Draws the boundary edges of the committed selection.
        """
        self._clear_artists(self.renderer.selected_artists)
        edges = self.grid.exposed_edges(self.selection.get_bins())
        self.renderer.selected_artists.extend(
            self.renderer.draw_edges(edges, self.renderer.plot_kwargs)
        )

    def _clear_artists(self, artists):
        """
        Remove a list of Matplotlib artists from the Axes.

        Parameters
        ----------
        artists : list
            Artists previously added to the Axes.

        Side Effects
        ------------
        - Calls ``artist.remove()`` on each artist.
        - Empties the input list in-place.
        """
        for artist in artists:
            artist.remove()
        artists.clear()

    # ------------------------------------------------------------------
    # Undo hook
    # ------------------------------------------------------------------

    def save_to_selection_history(self, last_bins):
        """
        Hold this space for a future implementation of 'undo' logic.

        Parameters
        ----------
        last_bins : frozenset[BinAddress]
            The committed selection state immediately before the most-recent
            gesture.
        """
