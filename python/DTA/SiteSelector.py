import numpy as np
import matplotlib.pyplot as plt


class SiteSelector:
    """
    Interactive bin selector for polar plots.

    Usage:
        fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
        ax.pcolormesh(theta_edges, r_edges, data)
        selector = PolarBinSelector(ax, theta_edges, r_edges)
        plt.show()

        # Access selected bins
        selected = selector.get_selected_bins()
    """

    def __init__(self, ax, theta_edges, r_edges,
                 hover_color='orange', selected_color='red',
                 hover_lw=1.5, selected_lw=2):
        """
        Initialize the polar bin selector.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Polar axes to attach selector to
        theta_edges : array-like
            Angular bin edges (radians)
        r_edges : array-like
            Radial bin edges
        hover_color : str, optional
            Color for hover preview
        selected_color : str, optional
            Color for confirmed selection
        hover_lw : float, optional
            Line width for hover preview
        selected_lw : float, optional
            Line width for confirmed selection
        """
        self.ax = ax
        self.theta_edges = np.asarray(theta_edges)
        self.r_edges = np.asarray(r_edges)
        self.hover_color = hover_color
        self.selected_color = selected_color
        self.hover_lw = hover_lw
        self.selected_lw = selected_lw

        # State
        self.selected_bins = set()
        self.hover_artists = []
        self.selected_artists = []
        self.drag_start = None
        self.mode = "replace"

        # Connect events
        self.fig = ax.figure
        self._connections = [
            self.fig.canvas.mpl_connect("button_press_event", self._on_press),
            self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion),
            self.fig.canvas.mpl_connect("button_release_event", self._on_release)
        ]

    def disconnect(self):
        """Disconnect all event handlers."""
        for cid in self._connections:
            self.fig.canvas.mpl_disconnect(cid)
        self.clear_selection()

    def get_selected_bins(self):
        """
        Get the currently selected bins.

        Returns
        -------
        list of tuple
            List of (theta_idx, r_idx) tuples
        """
        return list(self.selected_bins)

    def clear_selection(self):
        """Clear all selected bins and remove artists."""
        self.selected_bins.clear()
        self._clear_artists(self.selected_artists)
        self._clear_artists(self.hover_artists)
        self.fig.canvas.draw_idle()

    def set_selection(self, bins):
        """
        Programmatically set the selected bins.

        Parameters
        ----------
        bins : iterable of tuple
            Bins to select as (theta_idx, r_idx) tuples
        """
        self.selected_bins = set(bins)
        self._update_selected_display()

    def _clear_artists(self, artists):
        """Remove artists from plot and clear list."""
        for a in artists:
            a.remove()
        artists.clear()

    def _theta_in_bin(self, t0, t1, t_low, t_high):
        dt = (t1 - t0) % (2 * np.pi)
        if dt <= np.pi:
            theta_val = t0
        else:
            dt = (t0 - t1) % (2 * np.pi)
            theta_val = t1
        condition1 = (t_low - theta_val) % (2 * np.pi)
        condition2 = (t_high - theta_val) % (2 * np.pi)
        condition3 = (theta_val - t_low) % (2 * np.pi)
        return ((condition1 < dt) or (condition2 < dt) or (condition3 < (t_high - t_low)))

    def _bins_in_region(self, t0, t1, r0, r1):
        """Find all bins within a rectangular selection region."""
        r0, r1 = sorted([r0, r1])
        t0 = t0 % (2 * np.pi)
        t1 = t1 % (2 * np.pi)

        bins = []
        for ti in range(len(self.theta_edges) - 1):
            t_low = self.theta_edges[ti]
            t_high = self.theta_edges[ti + 1]
            if self._theta_in_bin(t0, t1, t_low, t_high):
                for ri in range(len(self.r_edges) - 1):
                    r_low = self.r_edges[ri]
                    r_high = self.r_edges[ri + 1]
                    if r_high >= r0 and r_low <= r1:
                        bins.append((ti, ri))
        return bins

    def _draw_outer_edges(self, mask_bins, color, lw):
        """Draw outline around selected bins."""
        if not mask_bins:
            return []

        n_r = len(self.r_edges) - 1
        n_t = len(self.theta_edges) - 1
        mask = np.zeros((n_r, n_t), dtype=bool)
        for ti, ri in mask_bins:
            mask[ri, ti] = True

        artists = []
        for ri in range(n_r):
            for ti in range(n_t):
                if not mask[ri, ti]:
                    continue

                # Check neighbors
                top = ri == n_r - 1 or not mask[ri + 1, ti]
                bottom = ri == 0 or not mask[ri - 1, ti]
                left = not mask[ri, (ti - 1) % n_t]
                right = not mask[ri, (ti + 1) % n_t]

                t0 = self.theta_edges[ti]
                t1 = self.theta_edges[(ti + 1) % n_t]
                r0 = self.r_edges[ri]
                r1 = self.r_edges[ri + 1]

                if top:
                    artists.append(self.ax.plot([t0, t1], [r1, r1],
                                               color=color, lw=lw, zorder=10)[0])
                if bottom:
                    artists.append(self.ax.plot([t0, t1], [r0, r0],
                                               color=color, lw=lw, zorder=10)[0])
                if left:
                    artists.append(self.ax.plot([t0, t0], [r0, r1],
                                               color=color, lw=lw, zorder=10)[0])
                if right:
                    artists.append(self.ax.plot([t1, t1], [r0, r1],
                                               color=color, lw=lw, zorder=10)[0])
        return artists

    def _get_bins_from_drag(self, x0, y0, x1, y1):
        """Get bins from drag coordinates."""
        dx = x1 - x0
        dy = y1 - y0
        threshold = 1e-8

        if abs(dx) < threshold and abs(dy) < threshold:
            # Single bin click
            ti = np.searchsorted(self.theta_edges, x1 % (2 * np.pi), side="right") - 1
            ri = np.searchsorted(self.r_edges, y1, side="right") - 1
            if 0 <= ti < len(self.theta_edges) - 1 and 0 <= ri < len(self.r_edges) - 1:
                return [(ti, ri)]
            return []
        else:
            # Rectangle selection
            return self._bins_in_region(x0, x1, y0, y1)

    def _update_hover_display(self, bins):
        """Update hover preview."""
        self._clear_artists(self.hover_artists)

        if self.mode == "replace":
            temp_mask = set(bins)
        elif self.mode == "add":
            temp_mask = self.selected_bins.union(bins)
        elif self.mode == "subtract":
            temp_mask = self.selected_bins.difference(bins)
        else:
            temp_mask = self.selected_bins.union(bins)

        self.hover_artists.extend(
            self._draw_outer_edges(temp_mask, self.hover_color, self.hover_lw)
        )
        self.fig.canvas.draw_idle()

    def _update_selected_display(self):
        """Update confirmed selection display."""
        self._clear_artists(self.selected_artists)
        self.selected_artists.extend(
            self._draw_outer_edges(self.selected_bins,
                                  self.selected_color, self.selected_lw)
        )
        self.fig.canvas.draw_idle()

    def _on_press(self, event):
        """Handle mouse press event."""
        if event.inaxes != self.ax or event.button != 1:
            return

        self.drag_start = (event.xdata, event.ydata)

        if event.key == "shift":
            self.mode = "add"
        elif event.key == "control":
            self.mode = "subtract"
        else:
            self.mode = "replace"
            self.selected_bins.clear()

    def _on_motion(self, event):
        """Handle mouse motion event."""
        if self.drag_start is None or event.inaxes != self.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        bins = self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            event.xdata, event.ydata
        )
        self._update_hover_display(bins)

    def _on_release(self, event):
        """Handle mouse release event."""
        if self.drag_start is None:
            return
        if event.xdata is None or event.ydata is None:
            self.drag_start = None
            return

        bins = self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            event.xdata, event.ydata
        )

        if self.mode == "replace":
            self.selected_bins.update(bins)
        elif self.mode == "add":
            self.selected_bins.update(bins)
        elif self.mode == "subtract":
            self.selected_bins.difference_update(bins)

        self._update_selected_display()
        self._clear_artists(self.hover_artists)
        self.drag_start = None


# Example usage
def create_example_plot():
    """Create an example plot with the polar bin selector."""
    # Data setup
    theta_edges = np.linspace(0, 2 * np.pi, 73)
    r_edges = np.linspace(0, 5, 41)

    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    Theta, R = np.meshgrid(theta_centers, r_centers)

    Z = np.sin(2 * Theta) * np.exp(-R / 2)

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    ax.pcolormesh(theta_edges, r_edges, Z, shading="auto")
    ax.set_title("Hover = preview | Release = confirm\nShift=add, Ctrl=subtract")
    ax.set_axis_off()

    # Create selector
    selector = SiteSelector(ax, theta_edges, r_edges)

    return fig, ax, selector
