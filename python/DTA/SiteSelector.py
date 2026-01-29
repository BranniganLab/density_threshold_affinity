from DTA.utils import unwrap_theta
from DTA.polar_bin_classes import PolarBinGrid, PolarBinRenderer
import numpy as np
import matplotlib.pyplot as plt


class SiteSelector:
    """
    Interactive polar bin selector.

    Conventions:
        - Bin indices: (ri, ti)
        - Coordinates: (r, theta)
    """

    def __init__(self, ax, theta_edges, r_edges):
        self.ax = ax
        self.fig = ax.figure

        self.grid = PolarBinGrid(theta_edges, r_edges)
        self.renderer = PolarBinRenderer(ax)

        # Selection state
        self.selected_bins = set()
        self.hover_artists = []
        self.selected_artists = []

        # Drag state
        self.drag_start = None   # (r, theta)
        self._last_theta = None
        self.mode = "replace"

        # Event wiring
        self._connections = [
            self.fig.canvas.mpl_connect("button_press_event", self._on_press),
            self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion),
            self.fig.canvas.mpl_connect("button_release_event", self._on_release),
        ]

        # UI cleanup
        self.fig.canvas.toolbar_visible = False
        self.fig.canvas.header_visible = False
        self.fig.canvas.footer_visible = False
        self.fig.canvas.resizable = False

    # ---------- public API ----------

    def get_selected_bins(self):
        return sorted(self.selected_bins)

    def clear_selection(self):
        self.selected_bins.clear()
        self._clear_artists(self.selected_artists)
        self._clear_artists(self.hover_artists)
        self.fig.canvas.draw_idle()

    # ---------- drawing helpers ----------

    def _clear_artists(self, artists):
        for a in artists:
            a.remove()
        artists.clear()

    def _draw_selection(self, bins, color, lw):
        edges = self.grid.exposed_edges(bins)
        return self.renderer.draw_edges(edges, color=color, lw=lw)

    # ---------- selection logic ----------

    def _get_bins_from_drag(self, r0, theta0, r1, theta1):
        dt = theta1 - theta0
        dr = r1 - r0
        threshold = 1e-8

        if abs(dt) < threshold and abs(dr) < threshold:
            bin_idx = self.grid.bin_at(r1, theta1)
            return [bin_idx] if bin_idx else []

        return self.grid.bins_in_region(r0, theta0, r1, theta1)

    def _update_hover(self, bins):
        self._clear_artists(self.hover_artists)
        bins = set(bins)

        if self.mode == "replace":
            temp = bins
        elif self.mode == "add":
            temp = self.selected_bins | bins
        else:  # subtract
            temp = self.selected_bins - bins

        self.hover_artists.extend(
            self._draw_selection(temp, color="orange", lw=1.5)
        )
        self.fig.canvas.draw_idle()

    def _update_selected(self):
        self._clear_artists(self.selected_artists)
        self.selected_artists.extend(
            self._draw_selection(self.selected_bins, color="red", lw=2)
        )
        self.fig.canvas.draw_idle()

    # ---------- event handlers ----------

    def _on_press(self, event):
        if event.inaxes != self.ax or event.button != 1:
            return

        # matplotlib gives (theta, r)
        self.drag_start = (event.ydata, event.xdata)
        self._last_theta = event.xdata

        if event.key == "shift":
            self.mode = "add"
        elif event.key == "control":
            self.mode = "subtract"
        else:
            self.mode = "replace"
            self.selected_bins.clear()
            self._clear_artists(self.selected_artists)

    def _on_motion(self, event):
        if (
            event.button != 1
            or self.drag_start is None
            or event.inaxes != self.ax
            or event.xdata is None
            or event.ydata is None
        ):
            return

        theta = unwrap_theta(self._last_theta, event.xdata)
        self._last_theta = theta

        bins = self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            event.ydata, theta
        )
        self._update_hover(bins)

    def _on_release(self, event):
        if self.drag_start is None:
            return
        if event.xdata is None or event.ydata is None:
            self.drag_start = None
            self._last_theta = None
            return

        theta = unwrap_theta(self._last_theta, event.xdata)

        bins = set(self._get_bins_from_drag(
            self.drag_start[0], self.drag_start[1],
            event.ydata, theta
        ))

        if self.mode in ("replace", "add"):
            self.selected_bins.update(bins)
        else:
            self.selected_bins.difference_update(bins)

        self._update_selected()
        self._clear_artists(self.hover_artists)

        self.drag_start = None
        self._last_theta = None

    def disconnect(self):
        for cid in self._connections:
            self.fig.canvas.mpl_disconnect(cid)
        self.clear_selection()


# Example usage
def create_example_plot():
    """Create an example plot with the polar bin selector."""
    # Data setup
    theta_edges = np.linspace(0, 2 * np.pi, 73)
    r_edges = np.linspace(0, 5, 41)

    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    theta_vals, r_vals = np.meshgrid(theta_centers, r_centers)
    example = np.sin(2 * theta_vals) * np.exp(-r_vals / 2)

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, layout='constrained')
    ax.pcolormesh(theta_edges, r_edges, example, shading="auto")
    ax.set_title("Click + Drag = preview | Release = confirm\nShift = add | Ctrl = subtract")
    ax.set_axis_off()

    # Create selector
    selector = SiteSelector(ax, theta_edges, r_edges)

    return fig, ax, selector
