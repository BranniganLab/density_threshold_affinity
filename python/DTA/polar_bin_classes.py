from dataclasses import dataclass
import numpy as np
from DTA.utils import bin_in_theta_arc


@dataclass(frozen=True)
class BinEdge:
    r: tuple
    theta: tuple


class PolarBinGrid:
    """
    Polar bin topology and geometry.

    Conventions:
        - Bin indices: (ri, ti)
        - Coordinates: (r, theta)
    """

    def __init__(self, theta_edges, r_edges):
        self.theta_edges = np.asarray(theta_edges)
        self.r_edges = np.asarray(r_edges)
        self.n_t = len(theta_edges) - 1
        self.n_r = len(r_edges) - 1

    # ---------- bin lookup ----------

    def bin_at(self, r, theta):
        ti = np.searchsorted(
            self.theta_edges, theta % (2 * np.pi), side="right"
        ) - 1
        ri = np.searchsorted(self.r_edges, r, side="right") - 1

        if 0 <= ri < self.n_r and 0 <= ti < self.n_t:
            return ri, ti
        return None

    def bins_in_region(self, r_start, theta_start, r_end, theta_end):
        r0, r1 = sorted([r_start, r_end])
        bins = []

        for ti in range(self.n_t):
            theta_bin_low = self.theta_edges[ti]
            theta_bin_high = self.theta_edges[ti + 1]
            if bin_in_theta_arc(theta_start, theta_end, theta_bin_low, theta_bin_high):
                for ri in range(self.n_r):
                    r_low = self.r_edges[ri]
                    r_high = self.r_edges[ri + 1]
                    if r_high >= r0 and r_low <= r1:
                        bins.append((ri, ti))

        return bins

    # ---------- masks & adjacency ----------

    def build_mask(self, bins):
        mask = np.zeros((self.n_r, self.n_t), dtype=bool)
        for ri, ti in bins:
            mask[ri, ti] = True
        return mask

    def exposed_edges(self, selected_bins):
        if not selected_bins:
            return []

        mask = self.build_mask(selected_bins)
        edges = []

        for ri, ti in zip(*np.where(mask)):
            edges.extend(self._exposed_edges_for_bin(mask, ri, ti))

        return edges

    def _exposed_edges_for_bin(self, mask, ri, ti):
        edges = []
        n_r, n_t = mask.shape

        if ri == n_r - 1 or not mask[ri + 1, ti]:
            edges.append(self._edge_geometry(ri, ti, "top"))
        if ri == 0 or not mask[ri - 1, ti]:
            edges.append(self._edge_geometry(ri, ti, "bottom"))
        if not mask[ri, (ti - 1) % n_t]:
            edges.append(self._edge_geometry(ri, ti, "left"))
        if not mask[ri, (ti + 1) % n_t]:
            edges.append(self._edge_geometry(ri, ti, "right"))

        return edges

    # ---------- geometry ----------

    def _edge_geometry(self, ri, ti, side):
        t0 = self.theta_edges[ti]
        t1 = self.theta_edges[(ti + 1) % self.n_t]
        r0 = self.r_edges[ri]
        r1 = self.r_edges[ri + 1]

        if side == "top":
            return BinEdge((r1, r1), (t0, t1))
        if side == "bottom":
            return BinEdge((r0, r0), (t0, t1))
        if side == "left":
            return BinEdge((r0, r1), (t0, t0))
        if side == "right":
            return BinEdge((r0, r1), (t1, t1))

        raise ValueError(f"Unknown edge side: {side}")


class PolarBinRenderer:
    """
    Draws bin edge geometry on a polar axis.
    """

    def __init__(self, ax):
        self.ax = ax

    def draw_edges(self, edges, **plot_args):
        artists = []
        for edge in edges:
            artists.append(
                self.ax.plot(
                    edge.theta, edge.r, **plot_args, zorder=10
                )[0]
            )
        return artists
