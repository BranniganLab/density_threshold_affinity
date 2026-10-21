"""Regression checks for heatmap colors and their colorbar."""

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import ListedColormap

from dta.bin_logic import PolarBinGrid
from dta.plotting import (
    HeatmapSettings,
    make_custom_colormap,
    make_density_enrichment_heatmap,
    plot_heatmap,
)


@pytest.mark.parametrize("maximum", [2, 5, 10])
@pytest.mark.parametrize("custom_colormap", [False, True])
def test_zero_is_dark_gray_and_positive_colors_are_unchanged(maximum, custom_colormap):
    """Even tiny positive values retain the previous palette and scaling."""
    data = np.array([[0, 1e-12, 0.5 / maximum, 1 / maximum,
                      0.75, 1, 2, maximum, 2 * maximum]])
    original = data.copy()
    settings = HeatmapSettings(["lipid"], ["leaflet"], (4, 4), maximum,
                               PolarBinGrid(0, 1, 1, data.size))
    settings.colormap = plt.colormaps['viridis'] if custom_colormap else make_custom_colormap()
    original_bad_color = settings.colormap.get_bad().copy()

    # Reference the palette and normalization used before issue 338.
    previous_palette = ListedColormap(np.concatenate([
        plt.colormaps['RdGy_r'](np.linspace(0.35, 0.5, 128)),
        plt.colormaps['bwr'](np.linspace(0.5, 1, 128)),
    ]))
    if custom_colormap:
        previous_palette = settings.colormap
    previous_colors = previous_palette(
        np.interp(data.ravel(), [1 / maximum, 1, maximum], [0, 0.5, 1])
    )

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    try:
        mesh = plot_heatmap(ax, data, settings)
        fig.canvas.draw()
        colors = mesh.get_facecolors()
        np.testing.assert_array_equal(colors[0], [0.4, 0.4, 0.4, 1])
        np.testing.assert_array_equal(colors[1:], previous_colors[1:])
        np.testing.assert_array_equal(data, original)
        np.testing.assert_array_equal(settings.colormap.get_bad(), original_bad_color)
    finally:
        plt.close(fig)


@pytest.mark.parametrize("data", [
    np.zeros((1, 4)),
    np.ma.array([[0, 0.1, np.nan, 1]], mask=[[False, True, False, True]]),
])
def test_zero_and_missing_bins_are_opaque_dark_gray(data):
    """Empty panels and pre-existing masks render without losing their masks."""
    settings = HeatmapSettings(["lipid"], ["leaflet"], (4, 4), 5,
                               PolarBinGrid(0, 1, 1, 4))
    settings.colormap = make_custom_colormap()
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    try:
        mesh = plot_heatmap(ax, data, settings)
        fig.canvas.draw()
        np.testing.assert_array_equal(mesh.get_facecolors(), np.tile([0.4, 0.4, 0.4, 1], (4, 1)))
    finally:
        plt.close(fig)


def test_multi_panel_colorbar_uses_heatmap_scale_with_helix_overlays():
    """The colorbar describes the density mesh rather than the helix scatter."""
    maximum = 3
    settings = HeatmapSettings(["lipid"], ["outer", "inner"], (4, 8), maximum,
                               PolarBinGrid(0, 1, 1, 4))
    settings.colormap = make_custom_colormap()
    data = np.array([[0, 0.1, 1, maximum]])
    helices = [np.array([0.5, 30]), np.array([0.5, 30])]
    fig, axes = make_density_enrichment_heatmap([data, data], helices, settings)
    try:
        fig.canvas.draw()
        mesh = axes[-1].collections[-1]
        colorbar = mesh.colorbar
        assert len(axes) == 2
        assert colorbar.mappable is mesh
        assert colorbar.norm is mesh.norm
        assert colorbar.cmap is mesh.cmap
        np.testing.assert_allclose(colorbar.get_ticks(), [1 / maximum, 1, maximum])
        np.testing.assert_allclose(colorbar.norm(colorbar.get_ticks()), [0, 0.5, 1])
        zero_key = next(ax for ax in colorbar.ax.child_axes if ax.get_label() == "zero_density_key")
        np.testing.assert_array_equal(zero_key.get_facecolor(), mesh.cmap.get_bad())
    finally:
        plt.close(fig)
