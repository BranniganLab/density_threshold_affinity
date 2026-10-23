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


@pytest.mark.parametrize("first_bin", [0, 1, 3])
@pytest.mark.parametrize("original_r_min", [0, 10])
def test_first_radial_bin_preserves_grid_data_and_retained_colors(first_bin, original_r_min):
    """Only the plotting copy is excluded; retained zero/positive colors agree."""
    grid = PolarBinGrid(original_r_min, original_r_min + 40, 4, 4)
    data = np.tile([0, 1, 2, 5], (4, 1))
    original = data.copy()
    original_edges = grid.r.bin_edges.copy()
    settings = HeatmapSettings(["lipid"], ["leaflet"], (4, 4), 5, grid)
    settings.colormap = make_custom_colormap()
    reference_fig, reference_ax = plt.subplots(subplot_kw={"projection": "polar"})
    fig = None
    try:
        reference_mesh = plot_heatmap(reference_ax, data, settings)
        reference_fig.canvas.draw()
        settings.first_radial_bin = first_bin
        fig, axes = make_density_enrichment_heatmap([data], [np.array([5, 30])], settings)
        fig.canvas.draw()
        mesh = axes[0].collections[-1]
        cutoff = grid.r.bin_edges[first_bin]
        assert settings.r_min == cutoff
        np.testing.assert_array_equal(data, original)
        np.testing.assert_array_equal(grid.r.bin_edges, original_edges)
        np.testing.assert_array_equal(mesh.get_coordinates(), reference_mesh.get_coordinates())
        np.testing.assert_array_equal(mesh.get_facecolors()[first_bin * 4:],
                                      reference_mesh.get_facecolors()[first_bin * 4:])
        plotted = mesh.get_array().reshape(data.shape)
        assert np.all(np.ma.getmaskarray(plotted)[:first_bin])
        assert np.all(np.isnan(plotted.data[:first_bin]))
        np.testing.assert_array_equal(plotted.data[first_bin:], data[first_bin:])
        keys = [ax.get_label() for ax in mesh.colorbar.ax.child_axes]
        assert keys == ["zero_density_key"]
        if cutoff > 0:
            disk = axes[0].patches[0]
            assert disk.get_height() == cutoff
            assert disk.get_hatch() is None
            np.testing.assert_array_equal(disk.get_facecolor(), [1, 1, 1, 1])
            boundary = axes[0].lines[-1]
            np.testing.assert_array_equal(boundary.get_ydata(), cutoff)
            np.testing.assert_allclose(boundary.get_xdata()[[0, -1]], [0, 2 * np.pi])
            outer_border = axes[0].spines["polar"]
            assert boundary.get_color() == outer_border.get_edgecolor()
            assert boundary.get_linewidth() == outer_border.get_linewidth()
        else:
            assert not axes[0].patches
    finally:
        plt.close(reference_fig)
        if fig is not None:
            plt.close(fig)


@pytest.mark.parametrize("index, error", [
    (-1, ValueError), (4, ValueError), (1.5, TypeError),
    (True, TypeError), ("1", TypeError),
])
def test_first_radial_bin_rejects_invalid_indices(index, error):
    settings = HeatmapSettings(["lipid"], ["leaflet"], (4, 4), 5,
                               PolarBinGrid(0, 40, 4, 4), first_radial_bin=index)
    settings.colormap = make_custom_colormap()
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    try:
        with pytest.raises(error, match="first_radial_bin"):
            plot_heatmap(ax, np.ones((4, 4)), settings)
    finally:
        plt.close(fig)


@pytest.mark.parametrize("colorbychain", [False, True])
@pytest.mark.parametrize("r_min", [0, 20, 40])
def test_helix_cutoff_preserves_colors_coordinates_and_boundary(colorbychain, r_min):
    from matplotlib.colors import to_rgba
    from dta.plotting import plot_helices

    helices = np.array([[5., 0., 20., 90., 30., 180.],
                        [30., 45., 5., 135., 20., 225.]])
    original = helices.copy()
    palette = ["purple", "blue", "green"]
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    try:
        plot_helices(helices, colorbychain, ax, colorlist=palette, r_min=r_min)
        for chain, (coordinates, markers) in enumerate(zip(helices, ax.collections)):
            visible = coordinates[::2] >= r_min
            expected = np.column_stack((np.deg2rad(coordinates[1::2][visible]),
                                        coordinates[::2][visible]))
            np.testing.assert_array_equal(markers.get_offsets(), expected)
            if np.any(visible):
                expected_colors = ([to_rgba(palette[chain])] if colorbychain else
                                   [to_rgba(color) for color, keep in zip(palette, visible) if keep])
                np.testing.assert_array_equal(markers.get_facecolors(), expected_colors)
                expected_edges = np.array(expected_colors)
                expected_edges[:, :3] *= 0.5
                np.testing.assert_array_equal(markers.get_edgecolors(), expected_edges)
        np.testing.assert_array_equal(helices, original)
    finally:
        plt.close(fig)


def test_heatmap_uses_first_radial_bin_to_hide_helix_markers():
    settings = HeatmapSettings(["lipid"], ["leaflet"], (4, 4), 5,
                               PolarBinGrid(0, 40, 4, 4), first_radial_bin=2)
    settings.colormap = make_custom_colormap()
    helices = np.array([5., 0., 20., 90., 30., 180.])
    original = helices.copy()
    fig, axes = make_density_enrichment_heatmap([np.ones((4, 4))], [helices], settings)
    try:
        markers = axes[0].collections[0]
        np.testing.assert_array_equal(markers.get_offsets()[:, 1], [20, 30])
        np.testing.assert_array_equal(helices, original)
    finally:
        plt.close(fig)
