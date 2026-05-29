"""Unit tests for matplotlib selection renderers."""

from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import pytest
from dta.gui import SelectionRenderer


DummyEdge = namedtuple("DummyEdge", ["endpoint1", "endpoint2"])


@pytest.fixture
def polar_ax():
    """Create a disposable polar Axes for renderer tests."""
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
    yield ax
    plt.close(fig)


def test_renderer_initializes_with_default_kwargs(polar_ax):
    """SelectionRenderer stores its Axes and default plotting options."""
    renderer = SelectionRenderer(polar_ax)

    assert renderer.ax is polar_ax
    assert renderer.preview_artists == []
    assert renderer.selection_artists == []

    assert renderer.default_preview_kwargs == {
        "color": "orange",
        "lw": 2.0,
        "zorder": 1000,
    }
    assert renderer.selection_kwargs == {
        "color": "red",
        "lw": 2.0,
        "zorder": 20,
    }


def test_renderer_merges_selection_kwargs_with_defaults(polar_ax):
    """Custom plot kwargs override selection defaults without affecting preview defaults."""
    renderer = SelectionRenderer(
        polar_ax,
        plot_kwargs={
            "color": "blue",
            "alpha": 0.5,
        },
    )

    assert renderer.selection_kwargs == {
        "color": "blue",
        "lw": 2.0,
        "zorder": 20,
        "alpha": 0.5,
    }
    assert renderer.default_preview_kwargs == {
        "color": "orange",
        "lw": 2.0,
        "zorder": 1000,
    }


def test_draw_bin_edges_draws_preview_edges_by_default(polar_ax):
    """draw_bin_edges stores preview artists when preview is left as True."""
    renderer = SelectionRenderer(polar_ax)

    edges = [
        DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi / 2)),
        DummyEdge(endpoint1=(2.0, np.pi), endpoint2=(3.0, np.pi)),
    ]

    renderer.draw_bin_edges(edges)

    assert len(renderer.preview_artists) == 2
    assert renderer.selection_artists == []
    assert len(polar_ax.lines) == 2
    assert renderer.preview_artists == list(polar_ax.lines)


def test_draw_bin_edges_draws_selection_edges_when_preview_false(polar_ax):
    """draw_bin_edges stores selection artists when preview is False."""
    renderer = SelectionRenderer(polar_ax)

    edges = [
        DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi / 2)),
        DummyEdge(endpoint1=(2.0, np.pi), endpoint2=(3.0, np.pi)),
    ]

    renderer.draw_bin_edges(edges, preview=False)

    assert renderer.preview_artists == []
    assert len(renderer.selection_artists) == 2
    assert len(polar_ax.lines) == 2
    assert renderer.selection_artists == list(polar_ax.lines)


def test_draw_bin_edges_uses_theta_as_x_and_radius_as_y(polar_ax):
    """draw_bin_edges converts BinEdge-style endpoints into polar plot coordinates."""
    renderer = SelectionRenderer(polar_ax)

    edge = DummyEdge(
        endpoint1=(1.0, 0.25 * np.pi),
        endpoint2=(2.0, 0.75 * np.pi),
    )

    renderer.draw_bin_edges([edge])
    artist = renderer.preview_artists[0]

    np.testing.assert_allclose(
        artist.get_xdata(),
        [0.25 * np.pi, 0.75 * np.pi],
    )
    np.testing.assert_allclose(
        artist.get_ydata(),
        [1.0, 2.0],
    )


def test_draw_bin_edges_uses_preview_kwargs_for_preview_edges(polar_ax):
    """Preview edges use the renderer's preview plotting kwargs."""
    renderer = SelectionRenderer(polar_ax)

    edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))

    renderer.draw_bin_edges([edge])
    artist = renderer.preview_artists[0]

    assert artist.get_color() == "orange"
    assert artist.get_linewidth() == 2.0
    assert artist.get_zorder() == 1000


def test_draw_bin_edges_uses_selection_kwargs_for_selection_edges(polar_ax):
    """Selection edges use the renderer's selection plotting kwargs."""
    renderer = SelectionRenderer(
        polar_ax,
        plot_kwargs={
            "color": "green",
            "lw": 4.0,
            "zorder": 30,
        },
    )

    edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))

    renderer.draw_bin_edges([edge], preview=False)
    artist = renderer.selection_artists[0]

    assert artist.get_color() == "green"
    assert artist.get_linewidth() == 4.0
    assert artist.get_zorder() == 30


def test_draw_bin_edges_clears_existing_preview_before_redrawing(polar_ax):
    """Drawing new edges removes stale preview artists before storing new ones."""
    renderer = SelectionRenderer(polar_ax)

    first_edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))
    second_edge = DummyEdge(endpoint1=(2.0, 0.0), endpoint2=(2.0, np.pi))

    renderer.draw_bin_edges([first_edge])
    old_artist = renderer.preview_artists[0]

    renderer.draw_bin_edges([second_edge])

    assert old_artist not in polar_ax.lines
    assert len(renderer.preview_artists) == 1
    assert len(polar_ax.lines) == 1

    np.testing.assert_allclose(
        renderer.preview_artists[0].get_ydata(),
        [2.0, 2.0],
    )


def test_draw_bin_edges_clears_existing_selection_before_redrawing(polar_ax):
    """Drawing new selection edges removes stale selection artists."""
    renderer = SelectionRenderer(polar_ax)

    first_edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))
    second_edge = DummyEdge(endpoint1=(2.0, 0.0), endpoint2=(2.0, np.pi))

    renderer.draw_bin_edges([first_edge], preview=False)
    old_artist = renderer.selection_artists[0]

    renderer.draw_bin_edges([second_edge], preview=False)

    assert old_artist not in polar_ax.lines
    assert len(renderer.selection_artists) == 1
    assert len(polar_ax.lines) == 1

    np.testing.assert_allclose(
        renderer.selection_artists[0].get_ydata(),
        [2.0, 2.0],
    )


def test_draw_bin_edges_clears_both_preview_and_selection_artists(polar_ax):
    """Drawing either preview or selection clears both existing artist groups."""
    renderer = SelectionRenderer(polar_ax)

    preview_edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))
    selection_edge = DummyEdge(endpoint1=(2.0, 0.0), endpoint2=(2.0, np.pi))

    renderer.draw_bin_edges([preview_edge])
    renderer.draw_bin_edges([selection_edge], preview=False)

    assert renderer.preview_artists == []
    assert len(renderer.selection_artists) == 1
    assert len(polar_ax.lines) == 1


def test_draw_bin_edges_with_empty_iterable_clears_existing_artists(polar_ax):
    """Drawing an empty edge list clears stale artists and draws nothing."""
    renderer = SelectionRenderer(polar_ax)

    edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))

    renderer.draw_bin_edges([edge])
    assert len(polar_ax.lines) == 1

    renderer.draw_bin_edges([])

    assert renderer.preview_artists == []
    assert renderer.selection_artists == []
    assert len(polar_ax.lines) == 0


def test_clear_artists_removes_preview_artists_only(polar_ax):
    """clear_artists removes preview artists when clear_preview is True."""
    renderer = SelectionRenderer(polar_ax)

    edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))
    renderer.draw_bin_edges([edge])

    renderer.clear_artists(clear_preview=True)

    assert renderer.preview_artists == []
    assert renderer.selection_artists == []
    assert len(polar_ax.lines) == 0


def test_clear_artists_removes_selection_artists_only(polar_ax):
    """clear_artists removes selection artists when clear_preview is False."""
    renderer = SelectionRenderer(polar_ax)

    edge = DummyEdge(endpoint1=(1.0, 0.0), endpoint2=(1.0, np.pi))
    renderer.draw_bin_edges([edge], preview=False)

    renderer.clear_artists(clear_preview=False)

    assert renderer.preview_artists == []
    assert renderer.selection_artists == []
    assert len(polar_ax.lines) == 0
